"""
    TTVCycle

Genuine multigrid V-cycle that lives ENTIRELY in tensor-train arithmetic, for the
backward-Euler CME step  (I − dt·L) x = b  on a nearest-neighbour reaction–diffusion
chain.  Graduated from `examples/tt_vcycle_skeleton.jl` (validated to machine
precision: MPO assembly 2.7e-13, V-cycle ρ ≈ 0.28–0.51 across K=2..4).

Three ingredients (Bolten/Kahl/Kressner/Macedo/Sokolović 2024, "MultigridAMEn"):

  (1) TT-GMRES SMOOTHER — matvec + inner product + axpy only, no diagonal inverse
      (a TT cannot form 1/diag).  `tt_gmres`.
  (2) FACTOR-WISE GALERKIN (their Prop. 1) — L = Σ Kronecker terms, so R·M·P with
      tensor-product transfers stays a sum of Kronecker terms: coarsen each MPO
      core's PHYSICAL (count) index independently, never touch the bond.  Coarse
      operator built in O(K), never assembled densely.  `galerkin_mpo`.
  (3) HEAVY SOLVE AT THE COARSEST LEVEL only — rank-adaptive two-site (DMRG/AMEn)
      ALS, or a many-iteration GMRES.  `tt_als` / coarse `tt_gmres`.

Coarsening is on the COUNT axis (mode size n_max+1 → halved per level), the
tensor-product transfer that keeps the MPO Kronecker structure.

High-level entry point: [`tt_be_solve`].  Build a fine MPO with [`build_M_mpo`]
(generic locals + optional pair generator) or [`bd_diffusion_chain_mpo`]
(immigration/death/diffusion convenience model).
"""
module TTVCycle

using LinearAlgebra

export build_M_mpo, bd_diffusion_chain_mpo, build_tt_hierarchy, tt_vcycle,
       tt_gmres, tt_als, tt_be_solve, tt_delta0, tt_zeros, tt_to_full,
       mpo_apply, mpo_to_dense, op_svd_decompose

# ── operator SVD: L_pair = Σ_s A_s ⊗ B_s  (small # of terms) ───────────────────
function op_svd_decompose(Lpair::Matrix{Float64}, nv::Int; tol = 1e-12)
    Q = zeros(nv*nv, nv*nv)
    for n1p in 0:nv-1, n2p in 0:nv-1, n1 in 0:nv-1, n2 in 0:nv-1
        Q[n1p*nv + n1 + 1, n2p*nv + n2 + 1] = Lpair[n1p*nv + n2p + 1, n1*nv + n2 + 1]
    end
    U, S, V = svd(Q)
    AB = Tuple{Matrix{Float64},Matrix{Float64}}[]
    for s in 1:length(S)
        S[s] > tol*S[1] || break
        # reshape fills column-major (minor index = count) ⇒ transpose to [n', n].
        A = Matrix(reshape(U[:, s] .* S[s], nv, nv)')   # A[n1', n1]
        B = Matrix(reshape(V[:, s],         nv, nv)')   # B[n2', n2]
        push!(AB, (A, B))
    end
    AB
end

# ── MPO assembly:  M = I − dt·L  as a nearest-neighbour matrix-product operator ─
# core shape [b_left, n', n, b_right]
function build_L_mpo(locals::Vector{Matrix{Float64}},
                     AB::Vector{Tuple{Matrix{Float64},Matrix{Float64}}}, nv::Int;
                     dscale::AbstractVector = ones(length(locals)))
    K = length(locals); P = length(AB); d = 2 + P
    Id = Matrix{Float64}(I, nv, nv)
    cores = Vector{Array{Float64,4}}(undef, K)
    for k in 1:K
        W = zeros(d, nv, nv, d)
        W[1, :, :, 1] = Id
        W[d, :, :, d] = Id
        W[1, :, :, d] = locals[k]                          # site-specific local term
        for s in 1:P
            W[1,   :, :, 1+s] = dscale[k] .* AB[s][1]       # emit A_s on bond k (scaled)
            W[1+s, :, :, d  ] = AB[s][2]                    # finish with B_s
        end
        if k == 1
            cores[k] = W[1:1, :, :, :]
        elseif k == K
            cores[k] = W[:, :, :, d:d]
        else
            cores[k] = W
        end
    end
    cores
end

identity_mpo(K, nv) = [reshape(Matrix{Float64}(I, nv, nv), 1, nv, nv, 1) for _ in 1:K]

mpo_scale_first(W, c) = (Z = deepcopy(W); Z[1] = c .* Z[1]; Z)

function mpo_add(A, B)                                       # block-diagonal in bond
    K = length(A); C = Vector{Array{Float64,4}}(undef, K)
    for k in 1:K
        Al, No, Ni, Ar = size(A[k]); Bl, _, _, Br = size(B[k])
        if k == 1
            W = zeros(1, No, Ni, Ar+Br)
            W[1,:,:,1:Ar] = A[k]; W[1,:,:,Ar+1:end] = B[k]
        elseif k == K
            W = zeros(Al+Bl, No, Ni, 1)
            W[1:Al,:,:,1] = A[k]; W[Al+1:end,:,:,1] = B[k]
        else
            W = zeros(Al+Bl, No, Ni, Ar+Br)
            W[1:Al,:,:,1:Ar] = A[k]; W[Al+1:end,:,:,Ar+1:end] = B[k]
        end
        C[k] = W
    end
    C
end

"""
    build_M_mpo(locals, pair, dt) -> Vector of MPO cores for M = I − dt·L

`locals` is a length-K vector of nv×nv single-site generators (column convention
dp/dt = L p, state n ↔ index n+1).  `pair` is an nv²×nv² nearest-neighbour
coupling generator applied on every edge (combined index n1*nv+n2+1, row-major),
or `nothing` for an uncoupled chain.  The pair is decomposed by operator-SVD into
a sum of Kronecker terms, so the resulting MPO has bond dim 2 + (#terms).

`dscale` gives a PER-BOND multiplier for the pair coupling (length K−1 for the
edges, or K with the last ignored).  Because diffusion is linear in its rate,
scaling the operator emitted on edge k scales exactly that edge's rate — so a
graded/non-uniform mesh (e.g. fine front voxels at rate D, coarse reservoir legs
at rate D/h²) is built by passing a UNIT pair and `dscale[k] = rate_k`.
"""
function build_M_mpo(locals::Vector{<:AbstractMatrix}, pair, dt::Real; dscale = nothing)
    K = length(locals); nv = size(locals[1], 1)
    locs = Matrix{Float64}.(locals)
    AB = pair === nothing ? Tuple{Matrix{Float64},Matrix{Float64}}[] :
                            op_svd_decompose(Matrix{Float64}(pair), nv)
    ds = dscale === nothing ? ones(K) : collect(float.(dscale))
    length(ds) == K - 1 && (ds = vcat(ds, 1.0))            # bond list → site list
    length(ds) == K || error("dscale must have length K or K−1 (got $(length(ds)))")
    Lmpo = build_L_mpo(locs, AB, nv; dscale = ds)
    mpo_add(identity_mpo(K, nv), mpo_scale_first(Lmpo, -float(dt)))   # I − dt·L
end

# ── immigration/death/diffusion convenience model (matches fig_tt_multigrid) ───
function _death_gen(nv::Int, k_d)
    nmax = nv - 1; G = zeros(nv, nv)
    for n in 1:nmax
        G[n,   n+1] += k_d * n
        G[n+1, n+1] -= k_d * n
    end
    G
end
function _imm_gen(nv::Int, λ)
    nmax = nv - 1; G = zeros(nv, nv)
    for n in 0:nmax-1
        G[n+2, n+1] += λ
        G[n+1, n+1] -= λ
    end
    G
end
# conservative-boundary two-site diffusion (no FSP leak at n=n_max)
function _diffusion_pair(nv::Int, D)
    nmax = nv - 1; N = nv*nv; L = zeros(N, N)
    for n1 in 0:nmax, n2 in 0:nmax
        col = n1*nv + n2 + 1
        if n1 > 0 && n2 < nmax
            L[(n1-1)*nv + (n2+1) + 1, col] += D * n1
            L[col, col]                    -= D * n1
        end
        if n2 > 0 && n1 < nmax
            L[(n1+1)*nv + (n2-1) + 1, col] += D * n2
            L[col, col]                    -= D * n2
        end
    end
    L
end

"""
    bd_diffusion_chain_mpo(K, nv, dt; λ=4.0, k_d=1.0, D=1.0)

Convenience fine MPO for the immigration(site 1)/death(everywhere)/diffusion
birth–death chain over K voxels with count axis size `nv` (n_max = nv-1).
"""
function bd_diffusion_chain_mpo(K::Int, nv::Int, dt::Real; λ = 4.0, k_d = 1.0, D = 1.0)
    locals = [k == 1 ? _death_gen(nv, k_d) + _imm_gen(nv, λ) : _death_gen(nv, k_d) for k in 1:K]
    build_M_mpo(locals, _diffusion_pair(nv, D), dt)
end

# ── TT vector arithmetic (cores [r_left, n, r_right]) ─────────────────────────
tt_delta0(K, nv) = (c = [zeros(1,nv,1) for _ in 1:K]; for k in 1:K; c[k][1,1,1]=1.0; end; c)
tt_zeros(K, nv)  = [zeros(1,nv,1) for _ in 1:K]
tt_scale(a, x)   = (y = deepcopy(x); y[1] = a .* y[1]; y)

function tt_dot(x, y)
    E = ones(1,1)
    for k in 1:length(x)
        xl,nv,xr = size(x[k]); yl,_,yr = size(y[k])
        Enew = zeros(xr, yr)
        for n in 1:nv
            Enew .+= (x[k][:,n,:]') * E * y[k][:,n,:]
        end
        E = Enew
    end
    E[1,1]
end
tt_norm(x) = sqrt(max(tt_dot(x,x), 0.0))

function tt_axpy(a, x, b, y)                                 # a·x + b·y (unrounded)
    K = length(x); z = Vector{Array{Float64,3}}(undef, K)
    for k in 1:K
        xl,nv,xr = size(x[k]); yl,_,yr = size(y[k])
        if k == 1
            W = zeros(1, nv, xr+yr)
            W[:,:,1:xr] = a .* x[k]; W[:,:,xr+1:end] = b .* y[k]
        elseif k == K
            W = zeros(xl+yl, nv, 1)
            W[1:xl,:,:] = x[k]; W[xl+1:end,:,:] = y[k]
        else
            W = zeros(xl+yl, nv, xr+yr)
            W[1:xl,:,1:xr] = x[k]; W[xl+1:end,:,xr+1:end] = y[k]
        end
        z[k] = W
    end
    z
end

# TT rounding: left-orthogonalise (QR) then right-to-left SVD truncation
function tt_round!(x, max_r; tol = 1e-10)
    K = length(x); K == 1 && return x
    for k in 1:K-1
        rl,nv,rr = size(x[k])
        Q, R = qr(reshape(x[k], rl*nv, rr)); Q = Matrix(Q)
        rnew = size(Q, 2)
        x[k] = reshape(Q, rl, nv, rnew)
        rl1,nv1,rr1 = size(x[k+1])
        x[k+1] = reshape(R * reshape(x[k+1], rl1, nv1*rr1), rnew, nv1, rr1)
    end
    for k in K:-1:2
        rl,nv,rr = size(x[k])
        Xk = reshape(x[k], rl, nv*rr)
        U,S,V = try
            svd(Xk)                                   # fast divide-and-conquer
        catch err
            err isa LinearAlgebra.LAPACKException || rethrow()
            svd(Xk; alg = LinearAlgebra.QRIteration())  # robust fallback
        end
        rkeep = max(1, min(max_r, count(>(tol*max(S[1],1e-300)), S)))
        x[k] = reshape(Matrix(V[:,1:rkeep]'), rkeep, nv, rr)
        US = U[:,1:rkeep] * Diagonal(S[1:rkeep])
        rl0,nv0,_ = size(x[k-1])
        x[k-1] = reshape(reshape(x[k-1], rl0*nv0, rl) * US, rl0, nv0, rkeep)
    end
    x
end

# ── MPO apply  y = M x  (bond dims multiply; round afterwards) ────────────────
function mpo_apply(W, x)
    K = length(x); y = Vector{Array{Float64,3}}(undef, K)
    for k in 1:K
        Bl,No,Ni,Br = size(W[k]); Rl,_,Rr = size(x[k])
        # contract physical input index Ni via one BLAS gemm, then reshuffle bonds:
        # y[(bl-1)*Rl+rl, no, (br-1)*Rr+rr] = Σ_ni W[bl,no,ni,br]·x[rl,ni,rr]
        Wm = reshape(permutedims(W[k], (1,4,2,3)), Bl*Br*No, Ni)   # (bl,br,no) × ni
        xm = reshape(permutedims(x[k], (2,1,3)),    Ni, Rl*Rr)     # ni × (rl,rr)
        T  = reshape(Wm * xm, Bl, Br, No, Rl, Rr)
        y[k] = reshape(permutedims(T, (4,1,3,5,2)), Bl*Rl, No, Br*Rr)
    end
    y
end

# ── factor-wise transfers (count-axis halving) and Galerkin ───────────────────
function build_RP(nv::Int)
    @assert iseven(nv) "count axis size $nv must be even to coarsen"
    m = nv ÷ 2
    R = zeros(m, nv); P = zeros(nv, m)
    for c in 0:m-1
        R[c+1, 2c+1] = 1.0; R[c+1, 2c+2] = 1.0      # coarse count = sum of fine pair
        P[2c+1, c+1] = 0.5; P[2c+2, c+1] = 0.5      # split evenly  ⇒  R·P = I
    end
    R, P
end

# coarse MPO: M_c[bl, n_c', n_c, br] = Σ R[n_c',n'] M[bl,n',n,br] P[n,n_c]
function galerkin_mpo(W, R, P)
    K = length(W); C = Vector{Array{Float64,4}}(undef, K)
    for k in 1:K
        Bl,No,Ni,Br = size(W[k])
        Wc = zeros(Bl, size(R,1), size(P,2), Br)
        for bl in 1:Bl, br in 1:Br
            Wc[bl,:,:,br] = R * W[k][bl,:,:,br] * P
        end
        C[k] = Wc
    end
    C
end

restrict_vec(x, R) = [reshape(R * reshape(permutedims(g,(2,1,3)), size(g,2), :),
                              size(R,1), size(g,1), size(g,3)) |> a -> permutedims(a,(2,1,3))
                      for g in x]
prolong_vec(x, P)  = [reshape(P * reshape(permutedims(g,(2,1,3)), size(g,2), :),
                              size(P,1), size(g,1), size(g,3)) |> a -> permutedims(a,(2,1,3))
                      for g in x]

# ── TT-GMRES(m): smoother (few iters) and coarsest solver (many iters) ────────
function tt_gmres(applyM, b, x0, max_r; m = 6, restarts = 1, tol = 1e-10)
    x = deepcopy(x0)
    for _ in 1:restarts
        r = tt_axpy(1.0, b, -1.0, applyM(x)); tt_round!(r, max_r)
        β = tt_norm(r); β < tol && return x
        V = [tt_scale(1/β, r)]
        H = zeros(m+1, m); jmax = m
        for j in 1:m
            w = applyM(V[j]); tt_round!(w, max_r)
            for i in 1:j
                H[i,j] = tt_dot(w, V[i])
                w = tt_axpy(1.0, w, -H[i,j], V[i]); tt_round!(w, max_r)
            end
            H[j+1,j] = tt_norm(w)
            if H[j+1,j] < 1e-12
                jmax = j; break
            end
            push!(V, tt_scale(1/H[j+1,j], w))
        end
        e1 = zeros(jmax+1); e1[1] = β
        yc = H[1:jmax+1, 1:jmax] \ e1
        for j in 1:jmax
            x = tt_axpy(1.0, x, yc[j], V[j]); tt_round!(x, max_r)
        end
    end
    x
end

# ── two-site ALS / AMEn coarsest solver ──────────────────────────────────────
# Minimises ‖M x − b‖ over the TT manifold by sweeping pairs of cores: at each
# bond we form the dense effective two-site system, solve it, and re-split with a
# truncated SVD — so the bond dimension ADAPTS (grows from rank 1 as needed).
function right_canonical!(x)                       # orthogonality centre → site 1
    for k in length(x):-1:2
        rl,nv,rr = size(x[k])
        F = qr(Matrix(reshape(x[k], rl, nv*rr)')); Q = Matrix(F.Q); R = F.R
        rnew = size(Q,2)
        x[k] = reshape(Matrix(Q'), rnew, nv, rr)
        rl0,nv0,_ = size(x[k-1])
        x[k-1] = reshape(reshape(x[k-1], rl0*nv0, rl) * Matrix(R'), rl0, nv0, rnew)
    end
    x
end

function _LM_update(LM, Xk, Wk)                    # grow left M-environment
    rl,nv,rr = size(Xk); bwl = size(Wk,1); bwr = size(Wk,4); new = zeros(rr,bwr,rr)
    for a0 in 1:rl, w0 in 1:bwl, a0p in 1:rl
        L = LM[a0,w0,a0p]; L == 0 && continue
        for np in 1:nv, a in 1:rr
            xb = Xk[a0,np,a]; xb == 0 && continue
            for n in 1:nv, w in 1:bwr
                wv = Wk[w0,np,n,w]; wv == 0 && continue
                Lxbw = L*xb*wv
                for ap in 1:rr; new[a,w,ap] += Lxbw*Xk[a0p,n,ap]; end
            end
        end
    end
    new
end
function _RM_update(RMn, Xk, Wk)                   # grow right M-environment
    rl,nv,rr = size(Xk); bwl = size(Wk,1); bwr = size(Wk,4); new = zeros(rl,bwl,rl)
    for a in 1:rl, np in 1:nv, b in 1:rr
        xb = Xk[a,np,b]; xb == 0 && continue
        for w in 1:bwl, n in 1:nv, w2 in 1:bwr
            wv = Wk[w,np,n,w2]; wv == 0 && continue
            xbw = xb*wv
            for ap in 1:rl, bp in 1:rr
                r = RMn[b,w2,bp]; r == 0 && continue
                new[a,w,ap] += xbw*Xk[ap,n,bp]*r
            end
        end
    end
    new
end
function _Lb_update(Lb, Xk, bk)
    rl,nv,rr = size(Xk); cbl = size(bk,1); cbr = size(bk,3); new = zeros(rr,cbr)
    for a0 in 1:rl, c0 in 1:cbl
        L = Lb[a0,c0]; L == 0 && continue
        for np in 1:nv, a in 1:rr, c in 1:cbr; new[a,c] += L*Xk[a0,np,a]*bk[c0,np,c]; end
    end
    new
end
function _Rb_update(Rbn, Xk, bk)
    rl,nv,rr = size(Xk); cbl = size(bk,1); cbr = size(bk,3); new = zeros(rl,cbl)
    for a in 1:rl, np in 1:nv, b in 1:rr
        xb = Xk[a,np,b]; xb == 0 && continue
        for c in 1:cbl, c2 in 1:cbr; new[a,c] += xb*bk[c,np,c2]*Rbn[b,c2]; end
    end
    new
end

function _solve_two(LMk, Wk, Wk1, RMk2, Lbk, bk, bk1, Rbk2)
    rl = size(LMk,1); rr = size(RMk2,1); nv = size(Wk,2)
    bwl = size(Wk,1); wm = size(Wk,4); bwr = size(Wk1,4)
    Wpair = zeros(bwl,nv,nv,nv,nv,bwr)                # contract shared MPO bond wm
    for w in 1:bwl, p1 in 1:nv, n1 in 1:nv, m in 1:wm
        a = Wk[w,p1,n1,m]; a == 0 && continue
        for p2 in 1:nv, n2 in 1:nv, w2 in 1:bwr
            c = Wk1[m,p2,n2,w2]; c == 0 && continue
            Wpair[w,p1,n1,p2,n2,w2] += a*c
        end
    end
    N = rl*nv*nv*rr
    lin(a,p1,p2,b) = (((a-1)*nv+(p1-1))*nv+(p2-1))*rr + b
    B = zeros(N,N)
    for a in 1:rl, ap in 1:rl, w in 1:bwl, w2 in 1:bwr
        L = LMk[a,w,ap]; L == 0 && continue
        for b_ in 1:rr, bp in 1:rr
            R = RMk2[b_,w2,bp]; R == 0 && continue
            LR = L*R
            for p1 in 1:nv, n1 in 1:nv, p2 in 1:nv, n2 in 1:nv
                wv = Wpair[w,p1,n1,p2,n2,w2]; wv == 0 && continue
                B[lin(a,p1,p2,b_), lin(ap,n1,n2,bp)] += LR*wv
            end
        end
    end
    c = zeros(N); cbl = size(bk,1); cm = size(bk,3); cbr = size(bk1,3)
    for a in 1:rl, c0 in 1:cbl
        L = Lbk[a,c0]; L == 0 && continue
        for b_ in 1:rr, c2 in 1:cbr
            R = Rbk2[b_,c2]; R == 0 && continue
            LR = L*R
            for p1 in 1:nv, mm in 1:cm, p2 in 1:nv
                v = bk[c0,p1,mm]*bk1[mm,p2,c2]; v == 0 && continue
                c[lin(a,p1,p2,b_)] += LR*v
            end
        end
    end
    y = B \ c
    Θ = zeros(rl,nv,nv,rr)
    for a in 1:rl, p1 in 1:nv, p2 in 1:nv, b_ in 1:rr; Θ[a,p1,p2,b_] = y[lin(a,p1,p2,b_)]; end
    Θ
end

function _split!(x, k, Θ, max_r, dir)
    rl = size(Θ,1); nv = size(Θ,2); rr = size(Θ,4)
    M2 = zeros(rl*nv, nv*rr)
    for a in 1:rl, p1 in 1:nv, p2 in 1:nv, b in 1:rr
        M2[(a-1)*nv+p1, (p2-1)*rr+b] = Θ[a,p1,p2,b]
    end
    U,S,V = svd(M2)
    rk = max(1, min(max_r, count(>(1e-12*max(S[1],1e-300)), S)))
    Xk = zeros(rl,nv,rk); Xk1 = zeros(rk,nv,rr)
    for a in 1:rl, p1 in 1:nv, s in 1:rk
        Xk[a,p1,s] = U[(a-1)*nv+p1, s] * (dir == :right ? S[s] : 1.0)
    end
    for s in 1:rk, p2 in 1:nv, b in 1:rr
        Xk1[s,p2,b] = V[(p2-1)*rr+b, s] * (dir == :left ? S[s] : 1.0)
    end
    x[k] = Xk; x[k+1] = Xk1
end

function tt_als(W, b, x0, max_r; sweeps = 5)
    K = length(W); x = deepcopy(x0)
    tt_norm(x) < 1e-14 && (x = deepcopy(b))            # rank-1 zero start is degenerate
    right_canonical!(x)
    RM = Vector{Any}(undef, K+1); RM[K+1] = ones(1,1,1)
    Rb = Vector{Any}(undef, K+1); Rb[K+1] = ones(1,1)
    for k in K:-1:1
        RM[k] = _RM_update(RM[k+1], x[k], W[k]); Rb[k] = _Rb_update(Rb[k+1], x[k], b[k])
    end
    LM = Vector{Any}(undef, K+1); LM[1] = ones(1,1,1)
    Lb = Vector{Any}(undef, K+1); Lb[1] = ones(1,1)
    for _ in 1:sweeps
        for k in 1:K-1                                  # left → right
            Θ = _solve_two(LM[k], W[k], W[k+1], RM[k+2], Lb[k], b[k], b[k+1], Rb[k+2])
            _split!(x, k, Θ, max_r, :left)
            LM[k+1] = _LM_update(LM[k], x[k], W[k]); Lb[k+1] = _Lb_update(Lb[k], x[k], b[k])
        end
        for k in K-1:-1:1                               # right → left
            Θ = _solve_two(LM[k], W[k], W[k+1], RM[k+2], Lb[k], b[k], b[k+1], Rb[k+2])
            _split!(x, k, Θ, max_r, :right)
            RM[k+1] = _RM_update(RM[k+2], x[k+1], W[k+1]); Rb[k+1] = _Rb_update(Rb[k+2], x[k+1], b[k+1])
        end
    end
    x
end

# ── hierarchy + recursive V-cycle ─────────────────────────────────────────────
"""
    build_tt_hierarchy(M_mpo) -> (Ms, RPs, nvs, K)

Build the count-axis multigrid hierarchy from a fine MPO by repeated factor-wise
Galerkin coarsening (halving the count axis while it is even and > 2).
"""
function build_tt_hierarchy(M_mpo)
    nv = size(M_mpo[1], 2); K = length(M_mpo)
    Ms = [M_mpo]; RPs = Tuple{Matrix{Float64},Matrix{Float64}}[]
    nvs = [nv]; cur = nv
    while iseven(cur) && cur > 2
        R, P = build_RP(cur)
        push!(Ms, galerkin_mpo(Ms[end], R, P)); push!(RPs, (R, P))
        cur ÷= 2; push!(nvs, cur)
    end
    (Ms = Ms, RPs = RPs, nvs = nvs, K = K)
end

"""
    tt_vcycle(H, b, x, level, max_r; ν=2, m_coarse=30, coarse_als=true, als_sweeps=6)

One multigrid V-cycle on the TT.  Pre/post-smooth with `ν` TT-GMRES iterations;
restrict the residual, recurse, prolong+correct.  At the coarsest level use
rank-adaptive two-site ALS (`coarse_als=true`) or many-iteration GMRES.
"""
function tt_vcycle(H, b, x, level, max_r; ν = 2, m_coarse = 30, coarse_als = true, als_sweeps = 6)
    applyM = z -> (y = mpo_apply(H.Ms[level], z); tt_round!(y, max_r); y)
    if level == length(H.Ms)                                       # coarsest (item 3)
        return coarse_als ? tt_als(H.Ms[level], b, x, max_r; sweeps = als_sweeps) :
                            tt_gmres(applyM, b, x, max_r; m = m_coarse, restarts = 2)
    end
    x = tt_gmres(applyM, b, x, max_r; m = ν)                       # pre-smooth  (item 1)
    r = tt_axpy(1.0, b, -1.0, applyM(x)); tt_round!(r, max_r)      # residual
    R, P = H.RPs[level]
    rc = restrict_vec(r, R)                                        # restrict   (item 2)
    ec = tt_vcycle(H, rc, tt_zeros(H.K, H.nvs[level+1]), level+1, max_r; ν, m_coarse, coarse_als, als_sweeps)
    x = tt_axpy(1.0, x, 1.0, prolong_vec(ec, P)); tt_round!(x, max_r)
    tt_gmres(applyM, b, x, max_r; m = ν)                           # post-smooth (item 1)
end

"""
    tt_be_solve(M_mpo, b; max_r=24, ν=3, coarse_als=true, tol=1e-8,
                maxcycles=40, m_coarse=30, verbose=false)
        -> (x, residual, residuals, cycles, hierarchy)

Solve the backward-Euler step `M_mpo x = b` (M = I − dt·L) by repeated TT
multigrid V-cycles until the TT residual norm falls below `tol`.  `b` and the
returned `x` are TT vectors (vectors of [r_left, n, r_right] cores); use
[`tt_delta0`] for a δ-at-zero RHS and [`tt_to_full`] to extract a dense
distribution for small chains.  `residuals` is the per-cycle TT residual-norm
history (for measuring the contraction factor when no dense reference exists).
"""
function tt_be_solve(M_mpo, b; max_r = 24, ν = 3, coarse_als = true,
                     tol = 1e-8, maxcycles = 40, m_coarse = 30, verbose = false)
    H = build_tt_hierarchy(M_mpo)
    x = tt_zeros(H.K, H.nvs[1])
    applyM = z -> (y = mpo_apply(M_mpo, z); tt_round!(y, max_r); y)
    res = Inf; cyc = 0; hist = Float64[]
    for c in 1:maxcycles
        cyc = c
        x = tt_vcycle(H, b, x, 1, max_r; ν, m_coarse, coarse_als)
        r = tt_axpy(1.0, b, -1.0, applyM(x)); tt_round!(r, max_r)
        res = tt_norm(r); push!(hist, res)
        verbose && println("  cycle ", lpad(c,2), "   ‖residual‖ = ", res)
        res < tol && break
    end
    (x = x, residual = res, residuals = hist, cycles = cyc, hierarchy = H)
end

# ── dense contractions (validation / small-chain extraction) ──────────────────
function tt_to_full(x)
    K = length(x); T = reshape(x[1][1,:,:], size(x[1],2), size(x[1],3))
    for k in 2:K
        rprev = size(T, 2); nv = size(x[k],2); rk = size(x[k],3)
        T = reshape(T * reshape(x[k], rprev, nv*rk), :, rk)
    end
    vec(T)
end

function mpo_to_dense(W)
    K = length(W); nv = size(W[1], 2); N = nv^K
    M = zeros(N, N)
    cols = Iterators.product(ntuple(_ -> 1:nv, K)...)
    for (cj, c) in enumerate(cols)
        e = tt_zeros(K, nv)
        for k in 1:K; e[k] = zeros(1,nv,1); e[k][1, c[k], 1] = 1.0; end
        M[:, cj] = tt_to_full(mpo_apply(W, e))
    end
    M
end

end # module TTVCycle
