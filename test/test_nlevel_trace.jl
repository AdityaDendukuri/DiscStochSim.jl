using DiscStochSim
using ExponentialUtilities

# ─── Tiny Schlögl: K=4, 1 level (4→2) ──────────────────────────────────────
K     = 4
model = SchloglModel1D(0.1)
fine_grid   = VoxelGrid(K, 1.0, 0)
coarse_grid = VoxelGrid(K÷2, 2.0, 0)

u0    = CartesianIndex(40, 0, 0, 0)
rates = Float64[]
dt    = 0.05

model_2h = DiscStochSim.coarsen_model(model, 2.0)

# ── warmup helper: direct fine FSP ───────────────────────────────────────────
function warmup_fine_fsp(u0, n_steps, dt, fine_sys, bc_fine)
    sp = DiscStochSim.StateSpace{CartesianIndex{K}, Float64}()
    DiscStochSim.add_state!(sp, u0, 1.0)
    DiscStochSim.expand!(sp, fine_sys, bc_fine; depth=3)
    for _ in 1:n_steps
        A, = DiscStochSim.build_generator(sp, fine_sys, rates, 0.0)
        sp.probs .= expv(dt, A, sp.probs; m=30)
        DiscStochSim.renormalize!(sp)
        DiscStochSim.prune_threshold!(sp, 1e-10)
        DiscStochSim.expand!(sp, fine_sys, bc_fine; depth=1)
    end
    sp
end

# ── helpers ──────────────────────────────────────────────────────────────────

function show_sp(label, sp)
    total = sum(sp.probs)
    ns    = length(sp)
    top5  = sortperm(sp.probs; rev=true)[1:min(5,ns)]
    println("  $label  |S|=$ns  Σp=$(round(total,digits=6))")
    for i in top5
        s = sp.states[i]; p = sp.probs[i]
        p > 0 && println("    $(Tuple(s))  p=$(round(p,digits=5))")
    end
end

function show_marginals(label, sp)
    μ = zeros(K)
    for (i,s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p*t[k]; end
    end
    println("  $label  μ=$(round.(μ, digits=2))")
    for j in 1:K÷2
        tot = μ[2j-1]+μ[2j]
        α   = tot > 0 ? abs(μ[2j-1]-μ[2j])/tot : 0.0
        println("    pair $j: μL=$(round(μ[2j-1],digits=2)) μR=$(round(μ[2j],digits=2))  α=$(round(α,digits=3))")
    end
end

# ── run ───────────────────────────────────────────────────────────────────────
function run_trace()
    fine_sys = DiscStochSim.build_schlogl_rdme_system(model, fine_grid)
    bc_fine  = s -> DiscStochSim.rdme_bc(s, 150)

    # Use concentrated IC directly — this tests whether the L1=2 fix
    # allows molecules to diffuse from pair 1 to pair 2
    sp_cur = let
        sp0 = DiscStochSim.StateSpace{CartesianIndex{K}, Float64}()
        DiscStochSim.add_state!(sp0, u0, 1.0)
        DiscStochSim.expand!(sp0, fine_sys, bc_fine; depth=1)
        sp0
    end

    println("═"^60)
    println("INITIAL (concentrated IC: $u0)")
    show_sp("fine", sp_cur)
    show_marginals("fine", sp_cur)

    for step in 1:3
        println()
        println("═"^60)
        println("STEP $step")

        # 1. Restrict — prune boundary (p=0) fine states first so sp_coarse_pre
        #    has no zero-coverage coarse states (which would be invisible to prolong step 3)
        sp_cur_pos = DiscStochSim.StateSpace{CartesianIndex{K}, Float64}()
        for i in eachindex(sp_cur.states)
            sp_cur.probs[i] > 1e-14 && DiscStochSim.add_state!(sp_cur_pos, sp_cur.states[i], sp_cur.probs[i])
        end
        sp_coarse     = DiscStochSim.restrict(sp_cur_pos)
        sp_coarse_pre = DiscStochSim._copy_sp(sp_coarse)
        show_sp("coarse (pre-expand)", sp_coarse)

        # 2. Coarse expand
        coarse_sys = DiscStochSim.build_schlogl_rdme_system(model_2h, coarse_grid)
        bc_coarse  = s -> DiscStochSim.rdme_bc(s, 300)
        DiscStochSim.expand!(sp_coarse, coarse_sys, bc_coarse; depth=1)
        println("  coarse after expand: |S|=$(length(sp_coarse))")

        # 3. Coarse solve
        A_c, = DiscStochSim.build_generator(sp_coarse, coarse_sys, rates, 0.0)
        sp_coarse.probs .= expv(dt, A_c, sp_coarse.probs; m=30)
        show_sp("coarse (post-solve)", sp_coarse)

        # 4. Prolongation — compare multiplicative vs conditional, and diagnose step 3
        n_pre_here = length(sp_coarse_pre)
        sp_mult = DiscStochSim.prolong_multiplicative(sp_cur_pos, sp_coarse_pre, sp_coarse;
                                                       prob_tol=1e-13)
        sp_cond = DiscStochSim.prolong_conditional(sp_cur_pos, sp_coarse_pre, sp_coarse;
                                                    prob_tol=1e-13)
        println("  mult-prolong:  |S|=$(length(sp_mult))  mass=$(round(sum(sp_mult.probs),digits=4))")
        println("  cond-prolong:  |S|=$(length(sp_cond))  mass=$(round(sum(sp_cond.probs),digits=4))")

        # Direct sanity check: what mass SHOULD mult-prolong give?
        expected_mass = 0.0
        n_no_pre = 0; n_no_post = 0; n_zero_pre = 0
        for i in eachindex(sp_cur_pos.states)
            xf = sp_cur_pos.states[i]; pf = sp_cur_pos.probs[i]
            tf = Tuple(xf)
            xc = CartesianIndex(tf[1]+tf[2], tf[3]+tf[4])
            i_pre  = get(sp_coarse_pre.index, xc, 0)
            i_post = get(sp_coarse.index,     xc, 0)
            p_pre  = i_pre  > 0 ? sp_coarse_pre.probs[i_pre]  : 0.0
            p_post = i_post > 0 ? sp_coarse.probs[i_post]     : 0.0
            if i_pre  == 0; n_no_pre  += 1; end
            if i_post == 0; n_no_post += 1; end
            if p_pre <= 1e-13; n_zero_pre += 1; continue; end
            expected_mass += pf * (p_post / p_pre)
        end
        println("  expected mult mass=$(round(expected_mass,digits=4))  n_no_pre=$n_no_pre  n_no_post=$n_no_post  n_zero_pre=$n_zero_pre")

        # Diagnose: how much mass is in expanded coarse states (indices n_pre_here+1..end)?
        mass_covered  = sum(sp_coarse.probs[i] for i in 1:n_pre_here)
        mass_expanded = sum(sp_coarse.probs[i] for i in n_pre_here+1:length(sp_coarse))
        println("  coarse mass: covered=$(round(mass_covered,digits=4))  expanded=$(round(mass_expanded,digits=4))")
        # Show a few expanded coarse states and whether they find a covered L1=1 neighbor
        println("  sample expanded states (checking L1=1 neighbor in sp_coarse_pre):")
        for i in (n_pre_here+1):min(n_pre_here+5, length(sp_coarse))
            xc = sp_coarse.states[i]; p  = sp_coarse.probs[i]
            tc = Tuple(xc)
            found = false
            for j in 1:K÷2, delta in (1, -1)
                t_try = ntuple(k -> k==j ? tc[k]-delta : tc[k], Val(K÷2))
                all(v->v>=0, t_try) || continue
                if get(sp_coarse_pre.index, CartesianIndex(t_try), 0) != 0
                    found = true; break
                end
            end
            println("    $tc  p=$(round(p,digits=5))  L1=1_covered=$found")
        end

        # Also try correction-based Binomial (original approach)
        sp_δ_pos = DiscStochSim.StateSpace{CartesianIndex{K÷2}, Float64}()
        sp_δ_neg = DiscStochSim.StateSpace{CartesianIndex{K÷2}, Float64}()
        n_pre = length(sp_coarse_pre)
        probs_pre_vec = sp_coarse_pre.probs
        for i in eachindex(sp_coarse.states)
            prob_before = i <= n_pre ? probs_pre_vec[i] : 0.0
            δ = sp_coarse.probs[i] - prob_before
            s = sp_coarse.states[i]
            if δ > 1e-13; DiscStochSim.add_state!(sp_δ_pos, s, δ)
            elseif δ < -1e-13; DiscStochSim.add_state!(sp_δ_neg, s, -δ); end
        end
        sp_δf_pos = DiscStochSim.prolong(sp_δ_pos, Val(K); weight_tol=1e-13, binom_tol=1e-4)
        sp_δf_neg = DiscStochSim.prolong(sp_δ_neg, Val(K); weight_tol=1e-13, binom_tol=1e-4)
        sp_binom = DiscStochSim._copy_sp(sp_cur)
        for i in eachindex(sp_δf_pos.states)
            s=sp_δf_pos.states[i]; δp=sp_δf_pos.probs[i]
            idx=get(sp_binom.index,s,0)
            if idx==0; DiscStochSim.add_state!(sp_binom,s,δp) else sp_binom.probs[idx]+=δp end
        end
        for i in eachindex(sp_δf_neg.states)
            s=sp_δf_neg.states[i]; δp=sp_δf_neg.probs[i]
            idx=get(sp_binom.index,s,0)
            if idx==0; DiscStochSim.add_state!(sp_binom,s,-δp) else sp_binom.probs[idx]-=δp end
        end
        for i in eachindex(sp_binom.probs); sp_binom.probs[i] < 0 && (sp_binom.probs[i]=0.0); end
        println("  binom-corr:    |S|=$(length(sp_binom))  mass=$(round(sum(sp_binom.probs),digits=4))")
        # show top-5 of the binom approach to see where the mass went
        top5b = sortperm(sp_binom.probs; rev=true)[1:min(5,length(sp_binom))]
        for i in top5b; p=sp_binom.probs[i]; p>0 && println("    $(Tuple(sp_binom.states[i]))  p=$(round(p,digits=4))"); end

        sp_new = sp_binom   # use binom-correction for now
        show_marginals("fine (binom-corr)", sp_new)

        # 5. Mass conservation check
        mass_before = sum(sp_cur.probs)
        mass_after  = sum(sp_new.probs)
        println("  mass: before=$(round(mass_before,digits=6))  after=$(round(mass_after,digits=6))")

        # 6. Renormalize, prune, expand
        DiscStochSim.renormalize!(sp_new)
        DiscStochSim.prune_threshold!(sp_new, 1e-12)
        DiscStochSim.expand!(sp_new, fine_sys, bc_fine; depth=1)
        println("  fine after expand: |S|=$(length(sp_new))")

        sp_cur = sp_new
    end

    println()
    println("═"^60)
    println("FINAL")
    show_sp("fine", sp_cur)
    show_marginals("fine", sp_cur)
end

run_trace()
