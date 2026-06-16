"""
    build_rdme_joint(voxel_rn, graph; D) -> DiscreteStochasticSystem{CartesianIndex{S*K}}

General RDME generator from a per-voxel Catalyst `ReactionSystem` and a
`VoxelGraph`. Each of the `S` base species is replicated across the `K` voxels:
component `(v-1)*S + s` of the joint state holds species `s` in voxel `v`.

`voxel_rn` is either

  * a single `ReactionSystem` — the same kinetics in every voxel (the genuine,
    spatially-uniform RDME); or
  * a function `v -> ReactionSystem` — voxel-dependent kinetics, e.g. a localized
    immigration source in voxel 1 only. Every voxel's network must declare the
    same `S` species in the same order (component `s` is the same physical species
    in every voxel, since species `s` is what diffuses between voxels).

Two reaction families:
  * local   – every base reaction applied independently inside each voxel;
  * diffusion – the linear hop `X_{s,i} → X_{s,j}` at rate `D · n_{s,i}` on every
    graph edge (both directions). Diffusion as a monomolecular reaction is exact,
    no closure.

The result plugs straight into `build_generator(sp, sys, P, t)`, where `P` is the
base parameter vector (same ordering as `Catalyst.get_ps`); the local propensities
reuse `P`, diffusion ignores it and uses the closed-over `D`. If the networks bake
numeric rate constants (no free parameters), `P` may be `Float64[]` / `nothing`.

Assembly is generic and the resulting sparse `A` feeds the same AMG / V-cycle /
coarsening unchanged.
"""
function build_rdme_joint(voxel_rn, graph::VoxelGraph; D::Real)
    K     = graph.n_voxels
    rn_of = voxel_rn isa ReactionSystem ? (v -> voxel_rn) : voxel_rn

    cache   = IdDict{Any,Any}()
    getbase(rn) = get!(() -> DiscreteStochasticSystem(rn), cache, rn)

    S = length(Catalyst.get_species(rn_of(1)))
    for v in 2:K
        length(Catalyst.get_species(rn_of(v))) == S ||
            error("every voxel network must declare the same $S species (voxel $v differs)")
    end

    N    = S * K
    valS = Val(S)
    Dval = float(D)

    stoichs = CartesianIndex{N}[]
    props   = Function[]

    # local reactions, replicated per voxel (kinetics may depend on the voxel)
    for v in 1:K
        off  = (v - 1) * S
        base = getbase(rn_of(v))
        for r in eachindex(base.stoichvecs)
            sv  = base.stoichvecs[r]                                   # CartesianIndex{S}
            big = CartesianIndex(ntuple(c -> (off < c <= off + S) ? sv[c - off] : 0, N))
            let off = off, fr = base.propensities[r], valS = valS
                push!(stoichs, big)
                push!(props, (x, rv, t) -> fr(CartesianIndex(ntuple(s -> x[off + s], valS)), rv, t))
            end
        end
    end

    # diffusion: linear hop per species, per edge, both directions
    for (i, j) in graph.edges
        offi = (i - 1) * S; offj = (j - 1) * S
        for s in 1:S
            ci = offi + s; cj = offj + s
            big_ij = CartesianIndex(ntuple(c -> c == ci ? -1 : (c == cj ? 1 : 0), N))
            big_ji = CartesianIndex(ntuple(c -> c == cj ? -1 : (c == ci ? 1 : 0), N))
            let ci = ci, Dval = Dval
                push!(stoichs, big_ij)
                push!(props, (x, rv, t) -> Dval * max(0, x[ci]))
            end
            let cj = cj, Dval = Dval
                push!(stoichs, big_ji)
                push!(props, (x, rv, t) -> Dval * max(0, x[cj]))
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, props)
end
