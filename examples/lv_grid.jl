# Algorithm visualization: expand / evolve / prune (lv_grid.pdf)
# 2×3 grid for toggle switch at t≈2:
#   Row 1: active states (truncation boundary) at each phase
#   Row 2: log₁₀(p) probability contours at each phase
#
# Run: julia --project examples/lv_grid.jl  (~2 min)
using DiscStochSim, Catalyst, Plots, Expokit, SparseArrays, LinearAlgebra

gr()
default(
    fontfamily    = "Computer Modern",
    titlefontsize  = 10,
    guidefontsize  = 9,
    tickfontsize   = 8,
    linewidth = 1.0,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

# ── Toggle switch model ──────────────────────────────────────────────────────
rn = @reaction_network begin
    @species U(t) V(t)
    @parameters η α1 β1 K1 d1 dU α2 β2 K2 d2
    η*(α1 + β1*K1^3/(K1^3 + V^3)),   0 --> U
    d1 + dU,                           U --> 0
    η*(α2 + β2*K2^3/(K2^3 + U^3)),   0 --> V
    d2,                                V --> 0
end
rates = [1.0, 20.0, 400.0, 100.0, 1.0, 0.1/1.1, 20.0, 400.0, 100.0, 1.0]
NGRID = 150

model = DiscreteStochasticSystem(rn)
bc    = x -> RectLatticeBoundaryCondition(x, (0, NGRID - 1))

# ── Advance to t≈2 ──────────────────────────────────────────────────────────
sp = StateSpace{CartesianIndex{2}, Float64}()
add_state!(sp, CartesianIndex(85, 5), 1.0)

A_old    = spzeros(Float64, 0, 0)
gids_old = Int[]
t_cur    = 0.0
t_target = 2.0

println("Advancing to t≈$t_target...")
while t_cur < t_target
    expand!(sp, model, bc; depth=1)
    A, _, _, _ = reconstruct_generator(sp, model, rates, t_cur, A_old, gids_old)
    A_old    = A
    gids_old = get_global_ids(sp)
    Φ = dot(-diag(A), sp.probs)
    dt = min(Φ > 0 ? 1.0/Φ : (t_target - t_cur), t_target - t_cur)
    sp.probs = expmv(dt, A, sp.probs)
    compress!(sp, model, rates, t_cur, 0.3; flux_tolerance=1.0)
    renormalize!(sp)
    t_cur += dt
end
println("t=$(round(t_cur; sigdigits=4)), |S|=$(length(sp))")

# ── Capture expand / evolve / prune phases ───────────────────────────────────
# Pre-expand state
pre_states = copy(sp.states)
pre_probs  = copy(sp.probs)

# After expand
expand!(sp, model, bc; depth=1)
exp_states = copy(sp.states)   # expanded set (new states have prob 0)

# After evolve
A, _, _, _ = reconstruct_generator(sp, model, rates, t_cur, A_old, gids_old)
Φ = dot(-diag(A), sp.probs)
dt = min(Φ > 0 ? 1.0/Φ : 0.05, 0.5)
sp.probs = expmv(dt, A, sp.probs)
evo_states = copy(sp.states)
evo_probs  = copy(sp.probs)

# After prune
compress!(sp, model, rates, t_cur, 0.3; flux_tolerance=1.0)
renormalize!(sp)
pru_states = copy(sp.states)
pru_probs  = copy(sp.probs)

println("Sizes: pre=$(length(pre_states)), expanded=$(length(exp_states)), " *
        "evolved=$(length(evo_states)), pruned=$(length(pru_states))")

# ── Build 2D grids ───────────────────────────────────────────────────────────
function make_grids(states, probs, N)
    active  = zeros(Bool, N, N)
    logprob = fill(NaN, N, N)
    for (s, p) in zip(states, probs)
        u, v = s[1]+1, s[2]+1
        1 <= u <= N && 1 <= v <= N || continue
        active[u, v] = true
        p > 0 && (logprob[u, v] = log10(p))
    end
    active, logprob
end

# Determine axis limits from combined active set
all_u = [s[1] for s in exp_states]
all_v = [s[2] for s in exp_states]
pad = 8
u_lo = max(minimum(all_u) - pad, 0);  u_hi = min(maximum(all_u) + pad, NGRID-1)
v_lo = max(minimum(all_v) - pad, 0);  v_hi = min(maximum(all_v) + pad, NGRID-1)

# Shared color range across all probability panels
all_lp = Float64[]
for (st, pr) in [(evo_states, evo_probs), (pru_states, pru_probs), (pre_states, pre_probs)]
    append!(all_lp, filter(p -> p > 0, pr) .|> log10)
end
vmin = isempty(all_lp) ? -10.0 : max(minimum(all_lp), -10.0)
vmax = isempty(all_lp) ?   0.0 : maximum(all_lp)

# ── Figure: 2 rows × 3 cols ─────────────────────────────────────────────────
titles = ["Expand", "Evolve", "Prune"]

# Row 1: boundary (scatter of active states)
function boundary_panel(active_states, added_states, title)
    au = [s[1] for s in active_states]
    av = [s[2] for s in active_states]
    p = scatter(au, av; markersize=2, markerstrokewidth=0,
                color=:black, label="", title=title,
                xlabel="U", ylabel="V",
                xlims=(u_lo, u_hi), ylims=(v_lo, v_hi),
                aspect_ratio=:equal)
    if !isnothing(added_states)
        nu = [s[1] for s in added_states]
        nv = [s[2] for s in added_states]
        scatter!(p, nu, nv; markersize=2, markerstrokewidth=0,
                 color=:gray, alpha=0.4, label="")
    end
    p
end

# Row 2: log probability heatmap
function prob_panel(states, probs, title)
    _, logprob = make_grids(states, probs, NGRID)
    u_range = (u_lo+1):(u_hi+1)
    v_range = (v_lo+1):(v_hi+1)
    sub = logprob[u_range, v_range]'  # transpose for (V, U) layout
    heatmap(u_lo:(u_hi), v_lo:(v_hi), sub;
            clims=(vmin, vmax), color=:viridis,
            xlabel="U", ylabel="V", title=title,
            xlims=(u_lo, u_hi), ylims=(v_lo, v_hi),
            aspect_ratio=:equal)
end

# Expand column: show original (black) + newly added states (gray) in boundary row
pre_set = Set(pre_states)
added   = [s for s in exp_states if s ∉ pre_set]

r1_expand = boundary_panel(pre_states, added, "Expand")
r1_evolve  = boundary_panel(evo_states, nothing, "Evolve")
r1_prune   = boundary_panel(pru_states, nothing, "Prune")

r2_expand = prob_panel(pre_states, pre_probs, "")
r2_evolve  = prob_panel(evo_states, evo_probs, "")
r2_prune   = prob_panel(pru_states, pru_probs, "")

fig = plot(r1_expand, r1_evolve, r1_prune,
           r2_expand, r2_evolve, r2_prune;
           layout=(2, 3), size=(900, 600), margin=4Plots.mm)

paper_dir = joinpath(@__DIR__, "..", "paper")
savefig(fig, joinpath(paper_dir, "lv_grid.pdf"))
savefig(fig, joinpath(paper_dir, "lv_grid.png"))
println("Saved lv_grid.{pdf,png} → paper/")
