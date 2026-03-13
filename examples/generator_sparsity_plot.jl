# Generator sparsity pattern: with vs. without bottleneck state
# Illustrates that pruning s_b makes A block-diagonal (disconnected).
# States: [s1, s2, sb, s3, s4] — s1,s2 in region A; s3,s4 in region B; sb bridges them.
using Plots
gr()

# ── Build generators ──────────────────────────────────────────────────────────
λ = 1.0   # internal metastable rate (s1↔s2, s3↔s4)
μ = 1.0   # bottleneck transition rate (sb ↔ {s1,s2,s3,s4})

# Full 5×5 generator (with sb at index 3)
A = zeros(5, 5)
A[2,1]=λ;  A[1,2]=λ     # s1 ↔ s2
A[3,1]=μ;  A[1,3]=μ     # s1 ↔ sb
A[3,2]=μ;  A[2,3]=μ     # s2 ↔ sb
A[4,3]=μ;  A[3,4]=μ     # sb ↔ s3
A[5,3]=μ;  A[3,5]=μ     # sb ↔ s4
A[5,4]=λ;  A[4,5]=λ     # s3 ↔ s4
for j in 1:5; A[j,j] = -sum(A[i,j] for i in 1:5 if i != j); end

# Pruned 4×4 generator (sb removed; states [s1,s2,s3,s4] = indices [1,2,4,5])
# Compressed generator: only transitions within the remaining set survive.
B = zeros(4, 4)
B[2,1]=λ;  B[1,2]=λ     # s1 ↔ s2
B[4,3]=λ;  B[3,4]=λ     # s3 ↔ s4  (sb↔s3/s4 connections are gone)
for j in 1:4; B[j,j] = -sum(B[i,j] for i in 1:4 if i != j); end

# ── Spy-plot helper ───────────────────────────────────────────────────────────
function make_spy(M, labels, title_str, block_ranges)
    n = size(M, 1)

    # Start with an empty plot (white background) — no heatmap to avoid interpolation
    p = plot(;
        title         = title_str,
        titlefontsize = 11,
        xticks        = (1:n, labels),
        yticks        = (1:n, labels),
        yflip         = true,
        aspect_ratio  = :equal,
        framestyle    = :box,
        grid          = false,
        tickfontsize  = 11,
        xlims         = (0.5, n + 0.5),
        ylims         = (0.5, n + 0.5),
        legend        = false,
        background_color_inside = :white,
    )

    # Draw a filled black square for each nonzero entry
    for i in 1:n, j in 1:n
        if M[i, j] != 0
            xs = [j-0.5, j+0.5, j+0.5, j-0.5, j-0.5]
            ys = [i-0.5, i-0.5, i+0.5, i+0.5, i-0.5]
            plot!(p, xs, ys;
                seriestype  = :shape,
                fillcolor   = :black,
                linecolor   = :black,
                linewidth   = 0,
                label       = false,
            )
        end
    end

    # Light gray cell-grid lines
    for k in 1.5:1.0:(n-0.5)
        plot!(p, [0.5, n+0.5], [k, k]; line = (0.5, :gray), label = false)
        plot!(p, [k, k], [0.5, n+0.5]; line = (0.5, :gray), label = false)
    end

    # Dashed boxes around each metastable block
    for (a, b) in block_ranges
        xs = [a-0.5, b+0.5, b+0.5, a-0.5, a-0.5]
        ys = [a-0.5, a-0.5, b+0.5, b+0.5, a-0.5]
        plot!(p, xs, ys; line = (:dash, 2.5, :white), label = false)
        plot!(p, xs, ys; line = (:dash, 0.8, :black), label = false)
    end
    p
end

labels_full  = ["\$s_1\$", "\$s_2\$", "\$s_b\$", "\$s_3\$", "\$s_4\$"]
labels_prune = ["\$s_1\$", "\$s_2\$", "\$s_3\$", "\$s_4\$"]

# Block ranges: region A = indices 1:2; region B = indices 4:5 (full) or 3:4 (pruned)
p1 = make_spy(A, labels_full,  "(a) with \$s_b\$",    [(1,2),(4,5)])
p2 = make_spy(B, labels_prune, "(b) without \$s_b\$", [(1,2),(3,4)])

fig = plot(p1, p2;
    layout = (1, 2),
    size   = (700, 340),
    left_margin   = 4Plots.mm,
    right_margin  = 2Plots.mm,
    bottom_margin = 2Plots.mm,
    top_margin    = 3Plots.mm,
)

mkpath("examples/output")
savefig(fig, "examples/output/generator_sparsity.pdf")
savefig(fig, "examples/output/generator_sparsity.png")
println("Saved → examples/output/generator_sparsity.{pdf,png}")
