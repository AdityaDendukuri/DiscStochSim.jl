using DiscStochSim
using Catalyst
using Random
using Plots

const OREGONATOR_FACE = RGB(0.92, 0.92, 0.92)

function parse_tfinal(args; default::Float64=0.11)::Float64
    isempty(args) && return default
    try
        parse(Float64, args[1])
    catch
        default
    end
end

function oregonator_problem(; tfinal::Float64=0.11, bounds=(0, 50_000))
    rn = @reaction_network begin
        @species X(t) Y(t) Z(t)
        @parameters k1 k2 k3 k4 k5
        k1, Y --> X
        k2, X + Y --> 0
        k3, X --> 2X + Z
        k4, 2X --> 0
        k5, Z --> Y
    end

    y1s, y2s, y3s = 500.0, 1000.0, 2000.0
    mu1s, mu2s = 2000.0, 50000.0
    rates = [
        mu1s / y2s,
        mu2s / (y1s * y2s),
        (mu1s + mu2s) / y1s,
        2 * mu1s / y1s^2,
        (mu1s + mu2s) / y3s,
    ]

    FSPProblem(rn, CartesianIndex(500, 1000, 2000), (0.0, tfinal), rates; bounds=bounds)
end

function solve_oregonator(; tfinal::Float64=0.11,
                          ε_dt::Float64=1.0,
                          prob_quantile::Float64=0.1,
                          flux_tolerance::Float64=1.0,
                          expansion_depth::Int=1,
                          save_interval::Int=1000,
                          seed::Int=1,
                          bounds=(0, 50_000))
    Random.seed!(seed)

    prob = oregonator_problem(; tfinal=tfinal, bounds=bounds)
    sol, _ = solve(prob, AdaptiveFSP(
        ε_dt=ε_dt,
        prob_quantile=prob_quantile,
        flux_tolerance=flux_tolerance,
        expansion_depth=expansion_depth,
        save_interval=save_interval,
    ))

    t = sol.t
    dt = vcat([0.0], t[2:end] .- t[1:end-1])
    traj = mean_trajectory(sol)

    (
        t = t,
        dt = dt,
        state_space_sizes = sol.state_space_sizes,
        traj = traj,
        sol = sol,
    )
end

function apply_oregonator_plot_style!()
    gr()
    default(
        fontfamily="Computer Modern",
        titlefontsize=18,
        guidefontsize=14,
        tickfontsize=10,
        legendfontsize=10,
        linewidth=2.3,
        framestyle=:box,
        grid=false,
        foreground_color=:black,
        background_color=:white,
        background_color_inside=OREGONATOR_FACE,
        dpi=300,
    )
end

function oregonator_paper_figure(data)
    apply_oregonator_plot_style!()

    t = data.t
    dt = data.dt
    s = data.state_space_sizes
    x = data.traj[:, 1]
    y = data.traj[:, 2]
    z = data.traj[:, 3]

    axis_kw = (
        tick_direction=:in,
        color=:black,
        grid=false,
        framestyle=:box,
    )

    p1 = plot(
        t[2:end], dt[2:end];
        title="(a)",
        xlabel="Time, t",
        ylabel="δt",
        yscale=:log10,
        color=:black,
        legend=false,
        axis_kw...,
    )

    p2 = plot(
        t, s;
        title="(b)",
        xlabel="Time, t",
        ylabel="|S|",
        color=:black,
        legend=false,
        axis_kw...,
    )

    p3 = plot(
        t, x;
        title="(c)",
        xlabel="Time, t",
        ylabel="Mean population",
        color=:black,
        linestyle=:solid,
        label="X (y1)",
        legend=:topright,
        axis_kw...,
    )
    plot!(p3, t, y; color=:black, linestyle=:dash, label="Y (y2)")
    plot!(p3, t, z; color=:black, linestyle=:dot, label="Z (y3)")

    p4 = plot(
        x, y;
        title="(d)",
        xlabel="X",
        ylabel="Y",
        color=:black,
        legend=false,
        axis_kw...,
    )

    p5 = plot(
        x, z;
        title="(e)",
        xlabel="X",
        ylabel="Z",
        color=:black,
        legend=false,
        axis_kw...,
    )

    p6 = plot(
        y, z;
        title="(f)",
        xlabel="Y",
        ylabel="Z",
        color=:black,
        legend=false,
        axis_kw...,
    )

    plot(
        p1, p2, p3, p4, p5, p6;
        layout=(2, 3),
        size=(1200, 800),
        margin=7Plots.mm,
    )
end

function save_oregonator_figure(; tfinal::Float64=0.11, outdir::String=joinpath(@__DIR__, "output"))
    data = solve_oregonator(; tfinal=tfinal, save_interval=1000)
    fig = oregonator_paper_figure(data)

    mkpath(outdir)
    suffix = replace(string(tfinal), "." => "_")
    png_out = joinpath(outdir, "oregonator_results_t$(suffix).png")
    pdf_out = joinpath(outdir, "oregonator_results_t$(suffix).pdf")

    Plots.savefig(fig, png_out)
    Plots.savefig(fig, pdf_out)

    traj = data.traj
    println("Final: X=$(traj[end,1]), Y=$(traj[end,2]), Z=$(traj[end,3])")
    println("Snapshots: $(length(data.t)), final |S| = $(data.state_space_sizes[end])")
    println("Saved $(png_out)")
    println("Saved $(pdf_out)")
end

save_oregonator_figure(; tfinal=parse_tfinal(ARGS; default=0.11))
