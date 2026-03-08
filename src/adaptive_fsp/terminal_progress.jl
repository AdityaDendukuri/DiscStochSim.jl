"""
    terminal_progress(species_names; tf=nothing, height=15, width=70, title="FSP progress",
                      xscale=:linear, transform=identity)

Return a progress callback for `AdaptiveFSP` that renders a live
`UnicodePlots` line plot of mean species trajectories in the terminal,
refreshing in-place on every snapshot.

# Usage
```julia
cb = terminal_progress(["A", "B", "C"]; tf=1e4)
sol, diag = solve(prob, AdaptiveFSP(save_interval=500, progress_callback=cb))
```

Set `xscale=:log10` for systems with large time horizons (e.g. Robertson, Bottleneck)
where dynamics span many orders of magnitude.

Use `transform` to pre-process the trajectory matrix (n_snaps × n_species) before
plotting — useful for rescaling species with very different magnitudes:

```julia
# Show 1000×B so the transient spike is visible alongside A and C
cb = terminal_progress(["A", "1000×B", "C"]; tf=1e4, xscale=:log10,
                       transform = m -> (m2 = copy(m); m2[:,2] .*= 1000; m2))
```

The plot is updated every `save_interval` iterations.  Press Ctrl-C to
abort early; the partial solution is not returned in that case.
"""
function terminal_progress(species_names::AbstractVector{<:AbstractString};
                            tf=nothing, height::Int=15, width::Int=70,
                            title::String="FSP progress",
                            xscale::Symbol=:linear,
                            transform=identity)
    prev_nlines = Ref(0)

    function callback(t_now, t_log, traj)
        # Erase the previous plot by moving cursor up and clearing to end
        if prev_nlines[] > 0
            print(stdout, "\e[$(prev_nlines[])A\e[0J")
        end

        n = length(t_log)
        n_sp = length(species_names)
        if n < 2
            # Not enough data yet — print a placeholder
            msg = "  Waiting for first snapshot…  (t=$(round(t_now; sigdigits=3)))\n"
            print(stdout, msg)
            prev_nlines[] = 1
            flush(stdout)
            return
        end

        # Optionally transform x-axis to log10, filtering t=0
        if xscale === :log10
            mask = t_log .> 0
            xs   = log10.(t_log[mask])
            ys   = traj[mask, :]
            xlbl_base = tf !== nothing ?
                "log₁₀(t)  ($(round(100*t_now/tf; digits=1))% of tf=$(round(tf; sigdigits=3)))" :
                "log₁₀(t)"
        else
            xs   = t_log
            ys   = traj
            xlbl_base = tf !== nothing ?
                "t  ($(round(100*t_now/tf; digits=1))% of tf=$(round(tf; sigdigits=3)))" :
                "t"
        end

        if length(xs) < 2
            msg = "  Waiting for first snapshot…  (t=$(round(t_now; sigdigits=3)))\n"
            print(stdout, msg)
            prev_nlines[] = 1
            flush(stdout)
            return
        end

        # Apply optional user transform (e.g. rescale one species)
        ys = transform(ys)

        # Compute y-range across all species so no curve is clipped.
        # Clamp lower bound to 0 — populations are non-negative.
        y_min = minimum(ys)
        y_max = maximum(ys)
        y_pad = max((y_max - y_min) * 0.05, 1.0)
        ylims = (max(0.0, y_min - y_pad), y_max + y_pad)

        # Auto-size plot width to fit terminal so lines don't wrap.
        # A UnicodePlots line occupies: ~14 (y-label/tick) + 1 (│) + w + 1 (│) + legend ≈ w+20
        term_cols = displaysize(stdout)[2]
        w = min(width, max(30, term_cols - 22))

        # First species initialises the plot
        plt = UnicodePlots.lineplot(
            xs, ys[:, 1];
            name    = n_sp > 1 ? species_names[1] : "",
            xlabel  = xlbl_base,
            ylabel  = "E[X]",
            title   = title * "  t=$(round(t_now; sigdigits=4))",
            height  = height,
            width   = w,
            ylim    = ylims,
        )
        for i in 2:min(n_sp, size(ys, 2))
            UnicodePlots.lineplot!(plt, xs, ys[:, i]; name=species_names[i])
        end

        # Render to string so we can count lines.
        # Split on \n to get actual line count (accounts for trailing newline).
        buf = IOBuffer()
        show(buf, plt)
        plot_str = String(take!(buf))
        lines = split(plot_str, '\n')
        # println adds one more \n after plot_str
        nlines = length(lines) + (isempty(last(lines)) ? 0 : 1)

        println(stdout, plot_str)
        flush(stdout)
        prev_nlines[] = nlines
    end

    callback
end
