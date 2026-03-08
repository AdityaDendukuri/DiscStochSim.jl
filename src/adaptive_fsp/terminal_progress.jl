# Helper: render two plot strings side by side in the terminal.
function _side_by_side(str1::String, str2::String, sep::String="  ")
    lines1 = split(str1, '\n')
    lines2 = split(str2, '\n')
    # Drop trailing blank lines
    while !isempty(lines1) && isempty(strip(last(lines1))); pop!(lines1); end
    while !isempty(lines2) && isempty(strip(last(lines2))); pop!(lines2); end
    n  = max(length(lines1), length(lines2))
    w1 = isempty(lines1) ? 0 : maximum(length.(lines1))
    while length(lines1) < n; push!(lines1, ""); end
    while length(lines2) < n; push!(lines2, ""); end
    join([rpad(l1, w1) * sep * l2 for (l1, l2) in zip(lines1, lines2)], '\n')
end

"""
    terminal_progress(species_names; tf=nothing, height=7, width=70,
                      title="FSP progress", xscale=:linear, transform=identity)

Return a progress callback for `AdaptiveFSP` that renders a live 2×2 grid of plots:

- **Top-left**: mean species trajectories E[X]
- **Top-right**: state space size |S|
- **Bottom-left**: time step dt
- **Bottom-right**: total flux Φ_total

All panels share the same x-axis (simulation time, optionally log-scaled).

# Usage
```julia
cb = terminal_progress(["A", "B", "C"]; tf=1e4, xscale=:log10)
sol, diag = solve(prob, AdaptiveFSP(save_interval=500, progress_callback=cb))
```

Set `xscale=:log10` for systems with large time horizons (Robertson, Bottleneck).

Use `transform` to pre-process the trajectory matrix before plotting — useful for
rescaling species with very different magnitudes:

```julia
cb = terminal_progress(["A", "1000×B", "C"]; tf=1e4, xscale=:log10,
                       transform = m -> (m2 = copy(m); m2[:,2] .*= 1000; m2))
```
"""
function terminal_progress(species_names::AbstractVector{<:AbstractString};
                            tf=nothing, height::Int=7, width::Int=70,
                            title::String="FSP progress",
                            xscale::Symbol=:linear,
                            transform=identity)
    prev_nlines = Ref(0)

    function callback(t_now, t_log, traj, extra=nothing)
        # Erase the previous frame
        if prev_nlines[] > 0
            print(stdout, "\e[$(prev_nlines[])A\e[0J")
        end

        n    = length(t_log)
        n_sp = length(species_names)

        if n < 2
            msg = "  Waiting for first snapshot…  (t=$(round(t_now; sigdigits=3)))\n"
            print(stdout, msg)
            prev_nlines[] = 1
            flush(stdout)
            return
        end

        # --- x-axis transform (shared by all panels) ---
        pct_str = tf !== nothing ? "  $(round(100*t_now/tf; digits=1))% of tf=$(round(tf; sigdigits=3))" : ""
        if xscale === :log10
            # Snapshot x (for traj panel)
            smask = t_log .> 0
            xs_snap = log10.(t_log[smask])
            ys_snap = traj[smask, :]
            xlbl    = "log₁₀(t)$pct_str"
            # Per-step x (for dt / size / flux panels)
            if extra !== nothing
                pmask   = extra.step_t .> 0
                xs_step = log10.(extra.step_t[pmask])
            end
        else
            xs_snap = t_log
            ys_snap = traj
            xlbl    = "t$pct_str"
            if extra !== nothing
                xs_step = extra.step_t
                pmask   = trues(length(xs_step))
            end
        end

        if length(xs_snap) < 2
            msg = "  Waiting for first snapshot…  (t=$(round(t_now; sigdigits=3)))\n"
            print(stdout, msg)
            prev_nlines[] = 1
            flush(stdout)
            return
        end

        # Apply user transform to trajectory
        ys_snap = transform(ys_snap)

        # --- auto-size plot width to avoid terminal line-wrapping ---
        term_cols = displaysize(stdout)[2]
        # Two plots side by side: each gets roughly half the terminal
        w = min(width, max(20, term_cols ÷ 2 - 24))

        # --- y-range for trajectory panel ---
        y_min = minimum(ys_snap)
        y_max = maximum(ys_snap)
        y_pad = max((y_max - y_min) * 0.05, 1.0)
        ylims_traj = (max(0.0, y_min - y_pad), y_max + y_pad)

        # ── Panel 1: mean trajectories ────────────────────────────────────────
        plt1 = UnicodePlots.lineplot(
            xs_snap, ys_snap[:, 1];
            name   = n_sp > 1 ? species_names[1] : "",
            xlabel = xlbl, ylabel = "E[X]",
            title  = title * "  t=$(round(t_now; sigdigits=4))",
            height = height, width = w,
            ylim   = ylims_traj,
        )
        for i in 2:min(n_sp, size(ys_snap, 2))
            UnicodePlots.lineplot!(plt1, xs_snap, ys_snap[:, i]; name=species_names[i])
        end

        buf1 = IOBuffer(); show(buf1, plt1); str1 = String(take!(buf1))

        if extra !== nothing && length(xs_step) >= 2
            # ── Panel 2: state space size |S| ────────────────────────────────
            ys_sz = Float64.(extra.size_log[pmask])
            plt2  = UnicodePlots.lineplot(
                xs_step, ys_sz;
                xlabel = xlbl, ylabel = "|S|",
                title  = "State space size",
                height = height, width = w,
            )
            buf2 = IOBuffer(); show(buf2, plt2); str2 = String(take!(buf2))

            # ── Panel 3: time step dt ─────────────────────────────────────────
            ys_dt = extra.dt_log[pmask]
            plt3  = UnicodePlots.lineplot(
                xs_step, ys_dt;
                xlabel = xlbl, ylabel = "dt",
                title  = "Time step",
                height = height, width = w,
            )
            buf3 = IOBuffer(); show(buf3, plt3); str3 = String(take!(buf3))

            # ── Panel 4: total flux Φ_total ───────────────────────────────────
            ys_fl = extra.flux_log[pmask]
            plt4  = UnicodePlots.lineplot(
                xs_step, ys_fl;
                xlabel = xlbl, ylabel = "Φ",
                title  = "Total flux Φ_total",
                height = height, width = w,
            )
            buf4 = IOBuffer(); show(buf4, plt4); str4 = String(take!(buf4))

            # ── Combine into 2×2 grid ─────────────────────────────────────────
            top = _side_by_side(str1, str2)
            bot = _side_by_side(str3, str4)
            plot_str = top * "\n" * bot
        else
            plot_str = str1
        end

        # Count lines for next erase (split is more reliable than counting \n)
        lines    = split(plot_str, '\n')
        nlines   = length(lines) + (isempty(last(lines)) ? 0 : 1)

        println(stdout, plot_str)
        flush(stdout)
        prev_nlines[] = nlines
    end

    callback
end
