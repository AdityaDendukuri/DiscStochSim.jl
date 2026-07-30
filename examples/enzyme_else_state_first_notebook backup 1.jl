### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 84414195-efbf-48a0-aef8-dd48216c62cf
begin
    using LinearAlgebra
    using Plots
    using Random
    using Statistics

    gr()
end

# ╔═╡ 87acfbdd-af24-4c63-a623-3cf52cb95de4
md"""
# State-first stochastic ELSE

Consider one enzyme moving rapidly through conformational states

\$\$
C_0 \rightleftharpoons C_1 \rightleftharpoons \cdots
\rightleftharpoons C_m,
\$\$

followed by rare product formation

\$\$
C_m \xrightarrow{\varepsilon} C_0+P.
\$\$

For a fixed product count \$P=p\$, the states \$C_0,\ldots,C_m\$ form one
transient subnetwork. ELSE replaces all fast internal reactions by one draw of:

1. the state from which the subnetwork escapes;
2. the escape time conditioned on that state.

After product formation, the enzyme resets to \$C_0\$ and the same precomputed
ELSE object is reused.
"""

# ╔═╡ 9f1b157f-482d-4b61-abaf-e049439a60cc
md"""
## 1. Killed subnetwork generator

We use the column convention

\$\$
\frac{d\boldsymbol p}{dt}=Q\boldsymbol p.
\$\$

The diagonal of \$Q\$ contains both internal reaction rates and rates leaving the
subnetwork. Therefore its column sums are nonpositive.
"""

# ╔═╡ b6b3f9a4-cd5a-436a-b693-c1460cf74ad8
function enzyme_generator(m, fast_rate, catalytic_rate)
    number_of_states = m + 1
    Q = zeros(number_of_states, number_of_states)
    exit_rate = zeros(number_of_states)
    internal_rate = zeros(number_of_states)

    for state in 0:m
        column = state + 1

        if state < m
            Q[column + 1, column] = fast_rate
            internal_rate[column] += fast_rate
        end

        if state > 0
            Q[column - 1, column] = fast_rate
            internal_rate[column] += fast_rate
        end

        if state == m
            exit_rate[column] = catalytic_rate
        end

        Q[column, column] =
            -(internal_rate[column] + exit_rate[column])
    end

    (Q=Q, exit_rate=exit_rate, internal_rate=internal_rate)
end

# ╔═╡ 3832d41b-fb73-44aa-a673-63771fab702c
md"""
## 2. Precompute the ELSE matrices

The fundamental matrix is

\$\$
Z=-Q^{-1}.
\$\$

Its entry \$Z_{ba}\$ is the expected time spent in state \$b\$ before escape when
the process starts in state \$a\$.

We also diagonalize \$Q\$ once:

\$\$
Q=U\operatorname{diag}(-\lambda)U^{-1}.
\$\$

This lets every subsequent time sample use scalar exponentials
\$e^{-\lambda_j t}\$ rather than recomputing \$e^{Qt}\$.
"""

# ╔═╡ 26c7e853-8475-4ac9-a88c-55c574cbb8b
function precompute_else(m, fast_rate, catalytic_rate)
    network = enzyme_generator(m, fast_rate, catalytic_rate)
    Q = network.Q
    Z = -inv(Q)

    decomposition = eigen(Q)
    U = ComplexF64.(decomposition.vectors)
    U_inverse = inv(U)
    λ = ComplexF64.(-decomposition.values)

    (
        m=m,
        Q=Q,
        Z=Z,
        U=U,
        U_inverse=U_inverse,
        λ=λ,
        exit_rate=network.exit_rate,
        internal_rate=network.internal_rate,
    )
end

# ╔═╡ 4417353c-076e-4a59-91ba-b8000b892c38
md"""
## 3. Sample the exit state first

Starting from state \$a\$, the probability of escaping from state \$b\$ is

\$\$
\Pr(B=b\mid X_0=a)=\Delta_b Z_{ba},
\$\$

where \$\Delta_b\$ is the total rate from \$b\$ to states outside the subnetwork.
"""

# ╔═╡ fbad24dd-8e91-40f0-9f4f-4ec08a016fb2
function sample_index(weights, rng)
    cutoff = rand(rng) * sum(weights)
    running_sum = 0.0

    for i in eachindex(weights)
        running_sum += weights[i]
        running_sum >= cutoff && return i
    end

    lastindex(weights)
end

# ╔═╡ 34e4a1b8-147a-4b73-ac71-c4c5317d8f21
function sample_exit_state(else_data, start_state, rng)
    weights =
        else_data.exit_rate .* else_data.Z[:, start_state]
    weights = max.(weights, 0.0)
    weights ./= sum(weights)

    sample_index(weights, rng)
end

# ╔═╡ 2dfca787-8311-4f05-81fb-5d1c5d89f408
md"""
## 4. Sample time conditional on the exit state

After selecting \$b\$, Equation 32 gives the conditional survival probability

\$\$
\Pr(\tau>T\mid B=b,X_0=a)
=
\frac{
\left[
U\operatorname{diag}
\left(e^{-\lambda T}/\lambda\right)
U^{-1}
\right]_{ba}
}{
Z_{ba}
}.
\$\$

We draw \$u\sim\operatorname{Uniform}(0,1)\$ and solve

\$\$
\Pr(\tau>T\mid B=b,X_0=a)=u
\$\$

by bisection.
"""

# ╔═╡ 8fc63f11-9a81-4e03-805d-325bb04186cf
function conditional_survival(else_data, start_state, exit_state, time)
    coefficients =
        else_data.U[exit_state, :] .*
        else_data.U_inverse[:, start_state] ./
        else_data.λ

    numerator =
        sum(coefficients .* exp.(-else_data.λ .* time))

    survival =
        real(numerator / else_data.Z[exit_state, start_state])

    clamp(survival, 0.0, 1.0)
end

# ╔═╡ 5122ec41-16c6-4d32-97a2-6fd755cf7db5
function sample_conditional_time(else_data, start_state, exit_state, rng)
    target_survival = max(rand(rng), eps(Float64))

    Z = else_data.Z
    conditional_mean =
        (Z * Z)[exit_state, start_state] /
        Z[exit_state, start_state]

    low = 0.0
    high = conditional_mean

    while conditional_survival(
        else_data, start_state, exit_state, high
    ) > target_survival
        high *= 2
    end

    for _ in 1:55
        middle = (low + high) / 2
        survival = conditional_survival(
            else_data, start_state, exit_state, middle
        )

        if survival > target_survival
            low = middle
        else
            high = middle
        end
    end

    (low + high) / 2
end

# ╔═╡ eb6f5fbf-b7d6-4298-a73e-36d6df11cb57
md"""
## 5. One stochastic ELSE product event

This is the complete state-first ordering:

1. sample \$B\$ using \$Z\$;
2. sample \$\tau\mid B\$ using Equation 32;
3. form one product and reset the enzyme to \$C_0\$.
"""

# ╔═╡ 1305a3bd-706a-4df9-a013-2883f411f624
function else_product_event(else_data, start_state, rng)
    exit_state =
        sample_exit_state(else_data, start_state, rng)

    waiting_time = sample_conditional_time(
        else_data, start_state, exit_state, rng
    )

    (
        waiting_time=waiting_time,
        exit_state=exit_state,
        next_state=1,  # C₀ after catalysis
    )
end

# ╔═╡ 929eb0f6-f746-499c-a9ec-96e0349ba8e4
md"""
## 6. Sampling loop replacing SSA

One loop iteration produces one product. It replaces all binding, unbinding, and
conformational reactions that SSA would have simulated before that product.
"""

# ╔═╡ aab809b2-b140-4ad6-b13f-ad6bf4cecc31
function else_product_trajectory(else_data, final_time; seed=1)
    rng = MersenneTwister(seed)
    time = 0.0
    product = 0
    enzyme_state = 1  # C₀

    times = Float64[time]
    products = Int[product]

    while time < final_time
        event =
            else_product_event(else_data, enzyme_state, rng)
        next_time = time + event.waiting_time

        next_time > final_time && break

        time = next_time
        product += 1
        enzyme_state = event.next_state

        push!(times, time)
        push!(products, product)
    end

    (times=times, products=products)
end

# ╔═╡ 7de2ee20-284a-4923-ae6d-c8dda3a25b06
md"""
## 7. Benchmark parameters

With \$m=10\$, \$k=50\,\mathrm{s}^{-1}\$, and
\$\varepsilon=10^{-2}\,\mathrm{s}^{-1}\$, SSA performs approximately

\$\$
\frac{2mk}{\varepsilon}=100{,}000
\$\$

fast internal reactions per product.
"""

# ╔═╡ 1d4f42a1-c774-4645-96a1-7cc2a6c24d79
begin
    m = 10
    fast_rate = 50.0
    catalytic_rate = 0.01

    else_data = precompute_else(m, fast_rate, catalytic_rate)
    start_state = 1  # C₀

    mean_escape_time = sum(else_data.Z[:, start_state])
    expected_fast_events = dot(
        else_data.internal_rate,
        else_data.Z[:, start_state],
    )

    (
        states=m + 1,
        mean_escape_time=mean_escape_time,
        expected_fast_events=expected_fast_events,
    )
end

# ╔═╡ b36c42db-5864-43e1-a6d7-aaac6e87ee66
md"""
## 8. Verify the sampled escape-time distribution
"""

# ╔═╡ c8b22f98-c35c-48b0-8e14-89c91a8a155c
begin
    number_of_waiting_times = 10_000
    waiting_time_rng = MersenneTwister(2026)

    sampled_waiting_times = [
        begin
            exit_state =
                sample_exit_state(else_data, start_state, waiting_time_rng)
            sample_conditional_time(
                else_data,
                start_state,
                exit_state,
                waiting_time_rng,
            )
        end
        for _ in 1:number_of_waiting_times
    ]
end

# ╔═╡ 3b2cf26c-bf0c-4cfc-a310-7ba6f048509e
begin
    exit_state = m + 1
    time_grid = range(
        0.0,
        5 * mean_escape_time;
        length=300,
    )

    exact_cdf = [
        1 - conditional_survival(
            else_data, start_state, exit_state, time
        )
        for time in time_grid
    ]

    empirical_cdf = [
        count(<=(time), sampled_waiting_times) /
        number_of_waiting_times
        for time in time_grid
    ]

    plot(
        time_grid,
        exact_cdf;
        label="Equation 32",
        xlabel="product waiting time",
        ylabel="CDF",
        title="State-first stochastic ELSE",
        lw=3,
    )
    plot!(
        time_grid,
        empirical_cdf;
        label="10,000 ELSE samples",
        lw=2,
        ls=:dash,
    )
end

# ╔═╡ a091e699-d59c-4f37-b33b-9ae5dd1c8e5a
(
    matrix_mean=mean_escape_time,
    sampled_mean=mean(sampled_waiting_times),
    sampled_standard_deviation=std(sampled_waiting_times),
)

# ╔═╡ 278bb0ca-38a5-44c0-b0fd-3a83927697b4
md"""
## 9. Visualize many product trajectories

Each trajectory is stochastic because Equation 32 supplies a new waiting-time
draw for every catalytic event.
"""

# ╔═╡ 6f96934f-4b69-4dff-b14d-568a91a27be9
begin
    number_of_trajectories = 500
    final_time = 5 * mean_escape_time

    product_trajectories = [
        else_product_trajectory(
            else_data,
            final_time;
            seed=10_000 + sample,
        )
        for sample in 1:number_of_trajectories
    ]
end

# ╔═╡ 4e566f18-c7d8-4e2f-b285-354d2b7a735b
let
    figure = plot(
        xlabel="time",
        ylabel="product count P",
        title="$(length(product_trajectories)) stochastic ELSE trajectories",
        legend=false,
        size=(850, 500),
    )

    for trajectory in product_trajectories
        plot!(
            figure,
            trajectory.times,
            trajectory.products;
            seriestype=:steppost,
            color=:steelblue,
            alpha=0.08,
            lw=0.8,
        )
    end

    figure
end

# ╔═╡ 830994b7-4f30-4091-96ea-c7c9b180ac01
md"""
## Interpretation

- The exit state is sampled first from \$\Delta_bZ_{ba}\$.
- The time is then sampled from Equation 32 conditional on that state.
- \$Z\$, \$U\$, \$U^{-1}\$, and \$\lambda\$ are computed once.
- Sampling does not repeatedly construct the matrix exponential \$e^{Qt}\$.
- The per-event bisection evaluates only scalar exponentials
  \$e^{-\lambda_jt}\$.
- The fast internal enzyme reactions are integrated out exactly.
"""

# ╔═╡ Cell order:
# ╠═84414195-efbf-48a0-aef8-dd48216c62cf
# ╟─87acfbdd-af24-4c63-a623-3cf52cb95de4
# ╟─9f1b157f-482d-4b61-abaf-e049439a60cc
# ╠═b6b3f9a4-cd5a-436a-b693-c1460cf74ad8
# ╟─3832d41b-fb73-44aa-a673-63771fab702c
# ╠═26c7e853-8475-4ac9-a88c-55c574cbb8b
# ╟─4417353c-076e-4a59-91ba-b8000b892c38
# ╠═fbad24dd-8e91-40f0-9f4f-4ec08a016fb2
# ╠═34e4a1b8-147a-4b73-ac71-c4c5317d8f21
# ╟─2dfca787-8311-4f05-81fb-5d1c5d89f408
# ╠═8fc63f11-9a81-4e03-805d-325bb04186cf
# ╠═5122ec41-16c6-4d32-97a2-6fd755cf7db5
# ╟─eb6f5fbf-b7d6-4298-a73e-36d6df11cb57
# ╠═1305a3bd-706a-4df9-a013-2883f411f624
# ╟─929eb0f6-f746-499c-a9ec-96e0349ba8e4
# ╠═aab809b2-b140-4ad6-b13f-ad6bf4cecc31
# ╟─7de2ee20-284a-4923-ae6d-c8dda3a25b06
# ╠═1d4f42a1-c774-4645-96a1-7cc2a6c24d79
# ╟─b36c42db-5864-43e1-a6d7-aaac6e87ee66
# ╠═c8b22f98-c35c-48b0-8e14-89c91a8a155c
# ╠═3b2cf26c-bf0c-4cfc-a310-7ba6f048509e
# ╠═a091e699-d59c-4f37-b33b-9ae5dd1c8e5a
# ╟─278bb0ca-38a5-44c0-b0fd-3a83927697b4
# ╠═6f96934f-4b69-4dff-b14d-568a91a27be9
# ╠═4e566f18-c7d8-4e2f-b285-354d2b7a735b
# ╟─830994b7-4f30-4091-96ea-c7c9b180ac01
