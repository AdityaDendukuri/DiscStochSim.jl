### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ a5d9129d-f389-48c0-8bb5-c0c110fc265c
begin
	using LinearOperators
	
	# --- toy active set + mapping (replace with yours) ---
	states = [CartesianIndex(0,0,0), CartesianIndex(1,0,0), CartesianIndex(0,1,0)]
	toidx  = Dict(s => i for (i,s) in enumerate(states))
	stoich = [CartesianIndex(1,0,0), CartesianIndex(0,1,0)]      # reactions S_k
	props  = [(x,p,t)->1.0, (x,p,t)->0.5]                        # α_k(x; rates,t)
	rates  = [1.0, 2.0];  t = 0.0
	m = length(states)
	
	# ---- fused in-place: y := α*A*x + β*y ----
	function cme_prod!(y::AbstractVector, x::AbstractVector, α::Number, β::Number)
	    if β == 0
	        fill!(y, 0.0)
	    elseif β != 1
	        @inbounds for i in eachindex(y); y[i] *= β; end
	    end
	    @inbounds for j in 1:m
	        xj = x[j]; xj == 0 && continue
	        sj = states[j]; outsum = 0.0
	        for k in eachindex(stoich)
	            i = get(toidx, sj + stoich[k], 0)
	            if i != 0
	                a = props[k](sj, rates, t)
	                if a != 0.0
	                    y[i] += α * a * xj   # off-diagonal
	                    outsum += a
	                end
	            end
	        end
	        if outsum != 0.0
	            y[j] += α * (-outsum) * xj   # diagonal
	        end
	    end
	    return y
	end
	
	# out-of-place matvec required by constructor (returns y = A*x)
	cme_prod(x) = (y = similar(x); cme_prod!(y, x, 1.0, 0.0))
	
	# Constructor: (T, m, n, symmetric, hermitian, prod, prod!)
	A = LinearOperator(Float64, m, m, false, false, cme_prod, cme_prod!)
	
	# --- use it ---
	p  = rand(m)
	y1 = A * p                    # out-of-place
	y2 = similar(p); mul!(y2, A, p)  # in-place
	
end

# ╔═╡ 686dd03b-3919-43c8-8587-a538a737a872


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearOperators = "5c8ed15e-5a4c-59e4-a42b-c7e8811fb125"

[compat]
LinearOperators = "~2.10.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "dc017c23f3299c0384317686d4b6f65a9f1504e0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearOperators]]
deps = ["FastClosures", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "TimerOutputs"]
git-tree-sha1 = "1894a798ed8887895c5ae70f1fe8331c0c1d8480"
uuid = "5c8ed15e-5a4c-59e4-a42b-c7e8811fb125"
version = "2.10.0"

    [deps.LinearOperators.extensions]
    LinearOperatorsAMDGPUExt = "AMDGPU"
    LinearOperatorsCUDAExt = "CUDA"
    LinearOperatorsChainRulesCoreExt = "ChainRulesCore"
    LinearOperatorsJLArraysExt = "JLArrays"
    LinearOperatorsLDLFactorizationsExt = "LDLFactorizations"
    LinearOperatorsMetalExt = "Metal"

    [deps.LinearOperators.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    JLArrays = "27aeb0d3-9eb9-45fb-866b-73c2ecf80fcb"
    LDLFactorizations = "40e66cde-538c-5869-a4ad-c39174c6795b"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═a5d9129d-f389-48c0-8bb5-c0c110fc265c
# ╠═686dd03b-3919-43c8-8587-a538a737a872
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
