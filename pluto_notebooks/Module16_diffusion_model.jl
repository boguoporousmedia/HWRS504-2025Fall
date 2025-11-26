### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 1094fba2-1de3-422d-b235-6b3756e0f9c3
using PlutoUI

# ╔═╡ bda33902-c9a2-11f0-b1af-713a5b2fa60a
md"""
# Module 16: Diffusion model
"""

# ╔═╡ 518efc97-8315-40ca-b12f-c18ff82ca468
md"""
Key references:

Song, Y., Sohl-Dickstein, J., Kingma, D.P., Kumar, A., Ermon, S. and Poole, B., 2020. Score-based generative modeling through stochastic differential equations. arXiv preprint arXiv:2011.13456.

Yang, L., Zhang, Z., Song, Y., Hong, S., Xu, R., Zhao, Y., Zhang, W., Cui, B. and Yang, M.H., 2023. Diffusion models: A comprehensive survey of methods and applications. ACM computing surveys, 56(4), pp.1-39.

Anderson, B.D., 1982. Reverse-time diffusion equation models. Stochastic Processes and their Applications, 12(3), pp.313-326.

Vincent, P., 2011. A connection between score matching and denoising autoencoders. Neural computation, 23(7), pp.1661-1674.
"""

# ╔═╡ 9af4c5c2-d3ff-4367-8e28-9628ba24f190
md"""
### What Is a Diffusion Model?

Diffusion models are a family of probabilistic generative models that progressively destruct data by injecting noise, then learn to reverse this process for sample generation.

"""

# ╔═╡ 2d152d61-fb5b-46fc-b75c-4e3b6043dcbb
local img = LocalResource("./figs/mod16_diffusion_model_concept.png",:width => "800px")

# ╔═╡ f35565f4-c9b4-43d2-b190-f137f3ecf9e5
md"""
### Forward Diffusion Process as an stochastic differential equation (SDE)

We start from a clean data sample

```math
\mathbf{x}_0 \sim p_0(\mathbf{x})
```

and gradually corrupt it with Gaussian noise until we reach pure noise

```math
\mathbf{x}_T \sim p_T(\mathbf{x}).
```

This corruption process can be written as a SDE:

```math
\mathrm{d}\mathbf{x} = \mathbf{f}(\mathbf{x},t)\,\mathrm{d}t + g(t)\,\mathrm{d}\mathbf{W}.
```

* Drift function:
  Controls how the sample moves deterministically through time.

* Diffusion function:
  Controls the strength of the noise injected at each time.

* Wiener process (``\mathbf{W}(t)``):
  A continuous-time stochastic process whose increments are Gaussian,
  representing “adding Gaussian noise at each step.”

"""


# ╔═╡ 21c42274-e025-44f4-8535-d839a0d61010
md"""
### Reverse SDE

This diffusion process can be reversed by solving another SDE backwards in time (Anderson, 1982), starting from a noise sample ``\mathbf{x}_T``.

The reverse-time SDE is:

```math
\mathrm{d}\mathbf{x}
=
\big[
  f(\mathbf{x},t)
  - g(t)^2 \nabla_\mathbf{x} \log p_t(\mathbf{x})
\big]\,\mathrm{d}t
+
g(t)\,\mathrm{d}\overline{\mathbf{W}},
```

where:

* ``p_t(\mathbf{x})`` is the probability density of ``\mathbf{x}_t``,
* ``\nabla_\mathbf{x} \log p_t(\mathbf{x})`` is the score,
* ``\overline{\mathbf{W}}`` is a Wiener process in reverse time.

The score term guides the stochastic dynamics back toward the data distribution,
allowing us to transform pure noise into realistic samples.

"""

# ╔═╡ ec943d20-724f-4d15-8bbb-be08264d989c
md"""
### Probability Flow ODE

Song et al. (2020) showed that the score-based models enable another numerical method for solving the reverse-time SDE. For all diffusion processes, there exists a corresponding deterministic process whose trajectories share the same marginal probability densities ``\{p_t(\mathbf{x})\}_{t=0}^T`` as the SDE.

This deterministic process satisfies the ODE:

```math
\mathrm{d}\mathbf{x}
=
\Big[
  f(\mathbf{x}, t)
  - \frac{1}{2} g(t)^2 \nabla_{\mathbf{x}} \log p_t(\mathbf{x})
\Big] \,\mathrm{d}t.
```

This is referred to as the probability flow ODE, which can be determined from the SDE once the score function is known. When the score function is approximated using a time-dependent score-based neural network, this becomes an example of a Neural ODE.
"""

# ╔═╡ cd9f5b0b-5d98-49de-8818-9d957e950353
md"""
To generate an image:

1. Sample ``\mathbf{x}_T \sim p_T(\mathbf{x})``(usually a Gaussian).

2. Solve the ODE (or reverse SDE) backwards from ``t = T`` to ``t = 0``:

```math
\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t}
=
f(\mathbf{x},t)
-
\frac{1}{2} g(t)^2 \nabla_\mathbf{x} \log p_t(\mathbf{x}).
```

However, we don’t know the true density ``p_t(\mathbf{x})``. **Idea**: Learn the "score function" using a neural network:

```math
\mathbf{s}(\mathbf{x}, t; \mathbf{\theta}) \;\approx\; \nabla_\mathbf{x} \log p_t(\mathbf{x}).
```

This leads to a neural ODE

```math
\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t}
=
f(\mathbf{x},t)
-
\frac{1}{2} g(t)^2 \mathbf{s}(\mathbf{x}, t; \mathbf{\theta}).
```

One can solve this neural ODE to generate an image.
"""

# ╔═╡ 97c0d982-438e-4648-bf74-af7e6cc87e8f
md"""
### Learning the Score Function (Diffusion Model)

The neural network needs to match the true score function, i.e.

```math
\mathcal{L}(\mathbf{\theta})
=
\mathbb{E}_{t,\, \mathbf{x}_t \sim p_t}
\Big[
\lVert s(\mathbf{x}_t, t; \mathbf{x}) - \nabla_{\mathbf{x}_t} \log p_t(\mathbf{x}_t) \rVert^2
\Big].
```

Vincent (2010) showed that this is equivalent to

```math
\mathcal{L}(\mathbf{\theta})
=
\mathbb{E}_{t,\, \mathbf{x}_0 \sim p_0,\, \mathbf{x}_t \sim p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0)}
\Big[
\lVert s(\mathbf{x}_t, t; \theta) - \nabla_{\mathbf{x}_t} \log p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0) \rVert^2
\Big]
+ C,
```

where ``p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0)`` is the transition probability and ``C`` is a constant.

Take the simple forward SDE ``\mathrm{d}x = \mathrm{d}\mathbf{W}`` as an example, the transition probability is

```math
p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0)
=
\frac{1}{\sqrt{2\pi t}}
\exp\!\left( -\frac{\lVert \mathbf{x}_t - \mathbf{x}_0 \rVert^2}{2t} \right),
```

which implies

```math
\nabla_{\mathbf{x}_t} \log p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0)
=
-\frac{1}{t}\,(\mathbf{x}_t - \mathbf{x}_0).
```

Thus the score-matching loss becomes

```math
\mathcal{L}(\theta)
=
\mathbb{E}_{t,\, \mathbf{x}_0 \sim p_0,\, \mathbf{x}_t \sim p_{0t}(\mathbf{x}_t \mid \mathbf{x}_0)}
\Big[
\lVert s(\mathbf{x}_t, t; \theta)
-
\frac{1}{t}(\mathbf{x}_0 - \mathbf{x}_t)
\rVert^2
\Big]
+ C.
```

**Interpretation:**
For this simple diffusion, the score function is just predicting the noise added to the image.

"""

# ╔═╡ f154c044-fd1d-4d53-8796-6dccb4762c23
md"""
### Applications of Diffusion Models in Earth Science

Diffusion models are not only useful for images and audio; they are a powerful way to model high-dimensional, non-Gaussian, stochastic fields, which are very common in Earth science.

Here are few potential application areas.

1. Meteorological Forcing & Hydroclimate Fields

Use diffusion models to generate or downscale:
- precipitation fields (storm structures, spatial rainfall patterns),
- temperature and evapotranspiration fields,
- stochastic meteorological scenarios for hydrologic models.

Rationale: These fields are highly non-Gaussian and spatio-temporally correlated. Diffusion models can learn their joint distribution and produce realistic ensembles for uncertainty analysis and climate-impact studies.

2. Hydrologic Time-Series Generation

Model and generate:
- streamflow hydrographs,
- groundwater head time series,
- snowpack or soil moisture evolution.

Rationale: Diffusion models can produce diverse, realistic realizations that respect observed statistics and extremes, useful for risk analysis and scenario generation.

3. Subsurface Property Fields

Generate 2-D/3-D subsurface structures:
- heterogeneous hydraulic conductivity and porosity fields,
- facies models and stratigraphy,
- pore-scale geometries from micro-CT data.

Rationale: Learn complex geological priors far beyond Gaussian random fields, improving stochastic groundwater modeling, uncertainty quantification, and pore-scale simulations.

4. Inverse Problems and Data Assimilation

Use diffusion models as powerful priors or conditional generators in:
- hydraulic conductivity inversion,
- contaminant source identification,
- groundwater model calibration,
- reconstructing missing or sparse observations (e.g., heads, soil moisture).

Rationale: Ill-posed inverse problems become better constrained when the unknown fields are sampled from a learned, physically realistic prior.

5. Remote Sensing Super-Resolution and Downscaling

Apply diffusion models to satellite products:
- precipitation and soil moisture downscaling,
- snow water equivalent and land-surface temperature,
- GRACE-type groundwater anomaly fields,
- SAR denoising and gap filling.

Rationale: Diffusion models are state-of-the-art in image super-resolution and inpainting, and the same ideas transfer directly to geophysical raster data.

6. Physics-Aware Diffusion Models for PDE Solutions

Combine diffusion models with PDEs describing:
- groundwater flow,
- advection–dispersion and reactive transport,
- multiphase flow in porous media

Potential ideas:
- Learn the distribution of PDE solutions under uncertain parameters.
- Use score-based models as surrogates or “generative solvers” for expensive simulations.
- Explore operator-learning versions of diffusion models for mapping inputs (parameters, forcings) to solution fields.

**Take-Home Message**

Diffusion models provide a flexible, high-dimensional probabilistic framework that is naturally compatible with:
- spatial fields (maps, 3-D volumes),
- temporal processes (time series),
- physics-based PDE models.

This makes them a promising tool for stochastic hydrology, subsurface characterization, and Earth system modeling, with potential roles in both forward simulation (generating fields) and inverse problems (inferring hidden states and parameters).
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "afb0c26499386d0dcd10f3ec249a5f3ceff09bd1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"
"""

# ╔═╡ Cell order:
# ╟─1094fba2-1de3-422d-b235-6b3756e0f9c3
# ╟─bda33902-c9a2-11f0-b1af-713a5b2fa60a
# ╟─518efc97-8315-40ca-b12f-c18ff82ca468
# ╟─9af4c5c2-d3ff-4367-8e28-9628ba24f190
# ╟─2d152d61-fb5b-46fc-b75c-4e3b6043dcbb
# ╟─f35565f4-c9b4-43d2-b190-f137f3ecf9e5
# ╟─21c42274-e025-44f4-8535-d839a0d61010
# ╟─ec943d20-724f-4d15-8bbb-be08264d989c
# ╟─cd9f5b0b-5d98-49de-8818-9d957e950353
# ╟─97c0d982-438e-4648-bf74-af7e6cc87e8f
# ╟─f154c044-fd1d-4d53-8796-6dccb4762c23
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
