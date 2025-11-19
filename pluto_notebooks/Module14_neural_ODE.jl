### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 6c6b08ad-4561-421f-949d-4c68d4b9e4a6
using PlutoUI

# ╔═╡ 5cd32a4e-c437-11f0-b519-193be7197e21
md"""
# Module 14: Neural ODE
"""

# ╔═╡ 23fc87b2-8e84-48cc-8f28-3b446320d1a1
md"""
### Brief Introduction to ResNets

Deep neural networks traditionally stack many layers of the form

```math
h_{k+1} = \sigma(W_k h_k + b_k),
```

where each layer applies a nonlinear transformation. However, as depth increases, these networks become hard to train — gradients vanish, optimization becomes unstable, and performance saturates.

**Residual Networks (ResNets)** solve this by changing the layer update to:

```math
h_{k+1} = h_k + F_\theta(h_k),
```

where ``F_\theta`` represents the mapping by a small block of neural networks (e.g., a two-layer block is used as an example). This structure adds a skip connection from ``h_k`` to ``h_{k+1}``.

"""

# ╔═╡ c84eda38-a798-4333-a2f8-6e547f4a8bf5
local img = LocalResource("./figs/mod14_resnetblock.png",:width => "500px")

# ╔═╡ c7dbd776-0473-4363-890b-eb548a306fdc
md"""
#### Why this works

The skip connection makes each layer learn only a *small change* (a residual) to the previous state, instead of the entire transformation. The residual connection resolves the problem of vanishing gradients in deep neural networks, and stabilizes the training and convergence of deep neural networks with hundreds of layers.

#### ResNet vs. forward Euler

The ResNet update

```math
h_{k+1} = h_k + F_\theta(h_k)
```

looks exactly like a forward Euler applied to an ODE:

```math
\frac{dh}{dt} = f_\theta(h,t).
```

If we view the layer index ``k`` as a discrete time step, the step size as ``Δt = 1``, and ``F_\theta`` as ``Δt · f_\theta``, then:

```math
h_{k+1} = h_k + Δt \, f_\theta(h_k, t_k).
```

This suggests that as we make the layers smaller and more numerous, ResNets approximate the solution of a continuous ODE.

#### This in part motivates the concept of neural ODEs

Neural ODEs replace the discrete residual update with a continuous-time differential equation

```math
\frac{dh}{dt} = f_\theta(h(t), t).
```

Instead of stacking many residual layers, we let a classical ODE solver integrate the dynamics defined by a neural network. This gives a model whose depth is continuous, adaptive, and mathematically interpretable as a dynamical system.

A ResNet is essentially forward Euler applied to a learned function (or vector field); Neural ODEs take the natural limit and learn the continuous-time dynamics directly.

"""

# ╔═╡ d02de953-3b02-4ffd-96ef-92ddee7cd935
md"""
### What is a Neural ODE?

#### Continuous-depth model

A **neural ODE** defines a continuous evolution of ``h(t)``
```math
\frac{dh}{dt} = f_\theta(h(t), t), \quad h(0) = h_0,
```
where ``f_\theta`` is a neural network in ``(h,t)``.

The output of the model at time ``T`` is
```math
h(T) = \text{ODESolve}\big(h_0, f_\theta, t \in [0,T]\big).
```

Conceptually:

* Replace **discrete layers** by a **learned vector field** ``f_\theta``.
* Use a **numerical ODE solver** (i.e., Euler, RK4, adaptive solvers) to get ``h(T)``.

#### Computing the forward pass

Given parameters ``\theta``, we integrate:

* Input ``x`` → initial state ``h_0 = g_\text{enc}(x)`` (optional encoder).
* Solve ODE: ``h(T) = \Phi_\theta(T,0) h_0``.
* Predict: ``\hat{y} = g_\text{dec}(h(T))``.

So the forward map is
```math
x \mapsto h_0 \mapsto h(T) \mapsto \hat{y},
```
where the middle step ``h_0 \mapsto h(T)``, i.e., ``h(T) = \Phi_\theta(T,0) h_0``, is a numerical solver for a learned ODE.

``\Phi_\theta(T,0)`` is the flow operator of the learned ODE (forcing represented by a neural network), i.e, a mapping that takes an initial state at time 0
and returns the state at time T. 

### How to train the neural ODE? The adjoint method 

To minimize a loss
```math
\mathcal{L}(\theta) = \sum_i \ell\big(h_i(T; \theta), y_i\big).
```

we need ``\frac{d\mathcal{L}}{d\theta}``.

Define the adjoint state as
```math
a(t) = \frac{\partial \mathcal{L}}{\partial h(t)}.
```

One can show (via calculus of variations / backprop through time) that ``a(t)`` satisfies another ODE:
```math
\frac{da}{dt} = - a^\top \frac{\partial f_\theta}{\partial h}(h(t), t), \quad a(T) = \frac{\partial \mathcal{L}}{\partial h(T)}.
```

NOTE: This can be shown by defining a new loss function with a(t) as a Lagrange multipler.

```math
J = \mathcal{L}(\theta) + \int_0^T a(t)^\top \left( \frac{dh(t)}{dt} - f_\theta (h(t), t) \right) \, dt.
```

* We integrate this *backwards in time* from ``t=T`` to ``t=0``.
* While doing so, we can accumulate the gradient w.r.t. parameters:
  ```math
  \frac{d\mathcal{L}}{d\theta}
  = - \int_T^0 a(t)^\top \frac{\partial f_\theta}{\partial \theta}(h(t), t) \, dt.
  ```

**Insights**

* Forward pass: solve the ODE for the state ``h(t)``.
* Backward pass: solve another ODE backward for the adjoint ``a(t)``.
* This is analogous to reverse-mode AD, but in **continuous** time.

"""

# ╔═╡ ae404fc0-d0e0-45ce-8c05-7eeefdb20fa8
md"""
### Neural ODEs in the context of PDEs

#### PDE can be converted to ODE via semi-discretization

Take a PDE
```math
\frac{\partial u}{\partial t} = \mathcal{N}(u, x, t),
```
discretize in space (finite difference / finite volume / finite element) with the degree of freedom ``n`` leads to an ODE system for ``u(t) \in \mathbb{R}^n``
```math
\frac{du}{dt} = F(u(t), t).
```

where ``F`` is derived from physics and discretization.

#### Where do neural ODEs fit?

1. **Learn unknown or partially known dynamics**

   Suppose the true dynamics are
   ```math
   \frac{du}{dt} = F_\text{phys}(u,t;\theta_\text{phys}) + F_\text{NN}(u,t;\theta_\text{NN}),
   ```
   where ``F_\text{phys}`` is a known PDE discretization, and ``F_\text{NN}`` represents closure terms, subgrid physics, or unmodeled physics.

   This can be viewed as a hybrid neural ODE:
   ```math
   \frac{du}{dt} = f_\Theta(u,t), \quad \Theta = (\theta_\text{phys}, \theta_\text{NN}),
   ```
   trained from data ``u(t)``.

2. **Latent-space neural ODE (model reduction)**

   * Use an encoder ``z = E(u)`` to reduce dimension.
   * Evolve latent state with a neural ODE:
     ```math
     \frac{dz}{dt} = f_\theta(z,t).
     ```
   * Decode: ``\hat{u}(t) = D(z(t))``.

   This yields a learned reduced-order model for PDE dynamics.

3. **Learned time integrator**

   Instead of learning ``du/dt`` (instantaneous rate of change or the velocity field), learn the time-advance operator.

   * A traditional solver advances ``u^n \to u^{n+1}`` using discretized physics.
   * A neural ODE learns a continuous-time flow such that
     ```math
     u^{n+1} \approx \Phi_\theta(\Delta t, 0) u^n.
     ```
   * For a fixed ``\Delta t``, this reduces to ResNet; for variable ``\Delta t``, this learns a continuous flow map.

"""

# ╔═╡ b272d6a2-6352-4b50-a8f2-186ef8092790
md"""
### Examples

#### Scalar ODE

```math
\frac{du}{dt} = \lambda u, \quad u(0) = u_0,
```
with solution
```math
u(t) = u_0 e^{\lambda t}.
```

Prepare data: generate (noisy) samples ``{(t_i, u(t_i))}`` from this system.

Neural ODE setup:

* State ``h(t) \in \mathbb{R}``.
* Vector field ``f_\theta(h,t) = \text{NN}_\theta(h,t)``.
* Loss:
  ```math
  \mathcal{L}(\theta) = \sum_i \big| h_\theta(t_i) - u(t_i)\big|^2,
  ```
  where ``h_\theta(t_i)`` is obtained by solving the neural ODE from ``t=0`` to ``t=t_i``.

#### Diffusion equation

```math
u_t = D u_{xx}
``` 
where ``x \in [0,1]``, with periodic BCs.

* Spatially discretize with finite differences ``→`` ODE
  ```math
  \frac{du}{dt} = F(u).
  ```
* Learn ``F(u)`` with a neural network.

"""

# ╔═╡ de08d098-067e-405e-a305-39bfa259c226
md"""

### Practical issues & trade-offs

#### Solver choice and stiffness

For PDEs, the semi-discrete ODEs are often stiff.

* Stiff systems require implicit or specialized solvers.
* Neural ODEs that model such dynamics must use stiff ODE solvers, which are more expensive.
* Neural ODEs inherit all the numerical challenges of ODE solvers for PDEs:

  * Stability.
  * Step-size control.
  * Error vs cost.

#### Computational cost

* Forward pass: integrate ODE ``→`` potentially many solver steps.
* Backward pass: integrate adjoint ODE ``→`` roughly similar cost.
* Memory vs speed:

  * Standard backprop through solver stores all intermediate states.
  * Adjoint method can save memory by recomputing the forward trajectory.

#### Regularization and smoothness

Because dynamics are continuous in time:

* We often want ``f_\theta`` to be reasonably smooth in (h) and (t).
* Regularization ideas:

  * Penalize large ``|f_\theta(h,t)|``.
  * Penalize large ``\partial f_\theta / \partial h`` (Lipschitz regularization).
  * Encourage stability (e.g., negative real parts of eigenvalues in linearized dynamics).

* A “bad” neural ODE vector field can produce trajectories that blow up in finite time, analogous to an unstable scheme for PDEs.
"""

# ╔═╡ 6dc48bff-c958-426e-8c84-418cb87a45ca
md"""
### Relationship to other ML-for-PDE methods

#### vs. PINNs

* PINNs: directly enforce the PDE residual ``\mathcal{N}(u,x,t)=0`` at collocation points.
* Neural ODEs: treat time evolution as an ODE defined by a neural net (or hybrid physics+NN), fit to trajectory data.

Insights:

* PINNs are equation-driven (physics is explicit).
* Neural ODEs are trajectory-driven (learn dynamics from data).
* Hybrid models combine both: neural ODE with physics-based structure in ``f_\theta``.

#### vs. DeepONet / FNO (operator learning)

* DeepONet / FNO: learn mappings like ``u_0(x) \mapsto u(x,t_1), ..., u(x,t_k)`` directly as operators.
* Neural ODE: learn a dynamical system in time:

  * Represent how the state changes infinitesimally rather than mapping initial condition directly to the entire trajectory.

Insights:

* Neural ODE ≈ continuous-time ResNet.
* DeepONet / FNO ≈ global-in-time operator mapping input functions to output functions.

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
# ╟─6c6b08ad-4561-421f-949d-4c68d4b9e4a6
# ╟─5cd32a4e-c437-11f0-b519-193be7197e21
# ╟─23fc87b2-8e84-48cc-8f28-3b446320d1a1
# ╟─c84eda38-a798-4333-a2f8-6e547f4a8bf5
# ╟─c7dbd776-0473-4363-890b-eb548a306fdc
# ╟─d02de953-3b02-4ffd-96ef-92ddee7cd935
# ╟─ae404fc0-d0e0-45ce-8c05-7eeefdb20fa8
# ╟─b272d6a2-6352-4b50-a8f2-186ef8092790
# ╟─de08d098-067e-405e-a305-39bfa259c226
# ╟─6dc48bff-c958-426e-8c84-418cb87a45ca
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
