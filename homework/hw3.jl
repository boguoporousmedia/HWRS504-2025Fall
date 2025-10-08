### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ ea254018-4624-49c5-9761-33215b35f6b5
md"""
### Homework 3

- **Assigned:** Tuesday, 7 October 2025  

- **Due:** Friday, 17 October 2025  

- **Instructor:** Bo Guo, University of Arizona  

- **Semester:** Fall 2025  

---
"""

# ╔═╡ 77601700-a3cb-11f0-22ea-0152abc5511e
md"""
#### **1. (15 points) Newton–Raphson method**

Use the Newton–Raphson method to find an approximation of  
```math
\sqrt{10}
```

with an error

```math
\epsilon < 10^{-8}.
```

You can choose your initial guess as 3.
Show (without using the square root function of a calculator) that your answer is indeed within ``10^{-8}`` of the true solution.

"""

# ╔═╡ 4f7967b5-f78e-469c-b092-96901bcfc170
md"""
---

#### **2. (55 points) One-dimensional transient advection–diffusion–reaction Equation**

Consider the PDE:

```math
\frac{\partial u}{\partial t}
+ V \frac{\partial u}{\partial x}
- D \frac{\partial^2 u}{\partial x^2}
+ K u = 0, \quad 0 < x < L, \; t > 0
```

Boundary and initial conditions:

```math
u(0,t) = 1, \quad \frac{\partial u}{\partial x}(L,t) = 0, \quad u(x,0) = 0.
```

**Tasks:**

1. Write a finite difference approximation using:

   * variable upstream weighting (with weight ``α``)
   * variably weighted Euler time-stepping (with weight ``θ``)
2. Write a code to solve the resulting system of algebraic equations.
3. The number of grid points ``N_x`` should be a user input.

**Use parameters:**

```math
L = 1, \quad V = 1, \quad D = 0.01, \quad K = 0.0001.
```

**Investigate:**

* Grid Peclet numbers:

  ```math
  Pe_G = \frac{V \, \Delta x}{D} \in \{0.1, 1, 2, 5, 10\}
  ```
* Different time steps ``Δt``:

  * For conditionally stable cases, find the threshold ``Δt`` above which the numerical solution becomes unstable.
* Upstream weighting:

  ```math
  α = 0, \; 0.5, \; 1
  ```
* Time weighting:

  ```math
  θ = 0, \; 0.5, \; 1
  ```

**Verification:**
Compare your numerical results with the analytical solution.

---
**Hints**

1. Use the semi-infinite analytical solution for comparison (``0 < x < ∞``).
   To ensure equivalence, choose a large enough domain ``[0, L]`` or small simulation time ``t`` such that the Neumann boundary at ``x = L`` is effectively zero flux.

2. Analytical solution (for semi-infinite domain):

```math
u(x, t) =
\frac{1}{2}
e^{\frac{Vx}{2D}}

```

```math
\left[
  e^{\frac{Vx}{2D} + x \sqrt{\frac{K}{D}}}
  \operatorname{erfc}\!
\left(
    \frac{x}{2\sqrt{Dt}} + \sqrt{\frac{V^2 t}{4D} + Kt}
  \right)
  + e^{-\frac{Vx}{2D} - x \sqrt{\frac{K}{D}}}
  \operatorname{erfc}\!\left(
    \frac{x}{2\sqrt{Dt}} - \sqrt{\frac{V^2 t}{4D} + Kt}
  \right)
\right]
```

where

```math
\operatorname{erfc}(x) =
\frac{2}{\sqrt{\pi}} \int_x^{+\infty} e^{-u^2} \, du
```

is the complementary error function (available in Python and MATLAB).

"""

# ╔═╡ cb289f54-6144-4b2a-b706-740c788a3646
md"""
#### **3. (10 points) Picard vs. Newton-Raphson**

Read the paper:
**Celia, M.A., Bouloutas, E.T., and Zarba, R.L. (1990)**
*A general mass-conservative numerical solution for the unsaturated flow equation.*
*Water Resources Research*, 26(7), 1483–1496.

Answer in your own words:

(a) What is Picard iteration? How is it different from the Newton–Raphson iteration? What are the advantages and disadvantages?

(b) Why is mass not conserved numerically when schemes other than the one by Celia et al. (1990) are used?

"""

# ╔═╡ ff977d2e-d27d-4470-8cc9-1b4e71869b27
md"""
#### **4. (10 points) Term project proposal**

Propose a term project problem that is:

* Challenging (several weeks of effort)
* Related to numerical modeling but not identical to your existing research

The final deliverables:

* Written report (5–10 pages) due in the last week of the semester
* Oral presentation (date TBD)

Examples of topics:

* Contaminant transport with nonlinear reactions (e.g., biodegradation)
* Modeling flow in unsaturated soils
* Employ deep learning methods to solve PDEs in your field (e.g., Richards equation)

For this homework, provide:

* Your proposed topic
* A short description of your intended approach or method

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╟─ea254018-4624-49c5-9761-33215b35f6b5
# ╟─77601700-a3cb-11f0-22ea-0152abc5511e
# ╟─4f7967b5-f78e-469c-b092-96901bcfc170
# ╟─cb289f54-6144-4b2a-b706-740c788a3646
# ╟─ff977d2e-d27d-4470-8cc9-1b4e71869b27
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
