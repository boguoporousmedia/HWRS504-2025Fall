### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ ca58c9b4-b993-4740-9882-fbbcd2ab50f8
using LaTeXStrings

# ╔═╡ 7c59dcf2-8e04-11f0-2712-49a020c07790
md"""
## Homework 2  
- Assigned: Thursday, 9 September 2025  
- Due: Tuesday, 23 September 2025  
- Instructor: Bo Guo, University of Arizona  
- Semester: Fall 2025
"""

# ╔═╡ e910408d-b9c6-45f8-ab32-1ffcc0dc5302
md"""
### Problem 1 (10 points) – Quadratic Spline Interpolation

Formulate an algorithm in words/equations (not in computer codes) for quadratic spline interpolation by following the steps below:

(a) What properties of the polynomial need to be continuous at the grid points?  

(b) Write the equations for the quadratic splines. Your equations and the variables therein must be functions of only the values at the grid points (e.g., the values at the rest of the domain are unknown; the values of the derivatives at the grid points are also unknown). Are boundary conditions required to solve the problem? If so, what boundary condition can you use? What is your argument?  

(c) What is the cost of this algorithm compared to cubic spline interpolation?  
(Hint: Solving a bidiagonal system is just as costly as solving a tridiagonal system; both scale linearly with the number of unknowns, i.e., the computational complexity ∼ O(N).)
"""

# ╔═╡ 32e95fc4-e4d6-4831-b09e-6dfcc7ff3ab6
md"""
### Problem 2 (15 points) – Volterra Integral Equation

Integral equations are very important in many areas of engineering and science. An inhomogeneous Volterra integral equation of the second kind can be (semi-)generically written as:

```math
f(x) + \int_0^x K(x,t) f(t)\,dt = g(x)
```

where the kernel ``K(x,t)`` and right-hand-side ``g(x)`` are given functions.

(a) Develop an algorithm to solve a generic Volterra integral equation based on the trapezoid method. What is the order of accuracy of this method (i.e., the order of the approximation for ``f(x)``)?

(b) Write a computer code (in Python or another language of your choice) to implement this algorithm to solve the integral equation with:

```math
K(x,t) = (x - t)^2, \quad g(x) = e^{-x^2}\cos(2\pi x), \quad x \in [0,1].
```

Demonstrate that your numerical implementation matches your predicted order of accuracy.

*Hint:*

1. Choose a very small grid size, compute the solution, and use that solution as the “true” solution.
2. Use L2-norm to measure the magnitude of the error between the numerical solutions and the “true” solution.
3. Plot the magnitude of the error vs. the grid size in log-log scale and compute the slope.

"""


# ╔═╡ 87582a91-3b66-484a-81ba-e4674011d907

md"""
### Problem 3 (10 points) – Finite Difference with Unequal Spacing

Consider three points denoted by ``x_{i-1}, x_i, x_{i+1}``, where the spacing is unequal, such that:

```math
\Delta x_i = x_{i+1} - x_i, \quad \nabla x_i = x_i - x_{i-1}, \quad \Delta x_i \neq \nabla x_i.
```

Consider a finite difference approximation for:

```math
\left.\frac{du}{dx}\right|_{x_i}
```

using these three points.

- Show that an infinite number of consistent approximations exists.  
- Choose the one that gives **second-order accuracy** for the derivative.  
- Can third-order accuracy ever be achieved using three points (with or without equal spacing)?
"""

# ╔═╡ 39a25c0b-74cd-46d8-9b62-6ed95a279585
md"""
### Problem 4 (15 points) - Error cancelation

For the variably upstream-weighted scheme with uniform grid spacing for a first-order derivative
```math
\left. \frac{du}{dx}\right|_{x_i}
= \alpha \frac{u_i - u_{i-1}}{\Delta x}
+ (1-\alpha)\frac{u_{i+1} - u_{i-1}}{2\Delta x}
```

where ``\alpha`` is an arbitrary parameter. Show that you can choose a value of ``\alpha`` so that the FDA to the steady-state advection–diffusion equation gives exact results.

**Hints:**

- Taylor expansion
- ``\left.\frac{d^m u}{dx^m}\right|_{x_i} = \left(\frac{V}{D}\right)^{m-1} \left.\frac{du}{dx}\right|_{x_i} ``

Take derivative for both sides:

```math
\frac{du}{dx} = \frac{V}{D} \frac{d^2 u}{dx^2}
\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\Rightarrow
\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\frac{d^2 u}{dx^2} = \frac{V}{D} \frac{d^3 u}{dx^3}
```

Keep taking the derivatives and you will quickly see that ``\left.\frac{d^m u}{dx^m}\right|_{x_i} = \left(\frac{V}{D}\right)^{m-1} \left.\frac{du}{dx}\right|_{x_i} ``.

- Apply Taylor expansion to the FDA and solve for ``\alpha`` by taking advantage of the above identity.


"""

# ╔═╡ 749c0c23-0dae-4462-b370-9a9a092b79e1
md"""
### Problem 5 (25 points) – Advection-Diffusion Transport Equation

The one-dimensional, steady-state advection–diffusion transport equation may be written as:

```math
V \frac{du}{dx} - D \frac{d^2 u}{dx^2} = 0, \quad 0 < x < L
```

For boundary conditions ``u(0) = 1`` and ``u(L) = 0``, the analytical solution is:

```math
u(x) = \frac{e^{Vx/D} - e^{VL/D}}{1 - e^{VL/D}}
```

(a) Write a centered-in-space, second-order finite difference approximation for the equation above, assuming constant grid spacing. Then, write a computer code (in Julia or Python or other programming languages that you like) to solve the resulting equations to obtain the finite difference approximation, using an arbitrary number of grid points (i.e., this should be a parameter that allows the user to set as input; you can assume constant spacing).

Perform calculations for ``L = 1, V = 1, D = 0.01``, and different values of ``\Delta x = 0.01, 0.02, 0.1, 0.2, 0.5``.
Comment on the behavior of the finite difference solution as a function of the **grid Peclet number**:

```math
Pe_G \equiv \frac{V \Delta x}{D}.
```

(b) Repeat part (a), but instead use a one-sided approximation for the advection term:

```math
\left.\frac{du}{dx}\right|_{x_i} \approx \frac{U_i - U_{i-1}}{x_i - x_{i-1}}
```

This approximation is often referred to as an **upstream-weighted** approximation for the advection term.

(c) Finally, write the approximation for the first derivative using a **variably upstream-weighted** approximation with weighting parameter $\alpha$:

```math
\left.\frac{du}{dx}\right|_{x_i} \approx \alpha \frac{U_i - U_{i-1}}{x_i - x_{i-1}} + (1-\alpha)\frac{U_{i+1} - U_{i-1}}{x_{i+1} - x_{i-1}}
```

Choose $\alpha = 0.5$ and repeat part (a).
Compare the solutions with the ones from parts (a) and (b), and comment on the behavior of the variably upstream-weighted approximation.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"

[compat]
LaTeXStrings = "~1.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "909cab9b0d3c33fad074c918f1e1b063e921183a"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"
"""

# ╔═╡ Cell order:
# ╟─ca58c9b4-b993-4740-9882-fbbcd2ab50f8
# ╟─7c59dcf2-8e04-11f0-2712-49a020c07790
# ╟─e910408d-b9c6-45f8-ab32-1ffcc0dc5302
# ╟─32e95fc4-e4d6-4831-b09e-6dfcc7ff3ab6
# ╟─87582a91-3b66-484a-81ba-e4674011d907
# ╟─39a25c0b-74cd-46d8-9b62-6ed95a279585
# ╟─749c0c23-0dae-4462-b370-9a9a092b79e1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
