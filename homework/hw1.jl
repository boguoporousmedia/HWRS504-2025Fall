### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ fd6d7bee-848b-11f0-1db7-c1791ffa18bf
md"""
# Homework 1
"""

# ╔═╡ a34e3683-fe10-4277-9ef3-9d0ce0c92717
md"""
- Assigned: Thursday, 28 August 2025

- Due: Thursday, 4 September 2025

- Instructor: Bo Guo, University of Arizona

- Semester: Fall 2025

- **NOTE**: Read Appendix A of CG and Appendix A of Moin (both uploaded to D2L) before working on the homework.
"""

# ╔═╡ 3b5eb3dc-b3d1-4bf9-9a83-1ba4e20b5f92
md"""
1. For the linear vector space `` P_n[\alpha, \beta] ``, derive a criterion by which the set of vectors  
``\{ p_i(x) \}_{i=1}^{n+1}``, with  

```math
p_i(x) = \alpha_{i0} + \alpha_{i1} x + \alpha_{i2} x^2 + \cdots + \alpha_{in} x^n,
```

can be tested for linear independence.

"""

# ╔═╡ 66ee76b7-26b4-48df-8226-82ed75c4a694
md"""

2. Lagrange polynomials are often used to interpolate data. Consider `` n+1 `` data points ``\{ U_i \}_{i=1}^{n+1} ``, measured at discrete locations `` \{ x_i\}_{i=1}^{n+1} ``. An `` n ``-th degree interpolation polynomial can then be defined by

```math
U(x) = \sum_{i=1}^{n+1} U_i \, \ell_{ni}(x)
```

Show that ``U(x)`` exactly matches the measurements at each node point. An alternative interpolation procedure is to use piecewise Lagrange polynomials. These are standard Lagrange polynomials that are of degree less than ``n`` and are defined over only “pieces” of the domain. For example, to “connect the data points by straight lines,” piecewise linear Lagrange polynomials are used as follows:

```math
U(x) = \sum_{i=1}^{n+1} U_i \, \phi_i(x)
```

where ``\phi_i(x)`` is a piecewise linear Lagrange polynomial defined as

```math
\phi_i(x) =
\begin{cases}
\dfrac{x - x_{i-1}}{x_i - x_{i-1}}, & x_{i-1} \leq x \leq x_i \\[1.2em]
\dfrac{x_{i+1} - x}{x_{i+1} - x_i}, & x_i \leq x \leq x_{i+1} \\[1.2em]
0, & \text{all other } x
\end{cases}
```

This interpolation retains the property that ``U(x_k) = U_k``, but differs from the previous one in that ``U(x)`` is only ``C^0[x_1, x_{n+1}]`` rather than ``C^\infty[x_1, x_{n+1}]``. Furthermore, the value of ``U(x)`` at any ``x_i \leq x \leq x_{i+1}`` is dependent only on the values ``U_i`` and ``U_{i+1}`` in the piecewise case, while in the previous case the dependence is on all of the ``U_j, \; j = 1, 2, \ldots, n+1``. For the data points below, compute ``U(x)`` using both of the interpolation methods above, and sketch the results.

| ``i`` | ``x_i`` | ``U_i`` |
|------|---------|---------|
| 1    | 0       | 1       |
| 2    | 2       | 1       |
| 3    | 3       | 2       |
| 4    | 4       | 0       |


"""

# ╔═╡ 4f1026e0-b2c0-4803-ae2b-5c84af5574de
md"""
3. Prove that the expression  

``\langle f(x), g(x) \rangle = \int_{\alpha}^{\beta} w(x) f(x) g(x) \, dx, \quad w(x) > 0``  

where ``w(x), f(x), g(x) \in C^0[\alpha, \beta]`` is a valid inner product for the space ``C^0[\alpha, \beta]``.

"""

# ╔═╡ e8d7347a-b1ed-405f-ac76-86bd7a6f98ae
md"""
4. Determine whether or not the differential operator  

```math
\mathcal{L} = \dfrac{d^2}{dx^2} + (x^2)\dfrac{d}{dx}
```  

is linear. State an appropriate domain and range for this operator.

"""

# ╔═╡ 0ea5b326-1bba-43dc-91df-6b19c4f0a067
md"""
5. Show that the expression  

```math
\|u(x)\| = \int_{\alpha}^{\beta} \left\{ [u(x)]^2 + \left[ \dfrac{du}{dx} \right]^2 \right\}^{1/2} dx
```  

is a valid norm for ``u(x) \in C^1[\alpha, \beta]``. Is this a valid norm for ``u(x) \in C^0[\alpha, \beta]``?

"""

# ╔═╡ 5d4d793f-4180-485d-97d4-8d698696f12d
md"""
6. (5 points) Watch this video on YouTube (https://youtu.be/ly4S0oi3Yz8), which provides a very nice introduction to partial differential equations using the Heat Equation as an example. The purpose of this video and the problem is to help you understand (or review in case you took partial differential equation classes before) what our primary target in this course (i.e., partial differential equations) is. 
    - Summarize the key points and insights. Use bullet points instead of paragraphs.

    - If you take the transient 1-D advection-diffusion equation presented in Lecture 1 and assume the fluid velocity V=0 and the source/sink term R=0, you obtain a so-called diffusion equation that describes the transport of a contaminant driven by only molecular diffusion. How is this diffusion equation related to the heat equation introduced in the video? Compare the variables and coefficients of the two equations directly. 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╟─fd6d7bee-848b-11f0-1db7-c1791ffa18bf
# ╟─a34e3683-fe10-4277-9ef3-9d0ce0c92717
# ╟─3b5eb3dc-b3d1-4bf9-9a83-1ba4e20b5f92
# ╟─66ee76b7-26b4-48df-8226-82ed75c4a694
# ╟─4f1026e0-b2c0-4803-ae2b-5c84af5574de
# ╟─e8d7347a-b1ed-405f-ac76-86bd7a6f98ae
# ╟─0ea5b326-1bba-43dc-91df-6b19c4f0a067
# ╟─5d4d793f-4180-485d-97d4-8d698696f12d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
