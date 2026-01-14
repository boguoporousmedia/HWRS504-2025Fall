### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ bb690cce-2489-43fb-a9f6-cc61b4701aa4
md"""

## Final Examination

### HWRS 504 – Numerical Methods for Environmental Transport Problems

- Date: December 16, 2025
- Time: 2 PM – 4 PM
- Location: Harshbarger 203
- Instructor: Bo Guo
- TA: Jianwen Du
"""

# ╔═╡ a0aa816c-ff0c-441d-b142-d7c603cd3fe3
md"""
### 1. (12 points)

Provide brief but precise explanations to the following questions.

(a) The method of characteristics (MOC) transforms certain PDEs into ODEs. Explain geometrically (i.e., by drawing graphs) and in words why MOC eliminates numerical dispersion for pure advection problems, e.g., ``
\frac{\partial{u}}{\partial{t}} + V \frac{\partial{u}}{\partial{x}} = 0.
``

(b) Explain the difference between cell-centered finite difference and finite volume methods in terms of how they enforce conservation and how they approximate spatial derivatives.

(c) Consider the diffusion equation ``\partial{u}/\partial{t} = D ∂^2 u / ∂ x^2``. Explain why forward Euler requires a severely restricted time step.

"""

# ╔═╡ df72dedd-dc90-4514-b998-65ef13d8b5a1
md"""
### 2. (28 points)

Consider the transient 1D advection–diffusion equation (``V>0``)

```math
\frac{\partial{u}}{\partial{t}} + V \frac{\partial{u}}{\partial{x}} = D \frac{\partial^2{u}}{\partial{x^2}}, \qquad 0 < x < L.
```

We use the following numerical scheme (first-order upwind for advection, central for diffusion):

```math
\frac{U_i^{n+1} - U_i^n}{\Delta t}
 + V \frac{U_i^n - U_{i-1}^n}{\Delta x}
  = D \frac{U_{i+1}^n - 2U_i^n + U_{i-1}^n}{\Delta x^2}.

```

(a) Derive the modified equation and identify the leading numerical diffusion term. Provide the expression for the numerical diffusion coefficient ``D_\text{num}``.

(b) Derive a stability condition using Fourier (von Neumann) stability analysis (you can assume an infinite or periodic domain for the von Neumann analysis). Express the condition in terms of the grid Peclet number and CFL number.

(c) Explain why increasing the grid Peclet number causes oscillations for central-difference advection schemes but not for upwind schemes.

"""

# ╔═╡ 74abad5c-dac8-4e6d-a62b-c431adeffe5b
md"""
### 3. (20 points)

Consider the pure advection equation (``V>0``)

```math
\frac{\partial{u}}{\partial{t}} + V \frac{\partial{u}}{\partial{x}} = 0.
```

For a uniform grid spacing, suppose we write the following three-point biased stencil:

```math
\frac{\partial{u}}{\partial{x}} \approx a U_{i-2} + b U_{i-1} + c U_i .
```

(a) Using Taylor expansion, determine conditions on ``a, b, c``, separately, for first-order and second-order accuracy.

(b) Using your result from part (a), explain why a linear three-point upwind stencil that is second-order accurate must include at least one coefficient that is negative.

(c) Give a qualitative explanation of why having a negative coefficient in such a stencil can lead to oscillations or wiggles in the numerical solution near sharp gradients. You can illustrate by constructing examples.

"""

# ╔═╡ ce70b986-d1ff-4604-8d25-d25f6755563a
md"""
### 4. (20 points) 

Consider the PINN formulation for the 1D steady-state diffusion equation:

```math
\frac{d^2 u}{dx^2} = f(x), \qquad u(0)=0, \quad u(1)=1.
```

The PINN loss is:

```math
\mathcal{L} =
\frac{1}{N_f}\sum_{j=1}^{N_f} \left( u''_\theta(x_j) - f(x_j) \right)^2

+ \left(u_\theta(0)\right)^2
+ \left(u_\theta(1)-1\right)^2 .
```

(a) Explain how automatic differentiation (AD) is used to obtain ``u''_\theta(x_j)`` and why AD avoids the discretization errors of finite differences.

(b) For a single-hidden-layer neural network
```math
u_\theta(x) = W_2 \sigma(W_1 x + b_1) + b_2,
```
where ``W₁`` is the weight matrix from the input layer to the hidden layer, write the backpropagation expression for ``\partial \mathcal{L} / \partial W_1.``

"""

# ╔═╡ 09ebf6e0-d629-11f0-b3b0-7940f0fce165
md"""
### 5. (20 points) 

We consider a family of PDEs with variable coefficients:

```math
\frac{d}{dx} \left( k(x) \frac{du}{dx} \right) = g(x), \qquad 0 < x < 1,
```
  with Dirichlet boundary conditions. Assume ``k(x)`` and ``g(x)`` vary across instances, and our goal is to learn the operator mapping:
```math
  (k, g) \mapsto u.
```

(a) Explain how a Neural Operator differs fundamentally from a standard neural network used for PINNs or supervised regression.

(b) For the Fourier Neural Operator (FNO), describe:

* how the Fourier transform is used in each layer,
* why global convolution in Fourier space can approximate PDE operators efficiently.

"""

# ╔═╡ 0c597f18-9917-4b5a-9ab7-8303e382ed41
md"""
### Appendix: 

Equation(s) that may be useful for answering the questions.

**Taylor series:**  
Let `` f(x) `` be `` C^n ``-continuous over the interval ``[a,b]``,  
`` f(x) \in C^n[a,b] ``. Then for any points ``x_0`` and ``x`` within ``[a,b]``, we have

```math
\begin{aligned}
f(x) =\; & f(x_0)
+ \left.\frac{df}{dx}\right|_{x_0} (x - x_0)
+ \left.\frac{d^2 f}{dx^2}\right|_{x_0} \frac{(x - x_0)^2}{2!}
+ \left.\frac{d^3 f}{dx^3}\right|_{x_0} \frac{(x - x_0)^3}{3!}
+ \cdots \\
&\qquad + \left.\frac{d^{n-1} f}{dx^{n-1}}\right|_{x_0}
  \frac{(x - x_0)^{n-1}}{(n-1)!}
  + R_n ,
\end{aligned}
\tag{5}
```

where the remainder term ``R_n`` is

```math
R_n \equiv 
\left.\frac{d^n f}{dx^n}\right|_{\xi}
\frac{(x - x_0)^n}{n!},
\qquad \xi \in [x_0, x].
\tag{6}
```

"""


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╟─bb690cce-2489-43fb-a9f6-cc61b4701aa4
# ╟─a0aa816c-ff0c-441d-b142-d7c603cd3fe3
# ╟─df72dedd-dc90-4514-b998-65ef13d8b5a1
# ╟─74abad5c-dac8-4e6d-a62b-c431adeffe5b
# ╟─ce70b986-d1ff-4604-8d25-d25f6755563a
# ╟─09ebf6e0-d629-11f0-b3b0-7940f0fce165
# ╟─0c597f18-9917-4b5a-9ab7-8303e382ed41
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
