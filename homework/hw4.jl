### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 58ede0c2-c045-11f0-a4dd-a7fcc1ce281d
md"""
# Homework 4: Neural networks and PINNs
"""

# ╔═╡ f4edd850-b737-490c-8c4b-443437eefeaf
md"""
### 1. Feedforward Networks + AD + Backprop

Use a simple feedforward neural network to approximate the solution to the 1D Poisson equation

```math
-u''(x) = f(x), \qquad x\in (0,1),
```
with Dirichlet boundary conditions ``u(0)=0, u(1)=0``.

* Let ``f(x) = \pi^2 \sin(\pi x)``, so the exact solution is ``u(x) = \sin(\pi x)``.
* Train a small network (e.g., 2–3 hidden layers, tanh activation).
* Implement the loss using automatic differentiation to compute ``u''(x)``.
* Do NOT use existing packages (e.g., PyTorch, Jax, or others) for training, rather
  * Manually code the forward pass.
  * Manually compute gradients using your own implementation of backpropagation (matrix calculus).
  * Update weights using mini-batch stochastic gradient descent

Provide the following *deliverables*

* Hand-written derivation of backpropagation for one of the network layers.
* Python or Julia code implementing the network and training loop.
* Plot comparing the NN solution and exact solution.
* Discuss the following questions: 
  * Why does AD give you derivatives for free? 
  * How is AD different from backprop?

*Rationale*: The problem forces you to implement a minimal neural network solver from scratch, so that you *feel* what AD and backprop are doing under the hood, rather than using the existing packages as a black box.

"""

# ╔═╡ b3697ce9-5597-481e-8a60-414f704cf2b5
md"""
### 2. Physics-Informed Neural Networks (PINNs)

Use a PINN to solve the time-dependent heat equation

```math
u_t = u_{xx}, \qquad  x\in(0,1),\ t\in(0,1),
```
with initial and boundary conditions
```math
u(x,0)=\sin(\pi x),\quad u(0,t)=u(1,t)=0.
```

* Construct a single network ``u_\theta(x,t)``.
* Train using a PINN loss:
  * PDE residual (interior points)
  * Boundary loss
  * Initial condition loss
* Use automatic differentiation to compute both ``u_t`` and ``u_{xx}``.
* Compare to analytical solution:
  ```math
  u(x,t)=e^{-\pi^2 t}\sin(\pi x).
  ```

Provide the following *deliverables*

* Plot PDE residual over training.
* Plots of the predicted solution vs analytical solution at various times.
* Where does the PINN struggle (e.g., stiffness, boundary layer, collocation density)? why?

*Rationale*: This problem showcases the core ideas of PINNs: Embedding physics in the loss function, balancing residual terms, and understanding why time-dependent PDEs require careful sampling.
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
# ╟─58ede0c2-c045-11f0-a4dd-a7fcc1ce281d
# ╟─f4edd850-b737-490c-8c4b-443437eefeaf
# ╟─b3697ce9-5597-481e-8a60-414f704cf2b5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
