### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ acd96492-cf2c-11f0-af29-73984b6f153e
md"""
# Module 17: Convolutional Neural Network and U-Net
"""

# ╔═╡ 022bad0c-9b52-4201-8ee9-f1a465ac6621
md"""
Introduction to CNN:

https://stanford.edu/~shervine/teaching/cs-230/cheatsheet-convolutional-neural-networks
"""

# ╔═╡ 170f50ca-f002-4bf4-a9a2-8b2ef9bacc74
md"""
### How to use CNN to solve PDEs?

Suppose:

* Example: 2D diffusion,
```math
u_t = \nabla \cdot (D \nabla u).
```
* Discretization: regular grid, finite-difference or finite-volume.
* Goal: Given input fields (e.g., ``D``, boundary conditions), predict the solution field ``\mathbf{u}(\mathbf{x},\mathbf{y},t)`` with a CNN.

**Idea**: Train the CNN to learn an update rule

```math
u^{n+1} = \mathcal{N}_\theta\big(u^n, \text{coefficients}; \text{BCs}\big)
```

Step 1 – Generate Training Data

For each trajectory:

1. Choose initial condition ``\mathbf{u}^0`` and parameters (diffusivity ``D``, etc.).
2. Use a classical time-stepping method to generate:
   ```math
   \mathbf{u}^0, \mathbf{u}^1, \mathbf{u}^2, \dots, \mathbf{u}^T.
   ```
3. Form training pairs:
   ```math
   X_n = (\mathbf{u}^n, D, \text{BCs}),\quad Y_n = \mathbf{u}^{n+1}.
   ```

Step 2 – CNN Architecture

Input channels:

* current state ``\mathbf{u}^n``
* coefficients / parameters
* boundary masks

Output:

* predicted next state ``\hat{\mathbf{u}}^{n+1} = \mathcal{N}_\theta(X_n)``.

Step 3 – Loss and Training

Loss:

```math
\mathcal{L}(\theta) = \frac{1}{N}\sum_n |\mathcal{N}_\theta(X_n) - Y_n|_2^2.
```

Train with standard gradient-based optimization (as before).

Step 4 – Using the CNN as a Solver

Algorithm (CNN Time-Stepper):

1. Initialize state ``\mathbf{u}^0`` from initial condition.
2. For ``n = 0,1,\dots,T-1``:

   1. Form input: ``X_n = (\mathbf{u}^n, D, \text{BCs})``.
   2. Predict: ``\hat{\mathbf{u}}^{n+1} = \mathcal{N}_{\theta^\star}(X_n)``.
   3. Set ``u^{n+1} \leftarrow \hat{\mathbf{u}}^{n+1}``.
3. Return ``\mathbf{u}^T`` as approximate solution.

This is conceptually similar to:

* explicit time-stepping, or
* a learned smoother/iteration in a nonlinear solver.

"""

# ╔═╡ 038caec1-4937-42df-8deb-3653dbda77aa
md"""
## CNN Architecture and convolution

#### What does the convolution layer do

A single 2D convolution layer takes an input tensor
```math
X \in \mathbb{R}^{C_{\text{in}} \times N_x \times N_y}
```
and produces an output tensor
```math
Y \in \mathbb{R}^{C_{\text{out}} \times N_x' \times N_y'}.
```

Each output channel ``c_{\text{out}}`` is computed as a weighted sum of local patches across all input channels.

For simplicity, take:

* kernel size = ``3 \times 3``
* stride = 1
* “same” padding (output size = input size)

Then for output channel ``m`` and spatial location ``(i,j)``:

```math
Y_{m,i,j} = b_m + \sum_{c=1}^{C_{\text{in}}} \sum_{p=-1}^{1} \sum_{q=-1}^{1}
W_{m,c,p,q}, X_{c,,i+p,,j+q}
```

Where:

* ``W_{m,c,p,q}`` are the learnable kernel weights
* ``b_m`` is a learnable bias
* ``(p,q)`` loop over the ``3\times 3`` neighborhood centered at ``(i,j)``

Interpretation:
For each output pixel, you take the 9 neighbors (in each input channel), multiply each by a learned coefficient, sum them up, add a bias.

This is exactly:

* a linear stencil over space,
* but with multiple input/output channels (like multiple PDE fields interacting).


#### Connection to PDE stencils

Consider a standard 5-point Laplacian stencil on a uniform grid:

```math
(\Delta u)*{i,j} \approx \frac{1}{h^2} \Big(-4u*{i,j} + u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} \Big)
```

You can see this as a fixed convolution with kernel:

```math
K_{\Delta} = \frac{1}{h^2}
\begin{bmatrix}
0 & 1 & 0 \newline
1 & -4 & 1 \newline
0 & 1 & 0
\end{bmatrix}
```

A convolution layer is exactly the same operation, except:

* The kernel is learned, not prescribed.
* You can have many kernels (one per output channel).
* You can act on multiple input fields simultaneously.

So one 2D conv layer with:

* ``C_{\text{in}}=1``, ``C_{\text{out}}=1``, kernel size ``3\times3``
  is mathematically identical to applying a 2D linear stencil like the Laplacian — just with arbitrary coefficients instead of the specific Laplace coefficients.

Intuition: A conv layer is a bunch of learned PDE-like stencils applied everywhere on the grid, sharing the same weights.

#### What “multiple channels” really mean

For example, for the input:

* Channel 1: ``k(x,y)`` (permeability)
* Channel 2: ``f(x,y)`` (forcing)
* Channel 3: boundary mask
* Channel 4: boundary values

``C_{\text{in}}=4``.

Now suppose we use:

* ``C_{\text{out}}=16``
* kernel size ``3\times3``

Then the kernel tensor has shape:
```math
W \in \mathbb{R}^{16 \times 4 \times 3 \times 3}.
```

For each output channel ``m``:

```math
Y_{m,i,j}
= b_m + \sum_{c=1}^{4} \sum_{p=-1}^{1} \sum_{q=-1}^{1}
  W_{m,c,p,q}, X_{c,i+p,j+q}.
```

So each of the 16 features at ``(i,j)`` is a linear combination of 3×3 patches from all 4 input fields.

We can think of it like this:

* One output channel might learn something like a local estimate of ``\nabla k`` (heterogeneity pattern).
* Another might learn a local approximation of ``\Delta f``.
* Another might encode something boundary-related.
* Deep layers combine these to form rich representations of the physics.

#### Adding nonlinearities: Stack convolution + activation

Not only applying one convolution, we apply conv ``→`` nonlinearity ``→`` conv ``→`` ...

A single layer is linear:
```math
Y = W * X + b
```

After applying, say, ReLU:

```math
Z = \sigma(Y) = \text{ReLU}(Y)
```

Stacked layers approximate highly nonlinear operators:

```math
u \approx \mathcal{N}_\theta(k,f,\text{BCs})
```

Even if each convolution looks like a linear stencil, the composition plus nonlinearities lets the network learn nonlinear PDE behavior, e.g.:

* Richards equation
* multiphase flow
* nonlinear diffusion
* advection-dominated flow (combined with other tricks)

#### Receptive field: how local stencils become nonlocal

Question: If the kernel only looks at ``3×3`` neighborhoods, how does the network know about far-away boundary conditions or global structure?

Answer: stacking conv layers grows the receptive field.

* 1 layer, ``3×3`` kernel ``→`` each output sees ``3×3`` patch of input
* 2 layers ``→`` each output sees 5×5 patch
* 3 layers ``→`` 7×7 patch
* in general, with ``L`` layers and kernel size ``k``, receptive field size ≈ ``1 + (k-1)L``.

So a deep CNN can approximate operators with long-range coupling, even though each layer is local.

This is completely analogous to:

* composing multiple local PDE updates
* or multigrid / diffusion over many steps.

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
# ╟─acd96492-cf2c-11f0-af29-73984b6f153e
# ╟─022bad0c-9b52-4201-8ee9-f1a465ac6621
# ╟─170f50ca-f002-4bf4-a9a2-8b2ef9bacc74
# ╟─038caec1-4937-42df-8deb-3653dbda77aa
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
