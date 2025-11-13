### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ fdc65eac-bace-11f0-9d77-6956731f466b
md"""
# Module 12: Introduction to Operator Learning for PDEs
"""

# ╔═╡ 08eab8a7-faed-4645-b43c-60a5ce5ee1b0
md"""
### What is Operator Learning?

In many applications, we solve the *same PDE form* repeatedly with different inputs. Examples:

* Changing forcing terms ``f(x)``
* Changing permeability fields ``k(x)``
* Changing boundary conditions
* Running uncertainty quantification with many realizations

Each solution requires a full PDE solve, often computationally heavy.

**Key Idea:**

The PDE solution defines a map:
```math
\mathcal{G}: \text{input function(s)} \rightarrow \text{solution function}
```

This map is called the solution operator. If we could learn a good approximation of this operator, then new PDE solutions become *instant lookup* rather than solving from scratch.

"""

# ╔═╡ 400657d0-f9d8-496c-b24e-2fd10d3c7d34
md"""
### What Is an Operator?

A regular neural network learns a function between finite-dimensional vectors:
```math
\mathbb{R}^n \rightarrow \mathbb{R}^m.
```

An operator takes a function as input and produces another function:
```math
\mathcal{G}: f(x) \mapsto u(x).
```

This is a map between infinite-dimensional spaces (in theory).

**Example**

Consider the 1D Poisson equation:
```math

u''(x) = f(x),\quad x \in (0,1), \quad u(0) = u(1) = 0.
```

For each input function ``f(x)``, the PDE returns a solution function ``u(x)``.

So the operator is:
```math
\mathcal{G}[f] = u.
```

**Once the operator is learned**:

* Fast inference
* Useful for inverse problems (find parameters from data)
* Useful for uncertainty quantification (solve many realizations)
* Can generalize across different spatial grids if structured correctly

"""

# ╔═╡ 2085ad02-80f3-4ac0-b39e-6ce522ac8657
md"""
### Two operator learning examples

**DeepONet**

Lu, L., Jin, P., Pang, G., Zhang, Z. and Karniadakis, G.E., 2021. Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators. Nature machine intelligence, 3(3), pp.218-229.

Idea: approximate the solution as a learned finite expansion:
```math
u(x) \approx \sum_{i=1}^m c_i \phi_i(x).
```

* Branch network learns coefficients ``c_i`` from the input function ``f``, which tell how strongly each basis function contributes to the solution.
* Trunk network learns the shape of the basis functions ``\phi_i(x)`` depending on the location ``x``.

Branch Network

* Input: sampled values of ``f``. We cannot feed a full function to a neural net, so we sample the input at a given resolution. **NOTE**: Because of this, DeepONet is resolution-dependent and cannot generalize to resolutions different from that used in the training. 
* Output: ``c = (c_1, c_2, ..., c_m) ``
* Architecture: typically fully-connected MLP (3–6 layers)

Trunk Network

* Input: spatial coordinate ``x``
* Output: ``m`` basis values at that location, i.e., basis vector ``\phi(x) = (\phi_1(x), ..., \phi_m(x))``
* Architecture: another MLP (often small)

Prediction
```math
u(x) = c^\top \phi(x).
```

Remarks:
- Many PDE solution spaces are low-dimensional even though the input function space is infinite-dimensional. DeepONet identifies the essential modes of variation.
- DeepONet *learns* the optimal (low-dimensional) basis set, not using predetermined bases such as polynomials or Fourier modes. It also learns how input functions activate those basis modes.
- DeepONet cannot generalize across different mesh resolutions.

"""

# ╔═╡ b3e14b4e-04f6-413f-9b47-8f2072495036
md"""
**Fourier Neural Operator (FNO)**

Li, Z., Kovachki, N., Azizzadenesheli, K., Liu, B., Bhattacharya, K., Stuart, A. and Anandkumar, A., 2020. Fourier neural operator for parametric partial differential equations. arXiv preprint arXiv:2010.08895.

Idea: Represent the operator using *learned convolutions in Fourier space*.

*Why Fourier?*

- Solution to many PDEs is an (possibly nonlinear) integral operator with a kernel ``K`` (e.g., Green's function). The kernel can be translation-invariant for many conditions, i.e., shifting the input shifts the output in exactly the same way. 
- Convolution in the spatial domain becomes multiplication in the Fourier domain: ``\widehat{K * f} = \hat{K}\cdot\hat{f}``. So if we learn ``\hat{K}``, we learn the operator.
- Nonlinearity and a local correction can be added separately.

Assumption: The operator is reasonably **smooth** in frequency space, i.e., It does not wildly amplify high-frequency components of the input. Instead, it tends to preserve or damp higher frequencies. This is true for most elliptic and diffusive systems, but not true for shock-forming hyperbolic PDEs.

Given an input field ``f(x)``, define a sequence of fields:
```math
v_{0}(x) = f(x), \quad v_{k+1} = \sigma\left(\text{IFFT}(R(k) \hat{v}(k)) + W * v_k + b\right).
```

One FNO layer does three steps:

1. *Fourier Transform*
   ```math
   \hat{v}(k) = \text{FFT}(v(x))
   ```

2. *Learn spectral mixing*

   Keep only low-frequency modes (say the first ``M``), and apply learned complex-valued weights:
   ```math
   \hat{v}_{\text{new}}(k) = R(k) \hat{v}(k), \quad \text{for } |k| \le M
   ```

3. *Inverse Fourier Transform*
   ```math
   v_{\text{new}}(x) = \text{IFFT}(\hat{v}_{\text{new}})
   ```

Then add nonlinearity and local correction:
```math
v_{k+1}(x) = \sigma \left( v_{\text{new}}(x) + W v(x) + b \right)
```

The neural network learns ``R(k)``, ``W``, and ``b``. Learning ``R(k)`` means learning how the operator scales each frequency component of the input — which fully determines the operator when it is translation-invariant.

Remarks:
- FNO learns PDE solution operators by learning how input fields interact in Fourier space. FNO handles global dependencies efficiently because Fourier modes span the entire domain, allowing the network to update the whole solution at once.
- Fourier modes correspond to physical frequencies, not grid indices. If you evaluate on a finer grid, you still evaluate the same low-frequency modes with the same ``R(k)``. Hence, after training on ``N``, you can infer on ``N' \gg N`` without retraining, as long as the target frequencies are within the retained band. Thus, the operator generalizes across different mesh resolutions, e.g., train at (``64 \times 64``), evaluate at (``256 \times 256``).
- The standard FNO assumes periodic or easily padded domains and struggles with sharp discontinuities (shocks, fractures). Modifications are required to handle these.
"""

# ╔═╡ 57d1c964-5802-41a0-bb43-f273af3d8c12
md"""
### Comparison

| Feature                  | Traditional Solver   | DeepONet            | FNO                                |
| ------------------------ | -------------------- | ------------------- | ---------------------------------- |
| Computes                 | One solution per run | Operator            | Operator                           |
| Cost per new solve       | High                 | Low                 | Low                                |
| Works across resolutions | Yes                  | No                  | Yes                                |
| Handles global couplings | Solver does          | Learns them         | Learns via Fourier modes |
| Data requirement         | None                 | Needs training data | Needs training data                |


In general, operator learning approaches require the following to be successful
1. The solution operator is *smooth* in the input function space (e.g., it would struggle with sharp jumps/discontinuities.)
2. The solution manifold is *low-dimensional* relative to full function space.
3. Representative training data covers the range of possible inputs.

"""

# ╔═╡ 160a1e3a-6fd7-4030-b40a-7e96f3bd2bf4
md"""

**Example:**
Consider the Poisson equation

```math
u''(x) = f(x),\quad u(0)=u(1)=0.
```

* Generate random smooth functions ``f(x)`` (e.g., sums of sine waves)
* Solve PDE with finite differences to get ``u(x)``
* Collect dataset of pairs ``(f, u)``

Then ask, does there seem to be a smooth relationship (i.e., mapping) between shapes of ``f`` and shapes of ``u``? DeepONet and FNO provide two different approaches to learn the mapping. 
"""

# ╔═╡ f1483064-efad-41a7-a3ac-34a005b2f98d
md"""

### Summary

* Operator learning approximates solution **mappings**, not individual solutions.
* DeepONet learns separable basis representations, while FNO learns convolution kernels in Fourier space.
* In general, their success relies on structural simplicity in solution manifolds. Any complexities will need to be handled by modifying the standard approaches.

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
# ╟─fdc65eac-bace-11f0-9d77-6956731f466b
# ╟─08eab8a7-faed-4645-b43c-60a5ce5ee1b0
# ╟─400657d0-f9d8-496c-b24e-2fd10d3c7d34
# ╟─2085ad02-80f3-4ac0-b39e-6ce522ac8657
# ╟─b3e14b4e-04f6-413f-9b47-8f2072495036
# ╟─57d1c964-5802-41a0-bb43-f273af3d8c12
# ╟─160a1e3a-6fd7-4030-b40a-7e96f3bd2bf4
# ╟─f1483064-efad-41a7-a3ac-34a005b2f98d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
