### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 858117aa-c048-11f0-9ed0-91cbbc7e6472
md"""
# Homework 5: Operator learning
"""

# ╔═╡ 83455544-e401-4a91-ae22-9e7ee972be3c
md"""
### Operator Learning (DeepONet + FNO)

Learn the solution operator of a parametric elliptic PDE
```math
-u''(x) = f(x), \qquad x\in(0,1),\quad u(0)=u(1)=0,
```
where ``f(x)`` varies across training samples.

1. Generate training data

   * Sample 100–1000 random forcing functions ``f(x)``:

     * e.g., random combinations of sine modes, or Gaussian random fields.
   * Solve the PDE for each sample using a finite-difference solver.

2. Train a DeepONet

   * Branch net takes the discretized ``f`` (e.g., 100-point vector).
   * Trunk net takes coordinate ``x``.
   * Output ``u(x)`` ≈ solution.

3. Train a Fourier Neural Operator (FNO)

   * Treat the PDE as mapping ``f(x)\to u(x)`` through spectral convolution layers.

4. Compare performance

   * Generalization to unseen ``f(x)``
   * Error vs training set size
   * Training time

Provide the following *deliverables*

* Code for DeepONet and FNO implementations.
* Error plots on training and test sets.
* Comparison and discussion:

  * Why does FNO implicitly capture translation invariance?
  * Why does DeepONet generalize well to functions not seen during training?

*Rationale*: This is the cleanest, simplest testbed for operator learning. The PDE is linear and the operator is smooth. You can see clearly how operator learning differs from PINNs.

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
# ╟─858117aa-c048-11f0-9ed0-91cbbc7e6472
# ╟─83455544-e401-4a91-ae22-9e7ee972be3c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
