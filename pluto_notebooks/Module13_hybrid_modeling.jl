### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ fec68b18-c0ac-11f0-8279-b533b0b6f33f
md"""
# Module 13: Hybrid modeling

Lecture notes adapted and expanded from ETH Zuerich *AI in the Sciences and Engineering (2024)* [Class notes](https://camlab.ethz.ch/teaching/ai-in-the-sciences-and-engineering-2024.html)
"""

# ╔═╡ 22077fa2-657d-470b-bea1-9e13179fd526
md"""
Hybrid modeling blends **physics-based models** with **machine learning components**.  

#### Why Hybrid Modeling?

**(a) Pure physics (PDE solvers)**
- Pros: interpretable, stable, physically grounded  
- Cons: slow for complex 3D problems, may require assumptions (e.g., constitutive laws), may not have sufficient information to parameterize.

**(b) Pure neural networks**
- Pros: very fast once trained, easy to run on GPUs  
- Cons: needs lots of data, struggles with extrapolation, hard to interpret  

**Limitations motivating hybrid methods**

From PINNs:
- Expensive to optimize  
- Hard to scale to multi-scale or high-frequency solutions 

From operator learning:
- Requires many training samples  
- Struggles outside training distribution
- Needs assumptions on encoded fields (Fourier, CNN, etc.)

Hybrid modeling uses ML **only where it helps**, not everywhere, for example:
- When physics is known but incomplete  
- When solvers are accurate but slow  
- When constitutive relations are uncertain

"""

# ╔═╡ f42be050-2ebf-4513-995a-1dfed6ff6f86
md"""
#### What Is a Hybrid Workflow?

Instead of replacing existing workflows with deep neural networks (DNNs) entirely, we insert DNNs *inside* the workflow of classical numerical methods algorithm to
1) Accelerate the workflow, or
2) Learn the parts we are unsure of / have incomplete knowledge

Examples:

- “residual modeling”

```math
y(x) \approx y_\text{physics}(x) + \text{NN}(x;\theta)
```

Or more generally:

```math
\text{HybridSolver}(u^n ; \theta) = \text{PhysicsStep}(u^n) \;\;+\;\; \text{DNNCorrection}(u^n ; \theta).
```

- The physics gives a *good prior*  
- The DNN only learns what physics misses  
- Smaller learning task → less data → better generalization  

Hybrid modeling is a **continuum**:
- Physics-only → PINNs → residual models → solver-in-the-loop → fully learned

"""

# ╔═╡ 05dc4602-db65-4565-aaa2-aa4ac6bb5320
md"""
#### Residual Modeling

Deep neural network learns residual correction to physics:

```math 
r(x) = y_\text{true}(x) - y_\text{physics}(x)
```

Use the DNN to approximate the residual correction

```math
NN(x;\theta) = \hat{r}(x;\theta) \approx r(x)
```

That is, minimize the loss function
```math
L(\theta) = \Sigma_i^N \left( NN(x;\theta) - r(x) \right)
```


Then prediction is:

```math
\hat{y}(x) = y_\text{physics}(x) + \hat{r}(x;\theta)
```

Advantages
- Easy to train  
- Good generalization  
- Physics still dictates the main structure  

Potential applications
- Constitutive laws (permeability–porosity, hydraulic conductivity curves)  
- Closure relations  
- Subgrid models  
- Soil/rock heterogeneity corrections  

"""

# ╔═╡ 562f69fc-ac0f-4e46-a908-9982b96da5ea
md"""
#### DNNs Inside PDE Solvers

Instead of post-processing the output, we can insert DNNs *inside* each step of a solver.

Example: Navier–Stokes fractional-step solver. 

Um et al, Solver-in-the-loop: Learning from differentiable physics to interact with iterative PDE-solvers, NeurIPS (2020)

A hybrid approach inserts:

```math
u^{n+1} = u^{n+1}_\text{low-res} + \text{NN}(u^{n+1}_\text{low-res}, p^{n+1}_\text{low-res}; \theta).
```

This gives:
- Low-resolution solver speed  
- High-resolution solver accuracy  
- End-to-end differentiability within PyTorch/JAX/Flux  

This is an **in-the-loop** model.

How to train the ``\text{NN}(u^{n+1}_\text{low-res}, p^{n+1}_\text{low-res}; \theta)``?

- Option 1: use pairs of low fidelity / high fidelity timesteps as training data
- Option 2: match outputs of hybrid solver to high-fidelity simulation directly

"""

# ╔═╡ 58a7c9e9-288e-44a9-95bb-d3fdb34be878
md"""
#### How to Train Hybrid Approaches: Autodifferentiation

Key insight:
Autodiff can differentiate **entire solvers**, not just neural networks.

Chains of computation:

```math
\text{HybridSolver}(u^0;\theta) \longrightarrow u^T
```

Loss:

```math
L(\theta) = \| u^T_{\text{hybrid}} - u^T_{\text{high-fidelity}} \|^2
```

Autodiff computes:

```math
\frac{dL}{d\theta} 
```

even when the solver includes:
- implicit steps  
- iterative solvers  
- PDE discretizations  

This is the foundation of *differentiable physics*.

"""


# ╔═╡ f7ad82aa-4dc4-4657-8792-3f894a2116f1
md"""
#### How to implement hybrid workflows in practice?

- Step 1: rewrite your traditional numerical methods algorithm in an autodifferentiation framework (e.g. PyTorch/JAX)
- Step 2: make parts of this algorithm learnable (either to accelerate it, or to improve accuracy)
- Step 3: get some training examples of what you want the input/output of the algorithm to be
- Step 4: train your algorithm by (auto)differentiating through it and using gradient descent

"""

# ╔═╡ dae17c47-05e3-44c9-b23f-484b422db9f3
md"""
#### When To Use Hybrid Models?

Use hybrid models if:
- Physics is mostly known but contains uncertain closure terms  
- High-fidelity simulation is expensive  
- Low-fidelity model is fast but inaccurate  
- You want physical interpretability + ML flexibility  
- You want faster simulations inside an optimization loop (inverse modeling, UQ)

Avoid hybrid models if:
- You have extremely abundant data → pure ML may be simpler  
- Physics is fully known and solver is fast → pure PDE solver may suffice  
- Solver is not differentiable and hard to rewrite  

Hybrid modeling sits between physics and ML and gets the best of both.

"""


# ╔═╡ e17bc021-7348-4fde-89ea-7f7606df70f5
md"""
#### Summary

1. Hybrid modeling adds ML only where physics is inaccurate or expensive.  
2. Residual modeling is the simplest hybrid approach.  
3. Solver-in-the-loop methods insert ML inside PDE solvers.  
4. Autodifferentiation makes the whole solver learnable end-to-end.  
5. Hybrid methods are powerful for large-scale PDE problems, multiphysics, and closure relations.

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
# ╟─fec68b18-c0ac-11f0-8279-b533b0b6f33f
# ╟─22077fa2-657d-470b-bea1-9e13179fd526
# ╟─f42be050-2ebf-4513-995a-1dfed6ff6f86
# ╟─05dc4602-db65-4565-aaa2-aa4ac6bb5320
# ╟─562f69fc-ac0f-4e46-a908-9982b96da5ea
# ╟─58a7c9e9-288e-44a9-95bb-d3fdb34be878
# ╟─f7ad82aa-4dc4-4657-8792-3f894a2116f1
# ╟─dae17c47-05e3-44c9-b23f-484b422db9f3
# ╟─e17bc021-7348-4fde-89ea-7f7606df70f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
