### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 1816242a-8054-11f0-152a-1d9edcfd46bc
md"""
# Module 3: Finite Difference Approximation
"""

# ╔═╡ f2cf33d1-948f-4d68-b52b-c9b750accbd7
md"""
## Numerical integration

- Rectangle rule (local: 3rd-order, global: 2nd-order)  
- Trapezoid rule (local: 3rd-order, global: 2nd-order)  
- Simpson’s rule (local: 5th-order, global: 4th-order)  
- Richardson extrapolation & Romberg integration (The idea of error cancelation)  
- Improve numerical integration with non-uniform grid spacing  
    - Adaptive quadrature  
    - Gauss quadrature  
"""


# ╔═╡ c1be78ce-9c9a-4f92-b12e-296679f4de8a
md"""
### Problem Statement

- Given a set of N+1 data points $x_i, f(x_i)$ with $i = 0,1,2, \dots, N$. Assume evenly spaced grid points with $x_{i+1} - x_i = h \; (\text{or } \Delta x)$  

- Approximate the derivative $f'(x_i)$  


### Two Methods

- Analytically differentiate interpolating function  
- Approximate the derivative with finite differences*  
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
# ╠═1816242a-8054-11f0-152a-1d9edcfd46bc
# ╠═f2cf33d1-948f-4d68-b52b-c9b750accbd7
# ╠═c1be78ce-9c9a-4f92-b12e-296679f4de8a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
