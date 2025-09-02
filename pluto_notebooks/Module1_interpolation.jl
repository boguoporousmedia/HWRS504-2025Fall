### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 79edeea7-beb8-4c8c-ad8a-e0ec889adf18
begin
    using PlutoUI
    using Plots
    using Polynomials
    using LaTeXStrings
    using Dierckx
end

# ╔═╡ 4d24e6e4-7bb2-11f0-3dd5-0ffe333a3418
md"""
# Module 1: Interpolation methods
"""

# ╔═╡ 33e65d3b-f855-496b-976b-1115575ceff1
md"""
## Review of lecture 1
1. **What do numerical methods do? Why are they useful and important?**

2. **Review**
   - ✔ The advection–diffusion equation that describes the transport of a contaminant  
   - ✔ Linear algebra (some of the concepts are listed below)  
       * Solution of `` {\bf{A}} \cdot {\bf{x}} = {\bf{b}} ``  
       * Norms ($L_1, L_2, L_\infty$); for both vector space and function space  
       * Adjoint operator; self-adjoint operator  
       * Eigenvalue problem; eigenvector; for both matrix and differential operators  
   - ✔ Notations from vector calculus  
   - ✔ Concepts of continuity; Taylor’s theorem
"""

# ╔═╡ 1b8b0759-9a2a-495f-9050-e15a118307f1
md"""
### General problem  

Given a finite set of values of some function ``f(x)``, say  
``f_1 \equiv f(x_1), f_2 \equiv f(x_2), \ldots, f_N \equiv f(x_N) ``,  
with the corresponding points $\{x_i\}_{i=1,2,\ldots,N}$.  

**Question:** How can we construct an approximation to $f(x)$ for all $x \in [a,b]$?

*Note:* Interpolation theory forms the basis for much of numerical analysis.  

---

### Theorem  

Let ``f(x) \in C^0[a,b]``, and let ``\delta`` be any positive number.  
Then, there is a polynomial ``P(x)`` such that for all ``x \in [a,b]``,  

```math
| f(x) - P(x) | < \delta
```

---

**This provides initial motivation to look to polynomials**

---

*Note:* We already have a polynomial representation, but this is based on function and derivative values at a single point.  

The **Taylor polynomial**:  

```math
\begin{align}
P_{T,n}(x) &= f(x_0) + \left.\frac{df}{dx}\right|_{x_0} (x-x_0) + \left.\frac{d^2 f}{dx^2}\right|_{x_0} \frac{(x-x_0)^2}{2!} + \left.\frac{d^3 f}{dx^3}\right|_{x_0} \frac{(x-x_0)^3}{3!} \\
&\quad + \left.\frac{d^{\,n-1} f}{dx^{\,n-1}}\right|_{x_0} 
      \frac{(x-x_0)^{n-1}}{(n-1)!} + R_n
\end{align}
```

However, $P_T$ is ineffective in the general interpolation problem defined above.
"""


# ╔═╡ 4127c185-3e75-4436-8d92-2054f1d24500
md"""
**Q:** Can we define an interpolating polynomial that passes through N interpolating points?  
**A:** Yes.  

---

### Theorem  
Given N distinct points ``x_1, x_2, …, x_2`` and N values ``f_1, f_2, …, f_n,``  
there exists a unique ``(N−1)^{\text{th}}`` degree polynomial ``P_{n-1}(x)`` for which  

```math
P(x_i) = f_i, \quad i=1,2,\ldots,N  
```

---

### How to form this polynomial?  

Let  

```math
p(x) = \alpha_0 + \alpha_1 x + \alpha_2 x^2 + \cdots + \alpha_{N-1} x^{N-1}  
```

with the interpolation conditions  

```math
p(x_i) = f_i, \quad (i=1,2,\ldots,N).  
```

This leads to the **Vandermonde system**:  

```math
\begin{bmatrix}
1 & x_1 & x_1^2 & \cdots & x_1^{N-1} \\
1 & x_2 & x_2^2 & \cdots & x_2^{N-1} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & x_N & x_N^2 & \cdots & x_N^{N-1}
\end{bmatrix}
\begin{bmatrix}
\alpha_0 \\ \alpha_1 \\ \vdots \\ \alpha_{N-1}
\end{bmatrix}
=
\begin{bmatrix}
f_1 \\ f_2 \\ \vdots \\ f_N
\end{bmatrix}  
```

or compactly,  

```math
{\bf{B}} \, {\bf{\alpha}} = {\bf{f}}  
```

where ``{\bf{B}}`` is the Vandermonde matrix.  

---

### Comments  
1. If the objective is **data interpolation**, this choice can be a bad one (ill-conditioning of Vandermonde).  
2. We can actually write the polynomial very easily **without solving the matrix equation** (e.g., Lagrange interpolation).                         
"""        


# ╔═╡ 78c3b4a3-e3e0-42e0-adcd-24ccfdbad999
md"""
## Langrange polynomials
"""

# ╔═╡ 0f15e3c2-f8aa-4f47-a4ee-6082227ee65e
md"""
The ``(N-1)^{\text{th}}`` degree polynomial ``P(x)`` has ``N`` undetermined coefficients ("degree of freedom").

We can write this equivalently as:

```math
P(x) = \sum_{j=1}^N C_j P_j(x)
```

Given polynomials ``P_j(x)`` of at most degree ``(N-1)``

Example: 
```math
C_j = \alpha_{j-1}, \quad P_j = x^{j-1} \quad (j = 1,2,\ldots,N)
```

---

Consider a second example, where $P_j(x)$ is a polynomial defined such that

```math
P_j(x_k) = \delta_{jk} =
\begin{cases}
1, & j = k \\
0, & j \neq k
\end{cases}
```

where $\delta_{jk}$ is the *Kronecker delta*.

---

For such a case,

```math
P(x_k) = \sum_{j=1}^N C_j P_j(x_k) = C_k = f_k
```

---

The set of $(N-1)$ degree polynomials that satisfy $P_j(x_k) = \delta_{jk}$ are called **Lagrange Polynomials**, and denoted by $l_j(x)$. Thus, the polynomial that passes through the $N$ points $\{ f(x_i) \}_{i=1}^N$ is written as

```math
P(x) = \sum_{j=1}^N f_j l_j(x) \; (= \hat{f}(x))
```

We use ``\hat{f}(x)`` as a general symbol for an approximation to ``f(x)``
"""

# ╔═╡ 0d507d1b-6ffe-4863-bae6-3699866684e7
md"""
**Q: How can we determine ``l_j(x)``?**

Require
```math
l_j(x_k) = 0 \quad \forall k \neq j 
```
This gives ``(N-1)`` constraints.

---

Thus, the points ``x_1, x_2, \ldots, x_{j-1}, x_{j+1}, \ldots, x_N`` are roots of the polynomial.  
So we can write

```math
l_j(x) = A_j (x - x_1)(x - x_2) \cdots (x - x_{j-1})(x - x_{j+1}) \cdots (x - x_N)
```

---

Final requirement is 
```math
l_j(x_j) = 1
```

Thus,
```math
1 = A_j (x_j - x_1)(x_j - x_2) \cdots (x_j - x_{j-1})(x_j - x_{j+1}) \cdots (x_j - x_N)
```

---

Therefore,
```math
A_j = \frac{1}{(x_j - x_1)(x_j - x_2)\cdots(x_j - x_{j-1})(x_j - x_{j+1})\cdots(x_j - x_N)}
```
"""

# ╔═╡ 8b82eb2d-79a0-4970-a075-9217c94092b4
md"""
**How do we determine the Lagrange basis ``l_j(x)``?**

Require the ``(N-1)`` constraints
```math
l_j(x_k)=0 \quad \text{for all } k \ne j .
```

Hence the points ``x_1,x_2,\ldots,x_{j-1},x_{j+1},\ldots,x_N `` are roots of ``l_j``, so
```math
l_j(x) = A_j (x-x_1)(x-x_2)\cdots(x-x_{j-1})(x-x_{j+1})\cdots(x-x_N).
```

Imposing ``l_j(x_j)=1`` gives
```math
1 = A_j (x_j-x_1)(x_j-x_2)\cdots(x_j-x_{j-1})(x_j-x_{j+1})\cdots(x_j-x_N),
```
so
```math
A_j = \frac{1}{(x_j-x_1)(x_j-x_2)\cdots(x_j-x_{j-1})(x_j-x_{j+1})\cdots(x_j-x_N)}.
```

Equivalently, the compact product form is
```math
l_j(x) = \prod_{\substack{i=1 \\ i\neq j}}^{N}\frac{x-x_i}{x_j-x_i}.
```
"""

# ╔═╡ 3cb4aabb-08fb-4f41-b02b-d234d71767a7
begin
	# 1) N = 5 interpolation nodes.
	# 2) Equally spaced nodes; change `xnodes` to any set you like.
	# 3) We plot l1(x) and l4(x) 
	local xnodes = collect(1.0:5.0)              # [x₁, …, x₅]
	local N = length(xnodes)

	# Lagrange basis constructor: returns a function ℓ_j(x)
	function lagrange_basis(j, xnodes)
	    function ℓ(x)
	        prod((x - xnodes[i])/(xnodes[j] - xnodes[i]) for i in eachindex(xnodes) if i != j)
	    end
	end

	ℓ1 = lagrange_basis(1, xnodes)
	ℓ4 = lagrange_basis(4, xnodes)

	# Plot range (a bit wider than the nodes)
	local xs = range(xnodes[1] - 0.5, xnodes[end] + 0.5, length = 600)

	# Plot two basis functions; zeros at all nodes; 1 at their own nodes.
	local p = plot(xs, ℓ1.(xs), lw=3, label="l₁(x)", legend=:outertopright)
	plot!(p, xs, ℓ4.(xs), lw=3, label="l₄(x)")

	# x-axis and node markers
	hline!(p, [0], lw=2, label="")
	scatter!(p, xnodes, zeros(N), ms=5, label="")           # nodes on the axis
	scatter!(p, [xnodes[1]],[1.0], ms=7, label="")          # l₁(x₁)=1
	scatter!(p, [xnodes[4]],[1.0], ms=7, label="")          # l₄(x₄)=1

	# Cosmetic axes settings
	xticks!(p, (xnodes, ["x₁","x₂","x₃","x₄","x₅"]))
	ylims!(p, -0.4, 1.2)
	xlabel!(p, "x")
	ylabel!(p, "value")
	title!(p, "Examples of Lagrange Basis Polynomials")

	p
end

# ╔═╡ ea9c422e-3353-4be9-a269-fc65b614410f
md"""
## Runge phenomenon
Oscillation occurs near the boundaries when using high-degree polynomials.

**Example:** Fit a cubic to 4 points ``(x_i, f_i)`` with ``x_1 < x_2 < x_3 < x_4``.  
Let ``\hat f(x)`` denote the cubic interpolant that satisfies ``\hat f(x_i)=f_i``.

Near the left boundary ``[x_1,x_2]``, an interior point ``x^*\in(x_1,x_2)`` can exhibit overshoot/undershoot:

**Case A (overshoot):**
```math
f(x^*) > f_1, \qquad f(x^*) > f_2 .
```

**Case B (undershoot):**
```math
f(x^*) < f_1, \qquad \hat f(x^*) < f_2 .
```

A useful diagnostic is that the **signs of the derivatives** of ``f`` and its interpolant ``\hat f`` can be opposite near the boundary:
```math
\operatorname{sign}\!\big(f'(x)\big) = -\,\operatorname{sign}\!\big(\hat f'(x)\big), \qquad x \in [x_1,x_2].
```

*Assumptions:* We consider distinct nodes ``\{x_i\}`` (e.g., equally spaced) and a smooth underlying function ``f``. The symbol ``\hat f`` denotes the unique polynomial interpolant through the data ``\{(x_i,f_i)\}``.
"""


# ╔═╡ ff9b8202-021f-4e6e-9b2f-1c46e8d767f0
begin
	# Four x–coordinates (equispaced for simplicity)
	local xnodes = [1.0, 2.0, 3.0, 4.0]
	
	# Sample values (you can change these)
	local ynodes_1 = [0.0, 0.0, 0.0, 1.0]
	local ynodes_2 = [0.0, 0.0, 0.0, -1.0]

	spl_1 = Spline1D(xnodes,ynodes_1;k=3,s=0.0)
	spl_2 = Spline1D(xnodes,ynodes_2;k=3,s=0.0)

	local xf = range(first(xnodes), last(xnodes), length=400)

	yf_1 = spl_1.(xf)
	yf_2 = spl_2.(xf)
	
	# Plot
	local plt_1=plot(xf, yf_1, label="cubic spline", lw=2)
	local plt_2=plot(xf, yf_2, label="cubic spline", lw=2)
	ylims!(plt_1, -1.2, 1.2)
	ylims!(plt_2, -1.2, 1.2)
	
	scatter!(plt_1, xnodes, ynodes_1, label="data", ms=5)
	scatter!(plt_2, xnodes, ynodes_2, label="data", ms=5)
		
	# Cosmetic touches: axis, ticks, labels
	hline!(plt_1, [0.0], lw = 1, color=:black, label="")
	xticks!(plt_1, (xnodes, ["x₁","x₂","x₃","x₄"]))
	xlabel!(plt_1, "x"); ylabel!(plt_1, "value")

	hline!(plt_2, [0.0], lw = 1, color=:black, label="")
	xticks!(plt_2, (xnodes, ["x₁","x₂","x₃","x₄"]))
	xlabel!(plt_2, "x"); ylabel!(plt_2, "value")
	
	plot(plt_1, plt_2; layout=(1,2), plot_title="Cubic spline interpolation of 4 points")
	
end


# ╔═╡ 3c390511-a3e9-4a7c-8a93-0964a11b96ad
md"""
	## Runge's Phenomenon (Julia + Pluto)

	Explore polynomial interpolation of the Runge function
	\[
	f(x) = \frac{1}{1+25x^2}, \quad x\in[-1,1],
	\]
	using either **equispaced** or **Chebyshev** nodes. Increase the degree to see oscillations appear near the interval ends for equispaced nodes.  
	Use the toggles to show **bars at node locations** and a **bar chart of node values**.
	"""

# ╔═╡ 320494de-6abb-468a-96d7-f249338bfba0
@bind n Slider(2:30; default=12, show_value=true)

# ╔═╡ db0eb000-49fc-4eb2-bf3f-6c3b82afbfee
@bind nodetype Select(["Equispaced","Chebyshev"]; default="Equispaced")

# ╔═╡ 467bfc49-f15a-4606-842d-d8499bf3be47
begin
    local f(x) = 1/(1 + 25x^2)
    local a, b = -1.0, 1.0
    local xs = range(a, b; length=2000)
    local fs = f.(xs)

    # Nodes
    if nodetype == "Equispaced"
        X = collect(range(a, b; length=n+1))
    else
        k = 0:n
        X = cos.((2k .+ 1) .* (pi/(2n+2)))
        X = (X .- (-1)) .* (b-a)/2 .+ a
    end
    Y = f.(X)

    # Interpolating polynomial
    poly = fit(X, Y)         
    ys = poly.(xs)
	nothing

	local plt = plot(xs, fs, label="f(x) = 1/(1+25x²)", legend=:outertopright)
	plot!(plt, xs, ys, label="degree $n interpolant ($nodetype)")
	scatter!(plt, X, Y, ms=5, label="nodes")
	plt
end

# ╔═╡ 824f0fe8-0d8d-4f7b-bc2e-52f9f78c860a
md"""
## Piecewise Interpolation

- **Chebyshev** points reduce oscillations of high–degree interpolants.
- However, if we’re **given** values on **equispaced** points, can we do better?
- **Idea:** use *localized*, lower–degree interpolants on subintervals and **piece** them together.

**Example:** Piecewise **linear** interpolation through 4 points.

*Locality:* Moving one node only changes the two segments adjacent to that node.  
*Trade-off:* we give up smoothness at the joints (continuous ``C^0``, not ``C^1``).

```math
\text{Piecewise linear } \Rightarrow \text{continuous (}C^0\text{), slopes can jump at nodes.}
```

"""

# ╔═╡ ea2d65f5-d7c8-4036-8873-eca549529a93
begin
	# Four x–coordinates (equispaced for simplicity)
	local xnodes = [1.0, 2.0, 3.0, 4.0]
	local ynodes = [0.0, 0.0, 0.0, 1.0]
	# Make the piecewise-linear interpolant: straight lines between nodes
	local plt = plot(xnodes, ynodes;
	           seriestype = :path,
	           marker = :circle,
	           lw = 3,
	           label = "piecewise linear", 
			   legend=:outertopright)
	
	# Add some “what-if” alternatives for the last point (dashed)
	alt_y4 = [0.5, -0.5]
	for y4 in alt_y4
	    plot!(plt,
	          [xnodes[end-1], xnodes[end]], [ynodes[end-1], y4];
	          ls = :dash, lw = 2, label = "")
	    scatter!(plt, [xnodes[end]], [y4], ms=5, label = "")
	end
	
	# Cosmetic touches: axis, ticks, labels
	hline!(plt, [0.0], lw = 1, color=:black, label="")
	xticks!(plt, (xnodes, ["x₁","x₂","x₃","x₄"]))
	xlabel!(plt, "x"); ylabel!(plt, "value")
	title!(plt, "Piecewise Linear Interpolation of 4 points")
	
	plt
end


# ╔═╡ 2706cd9b-4b4c-4b4d-9cd7-f05a99dffaec
md"""
We can make the following observations:

1. With **piecewise linear** polynomials, the interpolant ``\hat f(x)`` is continuous but has slope jumps at nodes:
```math
\hat f(x) \in C^{0}[a,b].
```
"""

# ╔═╡ 1e519a22-5f96-45f0-a982-9bb363741178
md"""
With a single global $(N-1)$-degree polynomial, the interpolant is smooth:

```math
\hat f(x) \in C^{\infty}[a,b].
```

2. **Piecewise linear segment** on $[x_j, x_{j+1}]$ depends only on $f_j, f_{j+1}$:

```math
\hat f(x) = \frac{f_{j+1}-f_j}{x_{j+1}-x_j}\,(x-x_j) + f_j
          = f_j\,\frac{x_{j+1}-x}{x_{j+1}-x_j} + f_{j+1}\,\frac{x-x_j}{x_{j+1}-x_j},
\qquad x_j \le x \le x_{j+1}.
```

Thus the value at node $j$ only influences the two pieces on $[x_{j-1},x_j]$ and $[x_j,x_{j+1}]$.

3. The **piecewise-linear Lagrange basis** (hat function) $\phi_j(x)$ is

```math
\phi_j(x) \equiv
\begin{cases}
\dfrac{x-x_{j-1}}{x_j-x_{j-1}}, & x_{j-1} \le x \le x_j, \\
\dfrac{x_{j+1}-x}{x_{j+1}-x_j}, & x_j \le x \le x_{j+1}, \\
0, & x < x_{j-1} \;\text{or}\; x > x_{j+1}.
\end{cases}
```

The interpolant can then be written as

```math
\hat f(x) = \sum_{j=1}^{N} f_j \,\phi_j(x).
```

"""

# ╔═╡ 4fefdcb1-a12b-4122-b6a5-2128dba7df5c
begin
    # Nodes (change as you like)
    local xnodes = [1.0, 2.0, 3.0, 4.0, 5.0]
    j = 3  # plot φ_j centered at x_j (here j=3 → neighbors are x₂ and x₄)

    @assert 2 ≤ j ≤ length(xnodes)-1 "Choose j with neighbors: 2 ≤ j ≤ N-1"

    # Factory for the piecewise-linear basis φ_j
    function phi_hat_factory(xnodes::AbstractVector, j::Int)
        xjm1, xj, xjp1 = xnodes[j-1], xnodes[j], xnodes[j+1]
        function ϕ(x)
            if x < xjm1 || x > xjp1
                return zero(float(x))
            elseif x <= xj
                return (x - xjm1) / (xj - xjm1)
            else
                return (xjp1 - x) / (xjp1 - xj)
            end
        end
    end

    ϕj = phi_hat_factory(xnodes, j)

    # Plot range around the support [x_{j-1}, x_{j+1}]
    local xs = range(xnodes[j-1] - 1, xnodes[j+1] + 1, length=600)

    local p = plot(xs, ϕj.(xs), lw=3, label=latexstring("φ_{", j, "}(x)"), legend=:outertopright)
    scatter!(p, [xnodes[j-1], xnodes[j], xnodes[j+1]], [0.0, 1.0, 0.0],
             ms=7, label="nodes")
    xt = xnodes[j-1:j+1]
    xticks!(p, (xt, [L"x_{j-1}",L"x_j",L"x_{j+1}"]))
    xlims!(p, minimum(xs), maximum(xs))
    ylims!(p, -0.05, 1.15)
    xlabel!(p, "x"); ylabel!(p, "value")
    title!(p, "Piecewise-linear Lagrange basis: " * string(L"φ_j(x)"))

    p
end

# ╔═╡ 8d1b8e82-5e3e-4627-aec9-a1ecf8723cdf
md"""
# Functional definition of piecewise polynomials

Let the nodes be ``x_1 < x_2 < \cdots < x_n`` and define the element (edge)
```math
e_i = [x_i,\, x_{i+1}],\quad i=1,\dots,n-1.
```

## Piecewise linear case (``C^0``-continuous)

Within each element $e_i$, approximate $f(x)$ by an affine function

```math
\hat f_i(x) = a_i + b_i x,\quad x\in e_i,
```

subject to nodal interpolation:

```math
\hat f_i(x_i) = f_i,\qquad \hat f_i(x_{i+1}) = f_{i+1}.
```

Solving gives the familiar barycentric (Lagrange) form on $e_i$:

```math
\hat f_i(x)
= f_i\,\frac{x_{i+1}-x}{x_{i+1}-x_i}
+ f_{i+1}\,\frac{x-x_i}{x_{i+1}-x_i},
\qquad x\in[x_i,x_{i+1}].
```

Equivalently, using local *hat* basis functions $\{\phi_j(x)\}$ with

```math
\phi_j(x_k)=\delta_{jk},\quad
\operatorname{supp}\phi_j=[x_{j-1},x_{j+1}],
```

the global interpolant is

```math
\hat f(x)=\sum_{j=1}^n f_j\,\phi_j(x).
```

On a uniform grid, each $\phi_j$ is piecewise linear, rising on $[x_{j-1},x_j]$ and falling on $[x_j,x_{j+1}]$.

"""

# ╔═╡ d1e04d74-8320-4b56-b895-f1bf5f41fe21
begin
	img1 = LocalResource("./figs/mod1_piecewise_linear_f.png");
	img1
end

# ╔═╡ 6a9b3e15-360e-49c4-98f4-dcbc2afc2524
md"""
## Representation with Lagrange pieces

For all ``x \in [a,b]``, a piecewise polynomial interpolant can be written as
```math
\hat f(x)=\sum_{i=1}^{N} f_i\,\phi_i(x),
```

where ``\phi_i(x)`` are the local *shape functions*.
For **piecewise linear** interpolation, each ``\phi_i`` is a hat function with
support ``[x_{i-1},x_{i+1}]``, continuous across nodes (``C^0`` on ``[a,b]``).

---

### Q: Can we use higher–degree Lagrange pieces?

**A:** Yes.

One can use **piecewise quadratic** Lagrange polynomials (or higher).  The
functions are still tied together to ensure $C^0$ continuity across element
interfaces.  With quadratic pieces, each basis function spans two elements and
the same global representation holds:

```math
\hat f(x)=\sum_{i=1}^{N} f_i\,\phi_i(x),
```

now with $\phi_i$ being **piecewise quadratic Lagrange polynomials**.

**Notes**

* Linear pieces ⇒ nodes at ``\{x_i\}_{i=1}^{N}`` and hat-shaped ``\phi_i``.
* Quadratic pieces ⇒ include midpoints (or other local nodes) in each element,
  giving curved segments while maintaining $C^0$ continuity.
* The coefficients ``f_i`` are typically the function values at the chosen
  interpolation nodes.
"""

# ╔═╡ 316fda26-172e-442c-999a-e05456d0cfa4
begin
	img2 = LocalResource("./figs/mod1_piecewise_quadratic_nonsmooth.png");
	img2
end

# ╔═╡ 4d58a662-5f99-498d-a17d-3acdebf30471
md"""
## How about increased smoothness?

```math
C^0[a,b] \;\rightarrow\; C^1[a,b]
```

* This leads to piecewise polynomials that satisfy a mixture of interpolation constraints and smoothness constraints. Such polynomials are called **splines**.

* **Note:** the additional constraints (smoothness) will require an increase in the degree of local polynomials.

"""

# ╔═╡ 56dba066-0d87-4c8c-9cf7-8837ff29d627
begin
	img3 = LocalResource("./figs/mod1_quadratic.png");
	img3
end

# ╔═╡ f3a4824b-e261-4f3b-8684-bbe7cf50e590
md"""
## The `` C^{1}[a,b] `` case

Let ``x_1< x_2<\cdots< x_N `` and elements `` e_i=[x_i,x_{i+1}] `` for `` i=1,\dots,N-1``.
On each element use a quadratic:
```math
\hat f_i(x)=a_i+b_i x+c_i x^{2}, \qquad x\in e_i .
```

**Unknowns**
```math
\text{coefficients per element}=3 \;\;\Rightarrow\;\; 3(N-1)\ \text{unknowns.}
```

**Constraints**

- *Interpolation (values)* for each element:
```math
\hat f_i(x_i)=f(x_i), \qquad \hat f_i(x_{i+1})=f(x_{i+1}),\quad i=1,\dots,N-1
```
giving
```math
2(N-1)\ \text{constraints}.
```

- *Smoothness (C¹ continuity)* at interior nodes:
```math
\frac{d\hat f_i}{dx}(x_{i+1})=\frac{d\hat f_{i+1}}{dx}(x_{i+1}),
\quad i=1,\dots,N-2,
```
giving
```math
(N-2)\ \text{constraints}.
```

**Total constraints**
```math
2(N-1)+(N-2)=3N-4.
```

Since ``3(N-1)=3N-3`` unknowns ``>`` ``3N-4`` constraints, **one more condition is needed**.
Common choices (one of the following):

- *Clamped* at the left end:
```math
\hat f'(x_1)=f'(x_1),
```
- or *Clamped* at the right end:
```math
\hat f'(x_N)=f'(x_N).
```

*(That additional boundary condition is either provided or assumed.)*

"""

# ╔═╡ 1b9f799f-6709-479a-bfc1-f368c30d1d6a
begin
	img4 = LocalResource("./figs/mod1_cubic.png");
	img4
end

# ╔═╡ a2da4b6e-ca17-4abd-a08d-d4291448c0ea
md"""
## The ``C^2[a, b]`` case

- Interpolation defined by cubic pieces between each \([xᵢ, xᵢ₊₁]\)  
  ⇒ **Cubic Splines**

- “Ultimate” interpolation
  - Sufficiently smooth (C²-continuous)
  - No oscillations

---

```math
\begin{aligned}
f_1(x) &= a_1 + b_1 x + c_1 x^2 + d_1 x^3, \quad x \in [x_1, x_2] \\
f_2(x) &= a_2 + b_2 x + c_2 x^2 + d_2 x^3, \quad x \in [x_2, x_3] \\
f_3(x) &= a_3 + b_3 x + c_3 x^2 + d_3 x^3, \quad x \in [x_3, x_4]
\end{aligned}
```

* Total number of undetermined coefficients:

  ```math
  4(N-1)
  ```

* Total number of constraints:

  * Interpolation:

    ```math
    2(N-1)
    ```
  * Smoothness:

    ```math
    2(N-2)
    ```
  * Total:

    ```math
    4N - 6
    ```

---

We need **two more constraints**.
Usually we impose:

```math
\frac{d^2 f_1}{dx^2}(x_1) = \frac{d^2 f}{dx^2}(x_1) \quad \text{or } 0
```

```math
\frac{d^2 f_N}{dx^2}(x_N) = \frac{d^2 f}{dx^2}(x_N) \quad \text{or } 0
```

* If imposed as $0$, it is called **natural splines**.

"""

# ╔═╡ 31d49e7a-5b7a-4e0c-b399-7ff7dcf84ee7
begin
	# target function
	f1(x) = exp(-x^2)

	# interpolation nodes
	local xnodes = collect(-2.5:0.5:2.5)
	local ynodes = f1.(xnodes)

	# finite-difference derivative at nodes (for markers on right plot)
	fd_deriv = similar(xnodes)
	local n = length(xnodes)
	for i in 1:n
		if i == 1
			fd_deriv[i] = (ynodes[i+1] - ynodes[i]) / (xnodes[i+1] - xnodes[i])
		elseif i == n
			fd_deriv[i] = (ynodes[i] - ynodes[i-1]) / (xnodes[i] - xnodes[i-1])
		else
			fd_deriv[i] = (ynodes[i+1] - ynodes[i-1]) / (xnodes[i+1] - xnodes[i-1])
		end
	end

	# fine grid for plotting
	local xf = range(-3, 3; length=800)
	(xnodes, ynodes, fd_deriv, xf)

	# degree n-1 polynomial that interpolates (xnodes, ynodes)
	pL = fit(Polynomial, xnodes, ynodes)  # Vandermonde fit → exact interpolation at nodes
	dpL = Polynomials.derivative(pL)                   # analytic derivative
	(pL, dpL)

	    # cubic spline, s=0 forces interpolation; default bc is not-a-knot
    sp = Spline1D(xnodes, ynodes; k=3, s=0.0)

    # spline derivative (scalar x)
    dsp(x) = Dierckx.derivative(sp, x)

    (sp, dsp)
	
	local p1 = plot(xf, f1.(xf), label="true f(x)", linestyle=:dash, alpha=0.6)
    scatter!(p1, xnodes, ynodes; label="Data", ms=5, mc=:white, msc=:blue, shape=:square)
    plot!(p1, xf, pL.(xf), label="Lagrange")
    plot!(p1, xf, sp.(xf), label="Cubic spline")
    xlabel!(p1, L"x"); ylabel!(p1, L"\hat{f}(x)"); title!(p1, L"f(x)=e^{-x^2}")
	ylims!(p1,(-1,2))

    local p2 = plot(xf, -2 .* xf .* exp.(-xf.^2), label="true f'(x)", linestyle=:dash, alpha=0.6)
    scatter!(p2, xnodes, fd_deriv; label="FD at nodes", ms=5, mc=:white, msc=:blue, shape=:square)
    plot!(p2, xf, dpL.(xf), label="Lagrange")     # use polynomial derivative
    plot!(p2, xf, dsp.(xf), label="Cubic spline") # use spline derivative
    xlabel!(p2, L"x"); ylabel!(p2, L"\hat{f}'(x)")
	ylims!(p2,(-2.5,2.5))
    plot(p1, p2; layout=(1,2))
	
end


# ╔═╡ e75a479d-cec4-40b0-8953-2fee3986e90a
md"""
### Summary of interpolation

- Smooth degree ``(N-1)^{\text{th}}`` Lagrange polynomial
  - Strong oscillations appear near the boundaries of the interval when using equally spaced interpolation points.
  - Oscillations **can** be remediated by smartly choosing the locations of points ``\{x_i\}_{i=1}^N``, though not a good option when data points are given.

- Piecewise lower degree Lagrange polynomials
  - No oscillations
  - The most popular one is the cubic splines, which is ``C^2[a,b]``.

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Dierckx = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"

[compat]
Dierckx = "~0.5.4"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.18"
PlutoUI = "~0.7.69"
Polynomials = "~4.1.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "fef1de61baf708620b744ef746e6fc1be9fd3c40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "76b3b7c3925d943edf158ddb7f693ba54eb297a5"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "7da4b14cc4c3443a1afc64abee17f4fcb45ad837"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.4"

[[deps.Dierckx_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3251f44b3cac6fec4cec8db45d3ab0bfed51c4d8"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.2.0+0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "52e1296ebbde0db845b356abbbe67fb82a0a116c"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.9"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ae7d4ddec2e13ad3bddf5c0796f7547cf682391"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.2+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "9a9216c0cf706cb2cc58fd194878180e3e51e8c0"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.18"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "2d7662f95eafd3b6c346acdbfc11a762a2256375"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.69"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "eb38d376097f47316fe089fc62cb7c6d85383a52"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "e1d5e16d0f65762396f9ca4644a5f4ddab8d452b"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "c5bf2dad6a03dfef57ea0a170a1fe493601603f2"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.5+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─79edeea7-beb8-4c8c-ad8a-e0ec889adf18
# ╟─4d24e6e4-7bb2-11f0-3dd5-0ffe333a3418
# ╟─33e65d3b-f855-496b-976b-1115575ceff1
# ╟─1b8b0759-9a2a-495f-9050-e15a118307f1
# ╟─4127c185-3e75-4436-8d92-2054f1d24500
# ╟─78c3b4a3-e3e0-42e0-adcd-24ccfdbad999
# ╟─0f15e3c2-f8aa-4f47-a4ee-6082227ee65e
# ╟─0d507d1b-6ffe-4863-bae6-3699866684e7
# ╟─8b82eb2d-79a0-4970-a075-9217c94092b4
# ╟─3cb4aabb-08fb-4f41-b02b-d234d71767a7
# ╟─ea9c422e-3353-4be9-a269-fc65b614410f
# ╟─ff9b8202-021f-4e6e-9b2f-1c46e8d767f0
# ╟─3c390511-a3e9-4a7c-8a93-0964a11b96ad
# ╟─320494de-6abb-468a-96d7-f249338bfba0
# ╟─db0eb000-49fc-4eb2-bf3f-6c3b82afbfee
# ╟─467bfc49-f15a-4606-842d-d8499bf3be47
# ╟─824f0fe8-0d8d-4f7b-bc2e-52f9f78c860a
# ╟─ea2d65f5-d7c8-4036-8873-eca549529a93
# ╟─2706cd9b-4b4c-4b4d-9cd7-f05a99dffaec
# ╟─1e519a22-5f96-45f0-a982-9bb363741178
# ╟─4fefdcb1-a12b-4122-b6a5-2128dba7df5c
# ╟─8d1b8e82-5e3e-4627-aec9-a1ecf8723cdf
# ╠═d1e04d74-8320-4b56-b895-f1bf5f41fe21
# ╟─6a9b3e15-360e-49c4-98f4-dcbc2afc2524
# ╟─316fda26-172e-442c-999a-e05456d0cfa4
# ╟─4d58a662-5f99-498d-a17d-3acdebf30471
# ╟─56dba066-0d87-4c8c-9cf7-8837ff29d627
# ╟─f3a4824b-e261-4f3b-8684-bbe7cf50e590
# ╟─1b9f799f-6709-479a-bfc1-f368c30d1d6a
# ╟─a2da4b6e-ca17-4abd-a08d-d4291448c0ea
# ╟─31d49e7a-5b7a-4e0c-b399-7ff7dcf84ee7
# ╟─e75a479d-cec4-40b0-8953-2fee3986e90a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
