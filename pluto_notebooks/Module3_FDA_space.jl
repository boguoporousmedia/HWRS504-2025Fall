### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 6b928336-5bad-48f7-b480-1bf337e1607f
using PlutoUI, Plots, LaTeXStrings, LinearAlgebra

# ╔═╡ 1816242a-8054-11f0-152a-1d9edcfd46bc
md"""
# Module 3: Finite Difference Approximation in Space (Boundary Value Problem)
"""

# ╔═╡ 02a95c05-026a-47fb-a325-f00eff3e0dd6
md"""
**Reading: CG Chapter 2**
"""

# ╔═╡ f2cf33d1-948f-4d68-b52b-c9b750accbd7
md"""
## Review of numerical integration

- Rectangle rule (local: 3rd-order, global: 2nd-order)  
- Trapezoid rule (local: 3rd-order, global: 2nd-order)  
- Simpson’s rule (local: 5th-order, global: 4th-order)  
- Richardson extrapolation & Romberg integration (The idea of *"error cancelation"*)  
- Improve numerical integration with non-uniform grid spacing  
    - Adaptive quadrature  
    - Gauss quadrature  
"""


# ╔═╡ c1be78ce-9c9a-4f92-b12e-296679f4de8a
md"""
## Problem Statement

- Given a set of ``N+1`` data points $x_i, f(x_i)$ with $i = 0,1,2, \dots, N$. Assume evenly spaced grid points with $x_{i+1} - x_i = h \; (\text{or } \Delta x)$  

- Approximate the derivative $f'(x_i)$  


### Two Methods

- Analytically differentiate interpolating function  
- Approximate the derivative with finite differences  
"""


# ╔═╡ cc63f48e-9e72-4b55-b3f2-414ce3bae01e
local img = LocalResource("./figs/mod3_FD.png",:width => "400px")

# ╔═╡ eb5ffd0f-4f84-4aa4-ad33-de9178c2218b
md"""
## Finite Difference Approximations

1. Backward difference
```math
\left.\frac{df}{dx}\right|_{x_i} \;=\; 
\frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}
```

2. Forward difference
```math
\left.\frac{df}{dx}\right|_{x_i} \;=\; 
\frac{f(x_{i+1}) - f(x_i)}{x_{i+1} - x_i}
```

3. Central difference
```math
\left.\frac{df}{dx}\right|_{x_i} \;=\; 
\frac{f(x_{i+1}) - f(x_{i-1})}{x_{i+1} - x_{i-1}}
```

- Which one to use? Which one is better? Are there others?

"""


# ╔═╡ b5933bf7-3f82-4e5f-95a7-829c0c0715c6
md"""
## Forward Differencing

- Consider Taylor expansion for $f(x_{i+1})$ about $x_i$:

```math
f(x_{i+1}) = f(x_i) + f'(x_i)(x_{i+1}-x_i) + 
\frac{f''(x_i)}{2!}(x_{i+1}-x_i)^2 + \cdots
```

If we let $\Delta x = x_{i+1} - x_i$:

```math
f(x_{i+1}) = f(x_i) + f'(x_i)\,\Delta x + \frac{f''(x_i)}{2}\,\Delta x^2 + \cdots
```

---

- Solve for $f'(x_i)$:

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{\Delta x} 
- \frac{f''(x_i)}{2}\,\Delta x + \cdots
```

This is the **forward finite difference approximation** with a **truncation error** of $\mathcal{O}(\Delta x)$.
"""


# ╔═╡ a9ed7ccd-4a73-4d04-8d1c-56c795c59f54
md"""
## Example — Forward Difference

- ``f(x)=e^{-x^2}``, ``f'(x)=-2x\,e^{-x^2}``
- Determine derivative at ``x=1``. Exact solution: ``f'(1) = -2e^{-1} = -0.73575888``

We illustrate the order of the truncation error of the *forward difference*

```math
f'(x_0) \approx \frac{f(x_0+h)-f(x_0)}{h}
```

by varying ``h`` and plotting ``|\varepsilon|=\big|\frac{f(x_0+h)-f(x_0)}{h}-f'(x_0)\big|`` on log–log axes.

"""


# ╔═╡ 702a4721-c5ee-40ae-bbf6-0eae201e9d61
begin
	# Function and exact derivative
	f(x)  = exp(-x^2)
	fp(x) = -2x*exp(-x^2)

	local x0 = 1.0
	local hs = 10 .^ range(-6, 0; length=200)   # h from 1e-6 to 1
	local fd(h) = (f(x0 + h) - f(x0)) / h
	local err = abs.(fd.(hs) .- fp(x0))

	# Reference ~ h^1 line (scaled to sit near the data)
	local c = err[1] / hs[1]
	local ref = 10c .* hs

	local plt = plot(hs, err; xscale=:log10, yscale=:log10, lw=2,
		label="Forward Difference", xlabel="h", ylabel="|ε|",
		legend=:bottomright, grid=true)
	plot!(plt, hs, ref; ls=:dash, lw=2, label="~h¹")

	plt
end


# ╔═╡ 45fbb6d6-dd13-4e1c-a2c8-702676aa20a5
md"""
## Backward Differencing

- Similar idea: expand ``f(x_{i-1})`` about ``x_i`` and solve for ``f'(x_i)``:

```math
f(x_{i-1}) = f(x_i) + f'(x_i)(x_{i-1}-x_i) 
+ \frac{f''(x_i)}{2!}(x_{i-1}-x_i)^2 + \cdots
```

If we let ``\Delta x = x_i - x_{i-1}``:

```math
f(x_{i-1}) = f(x_i) - f'(x_i)\,\Delta x + \frac{f''(x_i)}{2}\,\Delta x^2 - \cdots
```

---

- Solve for ``f'(x_i)``:

```math
f'(x_i) 
= \frac{f(x_i)-f(x_{i-1})}{\Delta x} 
+ \frac{\Delta x}{2}\,f''(x_i) + \cdots
```

---

- Same *truncation error* order (``\mathcal{O}(\Delta x)``) as forward differencing, but with the *opposite sign*.  
- Forward vs backward differencing will be important for some PDEs (this is where *upwinding / downwinding* ideas come in).

"""


# ╔═╡ bd6b2079-37e9-45c2-a451-adc504619922
md"""
## Central Differencing

- Can we do better by using both ``f(x_{i-1})`` and ``f(x_{i+1})``?  

---

Subtract the Taylor expansion of ``f(x_{i-1})`` from that of ``f(x_{i+1})``:

```math
f(x_{i+1}) - f(x_{i-1}) 
= 2 f'(x_i)\,\Delta x + \tfrac{1}{3} f^{(3)}(x_i)\,\Delta x^3 + \cdots
```

---

- Solve for ``f'(x_i)``:

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_{i-1})}{2\Delta x}
- \frac{\Delta x^2}{6}\,f^{(3)}(x_i) + \cdots
```

---

- Higher order of accuracy: second order (``\mathcal{O}(\Delta x^2)``).  
- However, requires a wider *stencil* (three-point vs two-point).

"""


# ╔═╡ 3c32791a-dee1-4873-bcf4-73abde73c57d
md"""
## Example

- ``f(x)=e^{-x^2}``, ``\quad f'(x)=-2(2x^2 - 1)\,e^{-x^2}``
- Determine derivative at ``x=1`` (Exact solution: $−0.73575888$)

We compare *forward*, *backward*, and *central* differences and plot
``|\varepsilon|=\big|\widehat{f'}(x_0)-f'(x_0)\big|`` versus ``h`` on log–log axes.
"""


# ╔═╡ ca7c7180-a026-40f4-b8bb-4bb9a8a3df1d
begin
	# Function and exact derivative
	local x0 = 1.0

	# Step sizes
	local hs = 10 .^ range(-6, -0.3; length=250)

	# Discrete derivatives
	forward(h)  = (f(x0 + h) - f(x0)) / h
	backward(h) = (f(x0) - f(x0 - h)) / h
	central(h)  = (f(x0 + h) - f(x0 - h)) / (2h)

	# Errors
	err_f = abs.(forward.(hs)  .- fp(x0))
	err_b = abs.(backward.(hs) .- fp(x0))
	err_c = abs.(central.(hs)  .- fp(x0))

	# Reference lines (~h^1 and ~h^2), offset so they don't overlap with curves
	c1 = err_f[1] / hs[1]              # slope-1 scaling
	ref1 = 10c1 .* hs                # move it down a bit

	c2 = err_c[60] / hs[60]^2          # slope-2 scaling from a mid point
	ref2 = 10c2 .* hs.^2               # move it above the central curve

	local plt = plot(hs, err_f; xscale=:log10, yscale=:log10, lw=2, label="Forward Difference",
	           xlabel="h", ylabel="|ε|", legend=:bottomright, grid=true)
	plot!(plt, hs, err_b; lw=2, label="Backward Difference")
	plot!(plt, hs, err_c; lw=2, label="Central Difference")
	plot!(plt, hs, ref1; ls=:dash, lw=2, label="~h¹")
	plot!(plt, hs, ref2; ls=:dash, lw=2, label="~h²")

	xlims!(plt, 1e-4, 1)

	plt
end


# ╔═╡ 6bafe47c-316c-4a7a-96de-1842ef0f3a0f
md"""
## Higher-Order Derivatives

- The same process applies for higher-order derivatives  

- Example: Second derivative $f''(x_i)$

```math
f(x_{i+1}) = f(x_i) 
+ f'(x_i)h 
+ \tfrac{1}{2} f''(x_i) h^2 
+ \tfrac{1}{6} f^{(3)}(x_i) h^3 
+ \tfrac{1}{24} f^{(4)}(x_i) h^4 
+ \cdots
```

```math
f(x_{i-1}) = f(x_i) 
- f'(x_i)h 
+ \tfrac{1}{2} f''(x_i) h^2 
- \tfrac{1}{6} f^{(3)}(x_i) h^3 
+ \tfrac{1}{24} f^{(4)}(x_i) h^4 
- \cdots
```

From Taylor expansions of $f(x_{i+1})$ and $f(x_{i-1})$:

```math
f(x_{i+1}) + f(x_{i-1})
= 2f(x_i) + f''(x_i) \, h^2 + \tfrac{1}{12} f^{(4)}(x_i) \, h^4 + \cdots
```

---

- Solve for $f''(x_i)$:

```math
f''(x_i) \;=\; \frac{f_{i-1} - 2f_i + f_{i+1}}{h^2}
- \tfrac{1}{12} h^2 f^{(4)}(x_i) + \cdots
```

- This is *second-order* accurate.

"""


# ╔═╡ a13d06f1-7383-4630-a4a5-73c9a54f92ed
md"""
## Nomenclature

- *Stencil*: Extent of grid points used for approximation  

    - First derivative central difference is a three-point stencil even though the coefficient of $f_i$ is zero  
    - Centered Scheme: Stencil is from $f_{i-k}$ to $f_{i+k}$, where $k$ is some integer  
    - One-Sided Scheme: Stencil is from $f_i$ to $f_{i+k}$ (or $f_{i-k}$ to $f_i$), where $k$ is some integer  
    - Forward (Backward)-Biased Scheme: Stencil is from $f_{i-k}$ to $f_{i+l}$, where $l > k$ (or $l < k$)  

---

- *Truncation error*: The portion truncated in a finite difference approximation  

```math
f'(x_i) \;=\; \frac{f(x_{i+1}) - f(x_i)}{h} + \tfrac{h}{2} f''(x_i) + \cdots
```

---

- *Order of approximation*: The order of the lowest-order term in the truncation error
"""


# ╔═╡ 4007ed66-f841-4c45-9a16-2c9d255529b8
md"""
## Requirement for all FDA's

- When the node spacing becomes arbitrarily small, the T.E. $\to 0$, so that the FDA $\to$ actual derivative.  

That is,  

```math
\lim_{\Delta x \to 0} \text{T.E.} = 0 \quad \text{(consistency condition)}
```

Any FDA for which $\lim_{\Delta x \to 0} \text{T.E.} = 0$ is called *consistent*.

---

- Order of approximation

   - The rate at which T.E. $\to 0$ as $\Delta x \downarrow$ is called the *order* of the approximation.  

   - The *leading term* (i.e., the lowest-order term) determines the order of approximation.

"""


# ╔═╡ 7243164e-8b07-4a60-aa74-3a0e650f6229
md"""
## Summary

- The wider the stencil, the greater the accuracy

    - Generally, for an $n$-th derivative computed with an $N$-point stencil, the order of accuracy is $N-n$.  

    - Sometimes, you get lucky and might gain an extra order with some error cancellation, either inherent to the method (e.g., centered second derivative on a uniform mesh) or due to the *particular function* being differentiated.  
    - Wider stencils require more computational cost and more bookkeeping.  

---

- Sometimes, for a given stencil, it might be beneficial to sacrifice accuracy for other beneficial properties. One example is to *avoid oscillation*. We will talk about this more later.

"""


# ╔═╡ 7c34e980-11cb-4744-98c1-4a01f23a9838
md"""
## General Formulation of FDA's (1D)

- All FDA's are of the general form  

```math
\left.\frac{d^m u}{dx^m}\right|_{x_*} \;\approx\; 
\gamma_1 u_1 + \gamma_2 u_2 + \cdots + \gamma_p u_p 
= \sum_{i=1}^p \gamma_i u_i,
```

where $p$ is the number of node points used.  

---

- Consistency Requirement:

```math
\sum_{i=1}^p \gamma_i u_i \;\approx\; 
\left.\frac{d^m u}{dx^m}\right|_{x_*} 
+ O((\Delta x)^r), \quad r > 0
```

---

- Taylor expansion:

```math
u_i = u_* + (x_i - x_*) \left.\frac{du}{dx}\right|_{x_*}
+ \frac{(x_i - x_*)^2}{2!} \left.\frac{d^2u}{dx^2}\right|_{x_*}
+ \cdots
```

"""


# ╔═╡ c1ed2624-2c4b-436f-b0b9-4089c802fcbb
md"""
## General Formulation of FDA's (1D)

- Substitute Taylor series into linear combination:

```math
\sum_{i=1}^p \gamma_i u_i
= F_0 u_* + F_1 \left.\frac{du}{dx}\right|_{x_*}
+ \cdots
+ F_m \left.\frac{d^m u}{dx^m}\right|_{x_*}
+ F_{m+1} \left.\frac{d^{m+1}u}{dx^{m+1}}\right|_{x_*}
+ \cdots
```

---

where:

```math
\begin{aligned}
F_0 &\equiv \gamma_1 + \gamma_2 + \cdots + \gamma_p \\
F_1 &\equiv \gamma_1(x_1 - x_*) + \gamma_2(x_2 - x_*) + \cdots + \gamma_p(x_p - x_*) \\
F_2 &\equiv \frac{\gamma_1(x_1 - x_*)^2}{2} + \frac{\gamma_2(x_2 - x_*)^2}{2} + \cdots + \frac{\gamma_p(x_p - x_*)^2}{2} \\
&\;\;\vdots \\
F_m &\equiv \frac{\gamma_1(x_1 - x_*)^m}{m!} + \frac{\gamma_2(x_2 - x_*)^m}{m!} + \cdots + \frac{\gamma_p(x_p - x_*)^m}{m!} \\
F_{m+1} &\equiv \frac{\gamma_1(x_1 - x_*)^{m+1}}{(m+1)!} + \frac{\gamma_2(x_2 - x_*)^{m+1}}{(m+1)!} + \cdots + \frac{\gamma_p(x_p - x_*)^{m+1}}{(m+1)!} \\
&\;\;\vdots
\end{aligned}
```

---
"""


# ╔═╡ 64687840-fa50-4282-b25c-a223f5ff4cd5
md"""
## General Formulation of FDA's (1D)

- From *Consistency Requirement*: $F_m = 1$

```math
\;\;\;\;\;\; \Rightarrow \text{ each } \gamma_i \sim \frac{1}{(\Delta x)^m}
```

---

- If $\gamma_i \sim \tfrac{1}{(\Delta x)^m}$, then  

```math
F_k \sim \frac{(\Delta x)^k}{(\Delta x)^m} = O\!\big((\Delta x)^{k-m}\big).
```

---

- For consistency, require $F_k$ for $k - m < 0$ to be zero.  

- Thus,  

```math
F_0 = F_1 = \cdots = F_{m-1} = 0, \quad F_m = 1
\;\;\;\;\;\;\; \text{(Required for consistency)}
```

---

- If $F_{m+1} \neq 0$, FDA is $O(\Delta x)$  
- If $F_{m+1} = 0$, and $F_{m+2} \neq 0$, FDA is $O(\Delta x^2)$  
$\cdots$

"""


# ╔═╡ e923d22e-01dd-4469-bb92-ff0a8721824b
md"""
- Consistency Requirements:

```math
\begin{bmatrix}
\frac{1}{(x_1 - x_*)} & \frac{1}{(x_2 - x_*)} & \cdots & \frac{1}{(x_p - x_*)} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{(x_1 - x_*)^{m-1}}{(m-1)!} & \frac{(x_2 - x_*)^{m-1}}{(m-1)!} & \cdots & \frac{(x_p - x_*)^{m-1}}{(m-1)!} \\
\frac{(x_1 - x_*)^m}{m!} & \frac{(x_2 - x_*)^m}{m!} & \cdots & \frac{(x_p - x_*)^m}{m!}
\end{bmatrix}
\begin{bmatrix}
\gamma_1 \\ \gamma_2 \\ \vdots \\ \gamma_{p-1} \\ \gamma_p
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 0 \\ \vdots \\ 0 \\ 1
\end{bmatrix}
```

---

* **Case 1:** ``m+1 = p \;\;\Rightarrow \;\; `` square matrix
  Then, $\det A$ is a *Vandermonde Determinant*, which is guaranteed to be nonzero.
  That is, `` \det A \neq 0 \;\;\Rightarrow \;\; `` unique solution. 

* **Case 2:** $p < m+1$
  ``\mathrm{rank}(A;b) \neq \mathrm{rank}(A) \;\;\Rightarrow\;\; `` No solution

* **Case 3:** $p > m+1$
  ``\mathrm{rank}(A;b) = \mathrm{rank}(A) < p \;\;\Rightarrow\;\; `` Infinite number of solutions

---

**Note:** Review Appendix A of CG
"""


# ╔═╡ 8fe1d36b-ada6-434d-a3d9-5b4ba8bc20e4
md"""
✔ **Summary:**

- A consistent FDA for an $m^{th}$-order derivative using $p$ points requires $p \geq m+1$.
- If ``p = m+1``, approximation is generally ``\mathcal{O}(\Delta x)``.
- If ``p > m+1``, can obtain higher order approximations by forcing ``F_{m+1} = 0``, ``F_{m+2} = 0``, etc.

---

✔ **Example:**

- Derive an approximation for 
```math
\left.\frac{du}{dx}\right|_{x_i}
```

in terms of the values of $u$ at points $x_{i-1}, x_i, x_{i+1}$.

---

For $m=1, \; p=3 \;\;\Rightarrow\;\; p > m+1 \;\;\Rightarrow\;\;$ we can find a consistent FDA:

```math
\begin{bmatrix}
1 & 1 & 1 \\
-\nabla x_i & 0 & \Delta x_i
\end{bmatrix}
\begin{bmatrix}
\gamma_{i-1} \\ \gamma_i \\ \gamma_{i+1}
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 1
\end{bmatrix}
```

---

**Solution (assuming $\nabla x_i = \Delta x_i$):**

```math
\gamma_{i-1} = \gamma_{i+1} - \frac{1}{\Delta x}, 
\quad \gamma_i = \frac{1}{\Delta x} - 2\gamma_{i+1}, 
\quad \gamma_{i+1} = \text{arbitrary}.
```

"""


# ╔═╡ 4cfb5136-d4a0-463e-9540-4758b3558ddd
md"""
Let $\gamma_{i+1} = \dfrac{\alpha}{\Delta x}$ ($\alpha$ is arbitrary). Then

```math
\left.\dfrac{du}{dx}\right|_{x_i} \;\approx\; 
\dfrac{(\alpha - 1) u_{i-1} + (1 - 2\alpha) u_i + \alpha u_{i+1}}{\Delta x}
```

1. ``\alpha = 1, \alpha = 0, \alpha = \tfrac{1}{2} `` correspond to forward, backward, and central differences.

2. Impose ``F_2 = 0 \Rightarrow \mathcal{O}(\Delta x^2) \Rightarrow ``

```math
\dfrac{\Delta x^2}{2} \cdot \dfrac{\alpha - 1}{\Delta x} + \dfrac{\Delta x^2}{2} \cdot \dfrac{\alpha}{\Delta x} = 0
```
Thus, ``\alpha = \tfrac{1}{2} ``

"""


# ╔═╡ 178467a9-b18b-4b51-914a-a1b42e65c8e2
md"""
What if $\nabla x_i \neq \Delta x_i$?

```math
\left.\frac{du}{dx}\right|_{x_i} \approx 
\frac{\alpha-1}{\nabla x_i} u_{i-1} +
\left[\frac{1}{\nabla x_i} - \alpha \frac{\nabla x_i + \Delta x_i}{\nabla x_i \Delta x_i}\right] u_i +
\frac{\alpha}{\Delta x_i} u_{i+1}
```

---

* ``\alpha = 0  \Rightarrow \dfrac{u_i - u_{i-1}}{\nabla x_i},  O(\Delta x)``
* ``\alpha = 1  \Rightarrow \dfrac{u_{i+1} - u_i}{\Delta x_i}, O(\Delta x)``
* ``\alpha = \tfrac{1}{2} \Rightarrow -\dfrac{1}{2\nabla x_i} u_{i-1} +
  \left(\dfrac{1}{2\nabla x_i} - \dfrac{1}{2\Delta x_i}\right) u_i +
  \dfrac{1}{2\Delta x_i} u_{i+1}, O(\Delta x)`` (except when ``\nabla x_i = \Delta x_i``)

---

- Can we get higher–order ``O(\Delta x^2)``?

```math
F_2 = \frac{\alpha - 1}{\nabla x_i} \nabla x_i^2 + \frac{\alpha}{\Delta x_i} \Delta x_i^2
= (\alpha - 1)\nabla x_i + \alpha \Delta x_i = 0
```

``\Rightarrow \alpha = \dfrac{\nabla x_i}{\nabla x_i + \Delta x_i} ``

(in particular ``\alpha = \tfrac{1}{2}`` when ``\nabla x_i = \Delta x_i``).
"""


# ╔═╡ 0ac21cf3-5e93-46f6-b1b5-0e9a948d3d4c
md"""
## Some standard FD Approximations

```math
\left.\frac{du}{dx}\right|_{x_i} \approx
```

```math
\frac{u_{i+1} - u_i}{\Delta x} \quad O(\Delta x)
```

```math
\frac{u_i - u_{i-1}}{\Delta x} \quad O(\Delta x)
```

```math
\frac{u_{i+1} - u_{i-1}}{2\Delta x} \quad O((\Delta x)^2)
```

```math
\frac{-u_{i+2} + 4u_{i+1} - 3u_i}{2\Delta x} \quad O((\Delta x)^2)
```

```math
\frac{-u_{i+2} + 8u_{i+1} - 8u_{i-1} + u_{i-2}}{12 \Delta x} \quad O((\Delta x)^4)
```

---

```math
\left.\frac{d^2 u}{dx^2}\right|_{x_i} \approx
```

```math
\frac{u_{i+1} - 2u_i + u_{i-1}}{(\Delta x)^2} \quad O((\Delta x)^2)
```

```math
\frac{-u_{i+2} + 16u_{i+1} - 30u_i + 16u_{i-1} - u_{i-2}}{12(\Delta x)^2} \quad O((\Delta x)^4)
```

"""


# ╔═╡ 99f3e744-f4d6-4481-b0a6-4a573a7db57d
md"""
## Example calculation 1

Solve

```math
D \frac{d^2 u}{dx^2} - k u = 0, \quad 0 \leq x \leq 1
```

with boundary conditions

```math
u(0) = 0, \quad u(1) = C_1
```

"""

# ╔═╡ a82a2ee5-c358-4c9f-899a-4a357bb58e9b
md"""
### Choose discretization

"""

# ╔═╡ 3325c9cf-c42f-45af-b26b-c09cbef5d094
local img = LocalResource("./figs/mod3_4_node_dm.png", :width => "400px")

# ╔═╡ 20b03877-84ac-4a85-810b-e5eb04a39bc5
md"""

```math
\text{4 nodes: 2 at boundaries, 2 interior}
```

```math
\Delta x = \tfrac{1}{3} \quad (\text{constant})
```

"""

# ╔═╡ 9da62fbc-742c-4801-90f3-e475d7e3fada
md"""

### Write algebraic equations

- At the boundary nodes, using boundary information:

```math
u_1 = 0, \quad u_4 = C_1
```

---

- At the interior nodes, write FDA equation:

```math
\left(D \frac{d^2 u}{dx^2} - k u \right)\bigg|_{x_2} = 0
\;\;\Rightarrow\;\;
D \left[\frac{u_1 - 2u_2 + u_3}{\Delta x^2} + O(\Delta x^2)\right] - k u_2 = 0
```

Neglecting the truncation error $O(\Delta x^2)$:

```math
D \frac{u_1 - 2u_2 + u_3}{\Delta x^2} - k u_2 = 0
```

So at **Node #2**:

```math
u_1 - \left(2 + \frac{k \Delta x^2}{D}\right) u_2 + u_3 = 0
```

Note that the coefficients are now *dimensionless* after the rearrangement.

---

Similarly, at **Node #3**:

```math
\left(D \frac{d^2 u}{dx^2} - k u \right)\bigg|_{x_3} = 0
\;\;\Rightarrow\;\;
u_2 - \left(2 + \frac{k \Delta x^2}{D}\right) u_3 + u_4 = 0
```

"""


# ╔═╡ 79430c03-b74e-4cff-90e1-31449e046db7
md"""

- Form matrix equation and solve

```math
\begin{bmatrix}
1 & 0 & 0 & 0 \\
1 & -(2 + \tfrac{k \Delta x^2}{D}) & 1 & 0 \\
0 & 1 & -(2 + \tfrac{k \Delta x^2}{D}) & 1 \\
0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
u_1 \\ u_2 \\ u_3 \\ u_4
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 0 \\ 0 \\ C_1
\end{bmatrix}
```

---

Or equivalently, focusing only on the interior nodes:

```math
\begin{bmatrix}
-(2 + \tfrac{k \Delta x^2}{D}) & 1 \\
1 & -(2 + \tfrac{k \Delta x^2}{D})
\end{bmatrix}
\begin{bmatrix}
u_2 \\ u_3
\end{bmatrix}
=
\begin{bmatrix}
0 \\ -C_1
\end{bmatrix}

```
"""

# ╔═╡ 5d7665df-c17e-44c0-9756-c53c2e80779d
md"""

``𝐷=0.01, 𝑘=0.1, 𝐶_1=1, Δ𝑥=1/3 \Rightarrow U_2=0.094, 𝑈_3=0.291``

```math
\hat u=(0, 0.094, 0.291, 1)
```

```math
u_{\text{true}}= (0, 0.083, 0.320, 1)
```


"""

# ╔═╡ d01b8353-3a58-47f1-88c6-05c072a14238
begin
	gr()

	# Parameters
	D  = 0.01
	k  = 0.1
	C1 = 1.0
	Δx = 1/3

	# Node locations
	x = [0.0, 1/3, 2/3, 1.0]

	# Discrete FD solution 
	u = [0.0, 0.094, 0.291, 1.0]

	# “True” (analytic) solution: u(x) = C1*sinh(λ x)/sinh(λ),  λ = √(k/D)
	λ = sqrt(k/D)
	u_true = C1 .* sinh.(λ .* x) ./ sinh(λ)

	# Plot
	plt = plot(x, u_true; lw=1.5, label=L"u_{\mathrm{true}}", xlim=(0,1), ylim=(0,1),
	           	xlabel="x", 
			   	legend=:topleft, 
			   	legendfontsize = 14, 
			   	guidefontsize = 16,   # axis labels
        	   	tickfontsize = 12,
			  	size = (400,300), )
	plot!(x, u; lw=1.5, label=L"\hat u")
	
	plt
end

# ╔═╡ 91a049da-181b-4995-adaa-bb944bf67383
md"""
## Example calculation 2

Solve

```math
V \frac{du}{dx} - D \frac{d^2 u}{dx^2} = 0, \quad 0 \leq x \leq L
```

with boundary conditions

```math
u(0) = 1, \quad u(L) = 0
```

---

### Choose discretization

"""

# ╔═╡ ca1a0c3a-c666-4534-abd2-6e578d4e660e
local img = LocalResource("./figs/mod3_n_node_dm.png", :width => "400px")

# ╔═╡ a9c3ec07-92b4-4f12-bf4e-155ae5f9ad28
md"""
### Central difference scheme for the first-order derivative

```math
\left.\frac{du}{dx}\right|_{x_i} = \frac{u_{i+1} - u_{i-1}}{2 \Delta x}
```
Only one option for the second derivative using a three-point stencil.

```math
\left.\frac{d^2 u}{dx^2}\right|_{x_i} = \frac{u_{i-1} - 2u_i + u_{i+1}}{(\Delta x)^2}
```

"""

# ╔═╡ e6cb7e94-6f0b-480c-a44a-2fedc062289c
md"""
### Write algebraic equations

```math
V \frac{u_{i+1} - u_{i-1}}{2 \Delta x}
- D \frac{u_{i-1} - 2u_i + u_{i+1}}{(\Delta x)^2} = 0
```

Rearrange:

```math
\left(-\frac{V \Delta x}{2D} - 1\right) u_{i-1}
+ 2u_i
+ \left(\frac{V \Delta x}{2D} - 1\right) u_{i+1} = 0
```

Note that the coefficients are *dimensionless* after the rearrangement.

---

Define grid Peclet number (*dimensionless*):

```math
\frac{V \Delta x}{D}
```

**Q:** Does Grid Peclet number matter for the solutions?

---

- If ``\frac{V \Delta x}{D} > 2`` → Oscillation!!!
- If ``\frac{V \Delta x}{D} = 2`` → Marginal stability
- If ``\frac{V \Delta x}{D} < 2`` → Stable solution

"""

# ╔═╡ 476097f8-3792-4d4a-8b23-571ad139b4e0
md"""
Let's illustrate the dependance of the solution behavior on the grid Peclet number assuming we have 3 nodes (``x_1=0``, ``x_2=L/2``, ``x_3=L``) for the entire domain.
"""

# ╔═╡ 73d4b8b5-db48-4e70-8a04-2b9da72f3215
local img = LocalResource("./figs/mod3_centralFD_oscillation.png", :width => "500px")

# ╔═╡ e5c1f03b-b604-478b-a771-8828b1f4e6d4
md"""
### Upstream-weighted difference scheme for the first-order derivative

```math
\left.\frac{du}{dx}\right|_{x_i} = \frac{u_{i} - u_{i-1}}{\Delta x}
```

Only one option for the second derivative using a three-point stencil.
```math
\left.\frac{d^2 u}{dx^2}\right|_{x_i} = \frac{u_{i-1} - 2u_i + u_{i+1}}{(\Delta x)^2}
```

---

### Write algebraic equations

```math
V \frac{u_i - u_{i-1}}{\Delta x}
- D \frac{u_{i-1} - 2u_i + u_{i+1}}{(\Delta x)^2} = 0
```

Rearrange:

```math
\left(-\frac{V \Delta x}{D} - 1\right) u_{i-1}
+ \left(\frac{V \Delta x}{D} + 2\right) u_i
- u_{i+1} = 0
```

---

**No oscillation for any** ``\frac{V \Delta x}{D}``!!!

"""


# ╔═╡ 4e3aa823-eabb-4a5d-be89-39d4269d1a4b
md"""
We can again illustrate the dependance of the solution behavior on the grid Peclet number assuming we have 3 nodes (``x_1=0``, ``x_2=L/2``, ``x_3=L``) for the entire domain.
"""

# ╔═╡ 22216780-22eb-4a92-b06d-0878d280a9d5
local img = LocalResource("./figs/mod3_forwardFD_no_oscillation.png", :width => "200px")

# ╔═╡ a3d0fb8d-61a5-492b-ad90-7223ceb73466
md"""
System of algebraic equations for ``N`` points

```math
\begin{bmatrix}
1 & 0 & 0 & \cdots & 0 \\
-\tfrac{V \Delta x}{D} - 1 & \tfrac{V \Delta x}{D} + 2 & -1 & \cdots & 0 \\
0 & \ddots & \ddots & \ddots & 0 \\
0 & \cdots & -\tfrac{V \Delta x}{D} - 1 & \tfrac{V \Delta x}{D} + 2 & -1 \\
0 & \cdots & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N
\end{bmatrix}
=
\begin{bmatrix}
1 \\ 0 \\ \vdots \\ 0 \\ 0
\end{bmatrix}
```

---

**Note:**

1. Main diagonals are greater than off-diagonals.
2. This is **not** the case if central difference is used for the advection term.

"""

# ╔═╡ a5e60b99-e675-402c-95a3-087e46a1611c
md"""
### Variably upstream-weighted approximation

```math
\left.\frac{du}{dx}\right|_{x_i}
= \alpha \frac{u_i - u_{i-1}}{\Delta x}
+ (1-\alpha)\frac{u_{i+1} - u_{i-1}}{2\Delta x}
```

where ``\alpha`` is an arbitrary parameter.

---

* Can you choose a value of ``\alpha`` so that the FDA to the steady-state advection–diffusion equation gives exact results?

**Hint:**

1. Taylor expansion
2. ``\left.\frac{d^m u}{dx^m}\right|_{x_i} = \left(\frac{V}{D}\right)^{m-1} \left.\frac{du}{dx}\right|_{x_i} ``

---

**Take derivative for both sides:**

```math
\frac{du}{dx} = \frac{V}{D} \frac{d^2 u}{dx^2}
\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\Rightarrow
\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\frac{d^2 u}{dx^2} = \frac{V}{D} \frac{d^3 u}{dx^3}
```

Keep taking the derivatives and you will quickly see that ``\left.\frac{d^m u}{dx^m}\right|_{x_i} = \left(\frac{V}{D}\right)^{m-1} \left.\frac{du}{dx}\right|_{x_i} ``.

Then, apply Taylor expansion to the FDA and solve for ``\alpha`` by taking advantage of the above identity.


"""


# ╔═╡ e8662f52-67bf-455e-9c8d-ef026f349595
md"""
### Evaluate solution behaviors (central vs. upstream-weighted)

Now, we evaluate the solutions using central difference vs. upstream-weighted difference for the advection term for different grid Peclet numbers. 

Set ``V=100, D=1, L=1``. Then, adjust ``\Delta x`` to achieve different grid Peclet numbers (``V \Delta x / D``).
"""

# ╔═╡ a8a76810-57c6-4fdc-a006-cae756dd7888
begin
    # --------- problem setup ----------
    local V  = 100.0          # advection speed
    local D  = 1.0            # diffusivity
    local L  = 1.0            # domain length
    local Pe_list = [0.5, 1.0, 2.0, 4.0, 8.0]   # grid Peclet numbers to test

    # Analytic solution: u(x) = (exp(V*L/D) - exp(V*x/D)) / (exp(V*L/D) - 1)
    E = exp(V*L/D)
    local u_true = x -> (E - exp(V*x/D)) / (E - 1)

    # Build & solve for a given alpha and grid Peclet number
    function solve_advecdiff(alpha, Pe_g; V=V, D=D, L=L)
        Δx_target = Pe_g * D / V
        # pick an integer grid size close to target and recompute Δx
        N = max(3, round(Int, 1 + L/Δx_target))
        Δx = L/(N-1)
        # coefficients for interior nodes (constant on uniform grid)
        a = V/Δx * (-(alpha + 1)/2) - D/Δx^2      # coeff of u_{i-1}
        b = V/Δx * alpha            + 2D/Δx^2      # coeff of u_i
        c = V/Δx * ((1 - alpha)/2)  - D/Δx^2      # coeff of u_{i+1}

        # Assemble tridiagonal system for i = 2..N-1
        main  = fill(b, N-2)
        lower = fill(a, N-3)
        upper = fill(c, N-3)
        A = Tridiagonal(lower, main, upper)

        # RHS (Dirichlet BCs: u1=1, uN=0)
        rhs = zeros(N-2)
        rhs[1]   -= a * 1.0       # u1
        rhs[end] -= c * 0.0       # uN

        u = zeros(N)
        u[1] = 1.0;  u[end] = 0.0
        u[2:end-1] = A \ rhs
        x = range(0.0, L, length=N)
        return x, u, Δx, N
    end

	local plt = plot(; size=(900,600),
	legend=:bottomleft, guidefontsize=14, tickfontsize=12, legendfontsize=12,
	xlabel=L"x", ylabel=L"u(x)")
	
	# analytic curve
	xs_true = range(0, L; length=400)
	plot!(plt, xs_true, u_true.(xs_true); lw=3, label=L"u_{\mathrm{true}}", color=:black)
	
	# generate a color palette with length(Pe_list) distinct colors
	colors = palette(:tab10, length(Pe_list))
	
	for (i, Pe) in enumerate(Pe_list)
	    col = colors[i]
	
	    # Upwind (alpha = 1) → solid line
	    xU, uU, ΔxU, _ = solve_advecdiff(1.0, Pe)
	    plot!(plt, xU, uU; lw=2, label="upwind, Pe=$(round(Pe,digits=2))", color=col, ls=:solid)
	
	    # Central (alpha = 0) → dashed line, same color
	    xC, uC, ΔxC, _ = solve_advecdiff(0.0, Pe)
	    plot!(plt, xC, uC; lw=2, label="central, Pe=$(round(Pe,digits=2))", color=col, ls=:dash)
	end
	
	# make legend box lighter but keep axes frame
	background_color_legend = :white
	foreground_color_legend = RGBA(0,0,0,0)
	
	plt
end


# ╔═╡ 0a29e6dd-7de8-4cf5-b555-5c12bad2debb
md"""
## Modified Equation Analysis

- How do the different discretization schemes change the behavior of the equation fundamentally?

    - Numerical scheme should not overly distort the critical physics that we are trying to model.

- How to analyze? ⇒ **Modified Equation Analysis**

- *Idea*: Using Taylor expansion to find another differential equation that is better approximated by our discretization scheme.

- We’ll discuss more about this when we look at transient problems.

"""

# ╔═╡ 83aa0392-70db-434c-acd3-28738e3217dd
md"""
### Upstream-weighted (or Upwind) Approximation

--- 

Taylor expansions:

```math
u_{i+1} = u_i 
+ \left.\frac{du}{dx}\right|_{x_i} \Delta x 
+ \frac{1}{2}\left.\frac{d^2 u}{dx^2}\right|_{x_i}\Delta x^2 
+ \frac{1}{6}\left.\frac{d^3 u}{dx^3}\right|_{x_i}\Delta x^3 
+ \frac{1}{24}\left.\frac{d^4 u}{dx^4}\right|_{x_i}\Delta x^4 + \cdots
```

```math
u_{i-1} = u_i 
- \left.\frac{du}{dx}\right|_{x_i} \Delta x 
+ \frac{1}{2}\left.\frac{d^2 u}{dx^2}\right|_{x_i}\Delta x^2 
- \frac{1}{6}\left.\frac{d^3 u}{dx^3}\right|_{x_i}\Delta x^3 
+ \frac{1}{24}\left.\frac{d^4 u}{dx^4}\right|_{x_i}\Delta x^4 - \cdots
```

---

```math
\frac{u_i - u_{i-1}}{\Delta x} 
= \left.\frac{du}{dx}\right|_{x_i} 
- \frac{\Delta x}{2} \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ O(\Delta x^2)
```

```math
\frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ \frac{1}{12}\left.\frac{d^4 u}{dx^4}\right|_{x_i}\Delta x^2 + \cdots
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} + O(\Delta x^2)
```

---

```math
\Rightarrow 
V \frac{u_i - u_{i-1}}{\Delta x} 
- D \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
= V \left.\frac{du}{dx}\right|_{x_i} 
- \frac{V \Delta x}{2} \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
- D \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ O(\Delta x^2)
```

- Upstream-weighted approximation modifies the original equation with an additional *diffusion* term 

```math
- \frac{V \Delta x}{2} \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
```

This additional term is caused by the numerical error (scales with $\Delta x$) and is often referred to as **numerical diffusion**. 

"""

# ╔═╡ 3b3f1e9d-ebba-4d69-9f83-c6500708c71c
md"""
### Central Difference Approximation

```math
\frac{u_{i+1} - u_{i-1}}{2 \Delta x} 
= \left.\frac{du}{dx}\right|_{x_i} 
+ \frac{1}{6}\left.\frac{d^3 u}{dx^3}\right|_{x_i} \Delta x^2 
+ O(\Delta x^4)
```

```math
\frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ \frac{1}{12}\left.\frac{d^4 u}{dx^4}\right|_{x_i} \Delta x^2 
+ \frac{1}{360}\left.\frac{d^6 u}{dx^6}\right|_{x_i} \Delta x^4 
+ \cdots
```

```math
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ \frac{1}{12}\left.\frac{d^4 u}{dx^4}\right|_{x_i} \Delta x^2 
+ O(\Delta x^4)
```

---

```math
\Rightarrow \;
V \frac{u_{i} - u_{i-1}}{\Delta x} 
- D \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
```

```math
= V \left[ \left.\frac{du}{dx}\right|_{x_i} 
+ \frac{1}{6}\left.\frac{d^3 u}{dx^3}\right|_{x_i} \Delta x^2 \right] 
```

```math
- D \left[ \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ \frac{1}{12}\left.\frac{d^4 u}{dx^4}\right|_{x_i} \Delta x^2 \right] 
+ O(\Delta x^3)
```

Central-difference approximation modifies the original equation with two additional terms

```math
\frac{1}{6}\left.\frac{d^3 u}{dx^3}\right|_{x_i} \Delta x^2  - \frac{1}{12}\left.\frac{d^4 u}{dx^4}\right|_{x_i} \Delta x^2
```

- The first third-order term shifts the phase speed of different Fourier modes. As a result, high-frequency (short-wavelength) waves travel at the wrong speed and introduces *oscillation*. Often referred to as **numerical dispersion**.

- The second fourth-order term damps out high-frequency oscillations more strongly than low-frequency ones (**numerical diffusion**). As a result, the scheme artificially smooths sharp gradients and oscillations, but less aggressively than the numerical diffusion term in an upstream-weighted scheme.

- The third-order term is more dominant than the fourth-order term, which leads to an overall oscillatory behavior (**numerical dispersion**).

**NOTE**: The concept of *numerical dispersion* here is different from the *macroscale dispersion* in groundwater solute transport.

"""



# ╔═╡ 6c004a80-5ada-45e9-8618-9995032244e5
md"""
### Downstream-weighted (or Downwind) Approximation

```math
\frac{u_{i+1} - u_i}{\Delta x} 
= \left.\frac{du}{dx}\right|_{x_i} 
+ \frac{\Delta x}{2} \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ O(\Delta x^2)
```

```math
\frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ \frac{1}{24}\left.\frac{d^4 u}{dx^4}\right|_{x_i}\Delta x^2 
+ \cdots
= \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ O(\Delta x^2)
```

---

```math
\Rightarrow \;
V \frac{u_{i+1} - u_i}{\Delta x} 
- D \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} 
```

```math
= V \left.\frac{du}{dx}\right|_{x_i} 
+ \frac{V \Delta x}{2} \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
- D \left.\frac{d^2 u}{dx^2}\right|_{x_i} 
+ O(\Delta x^2)
```

- The term ``\frac{V \Delta x}{2} \frac{d^2 u}{dx^2}`` acts like negative diffusion (anti-diffusion).

- If the effective diffusion coefficient ``(D - \frac{V \Delta x}{2})`` becomes negative, the solution has strong oscillations.

"""


# ╔═╡ dc3ad891-1e72-4ac7-9bc9-d7fe9c483d7b
md"""
## Boundary Conditions
"""

# ╔═╡ d5ac488d-247a-47aa-9600-6d11aa645ae9
md"""
### Dirichlet boundary condition
"""

# ╔═╡ 50b18840-ffcc-4132-9f22-a59d1d5d7425
local img = LocalResource("./figs/mod3_dirichlet_BC.png", :width => "400px")

# ╔═╡ 1dd43ab7-05c8-4773-91c2-14cea1fe133b
md"""
### Neumann boundary contion
"""

# ╔═╡ 6037e76b-bd82-4e18-825d-0a197d723a74
local img = LocalResource("./figs/mod3_neumann_BC.png", :width => "400px")

# ╔═╡ 85eaf88c-cfbc-4408-9c9c-0758790bc76d
md"""

```math
\left.\frac{du}{dx}\right|_{x_1} = C_1
```

- ``u_{-1}`` is a ghost node.
- ``u_1`` is unknown, we need to solve for it.

---

**Finite difference approximations at $x_1$:**

* Forward difference:

```math
\left.\frac{du}{dx}\right|_{x_1} = \frac{u_2 - u_1}{\Delta x}
```

* Backward difference:

```math
\left.\frac{du}{dx}\right|_{x_1} = \frac{u_1 - u_{-1}}{\Delta x}
```

* Central difference:

```math
\left.\frac{du}{dx}\right|_{x_1} = \frac{u_2 - u_{-1}}{2 \Delta x}
```

---

Often, we employ the following assumption

```math
u_1 = \frac{u_{-1} + u_2}{2}
```

"""


# ╔═╡ 5aaa13ab-e38c-4d3d-9699-c87fee236b07
md"""
## Multiple Dimensions
"""

# ╔═╡ f096d749-5ffe-4daa-b763-c6f99a8d73bf
md"""
- Rectangular grid
"""

# ╔═╡ 90f92c1e-746e-44d6-9298-a027ff95a1d8
local img = LocalResource("./figs/mod3_rectangle_grid.png", :width => "300px")

# ╔═╡ b786c0f4-65b1-4e1d-ad68-a6fa5adc2b3f
md"""
**One-dimensional Derivatives**

- Nothing new; add the index of the other direction.  
- Example (Second-order central difference):

```math
\frac{\partial u_{i,j}}{\partial x} 
= \frac{u_{i+1,j} - u_{i-1,j}}{2 \Delta x} 
- \frac{\Delta x^2}{6} \frac{\partial^3 u_{i,j}}{\partial x^3} 
+ \cdots
```

"""

# ╔═╡ b1e27e5d-c876-43ca-8288-5b5b309221c8
md"""
**Mixed derivatives**

- Use 1D derivative approximations in succession  
- Example:
```math
\frac{\partial^2 u}{\partial x \partial y}
```

* Note that

```math
\frac{\partial^2 u}{\partial x \partial y}
= \frac{\partial}{\partial x} \left( \frac{\partial u}{\partial y} \right)
= \frac{\partial}{\partial y} \left( \frac{\partial u}{\partial x} \right)
```

---

Let ``f = \frac{\partial u}{\partial y}``, then ``\frac{\partial^2 u}{\partial x \partial y} = \frac{\partial f}{\partial x}``

Using central difference for $f$:

```math
\frac{\partial f}{\partial x} = \frac{f_{i+1} - f_{i-1}}{2\Delta x}
- \frac{1}{6} f_{xxx}|_i \, \Delta x^2 + O(\Delta x^4)
```

where

```math
f_{i+1} = \left.\frac{\partial u}{\partial y}\right|_{i+1}
= \frac{u_{i+1,j+1} - u_{i+1,j-1}}{2\Delta y}
- \frac{1}{6} u_{yyy}|_{i+1}\, \Delta y^2 + O(\Delta y^4)
```

and

```math
f_{i-1} = \left.\frac{\partial u}{\partial y}\right|_{i-1}
= \frac{u_{i-1,j+1} - u_{i-1,j-1}}{2\Delta y}
- \frac{1}{6} u_{yyy}|_{i-1}\, \Delta y^2 + O(\Delta y^4)
```

---

Final approximation:

```math
\frac{\partial^2 u}{\partial x \partial y}
= \frac{u_{i+1,j+1} - u_{i+1,j-1} - u_{i-1,j+1} + u_{i-1,j-1}}{4 \Delta x \Delta y}
```

```math
- \frac{\Delta x^2}{6} u^{xxxy}_{i,j}
- \frac{\Delta y^2}{6} u^{xyyy}_{i,j}
+ O(\Delta x^4, \Delta y^4)
```

*(second-order accurate)*

"""

# ╔═╡ df3e1dfe-3347-46b6-9a62-f850ab0f5c68
local img = LocalResource("./figs/mod3_poisson_domain.png", :width => "400px")

# ╔═╡ 74627804-0ea5-44f6-8919-8992d41ca6a1
md"""
### Example: Poisson Equation

We want to solve:
```math
\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = f(x,u)
```

---

Approximations:

```math
\left.\frac{\partial^2 u}{\partial x^2}\right|_{x_i,y_i}
\approx \frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{\Delta x^2}
```

```math
\left.\frac{\partial^2 u}{\partial y^2}\right|_{x_i,y_i}
\approx \frac{u_{i,j+1} - 2u_{i,j} + u_{i,j-1}}{\Delta y^2}
```

---

Combining:

```math
\left[
\frac{\partial^2 u}{\partial x^2}
+ \frac{\partial^2 u}{\partial y^2}
- f(x,y)
\right]_{x_i,y_i}
\approx
\frac{1}{\Delta x^2} u_{i+1,j}
+ \frac{1}{\Delta x^2} u_{i-1,j}
+ \frac{1}{\Delta y^2} u_{i,j+1}
```

```math
+ \frac{1}{\Delta y^2} u_{i,j-1}
- \left( \frac{2}{\Delta x^2} + \frac{2}{\Delta y^2} \right) u_{i,j}
- f|_{x_i,y_i}
```

---

**Procedure:**

1. Number all the nodes.
2. Write the algebraic equation for each node.
3. Put together the linear matrix system.

**Notes:**

* The matrix is *5-diagonal*.
* Need to take care of *boundary conditions (BCs)*.

"""

# ╔═╡ 95ef9435-2502-41df-b14d-0cde79f1dee5
md"""
### Variable coefficients in the differential equation

- Constant coefficients:

```math
V \frac{du}{dx} - D \frac{d^2 u}{dx^2} = 0
```

- Variable coefficients:

```math
\frac{d(Vu)}{dx} - \frac{d}{dx} \left( D \frac{du}{dx} \right) = 0
```

---

- Advection term:

```math
\left. \frac{d(Vu)}{dx} \right|_i
= \frac{(Vu)_i - (Vu)_{i-1}}{\Delta x}
\quad \text{(upstream-weighted)}
```

```math
\left. \frac{d(Vu)}{dx} \right|_i
= \frac{(Vu)_{i+1} - (Vu)_{i-1}}{2\Delta x}
\quad \text{(central)}
```

- Diffusion term:

```math
\frac{d}{dx} \left( D \frac{du}{dx} \right)
= \frac{ \left(D \frac{du}{dx}\right)_{i+1/2}
- \left(D \frac{du}{dx}\right)_{i-1/2} }{\Delta x}
```

Expanding:

```math
= \frac{ D_{i+1/2} \frac{u_{i+1}-u_i}{\Delta x}
- D_{i-1/2} \frac{u_i-u_{i-1}}{\Delta x} }{\Delta x}
```

- Harmonic average:

```math
D_{i+1/2} = \frac{2}{ \tfrac{1}{D_i} + \tfrac{1}{D_{i+1}} }
```

- Final discretized form:

```math
\frac{ D_{i+1/2} u_{i+1} - (D_{i+1/2} + D_{i-1/2}) u_i + D_{i-1/2} u_{i-1} }{\Delta x^2}
```

"""


# ╔═╡ 345d7010-562f-4e42-a683-010ba03e74c5
md"""
## Point-wise grids vs. cell-centered grids for FDA
"""

# ╔═╡ 792a676e-a76e-4a69-8ae8-93f9119a3e51
md"""
**Point-wise grids**
"""

# ╔═╡ 8213a50b-ca98-456a-bd0a-cf592fd59d1e
local img = LocalResource("./figs/mod3_pointwise_1d.png", :width => "400px")

# ╔═╡ 35c8ab9e-43de-4e4d-94e8-9cd64ba3bce2
md"""
**Cell-centered grids**
"""

# ╔═╡ 2f025d7e-7598-4caa-8e0a-d265cadebdeb
local img = LocalResource("./figs/mod3_cell_centered_1d.png", :width => "400px")

# ╔═╡ ec58ba5d-76fd-4122-bbd8-32d054f6a808
md"""

* Each control volume corresponds to a **grid cell**
* Flux in at left face, flux out at right face

---

**Advection term:**

```math
\left. \frac{d(Vu)}{dx} \right|_i
= \frac{(Vu)_{i+1/2} - (Vu)_{i-1/2}}{\Delta x}
```

Approximation (central difference):

```math
\frac{(Vu)_i + (Vu)_{i+1}}{2\Delta x}
- \frac{(Vu)_{i-1} + (Vu)_i}{2\Delta x}
= \frac{(Vu)_{i+1} - (Vu)_{i-1}}{2\Delta x}
```

Approximation (upstream-weighted):

```math
\left. \frac{d(Vu)}{dx} \right|_i
= \frac{(Vu)_{i+1/2} - (Vu)_{i-1/2}}{\Delta x}
\approx \frac{(Vu)_i - (Vu)_{i-1}}{\Delta x}
```

---

**Notes:**

* Central difference: higher accuracy but can introduce oscillations.
* Upstream-weighted: adds numerical diffusion, but no oscillations.

"""


# ╔═╡ 3448eb07-97b0-48e5-ba5a-5c819e32c1b8
local img = LocalResource("./figs/mod3_cell_centered_1d_ghostcell.png", :width => "400px")

# ╔═╡ 7d0fcb65-fb46-4881-9b58-ba3992382745
md"""

**Boundary Conditions (BCs):**

- Dirichlet:
```math
\frac{u_1 + u_{-1}}{2} = u \bigg|_{x=0}
```

- Neumann:

```math
\frac{u_1 - u_{-1}}{2} = \left.\frac{du}{dx}\right|_{x=0}
```

"""

# ╔═╡ b683a35a-55de-4223-8741-1895facae7fe
md"""

**Diffusion term:**

```math
\left. \frac{d}{dx} \left( D \frac{du}{dx} \right) \right|_i
= \frac{ \left( D \frac{du}{dx} \right)_{i+1/2}
- \left( D \frac{du}{dx} \right)_{i-1/2} }{\Delta x}
```

Expanding:

```math
= \frac{ D_{i+1/2} \frac{u_{i+1} - u_i}{\Delta x}
- D_{i-1/2} \frac{u_i - u_{i-1}}{\Delta x} }{\Delta x}
```

Final form:

```math
= \frac{ D_{i+1/2} u_{i+1}
- (D_{i+1/2} + D_{i-1/2}) u_i
+ D_{i-1/2} u_{i-1} }{\Delta x^2}
```

* Here $D_{i+1/2}$ and $D_{i-1/2}$ are evaluated using the *harmonic mean*:

```math
D_{i+1/2} = \frac{2}{ \tfrac{1}{D_i} + \tfrac{1}{D_{i+1}} }
```

"""


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
LaTeXStrings = "~1.4.0"
Plots = "~1.40.18"
PlutoUI = "~0.7.69"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "a7ab5af3a6af871adac4f2e2b770be5270a576fa"

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
# ╟─6b928336-5bad-48f7-b480-1bf337e1607f
# ╟─1816242a-8054-11f0-152a-1d9edcfd46bc
# ╟─02a95c05-026a-47fb-a325-f00eff3e0dd6
# ╟─f2cf33d1-948f-4d68-b52b-c9b750accbd7
# ╟─c1be78ce-9c9a-4f92-b12e-296679f4de8a
# ╟─cc63f48e-9e72-4b55-b3f2-414ce3bae01e
# ╟─eb5ffd0f-4f84-4aa4-ad33-de9178c2218b
# ╟─b5933bf7-3f82-4e5f-95a7-829c0c0715c6
# ╟─a9ed7ccd-4a73-4d04-8d1c-56c795c59f54
# ╟─702a4721-c5ee-40ae-bbf6-0eae201e9d61
# ╟─45fbb6d6-dd13-4e1c-a2c8-702676aa20a5
# ╟─bd6b2079-37e9-45c2-a451-adc504619922
# ╟─3c32791a-dee1-4873-bcf4-73abde73c57d
# ╟─ca7c7180-a026-40f4-b8bb-4bb9a8a3df1d
# ╟─6bafe47c-316c-4a7a-96de-1842ef0f3a0f
# ╟─a13d06f1-7383-4630-a4a5-73c9a54f92ed
# ╟─4007ed66-f841-4c45-9a16-2c9d255529b8
# ╟─7243164e-8b07-4a60-aa74-3a0e650f6229
# ╟─7c34e980-11cb-4744-98c1-4a01f23a9838
# ╟─c1ed2624-2c4b-436f-b0b9-4089c802fcbb
# ╟─64687840-fa50-4282-b25c-a223f5ff4cd5
# ╟─e923d22e-01dd-4469-bb92-ff0a8721824b
# ╟─8fe1d36b-ada6-434d-a3d9-5b4ba8bc20e4
# ╟─4cfb5136-d4a0-463e-9540-4758b3558ddd
# ╟─178467a9-b18b-4b51-914a-a1b42e65c8e2
# ╟─0ac21cf3-5e93-46f6-b1b5-0e9a948d3d4c
# ╟─99f3e744-f4d6-4481-b0a6-4a573a7db57d
# ╟─a82a2ee5-c358-4c9f-899a-4a357bb58e9b
# ╟─3325c9cf-c42f-45af-b26b-c09cbef5d094
# ╟─20b03877-84ac-4a85-810b-e5eb04a39bc5
# ╟─9da62fbc-742c-4801-90f3-e475d7e3fada
# ╟─79430c03-b74e-4cff-90e1-31449e046db7
# ╟─5d7665df-c17e-44c0-9756-c53c2e80779d
# ╟─d01b8353-3a58-47f1-88c6-05c072a14238
# ╟─91a049da-181b-4995-adaa-bb944bf67383
# ╟─ca1a0c3a-c666-4534-abd2-6e578d4e660e
# ╟─a9c3ec07-92b4-4f12-bf4e-155ae5f9ad28
# ╟─e6cb7e94-6f0b-480c-a44a-2fedc062289c
# ╟─476097f8-3792-4d4a-8b23-571ad139b4e0
# ╟─73d4b8b5-db48-4e70-8a04-2b9da72f3215
# ╟─e5c1f03b-b604-478b-a771-8828b1f4e6d4
# ╟─4e3aa823-eabb-4a5d-be89-39d4269d1a4b
# ╟─22216780-22eb-4a92-b06d-0878d280a9d5
# ╟─a3d0fb8d-61a5-492b-ad90-7223ceb73466
# ╟─a5e60b99-e675-402c-95a3-087e46a1611c
# ╟─e8662f52-67bf-455e-9c8d-ef026f349595
# ╟─a8a76810-57c6-4fdc-a006-cae756dd7888
# ╟─0a29e6dd-7de8-4cf5-b555-5c12bad2debb
# ╟─83aa0392-70db-434c-acd3-28738e3217dd
# ╟─3b3f1e9d-ebba-4d69-9f83-c6500708c71c
# ╟─6c004a80-5ada-45e9-8618-9995032244e5
# ╟─dc3ad891-1e72-4ac7-9bc9-d7fe9c483d7b
# ╟─d5ac488d-247a-47aa-9600-6d11aa645ae9
# ╟─50b18840-ffcc-4132-9f22-a59d1d5d7425
# ╟─1dd43ab7-05c8-4773-91c2-14cea1fe133b
# ╟─6037e76b-bd82-4e18-825d-0a197d723a74
# ╟─85eaf88c-cfbc-4408-9c9c-0758790bc76d
# ╟─5aaa13ab-e38c-4d3d-9699-c87fee236b07
# ╟─f096d749-5ffe-4daa-b763-c6f99a8d73bf
# ╟─90f92c1e-746e-44d6-9298-a027ff95a1d8
# ╟─b786c0f4-65b1-4e1d-ad68-a6fa5adc2b3f
# ╟─b1e27e5d-c876-43ca-8288-5b5b309221c8
# ╟─df3e1dfe-3347-46b6-9a62-f850ab0f5c68
# ╟─74627804-0ea5-44f6-8919-8992d41ca6a1
# ╟─95ef9435-2502-41df-b14d-0cde79f1dee5
# ╟─345d7010-562f-4e42-a683-010ba03e74c5
# ╟─792a676e-a76e-4a69-8ae8-93f9119a3e51
# ╟─8213a50b-ca98-456a-bd0a-cf592fd59d1e
# ╟─35c8ab9e-43de-4e4d-94e8-9cd64ba3bce2
# ╟─2f025d7e-7598-4caa-8e0a-d265cadebdeb
# ╟─ec58ba5d-76fd-4122-bbd8-32d054f6a808
# ╟─3448eb07-97b0-48e5-ba5a-5c819e32c1b8
# ╟─7d0fcb65-fb46-4881-9b58-ba3992382745
# ╟─b683a35a-55de-4223-8741-1895facae7fe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
