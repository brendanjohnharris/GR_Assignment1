### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 09bee220-1249-11ec-3e2a-4fdde9976a59
begin
    import Pkg
    Pkg.activate(mktempdir())
	Pkg.add(url="https://github.com/brendanjohnharris/NonstationaryProcesses.jl")
    Pkg.add("DifferentialEquations")
	using NonstationaryProcesses, DifferentialEquations, Plots
	plotlyjs(); fourseas!(); # Plots initialisation
end

# ╔═╡ 81b414c9-2b23-4df4-9eaf-5549c83df269
md"# [Click here to view in a web browser](https://bl.ocks.org/brendanjohnharris/raw/4b08305f1b377d81d6c281a7d31c6227/?raw=true)"

# ╔═╡ 1cbf17f7-8322-4e3c-a895-163a48332205
md"# Question 3"

# ╔═╡ 08ca86b3-2c97-4486-bb4d-5c00d8239119
md"""
The equations of motion in the Schwarzschild metric (question 2) are:

$$a^t = -\frac{2m}{r} \left( 1 - \frac{2m}{r} \right)^{-1} u^t u^r,$$
$$a^r =  -\frac{m}{r^2} \left( 1 - \frac{2m}{r}\right) \left(u^t\right)^2 + \frac{m}{r^2} \left( 1 - \frac{2m}{r}\right)^{-1} \left( u^r\right)^2 + (r - 2m)\left( u^\phi \right)^2,$$
$$a^\phi = -\frac{2}{r} u^\phi u^r,$$

where $\theta$ has been set to $\frac{\pi}{2}$ for motion in the equatorial plane.

"""

# ╔═╡ 787bd797-e53a-4e56-878f-40bc8f1c9fa0
md"## a)"

# ╔═╡ 299da424-a0d6-43c5-a067-2faed5f691f9
md"""
Given:

$$x^\alpha_0 = (0, R, \frac{\pi}{2}, 0),$$

and that:

$$u^\alpha_0 = (u^t_0, 0, 0, C),$$

the initial condition for $u^t$ can be calculated from the normalisation of the four-velocity. .

The metric is:

$$g_{\alpha \beta} = \begin{bmatrix} -\left(1-\frac{2m}{r}\right) & 0 & 0 & 0 \\
		0 & \left( 1 - \frac{2m}{r} \right)^{-1} & 0 & 0\\
		0 & 0 & r^2 & 0 \\
		0 & 0 & 0 & r^2 \sin^2(\theta)
\end{bmatrix},$$

and the norm of the four-velocity is:

$$\mathbf{u}\cdot \mathbf{u} = g_{\alpha \beta} u^\alpha u^\beta,$$

Substituting the initial conditions:

$$\mathbf{u}\cdot \mathbf{u} = -\left(1 - \frac{2m}{R}\right) (u^t_0)^2 + R^2(u^\phi_0)^2$$
$$\implies \mathbf{u}\cdot \mathbf{u} = -\left(1 - \frac{2m}{R}\right) (u^t_0)^2 + R^2C^2$$
$$\implies \left(1 - \frac{2m}{R}\right) (u^t_0)^2 = R^2C^2 - \mathbf{u}\cdot \mathbf{u}$$
$$\implies (u^t_0)^2 = \frac{R^2C^2 - \mathbf{u}\cdot \mathbf{u}}{1 - \frac{2m}{R}}$$
$$\implies (u^t_0)^2 = \frac{(R^2C^2 - \mathbf{u}\cdot \mathbf{u})R}{R - 2m}$$
$$\implies u^t_0 = \sqrt{\frac{(R^2C^2 - \mathbf{u}\cdot \mathbf{u})R}{R - 2m}}$$
Taking the positive solution so that the coordinate time and the proper time pass in the same direction.
"""

# ╔═╡ 8d6c8823-f991-496b-8592-c7c82cdae82a
md"## c)"

# ╔═╡ d147f648-17cf-465a-b19f-887f0d76525f
md"First set the parameters and initial conditions:"

# ╔═╡ bd505786-290f-4a27-9f8e-d01c602fc43f
𝑚, 𝑅 = 1.0, 10.0;

# ╔═╡ 3963da44-b6f7-4f79-960f-3018b7190c18
𝑢ᵗ₀(𝐶, 𝑢𝑢) = sqrt(𝑅*(𝑅^2*𝐶^2 - 𝑢𝑢)/(𝑅 - 2*𝑚))

# ╔═╡ d4a4b5a5-77a4-49e2-94cf-1f6424e5e78c
md"Then the equations of motion, in a function for convenience::"

# ╔═╡ 8f553eaf-a008-4902-97a0-95a50eeb0e47
function f(𝐶, 𝑢𝑢, tmax)
	d𝐗d𝜏((𝑡, 𝑟, 𝜙, 𝑢ᵗ, 𝑢ʳ, 𝑢ᵠ), (𝑚, 𝐶, 𝑅), 𝜏) = [𝑢ᵗ, 𝑢ʳ, 𝑢ᵠ,
			(-2𝑚/𝑟^2)*(1-2𝑚/𝑟)^(-1)*𝑢ᵗ*𝑢ʳ, 									    # 𝑎ᵗ
			(-𝑚/𝑟^2)*(1-2𝑚/𝑟)*𝑢ᵗ^2 + (𝑚/𝑟^2)*(1-2𝑚/𝑟)^(-1)*𝑢ʳ^2 + (𝑟 - 2*𝑚)*𝑢ᵠ^2, # 𝑎ʳ
			(-2/𝑟)*𝑢ʳ*𝑢ᵠ ];													     # 𝑎ᵠ
	
	🚀 = Process(
      process = d𝐗d𝜏(P) = process2solution(P),
	  parameter_profile = [𝑚, 𝐶, 𝑅],
	  varnames = [:𝑡, :𝑟, :𝜙, :𝑢ᵗ, :𝑢ʳ, :𝑢ᵠ],
      X0 	   = [0.0, 𝑅, 0.0, 𝑢ᵗ₀(𝐶, 𝑢𝑢), 0, 𝐶],
      t0 = 0.0,
	  tmax = tmax,
      alg = Vern9(),
      solver_opts = Dict(:adaptive => true, :reltol => 1e-10, :abstol => 1e-10));
	
	return [times(🚀), (eachcol∘timeseries)(🚀)...] 
end;

# ╔═╡ dae3bcba-0db2-48ca-8135-f6495c367418
md"### i)$\quad C = (10\sqrt{7})^{-1}$"

# ╔═╡ 3ddbb46c-f02a-4137-90b1-1f8a7ce0b26f
𝜏₁, 𝑡₁, 𝑟₁, 𝜙₁, 𝑢ᵗ₁, 𝑢ʳ₁, 𝑢ᵠ₁ = f((10*sqrt(7))^-1, -1, 1000);

# ╔═╡ eb4112f9-9acb-4fea-a9a2-886035fded63
plot(𝜏₁, 𝑡₁, xguide="𝜏", yguide="𝑡", xlims=(-100, 1100), 
	title="Coordinate time")

# ╔═╡ 801391ab-f8a7-44ff-baac-2cb20f4448db
plot(𝜏₁, 𝑟₁, xguide="𝜏", yguide="𝑟", xlims=(-100, 1100), ylims=(0, 15), 
	title="Radius")

# ╔═╡ 5cb5cfb1-3374-49be-b01b-3c6e90841c84
plot(𝜏₁, 𝜙₁, xguide="𝜏", yguide="𝜙", xlims=(-100, 1100), ylims=(-1, 40), 
	title="Angle")

# ╔═╡ fc2a7f13-1d07-4238-bd47-54d719ac815b
plot(𝑟₁.*cos.(𝜙₁), 𝑟₁.*sin.(𝜙₁), xguide="𝑟cos(𝜙)", yguide="𝑟sin(𝜙)", 				 
			aspect_ratio=:equal, title="A circular orbit in the equatorial plane")

# ╔═╡ 64547853-09a1-45fc-9f31-823ec09d2733
md"""
Using the metric and equation for the norm of the four-velocity written in part a):
"""

# ╔═╡ da839751-1294-43a3-8eea-eeddf6c61835
𝐮𝐮(𝑟, 𝑢ᵗ, 𝑢ʳ, 𝑢ᵠ) = -(1-2𝑚/𝑟)*𝑢ᵗ^2 + (1-2𝑚/𝑟)^(-1)*𝑢ʳ^2 + 𝑟^2*𝑢ᵠ^2;

# ╔═╡ f3442ed7-b1cb-4782-85c6-9ef197ee34f4
plot(𝜏₁, 𝐮𝐮.(𝑟₁, 𝑢ᵗ₁, 𝑢ʳ₁, 𝑢ᵠ₁), xguide="𝜏", yguide="𝐮⋅𝐮", title="Norm of the four-velocity")

# ╔═╡ 0cd14138-57c0-4b40-bc61-1547c63e8c7a
md"Note that with `Vern9()` the error in the normof the four-velocity is to small to be appear in float values. For a poorer integrator, such as the midpoint method, it is on the order of $10^{-16}$."

# ╔═╡ 20303375-55f6-4d5f-abb8-4427631e960b
md"### ii)$\quad C = 1.1(10\sqrt{7})^{-1}$"

# ╔═╡ 5372aee7-a876-4fbb-8b66-2fd2fccc54f0
𝜏₂, 𝑡₂, 𝑟₂, 𝜙₂, 𝑢ᵗ₂, 𝑢ʳ₂, 𝑢ᵠ₂ = f(1.1*(10*sqrt(7))^-1, -1, 1000);

# ╔═╡ eae405cf-e409-43f6-b042-d2402f7f001e
plot(𝜏₂, 𝑡₂, xguide="𝜏", yguide="𝑡", xlims=(-100, 1100), 
	title="Coordinate time")

# ╔═╡ 97de977c-6c44-48d4-97ed-32da9e584300
plot(𝜏₂, 𝑟₂, xguide="𝜏", yguide="𝑟", xlims=(-100, 1100), ylims=(-1, 21), 
	title="Radius")

# ╔═╡ 29d74095-fdfb-4818-a177-fb815f768819
plot(𝜏₂, 𝜙₂, xguide="𝜏", yguide="𝜙", xlims=(-100, 1100), ylims=(-1, 22), 
	title="Angle")

# ╔═╡ 1635434d-cb88-4eb4-b15a-148686500da0
plot(𝑟₂.*cos.(𝜙₂), 𝑟₂.*sin.(𝜙₂), xguide="𝑟cos(𝜙)", yguide="𝑟sin(𝜙)", 				 
			aspect_ratio=:equal, title="A precessing orbit in the equatorial plane")

# ╔═╡ 391cbd14-cc90-4d37-a9fc-cba8d8d0b16b
plot(𝜏₂, 𝐮𝐮.(𝑟₂, 𝑢ᵗ₂, 𝑢ʳ₂, 𝑢ᵠ₂), xguide="𝜏", yguide="𝐮⋅𝐮", title="Norm of the four-velocity", ylims=(-1.01, 0.01))

# ╔═╡ 004a66f7-fbad-478e-8122-ebb0de7b8a63
md"This time the error is slightly greater, since the orbit is no longer a perfect circle:"

# ╔═╡ 187d5fc7-667c-4116-8c39-e47b1c9f85a2
plot(𝜏₂, 𝐮𝐮.(𝑟₂, 𝑢ᵗ₂, 𝑢ʳ₂, 𝑢ᵠ₂) .+ 1, xguide="𝜏", yguide="𝐮⋅𝐮 + 1", title="Error in the four-velocity norm")

# ╔═╡ ac443eb8-387c-425c-b878-510dd277d15e
md"""
This orbit precesses, unlike orbits under Netwonian gravity which are strictly conic sections. 
This is because general relativity, or the procedure to derive the path of an object from geodesic equaitons and the structure of spacetime, introduces an additional term to the effective potential of Newtonian gravity.
Both descriptions of gravity have stable minima at someradius that is determined by the mass of the central body and the angular momentum of the orbiting object.
Objects with apsides at this stable radius will follow a circular path, as for the first orbit shown in this question; their efective energy is minimised.
Objects with apsides displaced from this stable minimum will precess; their effective energy is greather than the potential at the stable minimum.
Their radius will oscillate between a point with an $r$ smaller than the stable minimum, and one with an $r$ greater (i.e. an elliptical orbit).
This occurs in both Newtonian gravity and general relativity.
However, Newtonian gravity has an effective potential with two terms: a radial force and an angular momentum. Hence elliptical orbits under Newtonian gravity are angularly stable and do not precess.
General relativity introduces a third term that is cubic in the radius and also depends on angular momentum; this term breaks the (angular) stability of elliptical orbits and causes them to precess (the angle $\phi$ of the apsides is not constant). 

(The effective potential described by general relativity also approaches negative infinity near the origin, so there are radially unstable states that do not appear in this assignment.)
"""

# ╔═╡ 79c379d0-e55b-4b78-aaa2-aa501e0b3a5f
md"# d)"

# ╔═╡ a9636721-42d8-471d-a22a-05fd9fa32bd7
md"""
For a massless particle, 𝐮⋅𝐮 is $0$.
Additionally, the variable $\tau$ is the equations of motion no longer represents proper time but an affine parameter ($0 \le \lambda \le 10$).

"""

# ╔═╡ 24090043-c604-42f4-8ef9-866ae0033c0d
md"### i) $\quad C = 0.5$"

# ╔═╡ 1b6ed9d6-fa3d-4d7d-9b24-0022cbfdbaa2
𝜏₃, 𝑡₃, 𝑟₃, 𝜙₃, 𝑢ᵗ₃, 𝑢ʳ₃, 𝑢ᵠ₃ = f(0.5, 0, 10);

# ╔═╡ 0d059280-6f6b-4c98-bcb1-11ecd17de395
𝜏₄, 𝑡₄, 𝑟₄, 𝜙₄, 𝑢ᵗ₄, 𝑢ʳ₄, 𝑢ᵠ₄ = f(1, 0, 10);

# ╔═╡ 4bc48a18-8d1d-4f41-8ff1-de690bc04d95
𝜏₅, 𝑡₅, 𝑟₅, 𝜙₅, 𝑢ᵗ₅, 𝑢ʳ₅, 𝑢ᵠ₅ = f(2, 0, 10);

# ╔═╡ 877d135e-6949-4a8c-8589-769f5d708b5a
plot(𝑟₃.*cos.(𝜙₃), 𝑟₃.*sin.(𝜙₃), xguide="𝑟cos(𝜙)", yguide="𝑟sin(𝜙)", ylims=(0, 200), 
			aspect_ratio=:equal, title="The path of a massless particle with 𝐶 = 0.5")

# ╔═╡ 3cd20c20-9500-44d8-ba0b-c043f9658749
plot(𝑟₄.*cos.(𝜙₄), 𝑟₄.*sin.(𝜙₄), xguide="𝑟cos(𝜙)", yguide="𝑟sin(𝜙)", ylims=(0, 200),
			aspect_ratio=:equal, title="The path of a massless particle with 𝐶 = 1")

# ╔═╡ 9093d939-2ce7-4ed4-8ca2-62527d726842
plot(𝑟₅.*cos.(𝜙₅), 𝑟₅.*sin.(𝜙₅), xguide="𝑟cos(𝜙)", yguide="𝑟sin(𝜙)", ylims=(0, 200), 
			aspect_ratio=:equal, title="The path of a massless particle with 𝐶 = 2")

# ╔═╡ 0bdd55ed-8451-46c2-8815-a40ea27d0bdc
md"""
In this question, the massless particle passes by a massive body and its trajectory (in Cartesian coordinates) is bent. This makes clear the existence of gravity as curvature in space-time, since it can influence even the paths of massless particles.
Increasing the intial component of the four-velocity, 𝐶, causes the total spatial length of the massless particle's path to increase.
However, the trajectory it ultimatley follows is the same; this is because the affine parameter is defined such that it can be substituted in place of (or, to generalise) $\tau$ in the equations of motion for a massive particle 
This gives a geodesic equation valid for a massless (0 four-velocity norm) particles.
As such, the affine parameter can be linearly rescaled to give the same trajectory in spacetime for a massless photon; the affine parameter does not have an absolute correspondence to a proper time. 
Then 𝐶 has no physical meaning as a velocity over a proper time for a massless partical, but does act to rescale the affine parameter. The result is that increasing 𝐶 produces a longer path (in space, or, Cartesian coordinates) for the same interval of affine parameters, as demonstrated above.
"""

# ╔═╡ Cell order:
# ╠═09bee220-1249-11ec-3e2a-4fdde9976a59
# ╟─81b414c9-2b23-4df4-9eaf-5549c83df269
# ╟─1cbf17f7-8322-4e3c-a895-163a48332205
# ╟─08ca86b3-2c97-4486-bb4d-5c00d8239119
# ╟─787bd797-e53a-4e56-878f-40bc8f1c9fa0
# ╟─299da424-a0d6-43c5-a067-2faed5f691f9
# ╟─8d6c8823-f991-496b-8592-c7c82cdae82a
# ╟─d147f648-17cf-465a-b19f-887f0d76525f
# ╠═bd505786-290f-4a27-9f8e-d01c602fc43f
# ╠═3963da44-b6f7-4f79-960f-3018b7190c18
# ╟─d4a4b5a5-77a4-49e2-94cf-1f6424e5e78c
# ╠═8f553eaf-a008-4902-97a0-95a50eeb0e47
# ╟─dae3bcba-0db2-48ca-8135-f6495c367418
# ╠═3ddbb46c-f02a-4137-90b1-1f8a7ce0b26f
# ╠═eb4112f9-9acb-4fea-a9a2-886035fded63
# ╠═801391ab-f8a7-44ff-baac-2cb20f4448db
# ╠═5cb5cfb1-3374-49be-b01b-3c6e90841c84
# ╠═fc2a7f13-1d07-4238-bd47-54d719ac815b
# ╟─64547853-09a1-45fc-9f31-823ec09d2733
# ╠═da839751-1294-43a3-8eea-eeddf6c61835
# ╠═f3442ed7-b1cb-4782-85c6-9ef197ee34f4
# ╟─0cd14138-57c0-4b40-bc61-1547c63e8c7a
# ╟─20303375-55f6-4d5f-abb8-4427631e960b
# ╠═5372aee7-a876-4fbb-8b66-2fd2fccc54f0
# ╠═eae405cf-e409-43f6-b042-d2402f7f001e
# ╠═97de977c-6c44-48d4-97ed-32da9e584300
# ╠═29d74095-fdfb-4818-a177-fb815f768819
# ╠═1635434d-cb88-4eb4-b15a-148686500da0
# ╠═391cbd14-cc90-4d37-a9fc-cba8d8d0b16b
# ╟─004a66f7-fbad-478e-8122-ebb0de7b8a63
# ╠═187d5fc7-667c-4116-8c39-e47b1c9f85a2
# ╟─ac443eb8-387c-425c-b878-510dd277d15e
# ╟─79c379d0-e55b-4b78-aaa2-aa501e0b3a5f
# ╟─a9636721-42d8-471d-a22a-05fd9fa32bd7
# ╟─24090043-c604-42f4-8ef9-866ae0033c0d
# ╠═1b6ed9d6-fa3d-4d7d-9b24-0022cbfdbaa2
# ╠═0d059280-6f6b-4c98-bcb1-11ecd17de395
# ╠═4bc48a18-8d1d-4f41-8ff1-de690bc04d95
# ╠═877d135e-6949-4a8c-8589-769f5d708b5a
# ╠═3cd20c20-9500-44d8-ba0b-c043f9658749
# ╠═9093d939-2ce7-4ed4-8ca2-62527d726842
# ╟─0bdd55ed-8451-46c2-8815-a40ea27d0bdc
