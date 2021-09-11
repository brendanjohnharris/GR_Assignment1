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

# ╔═╡ 1cbf17f7-8322-4e3c-a895-163a48332205
md"# Question 3"

# ╔═╡ 08ca86b3-2c97-4486-bb4d-5c00d8239119
md"""
The equations of motion in the Schwarzschild metric are:

$$\frac{dt}{d\tau} = e \left( 1 - \frac{2m}{r}\right)^{-1},$$
$$\left(\frac{dr}{d\tau} \right)^2 = e^2 - 1 + \frac{2m}{r} - \frac{l^2}{r^2} + \frac{2ml^2}{r^3},$$
$$\frac{d\phi}{d\tau} = \frac{l}{r^2},$$

where $\theta$ has been set to $\frac{\pi}{2}$ for motion in the equatorial plane.
Defining a potential $V$ as:

$$V(r) = -\frac{m}{r} + \frac{l^2}{2r^2} - \frac{ml^2}{r^3},$$

and using the initial conditions of $r(0) = R$, $u^{\phi}(0) = C$ and motion only in the $\phi$-direction ($dr/d\tau = 0$) gives values for the constants:

$$l = CR^2,$$
$$e^2 = 2V(R) + 1,$$

where positive $e$ is taken to be consistent with its interpretation as an energy density at large $r$. Then the equations of motion can be simplified to:

$$\frac{dt}{d\tau} = e \left( 1 - \frac{2m}{r}\right)^{-1},$$
$$\left(\frac{dr}{d\tau} \right)^2 = 2\left[V(R) - V(r)\right],$$
$$\frac{d\phi}{d\tau} = \frac{l}{r^2},$$
"""

# ╔═╡ 787bd797-e53a-4e56-878f-40bc8f1c9fa0
md"## a)"

# ╔═╡ 299da424-a0d6-43c5-a067-2faed5f691f9
md"""
................ find initial conditions..................
"""

# ╔═╡ 8d6c8823-f991-496b-8592-c7c82cdae82a
md"## b)"

# ╔═╡ d147f648-17cf-465a-b19f-887f0d76525f
md"First set the parameters and conserved quantities:"

# ╔═╡ bd505786-290f-4a27-9f8e-d01c602fc43f
𝑚, 𝐶, 𝑅 = 1.0, (10*√7)^-1, 10.0;

# ╔═╡ cf39fb2f-3f1c-4fe4-88fb-8bdec4afc2f6
𝑙 = 𝐶*𝑅^2;

# ╔═╡ f717d6ce-4a74-48f1-b1f1-fe4e16081dcc
md"Next define the potential"

# ╔═╡ 6a6c4ab1-378b-4055-bd13-53ff50ec1d77
𝑉(𝑟) = -𝑚/𝑟 + 𝑙^2/(2*𝑟^2) - 𝑚*𝑙^2/𝑟^3;

# ╔═╡ 00a0999b-6586-4aad-a334-fa83156a3ba6
𝑒 = sqrt(2*𝑉(𝑅) + 1);

# ╔═╡ d4a4b5a5-77a4-49e2-94cf-1f6424e5e78c
md"Then the equations of motion:"

# ╔═╡ f452c2ee-9e84-4abc-a5c0-9157922e1e0b
d𝐗d𝜏((𝑡, 𝑟, 𝜙), (𝑚, 𝐶, 𝑅), 𝜏) = [				     𝑒*(1-2*𝑚/𝑟)^(-1), 
										  		    sqrt(2*(𝑉(𝑅) - 𝑉(𝑟))), 
								  					𝑙/𝑟^2,
							    ];

# ╔═╡ 2a3c27b1-16ed-44a5-8ad8-33b201d14801
🚀 = Process(
      process = d𝐗d𝜏(P) = process2solution(P),
	  parameter_profile = [𝑚, 𝐶, 𝑅],
	  varnames = [:𝑡, :𝑟, :𝜙],
      X0 	   = [0.0, 𝑅, 0.0],
      t0 = 0.0,
	  tmax = 1000.0,
      alg = Vern9(), 
      solver_opts = Dict(:adaptive => true, :reltol => 1e-10, :abstol => 1e-10)); 

# ╔═╡ 35b961ec-8455-43c7-8779-1efe4884a3c8
𝜏, 𝑡, 𝑟, 𝜙 = times(🚀), (eachcol∘timeseries)(🚀)...; 

# ╔═╡ 376782f9-c3d3-4b64-a70f-8e88cff2c074
plot(𝑟, 𝑡, xguide="𝑟", yguide="𝑡")

# ╔═╡ 40b08363-94ed-4ee8-a779-9ed2412fd0c6
plot(𝑟, 𝜙, xguide="𝜙", yguide="𝑡")

# ╔═╡ fc2a7f13-1d07-4238-bd47-54d719ac815b
plot(𝑟.*cos.(𝜙), 𝑟.*sin.(𝜙), xguide="𝑟 cos(𝜙)", yguide="𝑟sin(𝜙)")

# ╔═╡ baff4215-c290-480d-91fe-1875bbaff5e6
norm(v) = -v[1]^2 + reduce(+, v[2:end].^2)

# ╔═╡ 8fa8907c-7fa4-4ef9-9afe-cced267494c1
𝑢 = [d𝐗d𝜏((𝑡[i], 𝑟[i], 𝜙[i]), (𝑚, 𝐶, 𝑅), 𝜏[i]) for i ∈ 1:length(𝜏)]

# ╔═╡ 5767121e-f50a-45fd-a4f5-7c141cd15c44
norm.(𝑢) 

# ╔═╡ 3c084e65-f147-4f2c-a437-dfd391cae54e
sqrt((𝑅^2*𝐶^2 - 1)*𝑅/(𝑅-2*𝑚))

# ╔═╡ 69a9f373-e300-467a-97c5-376cda8860ce
𝑅^2*𝐶^2 - 1

# ╔═╡ Cell order:
# ╠═09bee220-1249-11ec-3e2a-4fdde9976a59
# ╟─1cbf17f7-8322-4e3c-a895-163a48332205
# ╠═08ca86b3-2c97-4486-bb4d-5c00d8239119
# ╟─787bd797-e53a-4e56-878f-40bc8f1c9fa0
# ╠═299da424-a0d6-43c5-a067-2faed5f691f9
# ╟─8d6c8823-f991-496b-8592-c7c82cdae82a
# ╠═d147f648-17cf-465a-b19f-887f0d76525f
# ╠═bd505786-290f-4a27-9f8e-d01c602fc43f
# ╠═cf39fb2f-3f1c-4fe4-88fb-8bdec4afc2f6
# ╠═f717d6ce-4a74-48f1-b1f1-fe4e16081dcc
# ╠═6a6c4ab1-378b-4055-bd13-53ff50ec1d77
# ╠═00a0999b-6586-4aad-a334-fa83156a3ba6
# ╠═d4a4b5a5-77a4-49e2-94cf-1f6424e5e78c
# ╠═f452c2ee-9e84-4abc-a5c0-9157922e1e0b
# ╠═2a3c27b1-16ed-44a5-8ad8-33b201d14801
# ╠═35b961ec-8455-43c7-8779-1efe4884a3c8
# ╠═376782f9-c3d3-4b64-a70f-8e88cff2c074
# ╠═40b08363-94ed-4ee8-a779-9ed2412fd0c6
# ╠═fc2a7f13-1d07-4238-bd47-54d719ac815b
# ╠═baff4215-c290-480d-91fe-1875bbaff5e6
# ╠═8fa8907c-7fa4-4ef9-9afe-cced267494c1
# ╠═5767121e-f50a-45fd-a4f5-7c141cd15c44
# ╠═3c084e65-f147-4f2c-a437-dfd391cae54e
# ╠═69a9f373-e300-467a-97c5-376cda8860ce
