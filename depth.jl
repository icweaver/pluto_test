### A Pluto.jl notebook ###
# v0.11.4

using Markdown
using InteractiveUtils

# ╔═╡ 02bfa078-d62b-11ea-15df-d701431829b9
begin
	using Measurements, Unitful, UnitfulAstro, Markdown
	using PhysicalConstants.CODATA2018: G, k_B, m_u
	const amu, k = m_u, k_B
end;

# ╔═╡ 17302d74-d63b-11ea-3de3-49f0df0554ca
# Input params from study
studies = [
	(
		name = "Yea",
		μ    = 2.0*amu,
		α    = 0.0,
		K    = (346.0 ± 21)u"m/s", # newest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # newest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # newest transit data, from S&R16 
		Mₚ   = (1.34 ± 0.59)u"Mjup",
		RₚRₛ = 0.1113 ± 0.0010,
		aRₛ  = 4.26 ± 0.26,
		Tₛ   = (5734 ± 99.8735)u"K",
		Rₛ   = (1.1858169 ± 0.0424133)u"Rsun",
	),
	(
		name = "Ciceri et al. (2015)",
		μ    = 2.0*amu,
		α    = 0.0,
		K    = (368.5 ± 17.6)u"m/s",
		i    = (85.74 ± 0.95)u"°",
		P    = (1.21288287 ± 0.00000017)u"d",
		RₚRₛ = 0.11616 ± 0.00081,
		aRₛ  = 4.5459 ± 0.0919,
		Tₛ   = (5885 ± 72)u"K",
		Rₛ   = (1.089 ± 0.028)u"Rsun",
	),
	(
		name = "GAIA DR2",
		μ    = 2.0*amu,
		α    = 0.0,
		K    = (346.0 ± 21)u"m/s", # newest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # newest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # newest transit data, from S&R16
		RₚRₛ = 0.1113 ± 0.0010, # newest transit data, from S&R16
		aRₛ  = 4.16 ± 0.26,
		Tₛ   = (5734 ± 99.8735)u"K",
		Rₛ   = (1.1858169 ± 0.0424133)u"Rsun",
	),
];

# ╔═╡ 3f79c516-da77-11ea-1f6b-d3e7191a95d8
get_ρₛ(; P, aRₛ) = (3.0π / (G * P^2)) * aRₛ^3.0

# ╔═╡ f5b9abe4-da76-11ea-1e42-25ade960ed52
get_Mₛ(; ρₛ, Rₛ) = ρₛ * (4.0/3.0) * π * Rₛ^3.0

# ╔═╡ 480ed8b0-da72-11ea-18b0-07879c66ecfc
get_Mₚ(; K, i, P, Mₛ) = (K/sin(i)) * (P / (2.0π*G))^(1//3) * Mₛ^(2//3)

# ╔═╡ 2fa77fc8-da80-11ea-193f-911962ef9892
get_gₛ(; Mₛ, Rₛ) = G * Mₛ / Rₛ^2

# ╔═╡ ba7e96ce-d630-11ea-350d-cb961d23b482
get_g(; Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

# ╔═╡ 3f1ef4fe-d62b-11ea-3694-7bfea6c78d25
get_Tₚ(; Tₛ, aRₛ, α) = Tₛ * (1.0 - α)^(1//4) * (0.5/aRₛ)^(1//2)

# ╔═╡ 673a3e64-d635-11ea-2668-713470482653
get_H(; μ, Tₚ, g) = k * Tₚ / (μ * g)

# ╔═╡ 363e5c42-d639-11ea-24a5-31e2094480b9
get_Delta_D(; H, RₚRₛ, Rₛ) = 2.0 * H * RₚRₛ/Rₛ

# ╔═╡ 3833772c-d63f-11ea-09b5-f36d68e512ea
begin
	results = []
	for p in studies
		# Derived params
		ρₛ = (haskey(p, :ρₛ)) ? p.ρₛ : get_ρₛ(P=p.P, aRₛ=p.aRₛ)
		Mₛ = (haskey(p, :Mₛ)) ? p.Mₛ : get_Mₛ(ρₛ=ρₛ, Rₛ=p.Rₛ)
		Mₚ = (haskey(p, :Mₚ)) ? p.Mₚ : get_Mₚ(K=p.K, i=p.i, P=p.P, Mₛ=Mₛ)
		gₛ = get_gₛ(Mₛ=Mₛ, Rₛ=p.Rₛ)
		g  = get_g(Mₚ=Mₚ, RₚRₛ=p.RₚRₛ, Rₛ=p.Rₛ)
		Tₚ = get_Tₚ(Tₛ=p.Tₛ, aRₛ=p.aRₛ, α=p.α)
		H  = get_H(μ=p.μ, Tₚ=Tₚ, g=g)
		ΔD = get_Delta_D(H=H, RₚRₛ=p.RₚRₛ, Rₛ=p.Rₛ)

		m = md"""
		**Derived results for $(p.name)**:

		log gₛ (cm/s²) = $(log10(ustrip(uconvert(u"cm/s^2", gₛ))))

		ρₛ = $(uconvert(u"g/cm^3", ρₛ))

		Mₛ = $(uconvert(u"Msun", Mₛ))

		Mₚ = $(uconvert(u"Mjup", Mₚ))

		Tₚ = $Tₚ

		H = $(uconvert(u"km", H))

		ΔD = $(5 * upreferred(ΔD) * 1e6) ppm

		"""

		push!(results, m)
	end
end

# ╔═╡ 91767e6c-da98-11ea-3722-eff13b07d0d7
results

# ╔═╡ Cell order:
# ╠═17302d74-d63b-11ea-3de3-49f0df0554ca
# ╠═3833772c-d63f-11ea-09b5-f36d68e512ea
# ╠═91767e6c-da98-11ea-3722-eff13b07d0d7
# ╠═3f79c516-da77-11ea-1f6b-d3e7191a95d8
# ╠═f5b9abe4-da76-11ea-1e42-25ade960ed52
# ╠═480ed8b0-da72-11ea-18b0-07879c66ecfc
# ╠═2fa77fc8-da80-11ea-193f-911962ef9892
# ╠═ba7e96ce-d630-11ea-350d-cb961d23b482
# ╠═3f1ef4fe-d62b-11ea-3694-7bfea6c78d25
# ╠═673a3e64-d635-11ea-2668-713470482653
# ╠═363e5c42-d639-11ea-24a5-31e2094480b9
# ╠═02bfa078-d62b-11ea-15df-d701431829b9
