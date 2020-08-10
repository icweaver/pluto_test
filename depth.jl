### A Pluto.jl notebook ###
# v0.11.4

using Markdown
using InteractiveUtils

# ╔═╡ c9ac27ee-dac0-11ea-2a8c-2d144b034a82
md"""
# Exoplanet Calculator 🪐
"""

# ╔═╡ b2286b26-dac2-11ea-1ce0-c7da562aa641
md"Given exoplanet and host star parameters from the literature, calculate derived values relevant for detection of the planet's atmosphere"

# ╔═╡ 19b35ef4-dac3-11ea-2d25-97e5482ff6a0
md"### Literature values"

# ╔═╡ 75d6dcbe-db0a-11ea-2839-9542a238b679
md"Source: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=HAT-P-23%20b)"

# ╔═╡ 0b6821a4-dac3-11ea-27d7-911521f0d3c0
md"### Calculate all the things"

# ╔═╡ 02bfa078-d62b-11ea-15df-d701431829b9
begin
	using Measurements, Unitful, UnitfulAstro, Markdown
	using PhysicalConstants.CODATA2018: G, k_B, m_u
	const amu, k = m_u, k_B
end;

# ╔═╡ 17302d74-d63b-11ea-3de3-49f0df0554ca
# Input params from studies to compare
studies = [
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
		K    = (346.0 ± 21)u"m/s", # latest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		aRₛ  = 4.16 ± 0.26,
		Tₛ   = (5734 ± 99.8735)u"K",
		Rₛ   = (1.1858169 ± 0.0424133)u"Rsun",
	),
	(
		name = "GAIA DR2 w/ DR1 mass",
		μ    = 2.0*amu,
		α    = 0.0,
		K    = (346.0 ± 21)u"m/s", # latest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		Mₚ   = (1.34 ± 0.59)u"Mjup",
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		aRₛ  = 4.16 ± 0.26,
		Tₛ   = (5734 ± 99.8735)u"K",
		Rₛ   = (1.1858169 ± 0.0424133)u"Rsun",
	),
];

# ╔═╡ 3f79c516-da77-11ea-1f6b-d3e7191a95d8
begin
	# Star density
	get_ρₛ(; P, aRₛ) = (3.0π / (G * P^2)) * aRₛ^3.0

	# Star mass
	get_Mₛ(; ρₛ, Rₛ) = ρₛ * (4.0/3.0) * π * Rₛ^3.0

	#Planet mass
	get_Mₚ(; K, i, P, Mₛ) = (K/sin(i)) * (P / (2.0π*G))^(1//3) * Mₛ^(2//3)

	# Star surface gravity
	get_gₛ(; Mₛ, Rₛ) = G * Mₛ / Rₛ^2

	# Planet surface gravity
	get_gₚ(; Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

	# Planet equilibrium temperature
	get_Tₚ(; Tₛ, aRₛ, α) = Tₛ * (1.0 - α)^(1//4) * (0.5/aRₛ)^(1//2)

	# Planet scale height
	get_H(; μ, Tₚ, gₚ) = k * Tₚ / (μ * gₚ)

	# Estimated signal from planet atmosphere
	get_Delta_D(; H, RₚRₛ, Rₛ) = 2.0 * H * RₚRₛ/Rₛ
end;

# ╔═╡ 3833772c-d63f-11ea-09b5-f36d68e512ea
begin
	results = []
	for st in studies
		# Calculate secondary params if not given
		ρₛ = (haskey(st, :ρₛ)) ? st.ρₛ : get_ρₛ(P=st.P, aRₛ=st.aRₛ)
		Mₛ = (haskey(st, :Mₛ)) ? st.Mₛ : get_Mₛ(ρₛ=ρₛ, Rₛ=st.Rₛ)
		Mₚ = (haskey(st, :Mₚ)) ? st.Mₚ : get_Mₚ(K=st.K, i=st.i, P=st.P, Mₛ=Mₛ)
		
		# Calculate remaining params
		gₛ = get_gₛ(Mₛ=Mₛ, Rₛ=st.Rₛ)
		gₚ  = get_gₚ(Mₚ=Mₚ, RₚRₛ=st.RₚRₛ, Rₛ=st.Rₛ)
		Tₚ = get_Tₚ(Tₛ=st.Tₛ, aRₛ=st.aRₛ, α=st.α)
		H  = get_H(μ=st.μ, Tₚ=Tₚ, gₚ=gₚ)
		ΔD = get_Delta_D(H=H, RₚRₛ=st.RₚRₛ, Rₛ=st.Rₛ)

		# Collect results
		m = md"""
		**$(st.name):**

		log gₛ (cm/s²) = $(log10(ustrip(uconvert(u"cm/s^2", gₛ))))
		
		log gₚ (cm/s²) = $(log10(ustrip(uconvert(u"cm/s^2", gₚ))))

		ρₛ = $(uconvert(u"g/cm^3", ρₛ))

		Mₛ = $(uconvert(u"Msun", Mₛ))

		Mₚ = $(uconvert(u"Mjup", Mₚ))

		Tₚ = $Tₚ

		H = $(uconvert(u"km", H))

		ΔD = $(5 * upreferred(ΔD) * 1e6) ppm
		
		"""
		push!(results, m)
	end
	
	# Display results
	results
end

# ╔═╡ Cell order:
# ╟─c9ac27ee-dac0-11ea-2a8c-2d144b034a82
# ╟─b2286b26-dac2-11ea-1ce0-c7da562aa641
# ╟─19b35ef4-dac3-11ea-2d25-97e5482ff6a0
# ╟─75d6dcbe-db0a-11ea-2839-9542a238b679
# ╠═17302d74-d63b-11ea-3de3-49f0df0554ca
# ╟─0b6821a4-dac3-11ea-27d7-911521f0d3c0
# ╠═3833772c-d63f-11ea-09b5-f36d68e512ea
# ╠═3f79c516-da77-11ea-1f6b-d3e7191a95d8
# ╠═02bfa078-d62b-11ea-15df-d701431829b9
