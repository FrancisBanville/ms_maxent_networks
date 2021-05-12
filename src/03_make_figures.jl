## Read predicted networks 

# Simulated networks with simulated numbers of links
predicted_networks = load(joinpath("data", "sim", "predicted_networks.jld"))["data"]

# Simulated networks with empirical numbers of links and species
predicted_networks_empL = load(joinpath("data", "sim", "predicted_networks_empL.jld"))["data"]


## Compute measures of food webs of maximum entropy

# richness 
predicted_networks_S = richness.(predicted_networks)
# nestedness (spectral radius of the adjacency matrix)
predicted_networks_nestedness = ρ.(predicted_networks)
# maximum and average trophic levels
predicted_networks_tls = []
working = []
for i in 1:length(predicted_networks)
  try
  push!(predicted_networks_tls, values(trophic_level(predicted_networks[i])))
  push!(working, i)
  catch
  println("could not estimate trophic levels in network $(i)")
  end
end

predicted_networks_tl_max = maximum.(predicted_networks_tls)
predicted_networks_tl_avg = sum.(predicted_networks_tls)./length.(predicted_networks_tls)

## Measures of food webs archived on mangal.io

# Read food webs and convert them to unipartie networks 
Ns_mangal = network.(mangal_foodwebs.id)
Ns_mangal = convert.(UnipartiteNetwork, Ns_mangal)

# Compute food-web measures 
# richness 
Ns_mangal_S = richness.(Ns_mangal)
# nestedness (spectral radius of the adjacency matrix)
Ns_mangal_nestedness = ρ.(Ns_mangal)
# maximum and average trophic levels
Ns_mangal_tls = values.(trophic_level.(Ns_mangal))
Ns_mangal_tl_max = maximum.(Ns_mangal_tls)
Ns_mangal_tl_avg = sum.(Ns_mangal_tls)./length.(Ns_mangal_tls)



## Figures

# Nestedness (MaxEnt and Mangal food webs)
scatter(Ns_mangal_nestedness, predicted_networks_nestedness, alpha=0.5, lab="",
title="Nestedness", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(0.4, 0.9, 100), LinRange(0.4, 0.9, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs")
savefig(joinpath("figures", "nestedness.png"))


# Maximum trophic level (MaxEnt and Mangal food webs)
scatter(Ns_mangal_tl_max[working], predicted_networks_tl_max, alpha=0.5, lab="",
title="Maximum trophic level", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(2, 13, 100), LinRange(2, 13, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs")
savefig(joinpath("figures", "maximum_trophic_level.png"))

# Average trophic level (MaxEnt and Mangal food webs)
scatter(Ns_mangal_tl_avg[working], predicted_networks_tl_avg, alpha=0.5, lab="",
title="Average trophic level", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(1, 7, 100), LinRange(1, 7, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs")
savefig(joinpath("figures", "average_trophic_level.png"))

