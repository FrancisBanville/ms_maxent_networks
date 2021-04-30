## Read metadata of all food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl)

mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "mangal_foodwebs.csv")))


## Define flexible links model and infer model parameters

# S: number of species 
# R: number of flexible links that are realized (links above the minimum)

@model FL(S,R) = begin
  # number of trials
  N = length(S) 
  F = S.*S .- (S.-1)
  # parameters (priors)
  ϕ ~ Normal(3.0, 0.5)
  μ ~ Beta(3.0, 7.0)
  # flexible links model
  for i in 1:N
    R[i] ~ BetaBinomial(F[i], μ*exp(ϕ), (1-μ)*exp(ϕ))
  end
  return μ, ϕ
end

# Observations

# number of species
S = vec(mangal_foodwebs[:,:S]) 
# number of links
L = vec(mangal_foodwebs[:,:L]) 
# number of flexible links that are realized
R = L .- (S.-1) 

# Hamiltonian Monte Carlo sampler with static trajectory
chain = sample(FL(S,R), HMC(0.01,10), 3000)

# Diagnostic plot
plot(chain[200:end,:,:])


## Simulate counterfactuals

# Simulate the number of links for a range of species richness

# numbers of species
sp = collect(5:1000) 
# number of simulations
n = 1000 

predicted_links = zeros(Int64, n, length(sp))

p = Progress(length(sp))
Threads.@threads for s in sp
    i = findall(sp .== s)
    predicted_links[:, i] = predict_links(s, n, chain)
    next!(p)
end

# predicted number of links (median) for species richness in S
i = findall(in(S).(sp))
Sp = sort(unique(S))
L_med = convert.(Int64, round.(median.(eachcol(predicted_links[:, i]))))

## Predict the adjacency matrix of maximum entropy for a range of species richness and median predicted number of links

predicted_matrices = predict_adjacency_matrix.(S, L, 25)
predicted_networks = UnipartiteNetwork.(predicted_matrices)

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












## Example ##

# number of species
Sex = 8
# predicted number of links (take median)
Lex = predict_links(Sex, n, chain)
Lex = convert(Int64, round(median(Lex)))
# degree distribution of maximum entropy (probabilities)
dd_maxent_prob(Sex, Lex)
# degree distribution of maximum entropy (counts)
ddex = dd_maxent(Sex, Lex)
# joint-degree distribution of maximum entropy (probabilities)
jdd_maxent_prob(Sex, Lex)
# in-degree distribution of maximum entropy (counts)
inddex = indd_maxent(Sex, Lex)
# out-degree distribution of maximum entropy (non-adjusted counts)
outddex = ddex.- inddex
mex = length(ddex)
outddex[mex] = last(outddex) + (sum(inddex) - sum(outddex))
# first adjacency matrix constrained by in-degree and out-degree distributions
Aex = matrix_rowcolsum(inddex, outddex)
# adjacency matrix of maximum entropy
Aex = adjacency_matrix_maxent(Aex, 50)
# predict adjacency matrix of maximum entropy in one command
Bex = predict_adjacency_matrix(Sex, Lex, 50)

