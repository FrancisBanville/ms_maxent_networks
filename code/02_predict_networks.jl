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

save(joinpath("data", "sim", "predicted_links.jld"), "data", predicted_links)

# predicted number of links (median) for species richness in sp
L_med = convert.(Int64, round.(median.(eachcol(predicted_links))))

## Predict the adjacency matrix of maximum entropy for a range of species richness and median predicted number of links
i = findall(in(S).(sp))
Sp = sort(unique(S))

predicted_matrices = predict_adjacency_matrix.(Sp, L_med[i], 500)
predicted_networks = UnipartiteNetwork.(predicted_matrices)

save(joinpath("data", "sim", "predicted_networks.jld"), "data", predicted_networks)

## Predict the adjacency matrix of maximum entropy for networks with empirical numbers of links
predicted_matrices_empL = predict_adjacency_matrix.(S, L, 500)
predicted_networks_empL = UnipartiteNetwork.(predicted_matrices_empL)

save(joinpath("data", "sim", "predicted_networks_empL.jld"), "data", predicted_networks_empL)

## Predict the adjacency matrix for the entire range of worldwide species richness

# Minimum and maximum number of species (mammals, birds, and amphibians) in BiodiversityMapping maps
S_min = 5
S_max = 950

# Predict the adjacency matrix for this range of species richness (very long)
S_map = S_min:S_max
L_map = L_med[1:length(S_map)]

predicted_matrices_map = []
p = Progress(length(S_map))
Threads.@threads for i in 1:length(S_map)
    try
       N = predict_adjacency_matrix(S_map[i], L_map[i], 50)
       N = sparse(N)
       push!(predicted_matrices_map, N)
    catch
      println("could not simulate network $(i)")
    end
    next!(p)
end

# 5 networks were not predicted properly (i = 557, 573, 579, 603, 637)
i = vcat(557, 573, 579, 603, 637)
predicted_matrices_2ndtry = predict_adjacency_matrix.(i.+4, L_med[i.+4], 50)

for j in 1:length(i)
  insert!(predicted_matrices_map, i[j], predicted_matrices_2ndtry[j])
end

predicted_networks_map = UnipartiteNetwork.(predicted_matrices_map)

save(joinpath("data", "sim", "predicted_networks_map.jld"), "data", predicted_networks_map)

