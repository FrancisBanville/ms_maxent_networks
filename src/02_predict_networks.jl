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
L_med = convert.(Int64, round.(median.(eachcol(predicted_links))))

## Predict the adjacency matrix of maximum entropy for a range of species richness and median predicted number of links
i = findall(in(S).(sp))
Sp = sort(unique(S))

predicted_matrices = predict_adjacency_matrix.(Sp, L_med[i], 100)
predicted_networks = UnipartiteNetwork.(predicted_matrices)

save(joinpath("data", "sim", "predicted_networks.jld"), "data", predicted_networks)

## Predict the adjacency matrix of maximum entropy for networks with empirical numbers of links
predicted_matrices_empL = predict_adjacency_matrix.(S, L, 100)
predicted_networks_empL = UnipartiteNetwork.(predicted_matrices_empL)

save(joinpath("data", "sim", "predicted_networks_empL.jld"), "data", predicted_networks_empL)

