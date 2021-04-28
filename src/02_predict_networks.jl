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
nsp = length(sp)
# number of simulations
n = 1000 

predicted_links = zeros(Int64, n, length(sp))

p = Progress(length(sp))
Threads.@threads for s in sp
    i = findall(sp .== s)
    predicted_links[:, i] = predict_links(s, n, chain)
    next!(p)
end


## Predict the adjacency matrix of maximum entropy

# Example

# number of species
S = 50
# predicted number of links (take median)
L = predict_links(S, n, chain)
L = convert(Int64, round(median(L)))
# degree distribution of maximum entropy (probabilities)
dd_maxent_prob(S, L)
# degree distribution of maximum entropy (counts)
dd = dd_maxent(S, L)
# joint-degree distribution of maximum entropy (probabilities)
jdd_maxent_prob(S, L)
# in-degree distribution of maximum entropy (counts)
indd = indd_maxent(S, L)
# out-degree distribution of maximum entropy (non-adjusted counts)
outdd = dd.- indd
m = length(dd)
outdd[m] = last(outdd) + (sum(indd) - sum(outdd))
# first adjacency matrix constrained by in-degree and out-degree distributions
A = matrix_rowcolsum(indd, outdd)
# adjacency matrix of maximum entropy
A = adjacency_matrix_maxent(A, 50)
# predict adjacency matrix of maximum entropy in one command
B = predict_adjacency_matrix(S, L, 50)


