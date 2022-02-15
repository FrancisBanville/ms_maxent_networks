# Here we predict food webs using different models (flexible links, MaxEnt, and neutral models)

Random.seed!(123)

# Read processed data
N_all = load(joinpath("data", "proc", "N_all.jld"))["data"] # all empirical networks
N_abund = load(joinpath("data", "proc", "N_abund.jld"))["data"] # networks with abundance data
A_abund = load(joinpath("data", "proc", "A_abund.jld"))["data"] # abundance data

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

# number of species
S_all = richness.(N_all)
# number of links
L_all = links.(N_all)
# number of flexible links that are realized
R_all = L_all .- (S_all.-1) 

# Hamiltonian Monte Carlo sampler with static trajectory
chain = sample(FL(S_all, R_all), HMC(0.01,10), 3000)

# Diagnostic plot
plot(chain[200:end,:,:])


## Simulate numbers of links

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


## Species richness and numbers of links to be used in subsequent simulations
S_all = richness.(N_all)
L_all = links.(N_all)

## Simulate degree distributions of maximum entropy
dd_maxent_all = degree_dist_maxent.(S_all, L_all)

save(joinpath("data", "sim", "degree_dist_maxent", "dd_maxent_all.jld"), "data", dd_maxent_all)


## Simulate joint degree distributions of maximum entropy
jdd_maxent_all = joint_degree_dist_maxent.(S_all, L_all)

save(joinpath("data", "sim", "degree_dist_maxent", "jdd_maxent_all.jld"), "data", jdd_maxent_all)


## Simulated degree sequence from the joint degree distribution of maximum entropy
jds_maxent_all = simulate_degrees.(jdd_maxent_all)

save(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"), "data", jds_maxent_all)


## Simulate networks of maximum entropy (constrained by joint degree sequence)
nsteps = 2000 # number of steps
nchains = 4 # number of chains

# all networks 
N_maxentjds_all = []
p = Progress(length(N_all))
for i in 1:length(N_all)
  push!(N_maxentjds_all, network_maxent(N_all[i], "jds", nchains, nsteps))
  next!(p)
end

save(joinpath("data", "sim", "network_maxent", "N_maxentjds_all.jld"), "data", N_maxentjds_all)

# networks with abundance data
N_maxentjds_abund = []
p = Progress(length(N_abund))
for i in 1:length(N_abund)
  push!(N_maxentjds_abund, network_maxent(N_abund[i], "jds", nchains, nsteps))
  next!(p)
end

save(joinpath("data", "sim", "network_maxent", "N_maxentjds_abund.jld"), "data", N_maxentjds_abund)


## Simulate networks of maximum entropy (constrained by connectance)

# all networks 
N_maxentco_all = []
p = Progress(length(N_all))
for i in 1:length(N_all)
  push!(N_maxentco_all, network_maxent(N_all[i], "co", nchains, nsteps))
  next!(p)
end

save(joinpath("data", "sim", "network_maxent", "N_maxentco_all.jld"), "data", N_maxentco_all)

# networks with abundance data
N_maxentco_abund = []
p = Progress(length(N_abund))
for i in 1:length(N_abund)
  push!(N_maxentco_abund, network_maxent(N_abund[i], "co", nchains, nsteps))
  next!(p)
end

save(joinpath("data", "sim", "network_maxent", "N_maxentco_abund.jld"), "data", N_maxentco_abund)



## Simulate networks using the neutral model of relative abundances 

L_abund = links.(N_abund)
N_neutral_abund = neutral_model.(A_abund, L_abund)

save(joinpath("data", "sim", "neutral_models", "N_neutral_abund.jld"), "data", N_neutral_abund)



## Simulate networks using the type I null model of Fortuna & Bascompte, 2006

# all networks
N_nullco_all = null_model.(N_all, L_all, 1)

save(joinpath("data", "sim", "neutral_models", "N_nullco_all.jld"), "data", N_nullco_all)

# networks with abundance data
N_nullco_abund = null_model.(N_abund, L_abund, 1)

save(joinpath("data", "sim", "neutral_models", "N_nullco_abund.jld"), "data", N_nullco_abund)



## Simulate networks using the type I null model of Fortuna & Bascompte, 2006

# all networks
N_nulljds_all = null_model.(N_all, L_all, 2)

save(joinpath("data", "sim", "neutral_models", "N_nulljds_all.jld"), "data", N_nulljds_all)

# networks with abundance data
N_nulljds_abund = null_model.(N_abund, L_abund, 2)

save(joinpath("data", "sim", "neutral_models", "N_nulljds_abund.jld"), "data", N_nulljds_abund)