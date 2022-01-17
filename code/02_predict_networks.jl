# Here we predict food webs using different models (flexible links, MaxEnt, and neutral models)

## Read food webs and convert them to UnipartiteNetworks

# food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl) 
mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "raw", "mangal", "mangal_foodwebs.csv")))
N_mangal = network.(mangal_foodwebs.id)
N_mangal = convert.(UnipartiteNetwork, N)

save(joinpath("data", "raw", "mangal", "network_mangal.jld"), "data", N_mangal)

# New-Zealand food webs
NZ_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "new_zealand", "adjacency_matrices")),
                                        drop=[1])) # the first column is row names
N_NZ = convert.(Matrix{Bool}, NZ_foodwebs)
N_NZ = UnipartiteNetwork.(N_NZ)

save(joinpath("data", "raw", "new_zealand", "network_NZ.jld"), "data", N_NZ)

# Tuesday lake food webs
tuesday_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "tuesday_lake", "adjacency_matrices")),
                                        header=false)) # no column names here
N_tuesday = convert.(Matrix{Bool}, tuesday_foodwebs)
N_tuesday = UnipartiteNetwork.(N_tuesday)

save(joinpath("data", "raw", "tuesday_lake", "network_tuesday.jld"), "data", N_tuesday)

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
N = vcat(N_mangal, N_NZ, N_tuesday)
N = simplify.(N) # remove unconnected species

# number of species
S = richness.(N)
# number of links
L = links.(N)
# number of flexible links that are realized
R = L .- (S.-1) 

# Hamiltonian Monte Carlo sampler with static trajectory
chain = sample(FL(S,R), HMC(0.01,10), 3000)

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
N_mangal = simplify.(N_mangal)
S_mangal = richness.(N_mangal)
L_mangal = links.(N_mangal)

N_NZ = simplify.(N_NZ)
S_NZ = richness.(N_NZ)
L_NZ = links.(N_NZ)

N_tuesday = simplify.(N_tuesday)
S_tuesday = richness.(N_tuesday)
L_tuesday = links.(N_tuesday)


## Simulate degree distributions of maximum entropy
degree_dist_mangal = degree_dist_maxent.(S_mangal, L_mangal)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_mangal.jld"), "data", degree_dist_mangal)

degree_dist_NZ = degree_dist_maxent.(S_NZ, L_NZ)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_NZ.jld"), "data", degree_dist_NZ)

degree_dist_tuesday = degree_dist_maxent.(S_tuesday, L_tuesday)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_tuesday.jld"), "data", degree_dist_tuesday)


## Simulate joint degree distributions of maximum entropy
joint_degree_dist_mangal = joint_degree_dist_maxent.(S_mangal, L_mangal)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_mangal.jld"), "data", joint_degree_dist_mangal)

joint_degree_dist_NZ = joint_degree_dist_maxent.(S_NZ, L_NZ)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_NZ.jld"), "data", joint_degree_dist_NZ)

joint_degree_dist_tuesday = joint_degree_dist_maxent.(S_tuesday, L_tuesday)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_tuesday.jld"), "data", joint_degree_dist_tuesday)


## Simulate networks of maximum entropy
nsteps = 1000 # number of steps

network_maxent_mangal = network_maxent.(N_mangal, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_mangal.jld"), "data", network_maxent_mangal)

network_maxent_NZ = network_maxent.(N_NZ, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_NZ.jld"), "data", network_maxent_NZ)

network_maxent_tuesday = network_maxent.(N_tuesday, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_tuesday.jld"), "data", network_maxent_tuesday)



## Run neutral models

