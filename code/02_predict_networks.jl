# Here we predict food webs using different models (flexible links, MaxEnt, and neutral models)

Random.seed!(123)

## Read food webs and convert them to UnipartiteNetworks

# food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl) 
mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "raw", "mangal", "mangal_foodwebs.csv")))
N_mangal = network.(mangal_foodwebs.id)
N_mangal = convert.(UnipartiteNetwork, N_mangal)
N_mangal = UnipartiteNetwork.(N_mangal[i].edges for i in 1:length(N_mangal))

save(joinpath("data", "proc", "mangal", "networks_mangal.jld"), "data", N_mangal)

# New-Zealand food webs
NZ_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "new_zealand", "adjacency_matrices")),
                                        drop=[1])) # the first column is row names
N_NZ = convert.(Matrix{Bool}, NZ_foodwebs)
N_NZ = UnipartiteNetwork.(N_NZ, names.(NZ_foodwebs))

save(joinpath("data", "proc", "new_zealand", "networks_NZ.jld"), "data", N_NZ)

# Tuesday lake food webs
tuesday_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "tuesday_lake", "adjacency_matrices")),
                                        header=false)) # no column names here
N_tuesday = convert.(Matrix{Bool}, tuesday_foodwebs)
N_tuesday = UnipartiteNetwork.(N_tuesday)

save(joinpath("data", "proc", "tuesday_lake", "networks_tuesday.jld"), "data", N_tuesday)

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
N_all = vcat(N_mangal, N_NZ, N_tuesday)
N_all = simplify.(N_all) # remove unconnected species

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
dd_maxent_mangal = degree_dist_maxent.(S_mangal, L_mangal)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_mangal.jld"), "data", dd_maxent_mangal)

dd_maxent_NZ = degree_dist_maxent.(S_NZ, L_NZ)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_NZ.jld"), "data", dd_maxent_NZ)

dd_maxent_tuesday = degree_dist_maxent.(S_tuesday, L_tuesday)
save(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_tuesday.jld"), "data", dd_maxent_tuesday)


## Simulate joint degree distributions of maximum entropy
jdd_maxent_mangal = joint_degree_dist_maxent.(S_mangal, L_mangal)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_mangal.jld"), "data", jdd_maxent_mangal)

jdd_maxent_NZ = joint_degree_dist_maxent.(S_NZ, L_NZ)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_NZ.jld"), "data", jdd_maxent_NZ)

jdd_maxent_tuesday = joint_degree_dist_maxent.(S_tuesday, L_tuesday)
save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_tuesday.jld"), "data", jdd_maxent_tuesday)

## Simulated degree sequence from the joint degree distribution of maximum entropy
jdd_maxent_all = vcat(jdd_maxent_mangal, jdd_maxent_NZ, jdd_maxent_tuesday)
k_maxent_all = simulate_degrees.(jdd_maxent_all)

save(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_sequence_all.jld"), "data", k_maxent_all)


## Simulate networks of maximum entropy
nsteps = 2000 # number of steps
nchains = 4 # number of chains

N_maxent_mangal = network_maxent.(N_mangal, nchains, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_mangal.jld"), "data", N_maxent_mangal)

N_maxent_NZ = network_maxent.(N_NZ, nchains, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_NZ.jld"), "data", N_maxent_NZ)

N_maxent_tuesday = network_maxent.(N_tuesday, nchains, nsteps)
save(joinpath("data", "sim", "network_maxent", "network_maxent_tuesday.jld"), "data", N_maxent_tuesday)

## Run neutral models

## New Zealand data

# read taxonomic and abundance data in New Zealand (all networks)
abund_NZ = DataFrame(CSV.File(joinpath("data", "raw", "new_zealand", "taxa_dry_weight_abundance.csv")))
rename!(abund_NZ , 1 => :food_web, 4 => :no_m2)
fw_names = unique(abund_NZ[!, :food_web])

# get the name of all food webs (adjacency matrices) in New Zealand
matrix_names = readdir(joinpath("data", "raw", "new_zealand", "adjacency_matrices"))
matrix_names = replace.(matrix_names, ".csv" => "")

"""
abundance_data_NZ(fw_name::String) 
    fw_name: name of the food web (file name)
Returns the density (no/m2) of all species in the simplified food web with abundance data
"""
function abundance_data_NZ(fw_name::String)
  # filter abundance data for this food web and get species names
  abund_fw = abund_NZ[abund_NZ[!, :food_web] .== fw_name, :]
  abund_taxa = abund_fw[!, :taxa]
  # read food web data and get species names
  N_df = DataFrame.(CSV.File.(joinpath("data", "raw", "new_zealand", "adjacency_matrices", "$fw_name.csv"),
                    drop=[1])) # the first column is row names
  # convert to UnipartiteNetwork and make sure species have abundance data
  N = convert(Matrix{Bool}, N_df)
  N = UnipartiteNetwork(N, names(N_df))
  taxa_N = species(N)
  N.edges = N.edges[in(abund_taxa).(taxa_N), in(abund_taxa).(taxa_N)]
  N.S = taxa_N[in(abund_taxa).(taxa_N)]
  # simplify network and get species names
  simplify!(N)
  taxa_N = species(N)
  # get abundance data for all taxa in network
  abund = abund_fw[in(taxa_N).(abund_taxa), :no_m2]
  # put everything together
  network_abund = (name = fw_name, network = N, abundance = abund)
  return network_abund
end

# get the abundance data of all food webs in New Zealand
abund_NZ = abundance_data_NZ.(fw_names)
save(joinpath("data", "proc", "new_zealand", "abund_data_NZ.jld"), "data", abund_NZ)

# simulate networks in New Zealand using the neutral model of relative abundances
N_neutral_NZ = Any[]

for i in 1:length(abund_NZ)
  N_neutral = neutral_model(abund_NZ[i].abundance, chain)
  push!(N_neutral_NZ, N_neutral)
end

save(joinpath("data", "sim", "neutral_model", "neutral_networks_NZ.jld"), "data", N_neutral_NZ)


## Tuesday lake data

# read taxonomic and abundance data of Tuesday lake (all networks)
abund_tuesday = DataFrame.(CSV.File.(glob("*.csv", 
                                      joinpath("data", "raw", "tuesday_lake", "abundances"))))

# sum the numerical abundance of species in each trophic species cluster 
"""
abundance_data_tuesday(fw_name::String) 
    fw_name: name of the food web (file name)
Returns the total abundance of each trophic species in the food web 
"""
function abundance_data_tuesday(df::DataFrame, year::Int64, N::UnipartiteNetwork)
  df = df[df[!, :trophic_species] .!= " NA", :] # remove unknown trophic species
  
  # group each species by their trophic species and add their numerical abundances
  trophic_sp_tuesday = groupby(df, :trophic_species)
  trophic_sp_tuesday = combine(trophic_sp_tuesday, :numerical_abundance => sum)
  abund = vec(trophic_sp_tuesday[:, :numerical_abundance_sum])
  
  network_abund = (name = year, network = N, abundance = abund)
  return network_abund
end

# get the abundance data of all food webs of Tuesday lake
year = (1984, 1986)
abund_tuesday = abundance_data_tuesday.(abund_tuesday, year, N_tuesday)
save(joinpath("data", "proc", "tuesday_lake", "abund_data_tuesday.jld"), "data", abund_tuesday)

# simulate networks of Tuesday lake using the neutral model of relative abundances
N_neutral_tuesday = Any[]

for i in 1:length(abund_tuesday)
  neutral_network = neutral_model(abund_tuesday[i].abundance, chain)
  push!(N_neutral_tuesday, neutral_network)
end

save(joinpath("data", "sim", "neutral_model", "neutral_networks_tuesday.jld"), "data", N_neutral_tuesday)

