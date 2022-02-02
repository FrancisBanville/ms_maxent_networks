# Here we compute many measures of network structure, for all empirical and simulated networks 

#### Read data

## Empirical data 

# Mangal
N_mangal = load(joinpath("data", "proc", "mangal", "networks_mangal.jld"))["data"]

# New Zealand
N_NZ = load(joinpath("data", "proc", "new_zealand", "networks_NZ.jld"))["data"]
abund_NZ = load(joinpath("data", "proc", "new_zealand", "abund_data_NZ.jld"))["data"]

# Tuesday lake
N_tuesday = load(joinpath("data", "proc", "tuesday_lake", "networks_tuesday.jld"))["data"]
abund_tuesday = load(joinpath("data", "proc", "tuesday_lake", "abund_data_tuesday.jld"))["data"]

N_emp = simplify.(vcat(N_mangal, N_NZ, N_tuesday))

## Simulated data (maximum entropy)

# degree distribution
dd_maxent_mangal = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_mangal.jld"))["data"]
dd_maxent_NZ = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_NZ.jld"))["data"]
dd_maxent_tuesday = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_tuesday.jld"))["data"]

# joint degree distribution
jdd_maxent_mangal = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_mangal.jld"))["data"]
jdd_maxent_NZ = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_NZ.jld"))["data"]
jdd_maxent_tuesday = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_tuesday.jld"))["data"]

# network of maximum entropy
N_entropies_maxent_mangal = load(joinpath("data", "sim", "network_maxent", "network_maxent_mangal.jld"))["data"]
N_entropies_maxent_NZ = load(joinpath("data", "sim", "network_maxent", "network_maxent_NZ.jld"))["data"]
N_entropies_maxent_tuesday = load(joinpath("data", "sim", "network_maxent", "network_maxent_tuesday.jld"))["data"]

# access unipartite networks of maximum entropy

N_maxent_mangal = [] 
[push!(N_maxent_mangal, N_entropies_maxent_mangal[i].A) for i in 1:length(N_entropies_maxent_mangal)]

N_maxent_NZ = [] 
[push!(N_maxent_NZ, N_entropies_maxent_NZ[i].A) for i in 1:length(N_entropies_maxent_NZ)]

N_maxent_tuesday = [] 
[push!(N_maxent_tuesday, N_entropies_maxent_tuesday[i].A) for i in 1:length(N_entropies_maxent_tuesday)]

N_maxent = simplify.(vcat(N_maxent_mangal, N_maxent_NZ, N_maxent_tuesday))

## Simulated data (neutral models)

N_neutral_NZ = load(joinpath("data", "sim", "neutral_model", "neutral_networks_NZ.jld"))["data"]
N_neutral_tuesday = load(joinpath("data", "sim", "neutral_model", "neutral_networks_tuesday.jld"))["data"]


#### Join all empirical and simulated networks together and simplify them 

N_all = vcat(N_mangal, 
            N_maxent_mangal, 
            N_NZ, 
            N_maxent_NZ, 
            N_neutral_NZ, 
            N_tuesday, 
            N_maxent_tuesday, 
            N_neutral_tuesday)

N_all = simplify.(N_all)

# number of networks in each category (to be included in data frame)

n_mangal = length(N_mangal) 
n_NZ = length(N_NZ)
n_neutral_NZ = length(N_neutral_NZ) # not all networks in New Zealand have abundance data
n_tuesday = length(N_tuesday)


#### Table of network properties 

# create dataframe for all measures 
names_all = vcat(fill("N_mangal", n_mangal), 
                 fill("N_maxent_mangal", n_mangal),
                 fill("N_NZ", n_NZ),
                 fill("N_maxent_NZ", n_NZ),
                 fill("N_neutral_NZ", n_neutral_NZ),
                 fill("N_tuesday", n_tuesday),
                 fill("N_maxent_tuesday", n_tuesday),
                 fill("N_neutral_tuesday", n_tuesday))

metrics = DataFrame(network = names_all)


## species richness (S)
S_all = richness.(N_all)
insertcols!(metrics, :S => S_all)

## number of links (L)
L_all = links.(N_all)
insertcols!(metrics, :L => L_all)

## connectance (C)
C_all = connectance.(N_all)
insertcols!(metrics, :C => C_all)

## nestedness (rho)
rho_all = ρ.(N_all)
insertcols!(metrics, :rho => rho_all)

## maximum trophic level (maxtl)
maxtl_all = Union{Missing, Float64}[]
for i in 1:length(N_all)
      try 
      maxtl = maximum(values(trophic_level(N_all[i])))
      push!(maxtl_all, maxtl)
      catch
      push!(maxtl_all, missing) # no maximum trophic level found
      end
end

insertcols!(metrics, :maxtl => maxtl_all)

## network diameter - shortest distance between the two most distant nodes in the network (diam)
diam_all = maximum.(shortest_path.(N_all))
insertcols!(metrics, :diam => diam_all)

## SVD-entropy (entropy)
entropy_all = svd_entropy.(N_all)
insertcols!(metrics, :entropy => entropy_all)

## fraction of top species - species with no predators (T)
kin_all = values.(degree.(N_all, dims = 2))
T_all = sum.(x -> x == 0, kin_all) ./ S_all
insertcols!(metrics, :T => T_all)

## fraction of basal species - species with no preys (B)
kout_all = values.(degree.(N_all, dims = 1))
B_all = sum.(x -> x == 0, kout_all) ./ S_all
insertcols!(metrics, :B => B_all)

## fraction of intermediate species - species with preys and predators (I)
I_all = 1 .- T_all .- B_all
I_all[I_all .< 0] .= 0
insertcols!(metrics, :I => I_all)

## standard deviation of generality - number of preys normalized by link density (GenSD)
Gen_all = [kout_all[i] ./ (L_all[i] ./ S_all[i]) for i in 1:length(N_all)]
GenSD_all = std.(Gen_all)
insertcols!(metrics, :GenSD => GenSD_all)

## standard deviation of vulnerability (VulSD)
Vul_all = [kin_all[i] ./ (L_all[i] ./ S_all[i]) for i in 1:length(N_all)]
VulSD_all = std.(Vul_all)
insertcols!(metrics, :VulSD => VulSD_all)

## mean maximum similarity (MxSim)
MxSim_all = MxSim.(N_all)
insertcols!(metrics, :MxSim => MxSim_all)

## fraction of cannibal species (Cannib)
Cannib_all = [sum(diag(N_all[i].edges)) for i in 1:length(N_all)] ./ S_all
insertcols!(metrics, :Cannib => Cannib_all)

## fraction of omnivorous species - species that consume two or more species and have food chains of different lengths (Omniv)
Omniv_all = Union{Missing, Float64}[]
for i in 1:length(N_all)
    try 
        omniv = sum(values(omnivory(N_all[i])) .> 0) / richness(N_all[i])
        push!(Omniv_all, omniv)
    catch
        push!(Omniv_all, missing)
    end
end
insertcols!(metrics, :Omniv => Omniv_all)

## motifs (S1, S2, S3, S4, S5, D1, D2, D3, D4, D5, D6, D7, D8)
motifs = keys(unipartitemotifs()) # motifs names
motifs_all = count_motifs.(N_all) # count motifs
motifs_all = mapreduce(permutedims, vcat, motifs_all)

large_N = findall(S_all .> 500)

for i in 1:length(motifs)
    insertcols!(metrics, motifs[i] => motifs_all[:,i])
    allowmissing!(metrics)
    metrics[large_N, motifs[i]] .= missing # we didn't find the motifs of large networks 
end


## save data frame
CSV.write(joinpath("results", "metrics.csv"), metrics)


#### Average and standard deviation of all measures by type of networks (i.e. dataset and empirical / maxent / neutral)

# summarizing functions skipping missing values
avg(x) = all(ismissing, x) ? missing : mean(skipmissing(x))
sd(x) = all(ismissing, x) ? missing : std(skipmissing(x))

gmetrics = groupby(metrics, :network)
gmetrics = combine(gmetrics, valuecols(gmetrics) .=> avg,
                             valuecols(gmetrics) .=> sd)

CSV.write(joinpath("results", "gmetrics.csv"), gmetrics)



#### Difference between empirical networks and those of MaxEnt

## Divergence between the degree sequence of MaxEnt and empirical networks (MSD_ds)

# sorted degree sequence of empirical networks 
k_emp_all = joint_degree_seq.(N_emp)

k_emp_all_sorted = [sort(k_emp_all[i].kin .+ k_emp_all[i].kout, rev=true) for i in 1:length(k_emp_all)]

# sorted degree sequence of MaxEnt networks
k_maxent_all = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_sequence_all.jld"))["data"]

k_maxent_all_sorted = [sort(k_maxent_all[i].kin .+ k_maxent_all[i].kout, rev=true) for i in 1:length(k_maxent_all)]

# mean squared deviation between empirical and MaxEnt degree sequence 
MSD_ds_maxent = [mean(k_emp_all_sorted[i] .- k_maxent_all_sorted[i]).^2 for i in 1:length(k_maxent_all_sorted)]

## Create new data frame for all differences
metrics_diff = DataFrame(MSD_ds = MSD_maxent)



## Difference in SVD-entropy (entropy_diff)
entropy_emp = svd_entropy.(N_emp)
entropy_maxent = svd_entropy.(N_maxent)

entropy_diff = entropy_maxent .- entropy_emp

insertcols!(metrics_diff, :entropy_diff => entropy_diff)


CSV.write(joinpath("results", "metrics_diff.csv"), metrics_diff)




