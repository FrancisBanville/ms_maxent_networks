# Here we compute many measures of network structure, for all empirical and simulated networks 

#### Read data

## Empirical data 
N_all = load(joinpath("data", "proc", "N_all.jld"))["data"]
N_abund = load(joinpath("data", "proc", "N_abund.jld"))["data"]
A_abund = load(joinpath("data", "proc", "A_abund.jld"))["data"]

## Degree distributions of maximum entropy 
dd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "dd_maxent_all.jld"))["data"]
jdd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jdd_maxent_all.jld"))["data"]
jds_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"))["data"] 

## Networks of maximum entropy
N_ent_maxentco_all = load(joinpath("data", "sim", "network_maxent", "N_maxentco_all.jld"))["data"]
N_ent_maxentco_abund = load(joinpath("data", "sim", "network_maxent", "N_maxentco_abund.jld"))["data"]
N_ent_maxentjds_all = load(joinpath("data", "sim", "network_maxent", "N_maxentjds_all.jld"))["data"]
N_ent_maxentjds_abund = load(joinpath("data", "sim", "network_maxent", "N_maxentjds_abund.jld"))["data"]

# access unipartite networks of maximum entropy

N_maxentco_all = [] 
[push!(N_maxentco_all, N_ent_maxentco_all[i].A) for i in 1:length(N_ent_maxentco_all)]

N_maxentco_abund = [] 
[push!(N_maxentco_abund, N_ent_maxentco_abund[i].A) for i in 1:length(N_ent_maxentco_abund)]

N_maxentjds_all = [] 
[push!(N_maxentjds_all, N_ent_maxentjds_all[i].A) for i in 1:length(N_ent_maxentjds_all)]

N_maxentjds_abund = [] 
[push!(N_maxentjds_abund, N_ent_maxentjds_abund[i].A) for i in 1:length(N_ent_maxentjds_abund)]

## Neutral models 

N_neutral_abund = load(joinpath("data", "sim", "neutral_models", "N_neutral_abund.jld"))["data"]

N_nullco_all = load(joinpath("data", "sim", "neutral_models", "N_nullco_all.jld"))["data"]

N_nullco_abund = load(joinpath("data", "sim", "neutral_models", "N_nullco_abund.jld"))["data"]

N_nulljds_all = load(joinpath("data", "sim", "neutral_models", "N_nulljds_all.jld"))["data"]

N_nulljds_abund = load(joinpath("data", "sim", "neutral_models", "N_nulljds_abund.jld"))["data"]



# number of networks in each category (to be included in data frame)

n_all = length(N_all) 
n_abund = length(N_abund)

#### Table of network properties 
N_alls = vcat(N_all, 
            N_maxentco_all,
            N_maxentjds_all, 
            N_nullco_all, 
            N_nulljds_all)

N_abunds = vcat(N_abund, 
            N_maxentco_abund,
            N_maxentjds_abund, 
            N_nullco_abund, 
            N_nulljds_abund,
            N_neutral_abund)

Ns = vcat(N_alls, N_abunds)

# create dataframe for all measures 
names_all = vcat(fill("N_all", n_all), 
                fill("N_maxentco_all", n_all),
                fill("N_maxentjds_all", n_all),
                fill("N_nullco_all", n_all),
                fill("N_nulljds_all", n_all))

names_abund = vcat(fill("N_abund", n_abund), 
                fill("N_maxentco_abund", n_abund),
                fill("N_maxentjds_abund", n_abund),
                fill("N_nullco_abund", n_abund),
                fill("N_nulljds_abund", n_abund),
                fill("N_neutral_abund", n_abund))

namess = vcat(names_all, names_abund)

metrics = DataFrame(network = namess)


## species richness (S)
S_alls = richness.(Ns)
insertcols!(metrics, :S => S_alls)

## number of links (L)
L_alls = links.(Ns)
insertcols!(metrics, :L => L_alls)

## connectance (C)
C_alls = connectance.(Ns)
insertcols!(metrics, :C => C_alls)

## nestedness (rho)
rho_alls = ρ.(Ns)
insertcols!(metrics, :rho => rho_alls)

## maximum trophic level (maxtl)
maxtl_alls = Union{Missing, Float64}[]
for i in 1:length(Ns)
      try 
      maxtl = maximum(values(trophic_level(Ns[i])))
      push!(maxtl_alls, maxtl)
      catch
      push!(maxtl_alls, missing) # no maximum trophic level found
      end
end

insertcols!(metrics, :maxtl => maxtl_alls)

## network diameter - shortest distance between the two most distant nodes in the network (diam)
diam_alls = maximum.(shortest_path.(Ns))
insertcols!(metrics, :diam => diam_alls)

## SVD-entropy (entropy)
entropy_alls = svd_entropy.(Ns)
insertcols!(metrics, :entropy => entropy_alls)

## fraction of top species - species with no predators (T)
kin_alls = values.(degree.(Ns, dims = 2))
T_alls = sum.(x -> x == 0, kin_alls) ./ S_alls
insertcols!(metrics, :T => T_alls)

## fraction of basal species - species with no preys (B)
kout_alls = values.(degree.(Ns, dims = 1))
B_alls = sum.(x -> x == 0, kout_alls) ./ S_alls
insertcols!(metrics, :B => B_alls)

## fraction of intermediate species - species with preys and predators (I)
I_alls = 1 .- T_alls .- B_alls
I_alls[I_alls .< 0] .= 0
insertcols!(metrics, :I => I_alls)

## standard deviation of generality - number of preys normalized by link density (GenSD)
Gen_alls = [kout_alls[i] ./ (L_alls[i] ./ S_alls[i]) for i in 1:length(Ns)]
GenSD_alls = std.(Gen_alls)
insertcols!(metrics, :GenSD => GenSD_alls)

## standard deviation of vulnerability (VulSD)
Vul_alls = [kin_alls[i] ./ (L_alls[i] ./ S_alls[i]) for i in 1:length(Ns)]
VulSD_alls = std.(Vul_alls)
insertcols!(metrics, :VulSD => VulSD_alls)

## mean maximum similarity (MxSim)
MxSim_alls = MxSim.(Ns)
insertcols!(metrics, :MxSim => MxSim_alls)

## fraction of cannibal species (Cannib)
Cannib_alls = [sum(diag(Ns[i].edges)) for i in 1:length(Ns)] ./ S_alls
insertcols!(metrics, :Cannib => Cannib_alls)

## fraction of omnivorous species - species that consume two or more species and have food chains of different lengths (Omniv)
Omniv_alls = Union{Missing, Float64}[]
for i in 1:length(Ns)
    try 
        omniv = sum(values(omnivory(Ns[i])) .> 0) / richness(Ns[i])
        push!(Omniv_alls, omniv)
    catch
        push!(Omniv_alls, missing)
    end
end
insertcols!(metrics, :Omniv => Omniv_alls)

## motifs (S1, S2, S3, S4, S5, D1, D2, D3, D4, D5, D6, D7, D8)
motifs = keys(unipartitemotifs()) # motifs names
motifs_alls = count_motifs.(Ns) # count motifs
motifs_alls = mapreduce(permutedims, vcat, motifs_alls)

for i in 1:length(motifs)
    insertcols!(metrics, motifs[i] => motifs_alls[:,i])
    allowmissing!(metrics)
end


## save data frame
CSV.write(joinpath("results", "metrics.csv"), metrics)


#### Average and standard deviation of all measures by type of networks (i.e. dataset and empirical / maxent / neutral)

# summarizing functions skipping missing values
avg(x) = all(ismissing, x) ? missing : mean(skipmissing(x))
sd(x) = all(ismissing, x) ? missing : std(skipmissing(x))

metrics_all = metrics[in(names_all).(metrics.network),:]
metrics_abund = metrics[in(names_abund).(metrics.network),:]

gmetrics_all = groupby(metrics_all, :network)
gmetrics_all = combine(gmetrics_all, valuecols(gmetrics_all) .=> avg,
                                    valuecols(gmetrics_all) .=> sd)

CSV.write(joinpath("results", "gmetrics_all.csv"), gmetrics_all)


gmetrics_abund = groupby(metrics_abund, :network)
gmetrics_abund = combine(gmetrics_abund, valuecols(gmetrics_abund) .=> avg,
                                    valuecols(gmetrics_abund) .=> sd)

CSV.write(joinpath("results", "gmetrics_abund.csv"), gmetrics_abund)



#### Difference between empirical networks and those of MaxEnt

# empirical networks only
metrics_emp = metrics[metrics.network .== "N_all",:]

# MaxEnt networks only
metrics_maxentco = metrics[metrics.network .== "N_maxentco_all",:]
metrics_maxentjds = metrics[metrics.network .== "N_maxentjds_all",:]


## Divergence between the degree sequence of MaxEnt and empirical networks (MSD_ds)

# sorted degree sequence of empirical networks 
k_emp_all = joint_degree_seq.(N_all)

k_emp_all_sorted = [sort(k_emp_all[i].kin .+ k_emp_all[i].kout, rev=true) for i in 1:length(k_emp_all)]

# sorted degree sequence of MaxEnt networks
jds_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"))["data"]

k_maxent_all_sorted = [sort(jds_maxent_all[i].kin .+ jds_maxent_all[i].kout, rev=true) for i in 1:length(jds_maxent_all)]

# mean squared deviation between empirical and MaxEnt degree sequence 
MSD_ds_maxent = [mean(k_emp_all_sorted[i] .- k_maxent_all_sorted[i]).^2 for i in 1:length(k_maxent_all_sorted)]

## Create new data frame for all differences
metrics_diffco = DataFrame(MSD_ds = MSD_ds_maxent)
metrics_diffjds = DataFrame(MSD_ds = MSD_ds_maxent)

## Difference in SVD-entropy (entropy_diff)
entropy_all = svd_entropy.(N_all)
entropy_maxentco_all = svd_entropy.(N_maxentco_all)
entropy_maxentjds_all = svd_entropy.(N_maxentjds_all)

entropy_diffco= entropy_maxentco_all .- entropy_all
entropy_diffjds = entropy_maxentjds_all .- entropy_all

insertcols!(metrics_diffco, :entropy_diff => entropy_diffco)
insertcols!(metrics_diffjds, :entropy_diff => entropy_diffjds)


## Difference in nestedness (rho_diff)
rho_diffco = metrics_maxentco.rho .- metrics_emp.rho
rho_diffjds = metrics_maxentjds.rho .- metrics_emp.rho

insertcols!(metrics_diffco, :rho_diff => rho_diffco)
insertcols!(metrics_diffjds, :rho_diff => rho_diffjds)


## Jaccard distance (jaccard)

A_emp = [vec(Matrix(N_all[i].edges)) for i in 1:length(N_all)]

A_maxentco = [vec(Matrix(N_maxentco_all[i].edges)) for i in 1:length(N_maxentco_all)]

A_maxentjds = [vec(Matrix(N_maxentjds_all[i].edges)) for i in 1:length(N_maxentjds_all)]

jaccard_maxentco = jaccard.(A_emp, A_maxentco)
jaccard_maxentjds = jaccard.(A_emp, A_maxentjds)

insertcols!(metrics_diffco, :jaccard => jaccard_maxentco)
insertcols!(metrics_diffjds, :jaccard => jaccard_maxentjds)


CSV.write(joinpath("results", "metrics_diffco.csv"), metrics_diffco)
CSV.write(joinpath("results", "metrics_diffjds.csv"), metrics_diffjds)




