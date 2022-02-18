# Here we compute many measures of network structure, for all empirical and simulated networks 

#### Read data

## Empirical data 
N_all = load(joinpath("data", "proc", "N_all.jld"))["data"]
N_abund = load(joinpath("data", "proc", "N_abund.jld"))["data"]

## Networks of maximum entropy
N_ent_maxentco_all = load(joinpath("data", "sim", "network_maxent", "N_maxentco_all.jld"))["data"]
N_ent_maxentco_abund = load(joinpath("data", "sim", "network_maxent", "N_maxentco_abund.jld"))["data"]
N_ent_maxentjds_all = load(joinpath("data", "sim", "network_maxent", "N_maxentjds_all.jld"))["data"]
N_ent_maxentjds_abund = load(joinpath("data", "sim", "network_maxent", "N_maxentjds_abund.jld"))["data"]

# access unipartite networks of maximum entropy
N_maxentco_all = [N_ent_maxentco_all[i].A for i in 1:length(N_ent_maxentco_all)]

N_maxentco_abund = [N_ent_maxentco_abund[i].A for i in 1:length(N_ent_maxentco_abund)]

N_maxentjds_all = [N_ent_maxentjds_all[i].A for i in 1:length(N_ent_maxentjds_all)]

N_maxentjds_abund = [N_ent_maxentjds_abund[i].A for i in 1:length(N_ent_maxentjds_abund)]

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
Ns_all = vcat(N_all, 
            N_maxentco_all,
            N_maxentjds_all, 
            N_nullco_all, 
            N_nulljds_all)

Ns_abund = vcat(N_abund, 
            N_maxentco_abund,
            N_maxentjds_abund, 
            N_nullco_abund, 
            N_nulljds_abund,
            N_neutral_abund)

Ns = vcat(Ns_all, Ns_abund)

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

measures = DataFrame(network = vcat(names_all, names_abund))


## species richness (S)
S_Ns = richness.(Ns)
insertcols!(measures, :S => S_Ns)

## number of links (L)
L_Ns = links.(Ns)
insertcols!(measures, :L => L_Ns)

## connectance (C)
C_Ns = connectance.(Ns)
insertcols!(measures, :C => C_Ns)

## nestedness (rho)
rho_Ns = ρ.(Ns)
insertcols!(measures, :rho => rho_Ns)

## maximum trophic level (maxtl)
maxtl_Ns = Union{Missing, Float64}[]
for i in 1:length(Ns)
      try 
      maxtl = maximum(values(trophic_level(Ns[i])))
      push!(maxtl_Ns, maxtl)
      catch
      push!(maxtl_Ns, missing) # no maximum trophic level found
      end
end

insertcols!(measures, :maxtl => maxtl_Ns)

## network diameter - shortest distance between the two most distant nodes in the network (diam)
diam_Ns = maximum.(shortest_path.(Ns))
insertcols!(measures, :diam => diam_Ns)

## SVD-entropy (entropy)
entropy_Ns = svd_entropy.(Ns)
insertcols!(measures, :entropy => entropy_Ns)

## fraction of top species - species with no predators (T)
kin_Ns = values.(degree.(Ns, dims = 2))
T_Ns = sum.(x -> x == 0, kin_Ns) ./ S_Ns
insertcols!(measures, :T => T_Ns)

## fraction of basal species - species with no preys (B)
kout_Ns = values.(degree.(Ns, dims = 1))
B_Ns = sum.(x -> x == 0, kout_Ns) ./ S_Ns
insertcols!(measures, :B => B_Ns)

## fraction of intermediate species - species with preys and predators (I)
I_Ns = 1 .- T_Ns .- B_Ns
I_Ns[I_Ns .< 0] .= 0
insertcols!(measures, :I => I_Ns)

## standard deviation of generality - number of preys normalized by link density (GenSD)
Gen_Ns = [kout_Ns[i] ./ (L_Ns[i] ./ S_Ns[i]) for i in 1:length(Ns)]
GenSD_Ns = std.(Gen_Ns)
insertcols!(measures, :GenSD => GenSD_Ns)

## standard deviation of vulnerability (VulSD)
Vul_Ns = [kin_Ns[i] ./ (L_Ns[i] ./ S_Ns[i]) for i in 1:length(Ns)]
VulSD_Ns = std.(Vul_Ns)
insertcols!(measures, :VulSD => VulSD_Ns)

## mean maximum similarity (MxSim)
MxSim_Ns = MxSim.(Ns)
insertcols!(measures, :MxSim => MxSim_Ns)

## fraction of cannibal species (Cannib)
Cannib_Ns = [sum(diag(Ns[i].edges)) for i in 1:length(Ns)] ./ S_Ns
insertcols!(measures, :Cannib => Cannib_Ns)

## fraction of omnivorous species - species that consume two or more species and have food chains of different lengths (Omniv)
Omniv_Ns = Union{Missing, Float64}[]
for i in 1:length(Ns)
    try 
        omniv = sum(values(omnivory(Ns[i])) .> 0) / richness(Ns[i])
        push!(Omniv_Ns, omniv)
    catch
        push!(Omniv_Ns, missing)
    end
end
insertcols!(measures, :Omniv => Omniv_Ns)

## motifs (S1, S2, S3, S4, S5, D1, D2, D3, D4, D5, D6, D7, D8)
motifs = keys(unipartitemotifs()) # motifs names
motifs_Ns = count_motifs.(Ns) # count motifs
motifs_Ns = mapreduce(permutedims, vcat, motifs_Ns)

for i in 1:length(motifs)
    insertcols!(measures, motifs[i] => motifs_Ns[:,i])
    allowmissing!(measures)
end


## save data frame
CSV.write(joinpath("results", "measures.csv"), measures)




#### Average and standard deviation of all measures by type of networks (i.e. dataset and empirical / maxent / neutral)

# summarizing functions skipping missing values
replace!.(eachcol(measures), NaN => missing)

avg(x) = all(ismissing, x) ? missing : mean(skipmissing(x))
sd(x) = all(ismissing, x) ? missing : std(skipmissing(x))

measures_all = measures[in(names_all).(measures.network),:]
measures_abund = measures[in(names_abund).(measures.network),:]

# all networks 
gmeasures_all = groupby(measures_all, :network)
gmeasures_all = combine(gmeasures_all, valuecols(gmeasures_all) .=> avg,
                                    valuecols(gmeasures_all) .=> sd)

CSV.write(joinpath("results", "gmeasures_all.csv"), gmeasures_all)

# networks with abundance data
gmeasures_abund = groupby(measures_abund, :network)
gmeasures_abund = combine(gmeasures_abund, valuecols(gmeasures_abund) .=> avg,
                                    valuecols(gmeasures_abund) .=> sd)

CSV.write(joinpath("results", "gmeasures_abund.csv"), gmeasures_abund)



#### Standardized mean difference between models and empirical measures 

# all networks 

# select measures to be included in manuscript 
gmeasures_all_subset = gmeasures_all[:,vcat(
                                "rho_avg", 
                                "maxtl_avg", 
                                "diam_avg", 
                                "MxSim_avg", 
                                "Cannib_avg", 
                                "Omniv_avg", 
                                "entropy_avg")]

rename!(gmeasures_all_subset, [:rho, :maxtl, :diam, :MxSim, :Cannib, :Omniv, :entropy])

# compute standardized mean difference 
gmeasures_N_all = gmeasures_all_subset[gmeasures_all.network .== "N_all",:]

gmeasures_models_all = gmeasures_all_subset[Not(gmeasures_all.network .== "N_all"),:]

gmeasures_diff_all = (gmeasures_models_all .- gmeasures_N_all) ./ gmeasures_N_all

# add model names and reorder rows 
insertcols!(gmeasures_diff_all, 1, :model => vcat("MaxEnt-co", "MaxEnt-jds", "null 1", "null 2"))

gmeasures_diff_all = gmeasures_diff_all[vcat(3,1,4,2),:]

# format as Markdown table
table_all = latexify(gmeasures_diff_all, env=:mdtable, fmt="%.3f", latex=false, escape_underscores=true)

# export to file
table_all_path = joinpath("tables", "measures_all.md")
open(table_all_path, "w") do io
    print(io, table_all)
end


# networks with abundance data 

# select measures to be included in manuscript 
gmeasures_abund_subset = gmeasures_abund[:,vcat(
                                "rho_avg", 
                                "maxtl_avg", 
                                "diam_avg", 
                                "MxSim_avg", 
                                "Cannib_avg", 
                                "Omniv_avg", 
                                "entropy_avg")]

rename!(gmeasures_abund_subset, [:rho, :maxtl, :diam, :MxSim, :Cannib, :Omniv, :entropy])

# compute standardized mean difference 
gmeasures_N_abund = gmeasures_abund_subset[gmeasures_abund.network .== "N_abund",:]

gmeasures_models_abund = gmeasures_abund_subset[Not(gmeasures_abund.network .== "N_abund"),:]

gmeasures_diff_abund = (gmeasures_models_abund .- gmeasures_N_abund) ./ gmeasures_N_abund

# add model names and reorder rows 
insertcols!(gmeasures_diff_abund, 1, :model => vcat("MaxEnt-co", "MaxEnt-jds", "null 1", "null 2", "neutral"))

gmeasures_diff_abund = gmeasures_diff_abund[vcat(5,3,1,4,2),:]

# format as Markdown table
table_abund = latexify(gmeasures_diff_abund, env=:mdtable, fmt="%.3f", latex=false, escape_underscores=true)

# export to file
table_abund_path = joinpath("tables", "measures_abund.md")
open(table_abund_path, "w") do io
    print(io, table_abund)
end

