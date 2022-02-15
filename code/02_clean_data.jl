# Here we tidy food web data

## Read all food webs and convert them to UnipartiteNetworks

# food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl) 
mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "raw", "mangal", "mangal_foodwebs.csv")))
N_mangal = network.(mangal_foodwebs.id)
N_mangal = convert.(UnipartiteNetwork, N_mangal)
N_mangal = UnipartiteNetwork.(N_mangal[i].edges for i in 1:length(N_mangal))

# New-Zealand food webs
NZ_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "new_zealand", "adjacency_matrices")),
                                        drop=[1])) # the first column is row names
N_NZ = convert.(Matrix{Bool}, NZ_foodwebs)
N_NZ = UnipartiteNetwork.(N_NZ, names.(NZ_foodwebs))

# Tuesday lake food webs
tuesday_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "tuesday_lake", "adjacency_matrices")),
                                        header=false)) # no column names here
N_tuesday = convert.(Matrix{Bool}, tuesday_foodwebs)
N_tuesday = UnipartiteNetwork.(N_tuesday)

# simplify all food webs
N_all = vcat(N_mangal, N_NZ, N_tuesday)
N_all = simplify.(N_all)

# remove biggest food web
N_all = N_all[richness.(N_all) .!= maximum(richness.(N_all))]

save(joinpath("data", "proc", "N_all.jld"), "data", N_all)


## Read all food webs with abundance data and convert them to UnipartiteNetworks

## New Zealand food webs

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


## Tuesday lake food webs

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

# all food webs with abundance data (New Zealand and Tuesday lake)
N_abund = vcat([abund_NZ[i].network for i in 1:length(abund_NZ)], 
               [abund_tuesday[i].network for i in 1:length(abund_tuesday)])

N_abund = simplify.(N_abund)

save(joinpath("data", "proc", "N_abund.jld"), "data", N_abund)

# all abundance distributions 
A_abund = vcat([abund_NZ[i].abundance for i in 1:length(abund_NZ)], 
               [abund_tuesday[i].abundance for i in 1:length(abund_tuesday)])

save(joinpath("data", "proc", "A_abund.jld"), "data", A_abund)
