# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)

## Read empirical data

# Mangal
networks_mangal = load(joinpath("data", "proc", "mangal", "networks_mangal.jld"))["data"]

# New Zealand
networks_NZ = load(joinpath("data", "proc", "new_zealand", "networks_NZ.jld"))["data"]
abund_data_NZ = load(joinpath("data", "proc", "new_zealand", "abund_data_NZ.jld"))["data"]

# Tuesday lake
networks_tuesday = load(joinpath("data", "proc", "tuesday_lake", "networks_tuesday.jld"))["data"]
abund_data_tuesday = load(joinpath("data", "proc", "tuesday_lake", "abund_data_tuesday.jld"))["data"]


## Read simulated data

# predicted numbers of links
predicted_links = load(joinpath("data", "sim", "predicted_links.jld"))["data"]

# degree distribution
degree_dist_mangal_sim = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_mangal.jld"))["data"]
degree_dist_NZ_sim = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_NZ.jld"))["data"]
degree_dist_tuesday_sim = load(joinpath("data", "sim", "degree_dist_maxent", "degree_dist_tuesday.jld"))["data"]

# joint degree distribution
joint_degree_dist_mangal_sim = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_mangal.jld"))["data"]
joint_degree_dist_NZ_sim = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_NZ.jld"))["data"]
joint_degree_dist_tuesday_sim = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_tuesday.jld"))["data"]

# network of maximum svd-entropy
network_maxent_mangal = load(joinpath("data", "sim", "network_maxent", "network_maxent_mangal.jld"))["data"]
network_maxent_NZ = load(joinpath("data", "sim", "network_maxent", "network_maxent_NZ.jld"))["data"]
network_maxent_tuesday = load(joinpath("data", "sim", "network_maxent", "network_maxent_tuesday.jld"))["data"]

# access unipartite networks of maximum entropy

N_maxent_mangal = [] 
for i in 1:length(network_maxent_mangal)
      push!(N_maxent_mangal, network_maxent_mangal[i].A)
end

N_maxent_NZ = [] 
for i in 1:length(network_maxent_NZ)
      push!(N_maxent_NZ, network_maxent_NZ[i].A)
end

N_maxent_tuesday = [] 
for i in 1:length(network_maxent_tuesday)
      push!(N_maxent_tuesday, network_maxent_tuesday[i].A)
end

# neutral models
neutral_networks_NZ = load(joinpath("data", "sim", "neutral_model", "neutral_networks_NZ.jld"))["data"]
neutral_networks_tuesday = load(joinpath("data", "sim", "neutral_model", "neutral_networks_tuesday.jld"))["data"]


##### Figures 

### Density of mean degree constraints for different richness ###

# Different quantiles of species richness will be plotted
S_mangal = richness.(networks_mangal)
S_NZ = richness.(networks_NZ)
S_tuesday = richness.(networks_tuesday)
S_all = vcat(S_mangal, S_NZ, S_tuesday)

S015 = Int64(round(quantile(S_all, 0.015))) # 1.5% lower quantile
S500 = Int64(round(quantile(S_all, 0.5))) # median
S985 = Int64(round(quantile(S_all, 0.985))) # 1.5% upper quantile

"""
kavg_dist(S::Int64)
    S: number of species
Returns the predicted distribution of mean degrees of the flexible links model
"""
function kavg_dist(S::Int64)
  L = predicted_links[:, S-4]
  kavg = 2 .* L ./ S 
  return kavg
end

# Get predicted distribution of mean degrees for each level of species richness considered
kavg015 = kavg_dist(S015)
kavg500 = kavg_dist(S500)
kavg985 = kavg_dist(S985)

# Plot distributions of mean degrees for the 3 levels of species richness
plotA = density(kavg015, 
                linewidth=2, 
                label="$(S015) species",
                framestyle=:box, 
                grid=false,
                dpi=1000, 
                size=(800,500), 
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white, 
                legendfont=fonts)
density!(kavg500, 
        linewidth=2, 
        label="$(S500) species")
density!(kavg985,  
        linewidth=2, 
        label="$(S985) species")
xaxis!(xlabel="Mean degree (links/species)", 
      xlims=(0, maximum(kavg985)))
yaxis!(ylabel="Density", 
      ylims=(0,0.6))


### Degree distribution of MaxEnt for fixed richness and different numbers of links ###

# quantiles of the predicted numbers of links that will be plotted
q = [0.015, 0.055, 0.165, 0.500, 0.835, 0.945, 0.985]
Lquant = Int64.(round.(quantile(predicted_links[:, S500-4], q))) # median number of species
L015 = Lquant[1]
L055 = Lquant[2]
L165 = Lquant[3]
L500 = Lquant[4]
L835 = Lquant[5]
L945 = Lquant[6]
L985 = Lquant[7]

# get the degree distribution of MaxEnt for a network with a median number of species and different numbers of links
dd_maxent_Lquant = degree_dist_maxent.(S500, Lquant)

# plot different degree distributions of MaxEnt for a network with a median number of species
plotB = plot(0:S500,
            dd_maxent_Lquant[1], 
            color=:black, 
            alpha=0.7, 
            linewidth=2, 
            linestyle=:dot, 
            label="$L015 links (97% PI)", # 97% percentile interval
            framestyle=:box, 
            grid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts)
plot!(0:S500,
      dd_maxent_Lquant[2], 
      color=:black, 
      alpha=0.7, 
      linestyle=:dash, 
      label="$L055 links (89% PI)") # 89% PI
plot!(0:S500,
      dd_maxent_Lquant[3], 
      color=:black, 
      alpha=0.7, 
      linestyle=:solid, 
      label="$L165 links (67% PI)") # 67% PI
plot!(0:S500,
      dd_maxent_Lquant[5], 
      color=:grey, 
      alpha=0.7,
      linestyle=:solid, 
      label="$L835 links (67% PI)") # 67% PI
plot!(0:S500,
      dd_maxent_Lquant[6], 
      color=:grey, 
      alpha=0.7, 
      linestyle=:dash, 
      label="$L945 links (89% PI)") # 89% PI
plot!(0:S500,
      dd_maxent_Lquant[7], 
      color=:grey, 
      alpha=0.7, 
      linestyle=:dot, 
      label="$L985 links (97% PI)") # 97% PI
plot!(0:S500,
      dd_maxent_Lquant[4], 
      color=:darkblue, 
      linewidth=4, 
      label="$L500 links (median)")
xaxis!(xlabel="Degree k", 
      xlims=(0,S500))
yaxis!(ylabel="p(k)")

plot(plotA, plotB,
    title = ["(a)" "(b)"],
    titleloc=:right, 
    titlefont=fonts)

savefig(joinpath("figures","maxent_degree_dist_fl"))


### Heatmap of disconnected species ###

# Get the MaxEnt degree distribution for a range of species richness and different quantiles of the number of links
Ss = [5:5:100;] # range of species richness
qs = [0:0.05:1;] # quantiles of the number of links
p_disconnected = zeros(length(qs), length(Ss)) # create vector for the probability a species is disconnected given the number of links (rows) and species richness (columns)

for (i, S) in enumerate(Ss)
  Ls = round.(quantile(0:S^2, qs)) # get the number of links L associated to every quantile (proportion of S^2 links below L)
  Ls = convert(Vector{Int64}, Ls) 
  
  degree_dists = degree_dist_maxent.(S, Ls) # get the degree distributions of maximum entropy for a given S and all L
  degree_dists = reduce(hcat, degree_dists)
  p_disconnected[:,i] = degree_dists[1,:]
end

q_minS = Ss ./ (Ss.^2 .+ 1) # get the quantile associated to the minimum number of links S-1 that would allow all species to be connected

# plot the proportion of disconnected species for all combinations of S and L (L in quantile)
plotA = heatmap(Ss, qs, log.(p_disconnected), 
        c=:viridis,
        colorbar_title="log probability of a species being isolated",
        framestyle=:box, 
        dpi=1000, 
        size=(800,500), 
        margin=5Plots.mm, 
        guidefont=fonts, 
        xtickfont=fonts, 
        ytickfont=fonts)
# add the minimum number of links S-1 (quantile) on the plot 
scatter!(Ss, q_minS, 
      linewidth=1, 
      linecolor=:black,
      markerwidth=3, 
      markercolor=:black,
      lab="")  
plot!(Ss, q_minS,
      linewidth=1, 
      linecolor=:black,
      lab="") 
xaxis!(xlim=(minimum(Ss)-2.5, maximum(Ss)+2.5), "species richness")
yaxis!(ylim=(0,1), "quantile of the number of links")


### Correlation between kin and kout ###

## get in and out degrees for all networks (empirical and maxent)

# empirical data

# mangal
kout_emp_mangal = reduce(vcat, collect.(values.(degree.(networks_mangal, dims=1))))
kin_emp_mangal = reduce(vcat, collect.(values.(degree.(networks_mangal, dims=2))))

# New Zealand
kout_emp_NZ = reduce(vcat, collect.(values.(degree.(networks_NZ, dims=1))))
kin_emp_NZ = reduce(vcat, collect.(values.(degree.(networks_NZ, dims=2))))

# Tuesday lake
kout_emp_tuesday = reduce(vcat, collect.(values.(degree.(networks_tuesday, dims=1))))
kin_emp_tuesday = reduce(vcat, collect.(values.(degree.(networks_tuesday, dims=2))))

# all empirical networks 
kout_emp = vcat(kout_emp_mangal, kout_emp_NZ, kout_emp_tuesday)
kin_emp = vcat(kin_emp_mangal, kin_emp_NZ, kin_emp_tuesday)


## simulated networks (from the joint degree distribution of maximum entropy)
"""
simulate_degrees(JDD::Matrix{Float64})
    JDD: joint degree distribution 
Returns a simulated vector of in and out degrees using the joint degree distribution as weight
"""
function simulate_degrees(JDD::Matrix{Float64})
      S = size(JDD, 1) - 1 # number of species 
      deg = findall(JDD .>= 0) # get cartesian indices 
      deg_samp = sample(deg, Weights(vec(JDD)), S, replace=true) # select species degrees randomly 
      # get in and out degrees
      kin = zeros(Int64, S)
      kout = zeros(Int64, S)
      for i in 1:S
            kout[i] = deg_samp[i][1]
            kin[i] = deg_samp[i][2]
      end
      return (kin = kin, kout = kout)
end

degrees_maxent_mangal = simulate_degrees.(joint_degree_dist_mangal_sim)
degrees_maxent_NZ = simulate_degrees.(joint_degree_dist_NZ_sim)
degrees_maxent_tuesday = simulate_degrees.(joint_degree_dist_tuesday_sim)

function get_kin(degrees::Vector)
      kin_maxent = []
      for i in 1:length(degrees)
            push!(kin_maxent, degrees[i].kin)
      end
      return reduce(vcat, kin_maxent)
end

function get_kout(degrees::Vector)
      kout_maxent = []
      for i in 1:length(degrees)
            push!(kout_maxent, degrees[i].kout)
      end
      return reduce(vcat, kout_maxent)
end

# Mangal
kin_maxent_mangal = get_kin(degrees_maxent_mangal)
kout_maxent_mangal = get_kout(degrees_maxent_mangal)

# New Zealand
kin_maxent_NZ = get_kin(degrees_maxent_NZ)
kout_maxent_NZ = get_kout(degrees_maxent_NZ)

# Tuesday lake
kin_maxent_tuesday = get_kin(degrees_maxent_tuesday)
kout_maxent_tuesday = get_kout(degrees_maxent_tuesday)

# all simulated networks
kin_maxent = vcat(kin_maxent_mangal, kin_maxent_NZ, kin_maxent_tuesday)
kout_maxent = vcat(kout_maxent_mangal, kout_maxent_NZ, kout_maxent_tuesday)

# plot the association between in and out degrees for empirical and simulated data
plotB = scatter(kout_emp, 
                kin_emp, 
                alpha=0.2, 
                markersize=3, 
                label="",
                framestyle=:box, 
                grid=false,
                dpi=1000, 
                size=(800,500), 
                aspect_ratio=:equal,
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white,
                legendfont=fonts)
xaxis!(xlim=(-1, maximum(vcat(kout_emp_mangal, kout_emp_NZ, kout_emp_tuesday))), "Kout (number of preys)")
yaxis!(ylim=(-1, maximum(vcat(kin_emp_mangal, kin_emp_NZ, kin_emp_tuesday))), "Kin (number of predators)")
                
plotC = scatter(kout_maxent, 
                kin_maxent, 
                alpha=0.2, 
                markersize=3, 
                label="",
                framestyle=:box, 
                grid=false,
                dpi=1000, 
                size=(800,500), 
                aspect_ratio=:equal,
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white,
                legendfont=fonts)
xaxis!(xlim=(-1, maximum(vcat(kout_emp_mangal, kout_emp_NZ, kout_emp_tuesday))), "Kout (number of preys)")
yaxis!(ylim=(-1, maximum(vcat(kin_emp_mangal, kin_emp_NZ, kin_emp_tuesday))), "Kin (number of predators)")


l = @layout [a [b ; c]]

plot(plotA, plotB, plotC,
    layout = l,
    title = ["(a)" "(b) empirical data" "(c) simulated data"],
    titleloc=:right, 
    titlefont=fonts)

savefig(joinpath("figures","joint_degree_dist.png"))


################# TK TO DO #########################
# Figure: Measures of empirical and maximum entropy food webs 

# Convert food webs archieved on Mangal to UnipartiteNetworks
Ns_mangal = network.(mangal_foodwebs.id)
Ns_mangal = convert.(UnipartiteNetwork, Ns_mangal)

# Measure species richness
S_mangal = richness.(Ns_mangal)
S_maxent_empL = richness.(Ns_maxent_empL)
S_maxent_fl = richness.(Ns_maxent_fl)

# Measure nestedness (spectral radius of the adjacency matrix)
nestedness_mangal = ρ.(Ns_mangal)
nestedness_maxent_empL = ρ.(Ns_maxent_empL)
nestedness_maxent_fl = ρ.(Ns_maxent_fl)

# Measure maximum trophic levels
tl_mangal = values.(trophic_level.(Ns_mangal))
tl_max_mangal = maximum.(tl_mangal)

tl_maxent_empL = values.(trophic_level.(Ns_maxent_empL))
tl_max_maxent_empL = maximum.(tl_maxent_empL)

tl_maxent_fl = values.(trophic_level.(Ns_maxent_fl))
tl_max_maxent_fl = maximum.(tl_maxent_fl)

# Measure network diameter (shortest distance between the two most distant nodes in the network)
shortest_paths_mangal = shortest_path.(Ns_mangal)
diameter_mangal = maximum.(shortest_paths_mangal)

shortest_paths_maxent_empL = shortest_path.(Ns_maxent_empL)
diameter_maxent_empL = maximum.(shortest_paths_maxent_empL)

shortest_paths_maxent_fl = shortest_path.(Ns_maxent_fl)
diameter_maxent_fl = maximum.(shortest_paths_maxent_fl)

entropy_mangal = svd_entropy.(Ns_mangal)
entropy_maxent_empL = svd_entropy.(Ns_maxent_empL)
entropy_maxent_fl = svd_entropy.(Ns_maxent_fl)

# Make vector of measures of MaxEnt food webs (nb of links given by the flexible links model) the same length as the other vectors
n1 = length(Ns_mangal) 
n2 = length(Ns_maxent_fl)
i = [findall(in(mangal_foodwebs.S[i]).(S_maxent_fl)) for i in 1:n1]
i[n1-1] = [n2] # The second to last entry has no value. We'll change it manually.
i = reduce(vcat, i)

nestedness_maxent_fl_long = [nestedness_maxent_fl[j] for j in i]
tl_max_maxent_fl_long = [tl_max_maxent_fl[j] for j in i]
diameter_maxent_fl_long = [diameter_maxent_fl[j] for j in i]
entropy_maxent_fl_long = [entropy_maxent_fl[j] for j in i]

# Plot predicted as a function of empirical measures
plotA = scatter(nestedness_mangal, nestedness_maxent_empL, alpha=0.3, markersize=3, label="Empirical L",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:bottomright, legendfont=fonts,
      xlabel="Nestedness (empirical webs)",
      ylabel="Nestedness (MaxEnt webs)")
scatter!(nestedness_mangal, nestedness_maxent_fl_long, alpha=0.3, markersize=3, label="Median L")
plot!(LinRange(0.4, 0.9, 100), LinRange(0.4, 0.9, 100), lab="", color="grey", linestyle=:dot)

plotB = scatter(tl_max_mangal, tl_max_maxent_empL, alpha=0.3, markersize=3, label="Empirical L",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:bottomright, legendfont=fonts,
      ticks=[2,4,6,8,10,12],
      xlabel="Maximum trophic level (empirical webs)",
      ylabel="Maximum trophic level (MaxEnt webs)")
scatter!(tl_max_mangal, tl_max_maxent_fl_long, alpha=0.3, markersize=3, label="Median L")
plot!(LinRange(2, 12, 100), LinRange(2, 12, 100), lab="", color="grey", linestyle=:dot)

plotC = scatter(diameter_mangal, diameter_maxent_empL, alpha=0.3, markersize=3, label="Empirical L",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:bottomright, legendfont=fonts,
      ticks=[1,3,5,7,9,11,13],
      xlabel="Diameter (empirical webs)",
      ylabel="Diameter (MaxEnt webs)")
scatter!(diameter_mangal .+ 0.1, diameter_maxent_fl_long, alpha=0.3, markersize=3, label="Median L")
plot!(LinRange(1, 13, 100), LinRange(1, 13, 100), lab="", color="grey", linestyle=:dot)

plotD = scatter(entropy_mangal, entropy_maxent_empL, alpha=0.3, markersize=3, label="Empirical L",
framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
foreground_color_legend=nothing, background_color_legend=:white,
legend=:bottomright, legendfont=fonts,
xlabel="SVD-entropy (empirical webs)",
ylabel="SVD-entropy (MaxEnt webs)")
scatter!(entropy_mangal, entropy_maxent_fl_long, alpha=0.3, markersize=3, label="Median L")
plot!(LinRange(0.77, 1, 100), LinRange(0.77, 1, 100), lab="", color="grey", linestyle=:dot)

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "measures_mangal_maxent.png"))
#########################################################

################# TK TO DO #########################
# Figure: Plot measures as a function of species richness and of other measures
a = [5,10,20,50,100,400] # specified x-ticks 

plotA = scatter(S_mangal, nestedness_maxent_empL, alpha=0.3, label="Empirical L", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts,
      xlabel="Species richness",
      ylabel="Nestedness")
scatter!(S_mangal, nestedness_maxent_fl_long, alpha=0.3, label="Median L", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
scatter!(S_mangal, nestedness_mangal, alpha=0.3, label="Mangal", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
xaxis!(:log, xticks=(a,a))

plotB = scatter(S_mangal, tl_max_maxent_empL, alpha=0.3, label="Empirical L", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts,
      xlabel="Species richness",
      ylabel="Maximum trophic level",
      yticks=[2,4,6,8,10,12])
scatter!(S_mangal, tl_max_maxent_fl_long, alpha=0.3, label="Median L", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
scatter!(S_mangal, tl_max_mangal, alpha=0.3, label="Mangal", smooth=true,markersize=3, linestyle=:dot, linealpha=1)
xaxis!(:log, xticks=(a,a))

plotC = scatter(S_mangal, diameter_maxent_empL, alpha=0.3, label="Empirical L", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts, 
      xlabel="Species richness",
      ylabel="Diameter",
      yticks=[1,3,5,7,9,11,13])
scatter!(S_mangal, diameter_maxent_fl_long, alpha=0.3, label="Median L", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
scatter!(S_mangal, diameter_mangal, alpha=0.3, label="Mangal", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
xaxis!(:log, xticks=(a,a))

plotD = scatter(S_mangal, entropy_maxent_empL, alpha=0.3, label="Empirical L", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:bottomright, legendfont=fonts, 
      xlabel="Species richness",
      ylabel="SVD-entropy")
scatter!(S_mangal, entropy_maxent_fl_long, alpha=0.3, label="Median L", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
scatter!(S_mangal, entropy_mangal, alpha=0.3, label="Mangal", smooth=true, markersize=3, linestyle=:dot, linealpha=1)
xaxis!(:log, xticks=(a,a))

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "measures_richness.png"))

# Plot nestedness as a function of maximum entropy
scatter(tl_max_maxent_empL, nestedness_maxent_empL, alpha=0.5, label="Empirical L", smooth=true, markersize=4, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts,
      xticks=[2,4,6,8,10,12],
      xlabel="Maximum trophic level",
      ylabel="Nestedness")
scatter!(tl_max_maxent_fl, nestedness_maxent_fl, alpha=0.5, label="Median L", smooth=true, markersize=4, linestyle=:dot, linealpha=1)
scatter!(tl_max_mangal, nestedness_mangal, alpha=0.5, label="Mangal", smooth=true, markersize=4, linestyle=:dot, linealpha=1)

savefig(joinpath("figures", "maxtrophiclevel_nestedness.png"))
#############################################################


################# TK TO DO #########################
## Figure: Entropy and divergence of degree distributions
"""
get_divergence_dd(N1::UnipartiteNetwork{Bool, MangalNode}, N2::UnipartiteNetwork{Bool, String})  
    N1: Empirical food web
    N2: MaxEnt food web 
Returns the divergence in probabilistic degree distributions between the two networks
"""
function get_divergence_dd(N1::UnipartiteNetwork{Bool, MangalNode}, N2::UnipartiteNetwork{Bool, String})  
  # Maximum in species richness between the two networks (since some predicted networks didn't have the exact same number of species as their empirical counterpart)
  S = maximum(vcat(richness(N1), richness(N2)))

  # Get the probabilistic degree distributions (proportion of species with degree k) of the two networks
  dd_N1 = get_ranked_dd(N1)
  dd_N2 = get_ranked_dd(N2)

  dd_prob_N1 = zeros(Float64, S)
  dd_prob_N2 = zeros(Float64, S)

  for k in 1:S
    dd_prob_N1[k] = sum(dd_N1 .== k) / S
    dd_prob_N2[k] = sum(dd_N2 .== k) / S
  end

  # Compute the divergence in degree distributions between the two networks (sum of absolute differences)
  divergence = sum(abs.(dd_prob_N1 .- dd_prob_N2))

  return divergence
end

# Compute the divergence in degree distributions between food webs archived on Mangal and MaxEnt food webs (empirical L)
divergence_dd = get_divergence_dd.(Ns_mangal, Ns_maxent_empL)

# Plot the divergence in degree distributions as a function of entropy
plotA = scatter(entropy_maxent_empL, divergence_dd, alpha=0.3, label="", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts, 
      xlabel="SVD-entropy (MaxEnt food web)",
      ylabel="Divergence in degree distributions")

# Plot the divergence in degree distributions as a function of species richness
plotB = scatter(S_mangal, divergence_dd, alpha=0.3, label="", smooth=true, markersize=3, linestyle=:dot, linealpha=1,
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts, 
      xlabel="Species richness",
      ylabel="Divergence in degree distributions")
xaxis!(:log, xticks=(a,a))

plot(plotA, plotB,  
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)
     
savefig(joinpath("figures", "divergence_degree_distributions.png"))
################################################################################

################# TK TO DO #########################
## Figure: Distribution of entropy and z-scores of empirical food webs

plotA = density(entropy_maxent_empL, label="Empirical L",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts,
      ylims=(0,30),
      xlabel="SVD-entropy",
      ylabel="Density")
density!(entropy_maxent_fl, label="Median L")
density!(entropy_mangal, label="Mangal")

entropy_maxent_fl_avg = mean(entropy_maxent_fl)
entropy_maxent_fl_std = std(entropy_maxent_fl)
entropy_mangal_zscores = (entropy_mangal .- entropy_maxent_fl_avg) ./ entropy_maxent_fl_std

entropy_mangal_zscores_500 = quantile(entropy_mangal_zscores, 0.500)
entropy_mangal_zscores_015 = quantile(entropy_mangal_zscores, 0.015)
entropy_mangal_zscores_985 = quantile(entropy_mangal_zscores, 0.985)

plotB = density(entropy_mangal_zscores, label="",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white,
      legend=:topright, legendfont=fonts,
      color=:grey,
      ylims=(0,0.16),
      xlabel="z-score of SVD-entropy",
      ylabel="Density")
plot!([entropy_mangal_zscores_500], seriestype=:vline, color=:grey, ls=:dash, lab="")

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "entropy_distribution.png"))
#####################################################################