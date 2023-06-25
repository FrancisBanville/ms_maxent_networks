# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)


## Read empirical networks

## Empirical data 
N_all = load(joinpath("data", "proc", "N_all.jld"))["data"]
N_abund = load(joinpath("data", "proc", "N_abund.jld"))["data"]

## Simulated data
N_ent_maxentjds_all = load(joinpath("data", "sim", "network_maxent", "N_maxentjds_all.jld"))["data"]
N_maxentjds_all = [N_ent_maxentjds_all[i].A for i in 1:length(N_ent_maxentjds_all)]

## Degree distributions of maximum entropy 
dd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "dd_maxent_all.jld"))["data"]
jdd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jdd_maxent_all.jld"))["data"]
jds_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"))["data"] 

## Table of network properties
measures = DataFrame(CSV.File(joinpath("results", "measures.csv")))

# subset networks
measures_all = measures[measures.network .== "N_all",:]
measures_maxentco_all = measures[measures.network .== "N_maxentco_all",:]
measures_maxentjds_all = measures[measures.network .== "N_maxentjds_all",:]
measures_nullco_all = measures[measures.network .== "N_nullco_all",:]
measures_nulljds_all = measures[measures.network .== "N_nulljds_all",:]

measures_abund = measures[measures.network .== "N_abund",:]
measures_maxentco_abund = measures[measures.network .== "N_maxentco_abund",:]
measures_maxentjds_abund = measures[measures.network .== "N_maxentjds_abund",:]
measures_nullco_abund = measures[measures.network .== "N_nullco_abund",:]
measures_nulljds_abund = measures[measures.network .== "N_nulljds_abund",:]
measures_neutral_abund = measures[measures.network .== "N_neutral_abund",:]


## Read counterfactuals of the flexible links model
predicted_links = load(joinpath("data", "sim", "predicted_links.jld"))["data"]


##### Figures 

### Density of mean degree constraints for different richness ###

# Different quantiles of species richness will be plotted
S_all = measures_all.S

S015 = Int64(round(quantile(S_all, 0.015))) # 1.5% lower quantile
S500 = Int64(round(quantile(S_all, 0.5))) # median
S985 = Int64(round(quantile(S_all, 0.985))) # 1.5% upper quantile

function kavg_dist(S::Int64)
  # Returns the predicted distribution of mean degrees for a given level of species richness
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
                label="$(S015) species",
                linewidth=2, 
                framestyle=:box, 
                grid=false,
                minorgrid=false,
                dpi=5000, 
                size=(800,500), 
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white, 
                legendfont=fonts,
                legendfontpointsize=7,
                legendfontfamily="Arial")
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
            label="$L015 links", # 97% percentile interval
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=5000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial")
plot!(0:S500,
      dd_maxent_Lquant[2], 
      color=:black, 
      alpha=0.7, 
      linestyle=:dash, 
      label="$L055 links") # 89% PI
plot!(0:S500,
      dd_maxent_Lquant[3], 
      color=:black, 
      alpha=0.7, 
      linestyle=:solid, 
      label="$L165 links") # 67% PI
plot!(0:S500,
      dd_maxent_Lquant[4], 
      color=:darkblue, 
      linewidth=4, 
      label="$L500 links")
plot!(0:S500,
      dd_maxent_Lquant[5], 
      color=:grey, 
      alpha=0.7,
      linestyle=:solid, 
      label="$L835 links") # 67% PI
plot!(0:S500,
      dd_maxent_Lquant[6], 
      color=:grey, 
      alpha=0.7, 
      linestyle=:dash, 
      label="$L945 links") # 89% PI
plot!(0:S500,
      dd_maxent_Lquant[7], 
      color=:grey, 
      alpha=0.7, 
      linestyle=:dot, 
      label="$L985 links") # 97% PI
xaxis!(xlabel="Degree k", 
      xlims=(0,S500))
yaxis!(ylabel="p(k)",
      ylims=(0,0.23))

plot(plotA, plotB,
    title = ["(a)" "(b)"],
    titleloc=:right, 
    titlefont=fonts)

savefig(joinpath("figures","maxent_degree_dist_fl.svg"))


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
heatmap(Ss, qs, log.(p_disconnected), 
        c=:viridis,
        colorbar_title="log probability of a species being isolated",
        framestyle=:box, 
        dpi=1000, 
        size=(800,500), 
        margin=5Plots.mm, 
        guidefont=fonts,
        colorbar_titlefontsize=7,
        colorbar_titlefontfamily="Arial",
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

savefig(joinpath("figures","heatmap_disconnected.png"))


### Correlation between kin and kout ###

## get in and out degrees for all networks (empirical and maxent)

# empirical data
kout_emp = reduce(vcat, collect.(values.(degree.(N_all, dims=1))))
kin_emp = reduce(vcat, collect.(values.(degree.(N_all, dims=2))))


## simulated networks (from the joint degree distribution of maximum entropy)

# all simulated networks

kout_maxent = reduce(vcat,[jds_maxent_all[i].kout for i in 1:length(jds_maxent_all)])
kin_maxent = reduce(vcat,[jds_maxent_all[i].kin for i in 1:length(jds_maxent_all)])

# get degree sequence of empirical and MaxEnt networks
# species' out and in degree sequence will be sorted using their total degree
k_emp = kout_emp .+ kin_emp
k_maxent = kout_maxent .+ kin_maxent

# get network id (to be used when getting species rank)
N_id = reduce(vcat, [fill(i, measures_all.S[i]) for i in 1:length(N_all)])
S_id = reduce(vcat, [fill(measures_all.S[i], measures_all.S[i]) for i in 1:length(N_all)])

# put in data frame
k_emp_df = DataFrame(k_tot = k_emp,
                        kout = kout_emp,
                        kin = kin_emp,
                        kout_rel = kout_emp ./ S_id,
                        kin_rel = kin_emp ./ S_id,
                        N_id = N_id)
k_maxent_df = DataFrame(k_tot = k_maxent,
                        kout = kout_maxent,
                        kin = kin_maxent,
                        kout_rel = kout_maxent ./ S_id,
                        kin_rel = kin_maxent ./ S_id,
                        N_id = N_id)

# plot the association between in and out degrees for empirical and simulated data
plotA = scatter(k_emp_df.kout_rel, 
                k_emp_df.kin_rel, 
                alpha=0.2, 
                markersize=5, 
                label="",
                framestyle=:box, 
                grid=false,
                minorgrid=false,
                dpi=1000, 
                size=(800,500), 
                aspect_ratio=:equal,
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white,
                legendfont=fonts,
                legendfontpointsize=7,
                legendfontfamily="Arial")
xaxis!(xlim=(0, 1), "relative Kout")
yaxis!(ylim=(0, 1), "relative Kin")
                
plotB = scatter(k_maxent_df.kout_rel, 
                k_maxent_df.kin_rel, 
                alpha=0.2, 
                markersize=5, 
                label="",
                framestyle=:box, 
                grid=false,
                minorgrid=false,
                dpi=1000, 
                size=(800,500), 
                aspect_ratio=:equal,
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white,
                legendfont=fonts,
                legendfontpointsize=7,
                legendfontfamily="Arial")
xaxis!(xlim=(0, 1), "relative Kout")
yaxis!(ylim=(0, 1), "relative Kin")


### Difference of in and out degrees ###

# get the difference in out and in degree of all species in all networks (except the largest) using different sortings
function kin_kout_diff(k_emp_df::DataFrame, k_maxent_df::DataFrame, column_sort::String, output::String) 
      k_diff = []
      for i in 1:length(N_all) 
            # get degree sequence of network i and sort it
            k_emp_sorted = k_emp_df[k_emp_df.N_id .== i,:]
            sort!(k_emp_sorted, column_sort, rev=true)

            k_maxent_sorted = k_maxent_df[k_maxent_df.N_id .== i,:]
            sort!(k_maxent_sorted, column_sort, rev=true)

            # compute difference in out and in degree for network i
            k_diff_sorted = k_maxent_sorted[!, output] .- k_emp_sorted[!, output]

            push!(k_diff, k_diff_sorted)
      end

      k_diff = reduce(vcat, k_diff)
      
      return k_diff
end

kout_diff_sortedby_ktot = kin_kout_diff(k_emp_df, k_maxent_df, "k_tot", "kout_rel")
kin_diff_sortedby_ktot = kin_kout_diff(k_emp_df, k_maxent_df, "k_tot", "kin_rel")

kout_diff_sortedby_kout = kin_kout_diff(k_emp_df, k_maxent_df, "kout", "kout_rel")
kin_diff_sortedby_kout = kin_kout_diff(k_emp_df, k_maxent_df, "kout", "kin_rel")

kout_diff_sortedby_kin = kin_kout_diff(k_emp_df, k_maxent_df, "kin", "kout_rel")
kin_diff_sortedby_kin = kin_kout_diff(k_emp_df, k_maxent_df, "kin", "kin_rel")

# plot differences of in and out degrees between empirical and MaxEnt networks (sorted by total degree)
plotC = scatter(kout_diff_sortedby_ktot,
            kin_diff_sortedby_ktot,
            markersize=5,
            alpha=0.2,
            aspect_ratio=:equal,
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial",
            legend=:topleft,
            label="")
xaxis!(xlim=(-1, 1), "\\Delta relative Kout")
yaxis!(ylim=(-1, 1), "\\Delta relative Kin")
            
plot!([0], 
      seriestype=:vline, 
      linewidth=0.3,
      color=:black, 
      lab="")
plot!([0], 
      seriestype=:hline, 
      linewidth=0.3,
      color=:black, 
      lab="")

l = @layout [[a ; b] c]

plot(plotA, plotB, plotC,
      layout = l,
      title = ["(a) empirical" "(b) MaxEnt" "(c) difference"],
      titleloc=:right, 
      titlefont=fonts)
      
savefig(joinpath("figures","joint_degree_dist.png"))


# plot differences of in and out degrees between empirical and MaxEnt networks (sorted by out degree)
plotA = scatter(kout_diff_sortedby_kout,
            kin_diff_sortedby_kout,
            markersize=5,
            alpha=0.2,
            aspect_ratio=:equal,
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial",
            legend=:topleft,
            label="")
xaxis!(xlim=(-1, 1), "\\Delta relative Kout")
yaxis!(ylim=(-1, 1), "\\Delta relative Kin")
            
plot!([0], 
      seriestype=:vline, 
      linewidth=0.3,
      color=:black, 
      lab="")
plot!([0], 
      seriestype=:hline, 
      linewidth=0.3,
      color=:black, 
      lab="")

# plot differences of in and out degrees between empirical and MaxEnt networks (sorted by in degree)
plotB = scatter(kout_diff_sortedby_kin,
            kin_diff_sortedby_kin,
            markersize=5,
            alpha=0.2,
            aspect_ratio=:equal,
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial",
            legend=:topleft,
            label="")
xaxis!(xlim=(-1, 1), "\\Delta relative Kout")
yaxis!(ylim=(-1, 1), "\\Delta relative Kin")

plot!([0], 
      seriestype=:vline, 
      linewidth=0.3,
      color=:black, 
      lab="")
plot!([0], 
      seriestype=:hline, 
      linewidth=0.3,
      color=:black, 
      lab="")


plot(plotA, plotB,
      title = ["(a) sorted by Kout" "(b) sorted by Kin"],
      titleloc=:right, titlefont=fonts)
 
savefig(joinpath("figures", "kin_kout_difference.png"))
 
 
### KL divergence between in and out degree distributions ###

## get in and out degree sequences
# maximum entropy joint degree sequence
kin_maxent_all = [jds_maxent_all[i].kin for i in 1:length(jds_maxent_all)]
kout_maxent_all = [jds_maxent_all[i].kout for i in 1:length(jds_maxent_all)]

# empirical joint degree sequence
jds_all = joint_degree_seq.(N_all)
kin_all = [jds_all[i].kin .+ 1 for i in 1:length(jds_all)]
kout_all = [jds_all[i].kout .+ 1 for i in 1:length(jds_all)]

## get in and out degree distributions from in and out degree sequences
function prop_k(ks::Vector{Int64}) # ks: in or out degree sequence
      # number of species
      S = length(ks)
      # count the number of species having a in or out degree of k
      ks_count = zeros(S + 1)
      
      for k in 0:S
            ks_count[k+1] = sum(ks .== k) .+ 1 # add one for calculation of K-L divergence
      end

      return ks_count ./ S # in or out degree distribution
end

# maximum entropy in and out degree distributions
kin_dist_maxent_all = prop_k.(kin_maxent_all)
kout_dist_maxent_all = prop_k.(kout_maxent_all)

# empirical in and out degree distributions
kin_dist_all = prop_k.(kin_all)
kout_dist_all = prop_k.(kout_all)

## compute KL divergence between in and out degree distributions
# maximum entropy in and out degree distributions
kl_diverg_maxent_all = kl_divergence.(kin_dist_maxent_all, kout_dist_maxent_all)
# empirical in and out degree distributions
kl_diverg_all = kl_divergence.(kin_dist_all, kout_dist_all)

## plot density of both KL divergence distributions
plotA = density(kl_diverg_all,
                  linesize=3,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legend=:topright,
                  label="Empirical",
                  xlabel="KL divergence",
                  ylabel="Density")
density!(kl_diverg_maxent_all,
            linesize=3,
            label="MaxEnt",
            ylim=(0,3.6))

## plot difference in KL divergences as a function of connectance
plotB = scatter(measures_all.C,
                  kl_diverg_maxent_all.- kl_diverg_all,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legend=:topright,
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Connectance",
                  ylabel="\\Delta KL divergence")

plot(plotA, plotB,  
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "kl_divergence.png"))


### Measures of empirical and maximum entropy food webs ###

# Nestedness
plotA = scatter(measures_all.rho, 
                  measures_maxentjds_all.rho,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Nestedness (empirical)",
                  ylabel="Nestedness (MaxEnt)")
plot!(LinRange(0.4, 0.9, 100), LinRange(0.4, 0.9, 100), lab="", color="grey", linestyle=:dot)

# Maximum trophic level
plotB = scatter(measures_all.maxtl, 
                  measures_maxentjds_all.maxtl,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Maximum trophic level (empirical)",
                  ylabel="Maximum trophic level (MaxEnt)")
plot!(LinRange(2, 10, 100), LinRange(2, 10, 100), lab="", color="grey", linestyle=:dot)

plotC = scatter(measures_all.diam,
                  measures_maxentjds_all.diam,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Diameter (empirical)",
                  ylabel="Diameter (MaxEnt)")
plot!(LinRange(1, 10, 100), LinRange(1, 10, 100), lab="", color="grey", linestyle=:dot)

plotD = scatter(measures_all.entropy,
                  measures_maxentjds_all.entropy,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="SVD-entropy (empirical)",
                  ylabel="SVD-entropy (MaxEnt)")
plot!(LinRange(0.75, 1, 100), LinRange(0.75, 1, 100), lab="", color="grey", linestyle=:dot)

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "measures_emp_maxent.png"))



### Measures as a function of species richness ###

a = [5,10,20,50,100,400] # specified x-ticks 

# Nestedness and species richness
plotA = scatter(measures_all.S,
                  measures_all.rho,
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legend=:topright,
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Nestedness")
scatter!(measures_maxentjds_all.S, 
            measures_maxentjds_all.rho,
            alpha=0.3,
            label="MaxEnt", 
            smooth=true, 
            linealpha=0.9,
            linewidth=2,
            markersize=5)
xaxis!(:log, xticks=(a,a))

# Maximum trophic level and species richness

missing_tl = findall(ismissing, measures_maxentjds_all.maxtl) # remove missing values
measures_maxentjds_maxtl = measures_maxentjds_all.maxtl[Not(missing_tl),:] 
measures_maxentjds_S = measures_maxentjds_all.S[Not(missing_tl),:]

plotB = scatter(measures_all.S,
                  measures_all.maxtl,
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legend=:topleft,
                  legendfont=fonts,
                  legendfontpointsize=7,
                   legendfontfamily="Arial",
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Maximum trophic level")
scatter!(measures_maxentjds_S, 
            measures_maxentjds_maxtl,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            linealpha=0.9,
            linewidth=2,
            markersize=5)
xaxis!(:log, xticks=(a,a))

# Network diameter and species richness
plotC = scatter(measures_all.S,
                  measures_all.diam,
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legend=:topleft,
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Diameter")
scatter!(measures_maxentjds_all.S,
            measures_maxentjds_all.diam,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            linealpha=0.9,
            linewidth=2,
            markersize=5)
xaxis!(:log, xticks=(a,a))

# Entropy and species richness
plotD = scatter(measures_all.S,      
                  measures_all.entropy,
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  alpha=0.3,
                  markersize=5,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legend=:bottomleft,
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="SVD-entropy")
scatter!(measures_maxentjds_all.S,
            measures_maxentjds_all.entropy,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            linealpha=0.9,
            linewidth=2,
            markersize=5)
xaxis!(:log, xticks=(a,a))

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "measures_richness.png"))


### Plot nestedness as a function of maximum trophic level ###

missing_tl = findall(ismissing, measures_maxentjds_all.maxtl) # remove missing values
measures_maxentjds_maxtl = measures_maxentjds_all.maxtl[Not(missing_tl),:] 
measures_maxentjds_rho = measures_maxentjds_all.rho[Not(missing_tl),:]

scatter(measures_all.maxtl,
            measures_all.rho,
            smooth=true,
            linealpha=0.9,
            linewidth=2,
            alpha=0.3,
            markersize=6,
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legend=:topright,
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial",
            label="Empirical",
            xlabel="Maximum trophic level",
            ylabel="Nestedness")
scatter!(measures_maxentjds_maxtl,
            measures_maxentjds_rho,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            linealpha=0.9,
            linewidth=2,
            markersize=6)
      
savefig(joinpath("figures", "maxtrophiclevel_nestedness.png"))


### Divergence between degree sequences (MaxEnt vs empirical networks) ###

## Divergence between the degree sequence of MaxEnt and empirical networks (MSD_ds)

# sorted degree sequence of empirical networks 
k_emp_all = joint_degree_seq.(N_all)

k_emp_all_sorted = [sort(k_emp_all[i].kin .+ k_emp_all[i].kout, rev=true) for i in 1:length(k_emp_all)]

# sorted degree sequence of MaxEnt networks
jds_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"))["data"]

k_maxent_all_sorted = [sort(jds_maxent_all[i].kin .+ jds_maxent_all[i].kout, rev=true) for i in 1:length(jds_maxent_all)]

# mean squared deviation between empirical and MaxEnt degree sequence 
MSD_ds_maxent = [mean(k_emp_all_sorted[i] .- k_maxent_all_sorted[i]).^2 for i in 1:length(k_maxent_all_sorted)]

plotA = scatter(measures_all.S,
                  MSD_ds_maxent,
                  alpha=0.3,
                  markersize=6,
                  framestyle=:box, 
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Species richness",
                  ylabel="MSD of degree sequence",
                  ylims=(0,20.5))
xaxis!(:log, xticks=(a,a))


# Divergence in degree sequence and SVD-entropy
plotB = scatter(measures_all.entropy,
                  MSD_ds_maxent,
                  alpha=0.3,
                  markersize=6,
                  framestyle=:box, 
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="SVD-entropy (empirical)",
                  ylabel="MSD of degree sequence",
                  ylims=(0,20.5))

plot(plotA, plotB, 
      title = ["(a)" "(b)"],
      titleloc=:right, titlefont=fonts)
             
savefig(joinpath("figures", "divergence_degree_sequence.png"))
             



### Difference in SVD-entropy  ###

# Difference in SVD-entropy and species richness

entropy_all = svd_entropy.(N_all)
entropy_maxentjds_all = svd_entropy.(N_maxentjds_all)

entropy_diffjds = entropy_maxentjds_all .- entropy_all

plotA = scatter(measures_all.S,
                  entropy_diffjds,
                  alpha=0.3,
                  markersize=6,
                  framestyle=:box, 
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Species richness",
                  ylabel="\\Delta SVD-entropy")
xaxis!(:log, xticks=(a,a))

b = [10,100,1000,10000] # specified x-ticks 
# Difference in SVD-entropy and number of links
plotB = scatter(measures_all.L,
                  entropy_diffjds,
                  alpha=0.3,
                  markersize=6,
                  framestyle=:box, 
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Number of links",
                  ylabel="\\Delta SVD-entropy")
xaxis!(:log, xticks=(b,b))

# Difference in SVD-entropy and connectance
plotC = scatter(measures_all.C,
                  entropy_diffjds,
                  alpha=0.3,
                  markersize=6,
                  framestyle=:box, 
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="Connectance",
                  ylabel="\\Delta SVD-entropy")

plot(plotA, plotB, plotC,
      layout = grid(1,3),
      title = ["(a)" "(b)" "(c)"],
      titleloc=:right, titlefont=fonts)
             
savefig(joinpath("figures", "difference_entropy.png"))



### Distribution of entropy and z-scores of empirical food webs ###

# Distribution of entropies
plotA = density(measures_all.entropy,
                  linesize=3,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legend=:topleft,
                  label="Empirical",
                  xlabel="SVD-entropy",
                  ylabel="Density",
                  ylim=(0,21))
density!(measures_maxentjds_all.entropy,
            linesize=3,
            label="MaxEnt")

# compute z-scores
entropy_maxentjds_avg = mean(measures_maxentjds_all.entropy)
entropy_maxentjds_std = std(measures_maxentjds_all.entropy)

entropy_zscores = (measures_all.entropy .- entropy_maxentjds_avg) ./ entropy_maxentjds_std

entropy_zscores_500 = quantile(entropy_zscores, 0.500)
entropy_zscores_015 = quantile(entropy_zscores, 0.015)
entropy_zscores_985 = quantile(entropy_zscores, 0.985)

plotB = density(entropy_zscores, 
                  color=:grey,
                  linesize=3,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  label="",
                  xlabel="z-score of SVD-entropy",
                  ylabel="Density",
                  ylim=(0,0.27))
plot!([entropy_zscores_500], 
            seriestype=:vline, 
            color=:grey, 
            ls=:dash, 
            lab="")

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "entropy_distribution.png"))


### Jaccard dissimilarity and difference in nestedness and SVD-entropy ###

## Difference in nestedness (rho_diff)
rho_diffjds = measures_maxentjds_all.rho .- measures_all.rho


## Jaccard distance (jaccard)

A_all = [vec(Matrix(N_all[i].edges)) for i in 1:length(N_all)]

A_maxentjds_all = [vec(Matrix(N_maxentjds_all[i].edges)) for i in 1:length(N_maxentjds_all)]

jaccard_maxentjds = jaccard.(A_all, A_maxentjds_all)


plotA = scatter(entropy_diffjds,
            rho_diffjds,
            alpha=0.3,
            markersize=6,
            smooth=true,
            linealpha=0.9,
            linewidth=2,
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=7,
            legendfontfamily="Arial",
            legend=:topleft,
            label="",
            xlabel="\\Delta SVD-entropy",
            ylabel="\\Delta nestedness")

plotB = scatter(entropy_diffjds,
                  jaccard_maxentjds,
                  alpha=0.3,
                  markersize=6,
                  smooth=true,
                  linealpha=0.9,
                  linewidth=2,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legend=:topleft,
                  label="",
                  xlabel="\\Delta SVD-entropy",
                  ylabel="Jaccard distance")

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "difference_entropy_jaccard.png"))


### Motif distribution ###

# get motifs names
motifs = keys(unipartitemotifs())

# select proportion of each motifs in empirical networks and tidy data frame
motifs_all = select(measures_all, vcat([motifs[i] for i in 1:5]))
motifs_all = DataFrames.stack(motifs_all)
motifs_all = dropmissing(motifs_all)

# select proportion of each motifs in MaxEnt networks and tidy data frame
motifs_maxentjds_all = select(measures_maxentjds_all, vcat([motifs[i] for i in 1:5]))
motifs_maxentjds_all = DataFrames.stack(motifs_maxentjds_all)
motifs_maxentjds_all = dropmissing(motifs_maxentjds_all)

motifs_maxentco_all = select(measures_maxentco_all, vcat([motifs[i] for i in 1:5]))
motifs_maxentco_all = DataFrames.stack(motifs_maxentco_all)
motifs_maxentco_all = dropmissing(motifs_maxentco_all)

motifs_nulljds_all = select(measures_nulljds_all, vcat([motifs[i] for i in 1:5]))
motifs_nulljds_all = DataFrames.stack(motifs_nulljds_all)
motifs_nulljds_all = dropmissing(motifs_nulljds_all)

motifs_nullco_all = select(measures_nullco_all, vcat([motifs[i] for i in 1:5]))
motifs_nullco_all = DataFrames.stack(motifs_nullco_all)
motifs_nullco_all = dropmissing(motifs_nullco_all)

# join both datasets and label them
groups = vcat(fill("1-Empirical", nrow(motifs_all)),
              fill("2-Null 1", nrow(motifs_nullco_all)),
              fill("3-MaxEnt-co", nrow(motifs_maxentco_all)),
              fill("4-Null 2", nrow(motifs_nulljds_all)),
              fill("5-MaxEnt-jds", nrow(motifs_maxentjds_all)))

motifs_emp_maxent_neutral = [motifs_all; 
                              motifs_nullco_all;
                              motifs_maxentco_all; 
                              motifs_nulljds_all;
                              motifs_maxentjds_all]
insertcols!(motifs_emp_maxent_neutral, :group => groups)

# motif distribution
groupedboxplot(motifs_emp_maxent_neutral.variable, 
                  motifs_emp_maxent_neutral.value,
                  group=motifs_emp_maxent_neutral.group, 
                  alpha=0.7,
                  linewidth=1,
                  markersize=3,
                  framestyle=:box, 
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  grid=false,
                  minorgrid=false,
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legendposition=:outertopright,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  xminorgrid=false,
                  label=hcat("Empirical", 
                        "Null 1",
                        "MaxEnt 1",
                        "Null 2",
                        "MaxEnt 2"),
                  ylims=(0,1),
                  xaxis="Motifs", 
                  yaxis="Proportion")

savefig(joinpath("figures", "motifs_distribution.png"))



#### Motifs S1 vs S2 ####

plotA = scatter(measures_all.S1,
                  measures_all.S2,
                  alpha=0.2,
                  markersize=4,
                  linealpha=0.9,
                  linewidth=2,
                  smooth=true,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legend=:topleft,
                  label="Empirical",
                  xlabel="Proportion S1",
                  ylabel="Proportion S2",
                  xlims=(0,1),
                  ylims=(0,1))

scatter!(measures_nullco_all.S1,
                  measures_nullco_all.S2,
                  alpha=0.2, 
                  label="Null 1", 
                  smooth=true, 
                  markersize=4, 
                  linealpha=0.9,
                  linewidth=2)

scatter!(measures_maxentco_all.S1,
                  measures_maxentco_all.S2,
                  alpha=0.2, 
                  label="MaxEnt 1", 
                  smooth=true, 
                  markersize=4, 
                  linealpha=0.9,
                  linewidth=2)
                  
scatter!(measures_nulljds_all.S1,
                  measures_nulljds_all.S2,
                  alpha=0.2, 
                  label="Null 2", 
                  smooth=true, 
                  markersize=4, 
                  linealpha=0.9,
                  linewidth=2)

scatter!(measures_maxentjds_all.S1,
                  measures_maxentjds_all.S2,
                  alpha=0.2, 
                  label="MaxEnt 2", 
                  smooth=true, 
                  markersize=4, 
                  linealpha=0.9,
                  linewidth=2)
                  

              
plotB = scatter(measures_all.S4,
                  measures_all.S5,
                  alpha=0.2,
                  markersize=4,
                  linealpha=0.9,
                  linewidth=2,
                  smooth=true,
                  framestyle=:box, 
                  grid=false,
                  minorgrid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  legendfontpointsize=7,
                  legendfontfamily="Arial",
                  legend=:topright,
                  label="Empirical",
                  xlabel="Proportion S4",
                  ylabel="Proportion S5",
                  xlims=(0,1),
                  ylims=(0,1))

scatter!(measures_nullco_all.S4,
            measures_nullco_all.S5,
            alpha=0.2, 
            label="Null 1", 
            smooth=true, 
            markersize=4, 
            linealpha=0.9,
            linewidth=2)

scatter!(measures_maxentco_all.S4, 
            measures_maxentco_all.S5,
            alpha=0.2, 
            label="MaxEnt 1", 
            smooth=true, 
            markersize=4, 
            linealpha=0.9,
            linewidth=2)

scatter!(measures_nulljds_all.S4,
            measures_nulljds_all.S5,
            alpha=0.2, 
            label="Null 2", 
            smooth=true, 
            markersize=4, 
            linealpha=0.9,
            linewidth=2)

scatter!(measures_maxentjds_all.S4,
            measures_maxentjds_all.S5,
            alpha=0.2, 
            label="MaxEnt 2", 
            smooth=true, 
            markersize=4, 
            linealpha=0.9,
            linewidth=2)


plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "motifs_relations.png"))

