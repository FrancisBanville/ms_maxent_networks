# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)


## Read empirical networks

## Empirical data 
N_all = load(joinpath("data", "proc", "N_all.jld"))["data"]
N_abund = load(joinpath("data", "proc", "N_abund.jld"))["data"]
A_abund = load(joinpath("data", "proc", "A_abund.jld"))["data"]

## Degree distributions of maximum entropy 
dd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "dd_maxent_all.jld"))["data"]
jdd_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jdd_maxent_all.jld"))["data"]
jds_maxent_all = load(joinpath("data", "sim", "degree_dist_maxent", "jds_maxent_all.jld"))["data"] 

## Tables of network properties
metrics = DataFrame(CSV.File(joinpath("results", "metrics.csv")))
gmetrics_all = DataFrame(CSV.File(joinpath("results", "gmetrics_all.csv")))
gmetrics_abund = DataFrame(CSV.File(joinpath("results", "gmetrics_abund.csv")))

## Tables of differences in network properties
metrics_diffco = DataFrame(CSV.File(joinpath("results", "metrics_diffco.csv")))
metrics_diffjds = DataFrame(CSV.File(joinpath("results", "metrics_diffjds.csv")))

# subset networks
metrics_emp = metrics[metrics.network .== "N_all",:]
metrics_maxentco_all = metrics[metrics.network .== "N_maxentco_all",:]
metrics_maxentjds_all = metrics[metrics.network .== "N_maxentjds_all",:]
metrics_nullco_all = metrics[metrics.network .== "N_nullco_all",:]
metrics_nulljds_all = metrics[metrics.network .== "N_nulljds_all",:]


metrics_abund = metrics[metrics.network .== "N_abund",:]
metrics_maxentco_abund = metrics[metrics.network .== "N_maxentco_abund",:]
metrics_maxentjds_abund = metrics[metrics.network .== "N_maxentjds_abund",:]
metrics_nullco_abund = metrics[metrics.network .== "N_nullco_abund",:]
metrics_nulljds_abund = metrics[metrics.network .== "N_nulljds_abund",:]
metrics_neutral_abund = metrics[metrics.network .== "N_neutral_abund",:]


## Read counterfactuals of the flexible links model
predicted_links = load(joinpath("data", "sim", "predicted_links.jld"))["data"]


##### Figures 

### Density of mean degree constraints for different richness ###

# Different quantiles of species richness will be plotted
S_emp = metrics_emp.S

S015 = Int64(round(quantile(S_emp, 0.015))) # 1.5% lower quantile
S500 = Int64(round(quantile(S_emp, 0.5))) # median
S985 = Int64(round(quantile(S_emp, 0.985))) # 1.5% upper quantile

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
heatmap(Ss, qs, log.(p_disconnected), 
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


# plot the association between in and out degrees for empirical and simulated data
plotA = scatter(kout_emp, 
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
xaxis!(xlim=(-1, maximum(vcat(kout_emp, kout_maxent))), "Kout")
yaxis!(ylim=(-1, maximum(vcat(kin_emp, kin_maxent))), "Kin")
                
plotB = scatter(kout_maxent, 
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
xaxis!(xlim=(-1, maximum(vcat(kout_emp, kout_maxent))), "Kout")
yaxis!(ylim=(-1, maximum(vcat(kin_emp, kin_maxent))), "Kin")


### Difference of in and out degrees ###

# get degree sequence of empirical and MaxEnt networks
# species' out and in degree sequence will be sorted using their total degree
k_emp = kout_emp .+ kin_emp
k_maxent = kout_maxent .+ kin_maxent

# get network id (to be used when getting species rank)
N_id = reduce(vcat, [fill(i, metrics_emp.S[i]) for i in 1:length(N_all)])
S_id = reduce(vcat, [fill(metrics_emp.S[i], metrics_emp.S[i]) for i in 1:length(N_all)])

# put in data frame
k_emp_df = DataFrame(k_tot = k_emp,
                        kout = kout_emp,
                        kin = kin_emp,
                        N_id = N_id)
k_maxent_df = DataFrame(k_tot = k_maxent,
                        kout = kout_maxent,
                        kin = kin_maxent,
                        N_id = N_id)

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

      # remove the largest network (for graphical reasons)
      k_diff_nolarge = k_diff[S_id .!= maximum(S_id)]
      
      return k_diff_nolarge
end

kout_diff_sortedby_ktot = kin_kout_diff(k_emp_df, k_maxent_df, "k_tot", "kout")
kin_diff_sortedby_ktot = kin_kout_diff(k_emp_df, k_maxent_df, "k_tot", "kin")

kout_diff_sortedby_kout = kin_kout_diff(k_emp_df, k_maxent_df, "kout", "kout")
kin_diff_sortedby_kout = kin_kout_diff(k_emp_df, k_maxent_df, "kout", "kin")

kout_diff_sortedby_kin = kin_kout_diff(k_emp_df, k_maxent_df, "kin", "kout")
kin_diff_sortedby_kin = kin_kout_diff(k_emp_df, k_maxent_df, "kin", "kin")

S_id_nolarge = S_id[S_id .!= maximum(S_id)]

# plot differences of in and out degrees between empirical and MaxEnt networks (sorted by total degree)
plotC = scatter(kout_diff_sortedby_ktot,
            kin_diff_sortedby_ktot,
            markersize=S_id_nolarge ./ 10,
            alpha=0.2,
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
            legendfont=fonts,
            legend=:topleft,
            label="",
            xlabel="\\Delta Kout",
            ylabel="\\Delta Kin")
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
            markersize=S_id_nolarge ./ 10,
            alpha=0.2,
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
            legendfont=fonts,
            legend=:topleft,
            label="",
            xlabel="\\Delta Kout",
            ylabel="\\Delta Kin")
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
            markersize=S_id_nolarge ./ 10,
            alpha=0.2,
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
            legendfont=fonts,
            legend=:topleft,
            label="",
            xlabel="\\Delta Kout",
            ylabel="\\Delta Kin")
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
 
 

### Measures of empirical and maximum entropy food webs ###

# Nestedness
plotA = scatter(metrics_emp.rho, 
                  metrics_maxentjds.rho,
                  alpha=0.3,
                  markersize=3,
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
                  legendfont=fonts,
                  label="",
                  xlabel="Nestedness (empirical)",
                  ylabel="Nestedness (MaxEnt)")
plot!(LinRange(0.4, 0.9, 100), LinRange(0.4, 0.9, 100), lab="", color="grey", linestyle=:dot)

# Maximum trophic level
plotB = scatter(metrics_emp.maxtl, 
                  metrics_maxentjds.maxtl,
                  alpha=0.3,
                  markersize=3,
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
                  legendfont=fonts,
                  label="",
                  xlabel="Maximum trophic level (empirical)",
                  ylabel="Maximum trophic level (MaxEnt)")
plot!(LinRange(2, 10, 100), LinRange(2, 10, 100), lab="", color="grey", linestyle=:dot)

plotC = scatter(metrics_emp.diam,
                  metrics_maxentjds.diam,
                  alpha=0.3,
                  markersize=3,
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
                  legendfont=fonts,
                  label="",
                  xlabel="Diameter (empirical)",
                  ylabel="Diameter (MaxEnt)")
plot!(LinRange(1, 10, 100), LinRange(1, 10, 100), lab="", color="grey", linestyle=:dot)

plotD = scatter(metrics_emp.entropy,
                  metrics_maxentjds.entropy,
                  alpha=0.3,
                  markersize=3,
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
                  legendfont=fonts,
                  label="",
                  xlabel="SVD-entropy (empirical)",
                  ylabel="SVD-entropy (MaxEnt)")
plot!(LinRange(0.75, 1, 100), LinRange(0.75, 1, 100), lab="", color="grey", linestyle=:dot)

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "metrics_emp_maxent.png"))



### Metrics as a function of species richness ###

a = [5,10,20,50,100,400] # specified x-ticks 

# Nestedness and species richness
plotA = scatter(metrics_emp.S,
                  metrics_emp.rho,
                  smooth=true,
                  linestyle=:dot, 
                  linealpha=1,
                  alpha=0.3,
                  markersize=3,
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
                  legend=:topright,
                  legendfont=fonts,
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Nestedness")
scatter!(metrics_maxentjds.S, 
            metrics_maxentjds.rho,
            alpha=0.3,
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
xaxis!(:log, xticks=(a,a))

# Maximum trophic level and species richness

missing_tl = findall(ismissing, metrics_maxentjds.maxtl) # remove missing values
metrics_maxent_maxtl = metrics_maxentjds.maxtl[Not(missing_tl),:] 
metrics_maxent_S = metrics_maxentjds.S[Not(missing_tl),:]

plotB = scatter(metrics_emp.S,
                  metrics_emp.maxtl,
                  smooth=true,
                  linestyle=:dot, 
                  linealpha=1,
                  alpha=0.3,
                  markersize=3,
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
                  legend=:topright,
                  legendfont=fonts,
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Maximum trophic level")
scatter!(metrics_maxent_S, 
            metrics_maxent_maxtl,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
xaxis!(:log, xticks=(a,a))

# Network diameter and species richness
plotC = scatter(metrics_emp.S,
                  metrics_emp.diam,
                  smooth=true,
                  linestyle=:dot, 
                  linealpha=1,
                  alpha=0.3,
                  markersize=3,
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
                  legend=:topright,
                  legendfont=fonts,
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="Diameter")
scatter!(metrics_maxentjds.S,
            metrics_maxentjds.diam,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
xaxis!(:log, xticks=(a,a))

# Entropy and species richness
plotD = scatter(metrics_emp.S,      
                  metrics_emp.entropy,
                  smooth=true,
                  linestyle=:dot, 
                  linealpha=1,
                  alpha=0.3,
                  markersize=3,
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
                  legend=:topright,
                  legendfont=fonts,
                  label="Empirical",
                  xlabel="Species richness",
                  ylabel="SVD-entropy")
scatter!(metrics_maxentjds.S,
            metrics_maxentjds.entropy,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
xaxis!(:log, xticks=(a,a))

plot(plotA, plotB, plotC, plotD, 
     title = ["(a)" "(b)" "(c)" "(d)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "metrics_richness.png"))


### Plot nestedness as a function of maximum trophic level ###

missing_tl = findall(ismissing, metrics_maxentjds.maxtl) # remove missing values
metrics_maxent_maxtl = metrics_maxentjds.maxtl[Not(missing_tl),:] 
metrics_maxent_rho = metrics_maxentjds.rho[Not(missing_tl),:]

scatter(metrics_emp.maxtl,
            metrics_emp.rho,
            smooth=true,
            linestyle=:dot, 
            linealpha=1,
            alpha=0.3,
            markersize=3,
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
            legend=:topright,
            legendfont=fonts,
            label="Empirical",
            xlabel="Maximum trophic level",
            ylabel="Nestedness")
scatter!(metrics_maxent_maxtl,
            metrics_maxent_rho,
            alpha=0.3, 
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
      
savefig(joinpath("figures", "maxtrophiclevel_nestedness.png"))


### Divergence between degree sequences (MaxEnt vs empirical networks) ###

# Divergence in degree sequence and species richness
plotA = scatter(metrics_emp.S,
                  metrics_diffjds.MSD_ds,
                  alpha=0.3,
                  markersize=3,
                  framestyle=:box, 
                  reg=true,
                  grid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  label="",
                  xlabel="Species richness",
                  ylabel="MSD of degree sequence",
                  ylims=(0,20.5))
xaxis!(:log, xticks=(a,a))


# Divergence in degree sequence and SVD-entropy
plotB = scatter(metrics_emp.entropy,
                  metrics_diffjds.MSD_ds,
                  alpha=0.3,
                  markersize=3,
                  framestyle=:box, 
                  reg=true,
                  grid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
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
plotA = scatter(metrics_emp.S,
                  metrics_diffjds.entropy_diff,
                  alpha=0.3,
                  markersize=3,
                  framestyle=:box, 
                  reg=true,
                  grid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  label="",
                  xlabel="Species richness",
                  ylabel="\\Delta SVD-entropy")
xaxis!(:log, xticks=(a,a))

b = [10,100,1000,10000] # specified x-ticks 
# Difference in SVD-entropy and number of links
plotB = scatter(metrics_emp.L,
                  metrics_diffjds.entropy_diff,
                  alpha=0.3,
                  markersize=3,
                  framestyle=:box, 
                  reg=true,
                  grid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
                  label="",
                  xlabel="Number of links",
                  ylabel="\\Delta SVD-entropy")
xaxis!(:log, xticks=(b,b))

# Difference in SVD-entropy and connectance
plotC = scatter(metrics_emp.C,
                  metrics_diffjds.entropy_diff,
                  alpha=0.3,
                  markersize=3,
                  framestyle=:box, 
                  reg=true,
                  grid=false,
                  dpi=1000, 
                  size=(800,500), 
                  margin=5Plots.mm, 
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  legendfont=fonts,
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
plotA = density(metrics_emp.entropy,
                  linesize=3,
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
                  legendfont=fonts,
                  legend=:topleft,
                  label="Empirical",
                  xlabel="SVD-entropy",
                  ylabel="Density",
                  ylim=(0,20.5))
density!(metrics_maxentjds.entropy,
            linesize=3,
            label="MaxEnt")

# compute z-scores
entropy_maxent_avg = mean(metrics_maxentjds.entropy)
entropy_maxent_std = std(metrics_maxentjds.entropy)

entropy_emp_zscores = (metrics_emp.entropy .- entropy_maxent_avg) ./ entropy_maxent_std

entropy_emp_zscores_500 = quantile(entropy_emp_zscores, 0.500)
entropy_emp_zscores_015 = quantile(entropy_emp_zscores, 0.015)
entropy_emp_zscores_985 = quantile(entropy_emp_zscores, 0.985)

plotB = density(entropy_emp_zscores, 
                  color=:grey,
                  linesize=3,
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
                  legendfont=fonts,
                  label="",
                  xlabel="z-score of SVD-entropy",
                  ylabel="Density",
                  ylim=(0,0.27))
plot!([entropy_emp_zscores_500], 
            seriestype=:vline, 
            color=:grey, 
            ls=:dash, 
            lab="")

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "entropy_distribution.png"))


### Difference in nestedness and SVD-entropy ###

plotA = scatter(metrics_diffjds.entropy_diff,
            metrics_diffjds.rho_diff,
            alpha=0.3,
            markersize=3,
            smooth=true,
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
            legendfont=fonts,
            legend=:topleft,
            label="",
            xlabel="\\Delta SVD-entropy",
            ylabel="\\Delta nestedness")

plotB = scatter(metrics_diffjds.entropy_diff,
                  metrics_diffjds.jaccard,
                  alpha=0.3,
                  markersize=3,
                  smooth=true,
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
                  legendfont=fonts,
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
motifs_emp = select(metrics_emp, vcat([motifs[i] for i in 1:5]))
motifs_emp = DataFrames.stack(motifs_emp)
motifs_emp = dropmissing(motifs_emp)

# select proportion of each motifs in MaxEnt networks and tidy data frame
motifs_maxentjds = select(metrics_maxentjds, vcat([motifs[i] for i in 1:5]))
motifs_maxentjds = DataFrames.stack(motifs_maxentjds)
motifs_maxentjds = dropmissing(motifs_maxentjds)

motifs_maxentco = select(metrics_maxentco, vcat([motifs[i] for i in 1:5]))
motifs_maxentco = DataFrames.stack(motifs_maxentco)
motifs_maxentco = dropmissing(motifs_maxentco)

motifs_nulljds = select(metrics_nulljds_all, vcat([motifs[i] for i in 1:5]))
motifs_nulljds = DataFrames.stack(motifs_nulljds)
motifs_nulljds = dropmissing(motifs_nulljds)

motifs_nullco = select(metrics_nullco_all, vcat([motifs[i] for i in 1:5]))
motifs_nullco = DataFrames.stack(motifs_nullco)
motifs_nullco = dropmissing(motifs_nullco)

# join both datasets and label them
groups = vcat(fill("Empirical", nrow(motifs_emp)),
              fill("MaxEnt-co", nrow(motifs_maxentco)),
              fill("MaxEnt-jds", nrow(motifs_maxentjds)),
              fill("Null 1", nrow(motifs_nullco)),
              fill("Null 2", nrow(motifs_nulljds)))

motifs_emp_maxent_neutral = [motifs_emp; 
                              motifs_maxentco; 
                              motifs_maxentjds; 
                              motifs_nullco;
                              motifs_nulljds]
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
                  grid=:none,
                  guidefont=fonts, 
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  legendfont=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  xminorgrid=false,
                  ylims=(0,1),
                  xaxis="Motifs", 
                  yaxis="Proportion")

savefig(joinpath("figures", "motifs_distribution.png"))



#### Motifs S1 vs S2 ####

plotA = scatter(metrics_emp.S1,
                  metrics_emp.S2,
                  alpha=0.3,
                  markersize=3,
                  smooth=true,
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
                  legendfont=fonts,
                  legend=:topleft,
                  label="Empirical",
                  xlabel="Proportion S1",
                  ylabel="Proportion S2",
                  xlims=(0,1),
                  ylims=(0,1))

scatter!(metrics_maxentco.S1,
                  metrics_maxentco.S2,
                  alpha=0.3, 
                  label="MaxEnt-co", 
                  smooth=true, 
                  markersize=3, 
                  linestyle=:dot, 
                  linealpha=1)
                  
scatter!(metrics_maxentjds.S1,
                  metrics_maxentjds.S2,
                  alpha=0.3, 
                  label="MaxEnt-jds", 
                  smooth=true, 
                  markersize=3, 
                  linestyle=:dot, 
                  linealpha=1)

scatter!(metrics_nullco_all.S1,
                  metrics_nullco_all.S2,
                  alpha=0.3, 
                  label="Null 1", 
                  smooth=true, 
                  markersize=3, 
                  linestyle=:dot, 
                  linealpha=1)
                  
scatter!(metrics_nulljds_all.S1,
                  metrics_nulljds_all.S2,
                  alpha=0.3, 
                  label="Null 2", 
                  smooth=true, 
                  markersize=3, 
                  linestyle=:dot, 
                  linealpha=1)

              
plotB = scatter(metrics_emp.S4,
                  metrics_emp.S5,
                  alpha=0.3,
                  markersize=3,
                  smooth=true,
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
                  legendfont=fonts,
                  legend=:topright,
                  label="Empirical",
                  xlabel="Proportion S4",
                  ylabel="Proportion S5",
                  xlims=(0,1),
                  ylims=(0,1))

scatter!(metrics_maxentco.S4, 
            metrics_maxentco.S5,
            alpha=0.3, 
            label="MaxEnt-co", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)

scatter!(metrics_maxentjds.S4,
            metrics_maxentjds.S5,
            alpha=0.3, 
            label="MaxEnt-jds", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)

scatter!(metrics_nullco_all.S4,
            metrics_nullco_all.S5,
            alpha=0.3, 
            label="Null 1", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)

scatter!(metrics_nulljds_all.S4,
            metrics_nulljds_all.S5,
            alpha=0.3, 
            label="Null 2", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "motifs_relations.png"))
