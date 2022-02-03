# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)


## Read empirical networks

N_mangal = load(joinpath("data", "proc", "mangal", "networks_mangal.jld"))["data"]

N_NZ = load(joinpath("data", "proc", "new_zealand", "networks_NZ.jld"))["data"]

N_tuesday = load(joinpath("data", "proc", "tuesday_lake", "networks_tuesday.jld"))["data"]

## Read joint degree distributions (MaxEnt)

jdd_maxent_mangal = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_mangal.jld"))["data"]

jdd_maxent_NZ = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_NZ.jld"))["data"]

jdd_maxent_tuesday = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_dist_tuesday.jld"))["data"]

k_maxent_all = load(joinpath("data", "sim", "joint_degree_dist_maxent", "joint_degree_sequence_all.jld"))["data"]

## Read tables of network properties

metrics = DataFrame(CSV.File(joinpath("results", "metrics.csv")))
gmetrics = DataFrame(CSV.File(joinpath("results", "gmetrics.csv")))

## Read tables of differences in network properties

metrics_diff = DataFrame(CSV.File(joinpath("results", "metrics_diff.csv")))

# empirical networks only
metrics_emp = metrics[in(vcat("N_mangal", "N_NZ", "N_tuesday")).(metrics[!,:network]),:]

# MaxEnt networks only
metrics_maxent = metrics[in(vcat("N_maxent_mangal", "N_maxent_NZ", "N_maxent_tuesday")).(metrics[!,:network]),:]

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

N_emp = simplify.(vcat(N_mangal, N_NZ, N_tuesday))

kout_emp = reduce(vcat, collect.(values.(degree.(N_emp, dims=1))))
kin_emp = reduce(vcat, collect.(values.(degree.(N_emp, dims=2))))


## simulated networks (from the joint degree distribution of maximum entropy)

function get_kin(degrees::Vector)
      # retrieve kin from the simulated degree sequence
      kin_maxent = []
      for i in 1:length(degrees)
            push!(kin_maxent, degrees[i].kin)
      end
      return reduce(vcat, kin_maxent)
end

function get_kout(degrees::Vector)
      # retrieve kout from the simulated degree sequence
      kout_maxent = []
      for i in 1:length(degrees)
            push!(kout_maxent, degrees[i].kout)
      end
      return reduce(vcat, kout_maxent)
end

# all simulated networks
kin_maxent = get_kin(k_maxent_all)
kout_maxent = get_kout(k_maxent_all)

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
xaxis!(xlim=(-1, maximum(vcat(kout_emp, kout_maxent))), "Kout (number of preys)")
yaxis!(ylim=(-1, maximum(vcat(kin_emp, kin_maxent))), "Kin (number of predators)")
                
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
xaxis!(xlim=(-1, maximum(vcat(kout_emp, kout_maxent))), "Kout (number of preys)")
yaxis!(ylim=(-1, maximum(vcat(kin_emp, kin_maxent))), "Kin (number of predators)")


l = @layout [a [b ; c]]

plot(plotA, plotB, plotC,
    layout = l,
    title = ["(a)" "(b) empirical data" "(c) simulated data"],
    titleloc=:right, 
    titlefont=fonts)

savefig(joinpath("figures","joint_degree_dist.png"))


### Measures of empirical and maximum entropy food webs ###

# Nestedness
plotA = scatter(metrics_emp.rho, 
                  metrics_maxent.rho,
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
                  metrics_maxent.maxtl,
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
                  metrics_maxent.diam,
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
                  metrics_maxent.entropy,
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
scatter!(metrics_maxent.S, 
            metrics_maxent.rho,
            alpha=0.3,
            label="MaxEnt", 
            smooth=true, 
            markersize=3, 
            linestyle=:dot, 
            linealpha=1)
xaxis!(:log, xticks=(a,a))

# Maximum trophic level and species richness

missing_tl = findall(ismissing, metrics_maxent.maxtl) # remove missing values
metrics_maxent_maxtl = metrics_maxent.maxtl[Not(missing_tl),:] 
metrics_maxent_S = metrics_maxent.S[Not(missing_tl),:]

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
scatter!(metrics_maxent.S,
            metrics_maxent.diam,
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
scatter!(metrics_maxent.S,
            metrics_maxent.entropy,
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

missing_tl = findall(ismissing, metrics_maxent.maxtl) # remove missing values
metrics_maxent_maxtl = metrics_maxent.maxtl[Not(missing_tl),:] 
metrics_maxent_rho = metrics_maxent.rho[Not(missing_tl),:]

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
                  metrics_diff.MSD_ds,
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
                  metrics_diff.MSD_ds,
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
                  metrics_diff.entropy_diff,
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
                  ylabel="SVD-entropy (difference)")
xaxis!(:log, xticks=(a,a))

b = [10,100,1000,10000] # specified x-ticks 
# Difference in SVD-entropy and number of links
plotB = scatter(metrics_emp.L,
                  metrics_diff.entropy_diff,
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
                  ylabel="SVD-entropy (difference)")
xaxis!(:log, xticks=(b,b))

# Difference in SVD-entropy and connectance
plotC = scatter(metrics_emp.C,
                  metrics_diff.entropy_diff,
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
                  ylabel="SVD-entropy (difference)")

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
density!(metrics_maxent.entropy,
            linesize=3,
            label="MaxEnt")

# compute z-scores
entropy_maxent_avg = mean(metrics_maxent.entropy)
entropy_maxent_std = std(metrics_maxent.entropy)

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

plotA = scatter(metrics_diff.entropy_diff,
            metrics_diff.rho_diff,
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
            xlabel="SVD-entropy (difference)",
            ylabel="Nestedness (difference)")

plotB = scatter(metrics_diff.entropy_diff,
                  metrics_diff.jaccard,
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
                  xlabel="SVD-entropy (difference)",
                  ylabel="Jaccard distance")

plot(plotA, plotB,
     title = ["(a)" "(b)"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "difference_entropy_jaccard.png"))


### Motif distribution ###

# get motifs names
motifs = keys(unipartitemotifs())

# select proportion of each motifs in empirical networks and tidy data frame
motifs_emp = select(metrics_emp, vcat([motifs[i] for i in 1:length(motifs)]))
motifs_emp = DataFrames.stack(motifs_emp)
motifs_emp = dropmissing(motifs_emp)

# select proportion of each motifs in MaxEnt networks and tidy data frame
motifs_maxent = select(metrics_maxent, vcat([motifs[i] for i in 1:length(motifs)]))
motifs_maxent = DataFrames.stack(motifs_maxent)
motifs_maxent = dropmissing(motifs_maxent)

# join both datasets and label them
groups = vcat(fill("Empirical", nrow(motifs_emp)),
              fill("MaxEnt", nrow(motifs_emp)))

motifs_emp_maxent = [motifs_emp; motifs_maxent]
insertcols!(motifs_emp_maxent, :group => groups)

# motif distribution
groupedboxplot(motifs_emp_maxent.variable, 
                  motifs_emp_maxent.value,
                  group=motifs_emp_maxent.group, 
                  alpha=0.9,
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
                  legend_font=fonts,
                  foreground_color_legend=nothing, 
                  background_color_legend=:white, 
                  xminorgrid=false,
                  ylims=(0,1),
                  xaxis="Motifs", 
                  yaxis="Proportion")

savefig(joinpath("figures", "motifs_distribution.png"))

