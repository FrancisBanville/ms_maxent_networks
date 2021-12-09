# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)


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
        colorbar_title="log probability",
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

savefig(joinpath("figures", "prop_disconnected_species.png"))




################# TK TO DO #########################
## Figure: Mean degree constraints and MaxEnt degree distributions
species015 = Int64(round(species015))
species500 = Int64(round(species500))
species985 = Int64(round(species985))

"""
get_avgk_dist(S::Int64)
    S: number of species
Returns the distribution of mean degrees according to the predictions of the flexible links model
"""
function get_avgk_dist(S::Int64)
  L = predicted_links[:, S-4]
  kavg = 2 .* L ./ S 
  return kavg
end

# Get distribution of mean degrees from species richness
kavg015 = get_avgk_dist(species015)
kavg500 = get_avgk_dist(species500)
kavg985 = get_avgk_dist(species985)

# Plot distributions of mean degrees of 3 different levels of species richness
plotA = density(kavg015, linewidth=2, label="$(species015) species",
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white, legendfont=fonts)
density!(kavg500, linewidth=2, label="$(species500) species")
density!(kavg985, linewidth=2, label="$(species985) species")
xaxis!(xlabel="Mean degree (links/species)", xlims=(0,maximum(kavg985)))
yaxis!(ylabel="Density", ylims=(0,0.55))

"""
dd_maxent_prob(S::Int64, q::Float64)
    S: number of species
    q: quantile for the number of links using the predictions of the flexible links model 
Returns the degree distribution of MaxEnt (analytical solution) for a network of S species and qth quantile of the number of links from the predictions of the flexible links model 
"""
function dd_maxent_prob(S::Int64, q::Float64)
  L = Int64(round(quantile(predicted_links[:, S-4], q)))
  dd_maxent = dd_maxent_prob(S, L)
  return dd_maxent
end

# Get the degree distribution of MaxEnt (analytical solution) for a network with a median number of species and different numbers of links
dd015 = dd_maxent_prob(species500, 0.015)
dd055 = dd_maxent_prob(species500, 0.055)
dd165 = dd_maxent_prob(species500, 0.165)
dd500 = dd_maxent_prob(species500, 0.5)
dd835 = dd_maxent_prob(species500, 0.835)
dd945 = dd_maxent_prob(species500, 0.945)
dd985 = dd_maxent_prob(species500, 0.985)

# Plot different degree distributions of MaxEnt (analytical solution) for a network with a median number of species
a = [1,2,5,10,15,20,25] # specified x-ticks 

plotB = plot(dd015, color=:black, alpha=0.7, linewidth=2, linestyle=:dot, label="97% PI", # 97% PI
    framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
    guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
    foreground_color_legend=nothing, background_color_legend=:white, legendfont=fonts)
plot!(dd055, color=:black, alpha=0.7, linestyle=:dash, label="89% PI") # 89% PI
plot!(dd165, color=:black, alpha=0.7, linestyle=:solid, label="67% PI") # 67% PI
plot!(dd835, color=:grey, alpha=0.7, linestyle=:solid, label="") # 67% PI
plot!(dd945, color=:grey, alpha=0.7, linestyle=:dash, label="") # 89% PI
plot!(dd985, color=:grey, alpha=0.7, linestyle=:dot, label="") # 97% PI
plot!(dd500, color=:darkblue, linewidth=4, label="median")
xaxis!(:log, xlabel="Degree k", xticks=(a, a), xlims=(1,species500))
yaxis!(ylabel="p(k)")

plot(plotA, plotB,
    title = ["(a)" "(b)"],
    titleloc=:right, titlefont=fonts)

savefig(joinpath("figures","maxent_degree_distributions"))
###########################################################

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