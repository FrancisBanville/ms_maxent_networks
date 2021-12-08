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













################## OLD FIGURES ###############################

## Figure: Divergence between numerical and analytical solutions of the MaxEnt degree distributions

## Load metadata of all food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl)
mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "mangal_foodwebs.csv")))

# Load predicted links (generated from 02_predict_networks.jl)
# Columns: S from 5 to 1000, Rows: 1000 simulations
predicted_links = load(joinpath("data", "sim", "predicted_links.jld"))["data"]

# Compute average degree for all values of species richness and median predicted numbers of links
sp = collect(5:1000) 
L_med = convert.(Int64, round.(median.(eachcol(predicted_links))))
kavg = 2 .* L_med ./ sp

"""
estimate_lambda(S::Int64, kavg::Float64) 
    S: number of species
    kavg: mean degree constraint (2L/S)
Returns the numerical value of the Lagrange multiplier 
"""
function estimate_lambda(S::Int64, kavg::Float64) 
  model = Model(Ipopt.Optimizer)
  set_silent(model)

  @variable(model, x >= 0)  # lambdas (x) are positive real numbers 
  @NLobjective(model, Max, x) # could have been another objective function
  @NLconstraint(model, sum(k * exp(-x * k) for k = 1:S) / sum(exp(-x * k) for k = 1:S) == kavg)
  optimize!(model)
  return value.(x)
end

"""
compute_divergence(S::Int64, kavg::Float64)
    S: number of species
    kavg: mean degree constraint (2L/S)
Returns the divergence (sum of absolute differences) between the numerical and analytical solutions of the MaxEnt degree distributions
"""
function compute_divergence(S::Int64, kavg::Float64)
  lambda = estimate_lambda(S, kavg)
  pk_numerical = [exp(-lambda * k) for k in 1:S] ./ sum([exp(-lambda * k) for k = 1:S])

  c = 1 / (kavg - 1)
  r = (kavg - 1) / kavg
  pk_analytical = (c * r^k for k in 1:S)

  divergence = sum(abs.(pk_numerical .- pk_analytical))
  return divergence
end

# Compute divergence for all values of species richness
divergence = compute_divergence.(sp, kavg) 

# Compute quantiles of species richness for all networks archived on Mangal (to be plotted)
species500 = quantile(mangal_foodwebs.S, 0.5)
species015 = quantile(mangal_foodwebs.S, 0.015)
species985 = quantile(mangal_foodwebs.S, 0.985)

# Plot divergence between numerical and analytical solutions of the MaxEnt degree distributions
scatter(sp[1:150], divergence[1:150], alpha=0.4, lab="", 
  framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
  guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!([species500], seriestype=:vline, color=:grey, ls=:dash, lab="", ylim=(0,0.12))
plot!([species015, species985], [0.12, 0.12], fill=(0, :grey, 0.12), c=:transparent, lab="")
plot!([species015], seriestype=:vline, color=:grey, ls=:dot, lab="")
plot!([species985], seriestype=:vline, color=:grey, ls=:dot, lab="")
xaxis!(xlabel="Species richness")
yaxis!(ylabel="Divergence")

savefig(joinpath("figures", "divergence_numerical_analytical.png"))


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


## Figure: MaxEnt and empirical degree distributions

# Get Mangal food webs with given numbers of species (median and 97% PI)
species985 = 81 # no network with 85 species (but some with 81 and 89 species)

Ns_mangal500 = mangal_foodwebs[mangal_foodwebs.S .== species500, :]
Ns_mangal015 = mangal_foodwebs[mangal_foodwebs.S .== species015, :]
Ns_mangal985 = mangal_foodwebs[mangal_foodwebs.S .== species985, :] 

# Get food webs with minimum and maximum numbers of links (from food webs with median and 97% PI of species richness) (arbitrary criteria to select networks that will be plotted)
N_mangal500min = Ns_mangal500[Ns_mangal500.L .== minimum(Ns_mangal500.L), :][1,:]
N_mangal500max = Ns_mangal500[Ns_mangal500.L .== maximum(Ns_mangal500.L), :][1,:]
N_mangal015min = Ns_mangal015[Ns_mangal015.L .== minimum(Ns_mangal015.L), :][2,:] # 2 food webs with 7 species and 8 links (get the one with predatory links)
N_mangal015max = Ns_mangal015[Ns_mangal015.L .== maximum(Ns_mangal015.L), :][1,:]
N_mangal985min = Ns_mangal985[Ns_mangal985.L .== minimum(Ns_mangal985.L), :][1,:]
N_mangal985max = Ns_mangal985[Ns_mangal985.L .== maximum(Ns_mangal985.L), :][1,:]

"""
get_ranked_dd(N::DataFrameRow)
    N: Data frame row of Mangal metadata
Returns the ranked degree distribution of the input network
"""
function get_ranked_dd(N::DataFrameRow)
    # Convert to UnipartiteNetwork
    N = convert(UnipartiteNetwork, network(N.id))
    # Get sorted degree distribution 
    dd = sort(collect(values(degree(N))), rev=true)
    return dd
end

# Get degree distribution for all selcted food webs
dd_mangal500min = get_ranked_dd(N_mangal500min)
dd_mangal500max = get_ranked_dd(N_mangal500max)
dd_mangal015min = get_ranked_dd(N_mangal015min)
dd_mangal015max = get_ranked_dd(N_mangal015max)
dd_mangal985min = get_ranked_dd(N_mangal985min)


# Get simulated networks with simulated numbers of links (median of flexible links model)
Ns_maxent_fl = load(joinpath("data", "sim", "predicted_networks.jld"))["data"]

# Get simulated networks with same numbers of species as the food webs previously selected
S_maxent_fl = richness.(Ns_maxent_fl)

N_maxent_fl500 = Ns_maxent_fl[S_maxent_fl .== species500][1]
N_maxent_fl015 = Ns_maxent_fl[S_maxent_fl .== species015][1]
N_maxent_fl985 = Ns_maxent_fl[S_maxent_fl .== species985][1]

"""
get_ranked_dd(N::UnipartiteNetwork)
    N: Unipartite Network
Returns the ranked degree distribution of the input network
"""
function get_ranked_dd(N::UnipartiteNetwork)
  # Get sorted degree distribution 
  dd = sort(collect(values(degree(N))), rev=true)
  return dd
end

# Get degree distributions of all selected simulated networks
dd_maxent_fl500 = get_ranked_dd(N_maxent_fl500)
dd_maxent_fl015 = get_ranked_dd(N_maxent_fl015)
dd_maxent_fl985 = get_ranked_dd(N_maxent_fl985)


# Get simulated networks with the same numbers of links as food webs on Mangal
Ns_maxent_empL = load(joinpath("data", "sim", "predicted_networks_empL.jld"))["data"]

"""
get_network_empL(S::Int64, fun::String) 
    S: Species richness
    fun: Either "min" or "max"
Returns the simulated network of S species with the lowest (fun = "min") or higest (fun = "max") numbers of links
"""
function get_network_empL(S::Int64, fun::String) 
  # Get simulated networks with S species
  Ns_S = Ns_maxent_empL[richness.(Ns_maxent_empL) .== S]
  
  # Get the network with the lowest (min) or highest (max) numbers of links
  if fun == "min"
    N = Ns_S[links.(Ns_S) .== minimum(links.(Ns_S))][1]
  elseif fun == "max"
    N = Ns_S[links.(Ns_S) .== maximum(links.(Ns_S))][1]
  else 
    message("The function must be min or max")
  end
  return N
end

# Get simulated networks with same numbers of species and links as previously selected food webs
N_maxent_empL500min = get_network_empL(species500, "min")
N_maxent_empL500max = get_network_empL(species500, "max")
N_maxent_empL015min = get_network_empL(species015, "min")
N_maxent_empL015max = get_network_empL(species015, "max")
N_maxent_empL985min = get_network_empL(species985, "min")

# Get degree distributions of all selected simulated networks
dd_maxent_empL500min = get_ranked_dd(N_maxent_empL500min)
dd_maxent_empL500max = get_ranked_dd(N_maxent_empL500max)
dd_maxent_empL015min = get_ranked_dd(N_maxent_empL015min)
dd_maxent_empL015max = get_ranked_dd(N_maxent_empL015max)
dd_maxent_empL985min = get_ranked_dd(N_maxent_empL985min)


"""
plot_dd(mangalmin::T, mangalmax::T, empLmin::T, empLmax::T, fl::T) where {T <: Vector{Int64}}
    mangalmin: Degree distribution of Mangal food web of S species with lowest number of links
    mangalmax: Degree distribution of Mangal food web of S species and highest number of links
    empLmin: Degree distribution of simulated food web of S species and lowest number of links (same number of links as Mangal food webs)
    empLmax: Degree distribution of simulated food web of S species and highest number of links (same number of links as Mangal food webs)
    fl: Degree distribution of simulated food web of S and number of links given by the median of the flexible links model
Returns a plot of the degree distributions of Mangal and empirical food webs with S species
"""
function plot_dd(mangalmin::T, mangalmax::T, empLmin::T, empLmax::T, fl::T) where {T <: Vector{Int64}}
    plot(mangalmin, color=:black, alpha=0.7, linewidth=2, linestyle=:dot, label="Mangal", 
      framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
      guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
      foreground_color_legend=nothing, background_color_legend=:white, legendfont=fonts,
      xlab="Rank", ylab="Degree k")
  plot!(mangalmax, color=:darkblue, alpha=0.7, linewidth=2, linestyle=:dot, label="")
  plot!(empLmin, color=:black, alpha=0.7, linewidth=2, linestyle=:solid, label="MaxEnt")
  plot!(empLmax, color=:darkblue, alpha=0.7, linewidth=2, linestyle=:solid, label="")
  plot!(fl, color=:darkred, alpha=0.7, linewidth=2, linestyle=:solid, label="")
end

# Plot the degree distributions of Mangal and empirical food webs with median and 97% PI numbers of species
plotA = plot_dd(dd_mangal015min, dd_mangal015max, dd_maxent_empL015min, dd_maxent_empL015max, dd_maxent_fl015)
plotB = plot_dd(dd_mangal500min, dd_mangal500max, dd_maxent_empL500min, dd_maxent_empL500max, dd_maxent_fl500)
plotC = plot_dd(dd_mangal985min, dd_mangal985min, dd_maxent_empL985min, dd_maxent_empL985min, dd_maxent_fl985)

plot(plotA, plotB, plotC, 
    layout=(1,3),
    title = ["(a) S=$(species015) species" "(b) S=$(species500) species" "(c) S=$(species985) species"],
    titleloc=:right, titlefont=fonts)

savefig(joinpath("figures","mangal_maxent_degree_distributions"))



## Figures: Measures of empirical and maximum entropy food webs 

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

# Measure SVD-entropy
"""
svd_entropy(N::T) where {T <: UnipartiteNetwork}
    N: unipartite network
Returns the entropy of the adjacency matrix of a deterministic network using the singular values from a Singualar Value Decomposition
"""
function svd_entropy(N::T) where {T <: UnipartiteNetwork}
    A = convert(Matrix, N.edges)
    F = svd(A)
    Λ = F.S[1:rank(A)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

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


# Plot measures as a function of species richness
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


## Figure: Motif distribution

"""
count_motifs(N::T) where {T <: UnipartiteNetwork}
    N: unipartite network
Returns the number of times each motif was found in the unipartite network
"""
function count_motifs(N::T) where {T <: UnipartiteNetwork}
  ms = unipartitemotifs() # list motifs
  ms_count = [length(find_motif(N, ms[i])) for i in 1:13]
end

"""
count_motifs(Ns::T) where {T <: Union{Vector{UnipartiteNetwork{Bool, String}}, Vector{UnipartiteNetwork{Bool, MangalNode}}}}
    Ns: vector of unipartite networks
Returns the proportion of each motif in each unipartite network
Rows are networks and columns are different motifs
"""
function count_motifs(Ns::T) where {T <: Union{Vector{UnipartiteNetwork{Bool, String}}, Vector{UnipartiteNetwork{Bool, MangalNode}}}}
  # remove the biggest network (too long to find all motifs in networks that are too large)
  Ns_small = Ns[richness.(Ns) .!== maximum(richness.(Ns))]
  # count each of the 13 possible motifs for all networks in Ns
  motifs_count = zeros(Float64, length(Ns_small), length(unipartitemotifs()))
  
  for i in 1:length(Ns_small)
    motifs_count[i,:] = count_motifs(Ns_small[i]) 
  end

  # returns the proportion of each motif for all networks
  return motifs_count ./ sum.(eachrow(motifs_count))
end

# Get the proportion of motifs for all food webs (empirical and predicted)
motifs_prop_mangal = count_motifs(Ns_mangal)
motifs_prop_maxent_empL = count_motifs(Ns_maxent_empL)
motifs_prop_maxent_fl = count_motifs(Ns_maxent_fl)

"""
plot_motifs(ms::T) where {T <: Matrix{Float64}}
    ms: matrix of proportions of motifs in networks (rows = networks, columns = motifs)
Returns the plot of the distribution of all motifs in a set of food webs 
Gives a violon plot, boxplot and dotplot
"""
function plot_motifs(ms::T) where {T <: Matrix{Float64}}
  
  # change data format for plotting
  ms_df = convert(DataFrame, ms)
  rename!(ms_df, vcat(["S$i" for i in 1:5],["D$j" for j in 1:8]))
  ms_df = DataFrames.stack(ms_df)
  
  # plot proportions of all modules
  @df ms_df violin(string.(:variable), :value, linewidth=0, label="",
              framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
              guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
              ylims=(0,1),
              xaxis="Motifs", yaxis="Proportion")
  @df ms_df boxplot!(string.(:variable), :value, fillalpha=0.2, markersize=3, linewidth=1, color=:grey, label="")
  @df ms_df dotplot!(string.(:variable), :value, color=:black, markersize=1.5, alpha=0.15, label="")
end

# Plot the proportions of all motifs for empirical and predicted food webs 
# Food webs with median L are not plotted since their distribution is very similar to the one of MaxEnt food webs with empirical L
plotA = plot_motifs(motifs_prop_mangal)
plotB = plot_motifs(motifs_prop_maxent_empL)

plot(plotA, plotB,
     title = ["(a) Empirical food webs" "(b) MaxEnt food webs"],
     titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "motifs_distribution.png"))