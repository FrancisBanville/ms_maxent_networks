# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)


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
      foreground_color_legend=nothing, background_color_legend=:white)
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
dd015 = dd_maxent_prob(species50, 0.015)
dd055 = dd_maxent_prob(species50, 0.055)
dd165 = dd_maxent_prob(species50, 0.165)
dd500 = dd_maxent_prob(species50, 0.5)
dd835 = dd_maxent_prob(species50, 0.835)
dd945 = dd_maxent_prob(species50, 0.945)
dd985 = dd_maxent_prob(species50, 0.985)

# Plot different degree distributions of MaxEnt (analytical solution) for a network with a median number of species
a = [1,2,5,10,15,20,25] # specified x-ticks 

plotB = plot(dd015, color=:black, alpha=0.7, linewidth=2, linestyle=:dot, label="97% PI", # 97% PI
    framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
    guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
    foreground_color_legend=nothing, background_color_legend=:white)
plot!(dd055, color=:black, alpha=0.7, linestyle=:dash, label="89% PI") # 89% PI
plot!(dd165, color=:black, alpha=0.7, linestyle=:solid, label="67% PI") # 67% PI
plot!(dd835, color=:grey, alpha=0.7, linestyle=:solid, label="") # 67% PI
plot!(dd945, color=:grey, alpha=0.7, linestyle=:dash, label="") # 89% PI
plot!(dd_985, color=:grey, alpha=0.7, linestyle=:dot, label="") # 97% PI
plot!(dd500, color=:darkblue, linewidth=4, label="median")
xaxis!(:log, xlabel="Degree k", xticks=(a, a), xlims=(1,species50))
yaxis!(ylabel="p(k)")

plot(plotA, plotB)

savefig(joinpath("figures","maxent_degree_distributions"))




## Read metadata of all food webs archived on mangal.io (generated from 01_import_mangal_metadata.jl)
mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "mangal_foodwebs.csv")))

## Read predicted networks 

# Simulated networks with simulated numbers of links
predicted_networks = load(joinpath("data", "sim", "predicted_networks.jld"))["data"]

# Simulated networks with empirical numbers of links and species
predicted_networks_empL = load(joinpath("data", "sim", "predicted_networks_empL.jld"))["data"]


## Compute measures of food webs of maximum entropy (simulated numbers of links)

# richness 
predicted_networks_S = richness.(predicted_networks)
# nestedness (spectral radius of the adjacency matrix)
predicted_networks_nestedness = ρ.(predicted_networks)
# maximum and average trophic levels
predicted_networks_tls = []
working = []
for i in 1:length(predicted_networks)
  try
  push!(predicted_networks_tls, values(trophic_level(predicted_networks[i])))
  push!(working, i)
  catch
  println("could not estimate trophic levels in network $(i)")
  end
end

trophic_level(transpose(predicted_networks))
predicted_networks_tl_max = maximum.(predicted_networks_tls)
predicted_networks_tl_avg = sum.(predicted_networks_tls)./length.(predicted_networks_tls)


## Compute measures of food webs of maximum entropy (empirical numbers of links)

# richness 
predicted_networks_empL_S = richness.(predicted_networks_empL)
# nestedness (spectral radius of the adjacency matrix)
predicted_networks_empL_nestedness = ρ.(predicted_networks_empL)
# maximum and average trophic levels
predicted_networks_empL_tls = []
working_empL = []
for i in 1:length(predicted_networks_empL)
  try
  push!(predicted_networks_empL_tls, values(trophic_level(predicted_networks_empL[i])))
  push!(working_empL, i)
  catch
  println("could not estimate trophic levels in network $(i)")
  end
end

predicted_networks_empL_tl_max = maximum.(predicted_networks_empL_tls)
predicted_networks_empL_tl_avg = sum.(predicted_networks_empL_tls)./length.(predicted_networks_empL_tls)


## Measures of food webs archived on mangal.io

# Read food webs and convert them to unipartie networks 
Ns_mangal = network.(mangal_foodwebs.id)
Ns_mangal = convert.(UnipartiteNetwork, Ns_mangal)

# Compute food-web measures 
# richness 
Ns_mangal_S = richness.(Ns_mangal)
# nestedness (spectral radius of the adjacency matrix)
Ns_mangal_nestedness = ρ.(Ns_mangal)
# maximum and average trophic levels
Ns_mangal_tls = values.(trophic_level.(Ns_mangal))
Ns_mangal_tl_max = maximum.(Ns_mangal_tls)
Ns_mangal_tl_avg = sum.(Ns_mangal_tls)./length.(Ns_mangal_tls)

Ns_mangal_nestedness2 = zeros(length(Sp))
Ns_mangal_tl_max2 = zeros(length(Sp))
Ns_mangal_tl_avg = zeros(length(Sp))

for s in Sp, i in 1:length(Sp)
  j = findall(Ns_mangal_S .== s)
  Ns_mangal_nestedness2[i] = median(Ns_mangal_nestedness[j])
end 


## Figures

# Nestedness (MaxEnt and Mangal food webs)
plotA = scatter(Ns_mangal_nestedness, predicted_networks_empL_nestedness, alpha=0.5, lab="", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(0.4, 0.9, 100), LinRange(0.4, 0.9, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs (empirical number of links")

plot(plotA, plotB, layout=(1,2), legend=false,
    title=["flexible links" "empirical links"], titleloc=:right, titlefont=font(6))

savefig(joinpath("figures", "nestedness.png"))

# Maximum trophic level (MaxEnt and Mangal food webs)
scatter(Ns_mangal_tl_max[working], predicted_networks_tl_max, alpha=0.5, lab="",
title="Maximum trophic level", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(2, 13, 100), LinRange(2, 13, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs")
savefig(joinpath("figures", "maximum_trophic_level.png"))

# Average trophic level (MaxEnt and Mangal food webs)
scatter(Ns_mangal_tl_avg[working], predicted_networks_tl_avg, alpha=0.5, lab="",
title="Average trophic level", framestyle=:box, guidefont=fonts, xtickfont=fonts, ytickfont=fonts)
plot!(LinRange(1, 7, 100), LinRange(1, 7, 100), lab="", color="grey")
xaxis!("Empirical food webs")
yaxis!("MaxEnt food webs")
savefig(joinpath("figures", "average_trophic_level.png"))

