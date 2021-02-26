# Here we estimate the parameters of the flexible links model using Bayesian inference
# We also predict the number of links for a range of numbers of species

## Read metadata of all food webs archived on mangal.io

mangal_foodwebs = CSV.read(joinpath("data", "mangal_foodwebs.csv"), DataFrame)


## Define flexible links model and infer model parameters

# S: number of species in the network
# R: number of flexible links in the network (links above the minimum)

@model FL(S,R) = begin
  N = length(S)
  # Number of trials
  F = S.*S .- (S.-1)
  # Parameters (priors)
  ϕ ~ Normal(3.0, 0.5)
  μ ~ Beta(3.0, 7.0)
  # Flexible links model
  for i in 1:N
    R[i] ~ BetaBinomial(F[i], μ*exp(ϕ), (1-μ)*exp(ϕ))
  end
  return μ, ϕ
end

# Observations
S = vec(mangal_foodwebs[:,:S])
L = vec(mangal_foodwebs[:,:L])
R = L .- (S.-1)

# Hamiltonian Monte Carlo sampler with static trajectory
chain = sample(FL(S,R), HMC(0.01,10), 3000)

# Diagnostic plot
plot(chain[200:end,:,:])


## Simulate counterfactuals

# Predict the number of links from a given number of species
function prediction(chain, S)
  p = get_params(chain[200:end,:,:])
  i = rand(1:length(p.μ))
  μ, ϕ = p.μ[i], p.ϕ[i]
  return rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))) + (S-1)
end

# Simulate the number of links for a range of species richness
sp = collect(5:1000) # numbers of species
nsp = length(sp)
nsim = 1000 # number of simulations

predicted_links = DataFrame()

for i in 1:nsim
  colname = "Lhat$i"
  Lhat = [prediction(chain, S) for S in sp]
  predicted_links[!, colname] = Lhat
end

predicted_links[!, :S] = sp


## Write file

CSV.write(joinpath("data", "predicted_links.csv"), predicted_links)
