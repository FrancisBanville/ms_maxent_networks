# Here we predict the number of links in food webs from their number of species (flexible links model)
# We then predict their joint degree distributions from their number of species and number of links (MaxEnt)
# We then predict their adjacency matrix from their joint degree distribution and select the one with maximum entropy

## Read metadata of all food webs archived on mangal.io

mangal_foodwebs = DataFrame(CSV.File(joinpath("data", "mangal_foodwebs.csv")))


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

# Predict the number of links for a given number of species
function predict_links(chain, S::Int64, nsim::Int64)
  p = get_params(chain[200:end,:,:])
  i = rand(1:length(p.μ))
  μ, ϕ = p.μ[i], p.ϕ[i]
  return rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ)), nsim) .+ (S-1)
end

# Simulate the number of links for a range of species richness
sp = collect(5:1000) # numbers of species
nsp = length(sp)
nsim = 1000 # number of simulations

predicted_links = zeros(Int64, nsim, length(sp))

p = Progress(length(sp))
Threads.@threads for s in sp
    i = findall(sp .== s)
    predicted_links[:, i] = predict_links(chain, s, nsim)
    next!(p)
end



## Predict the degree distribution of maximum entropy (approximation of the analytical solution)

function p_k(S::Int64, L::Int64, k::Int64)
  # S: Number of species
  # L: Number of links
  # k: Species degree

  # Mean degree constraint
  kavg = 2*L/S 

  # Degree distribution of maximum entropy
  c = 1/(kavg - 1)
  r = (kavg - 1)/kavg
  p_k = c*(r^k)

  # Return the probability that a species has k links
  return p_k
end


function dd_maxent_prob(S::Int64, L::Int64)
  # S: Number of species
  # L: Number of interactions

  dd_maxent_prob = [p_k(S, L, k) for k in 1:S]

  # Return the degree distribution (probabilistic)
  return dd_maxent_prob
end


function dd_maxent(S::Int64, L::Int64)
  # S: Number of species
  # L: Number of interactions

  # Predict the degree distribution (probabilistic)
  dd_maxent_proba = dd_maxent_prob(S, L)

  # Cumulative number of species expected to have degree k or lower
  dd_maxent_cum = cumsum(dd_maxent_proba) .* S

  # Number of species expected to have degree k
  dd_maxent = zeros(Int64, S)

  K = 1:S

  for k in K
    dd_maxent[k] = convert(Int64, round(dd_maxent_cum[k] - sum(dd_maxent)))
  end
  
  # Return the expected degree for all species (ordered)
  dd_maxent = vcat(fill.(K, dd_maxent)...)

  return dd_maxent

end


# Example 1
S = 110
L = predict_links(chain, S, nsim)
L = convert(Int64, round(median(L)))
ddist = dd_maxent(S, L)
sum(ddist)
L*2


## Predict the joint degree distribution of maximum entropy (analytical approximation)

function p_kin_kout(S::Int64, L::Int64, kin::Int64, kout::Int64)
  # S: Number of species
  # L: Number of links
  # kin: Species in-degree
  # Kout: Species out-degree

  # A species cannot have a degree of k = 0
  if (kin + kout == 0 || kin + kout > S)
    p_kin_kout = 0
  else

  # Mean degree constraint
  kavg = 2L/S 

  # Joint degree distribution of maximum entropy
  k = kin + kout
  c = 1/(kavg - 1)
  r = (kavg - 1)/kavg
  p_kin_kout = binomial(BigInt(k), BigInt(kin))*c*(r/2)^k
 
   # Return the probability that a species has a kin predators and kout preys
   return p_kin_kout

  end
 end


 function jdd_maxent_prob(S::Int64, L::Int64)
  # S: Number of species
  # L: Number of interactions

  jdd_maxent_prob = zeros(Float64, S+1, S+1)
  
  for kin in 0:S, kout in 0:S
    jdd_maxent_prob[kin+1, kout+1] = p_kin_kout(S, L, kin, kout)
  end
  
  # Return the joint degree distribution (probabilistic)
  return jdd_maxent_prob
end



function jdd_maxent(S::Int64, L::Int64)
  # S: Number of species
  # L: Number of interactions

  # Predict the degree distribution (probabilistic)
  jdd_maxent_proba = jdd_maxent_prob(S, L)

  # Cumulative number of species expected to have an in-degree kin or lower
  jdd_maxent_cum = cumsum(sum.(eachcol(jdd_maxent_proba))) .* S

  # Number of species expected to have an in-degree k
  jdd_maxent = zeros(Int64, S+1)

  Kin = 0:S

  for kin in Kin
    jdd_maxent[kin+1] = convert(Int64, round(jdd_maxent_cum[kin+1] - sum(jdd_maxent)))
  end
  
  # Return the expected in-degree for all species (ordered)
  jdd_maxent = vcat(fill.(Kin, jdd_maxent)...)

  return jdd_maxent

end


# Example 2
S = 225
L = predict_links(chain, S, nsim)
L = convert(Int64, round(median(L)))
ddist = dd_maxent(S, L)
jddin = jdd_maxent(S, L)
jddout = ddist .- jddin




