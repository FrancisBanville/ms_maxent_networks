"""
neutral_model(n::Vector, L::Int64)
    n: Vector of abundances
    L: number of interactions
Returns a boolean unipartite network with sampling weights given by the product of relative abundances
"""
function neutral_model(n::Vector, L::Int64)
    rel_n = n ./ sum(n) # vector of relative abundances
    na_mat = rel_n * rel_n' # neutral abundance matrix (proportional to adjacency matrix)

    S = length(n) # species richness
    
    # draw L interactions using the neutral abundance matrix as weights
    # repeat multiple times and keep the L higher values
    N_sim = zeros(Float64, S, S)
    N_neutral = copy(N_sim)

    for i in 1:100
        N_samp = sample(eachindex(na_mat), Weights(vec(na_mat)), L, replace=false)
        N_sim[N_samp] .= N_sim[N_samp] .+ 1.0
    end

    thres = sort(vec(N_sim), rev=true)[L] # Lth maximum value
    N_neutral[N_sim .< thres] .= 0.0
    N_neutral[N_sim .> thres] .= 1.0

    i = sample(findall(N_sim .== thres), Int64(L-sum(N_neutral)), replace=false)
    N_neutral[i] .= 1.0

    # convert to unipartite network 
    N_neutral = convert(Matrix{Bool}, N_neutral)
    N_neutral = UnipartiteNetwork(N_neutral)
    return N_neutral
end


"""
nullco_model(N::UnipartiteNetwork, L::Int64, model::Int64)
    N: Unipartite directed simple network
    L: number of interactions
    model: 1 = type I null model; 2 = type II null model
Returns a boolean unipartite network with sampling probabilities given by the type I or II null model
"""
function null_model(N::UnipartiteNetwork, L::Int64, model::Int64)
    # run null model
    if model == 1
        N_null_prob = Matrix(null1(N).edges) 
    elseif model == 2
        N_null_prob = Matrix(null2(N).edges) 
    end
    
    S = richness(N) # species richness
    
    # draw L interactions using the type I null model matrix as weights
    # repeat multiple times and keep the L higher values
    N_sim = zeros(Float64, S, S)
    N_null = copy(N_sim)

    for i in 1:100
        N_samp = sample(eachindex(N_null_prob), Weights(vec(N_null_prob)), L, replace=false)
        N_sim[N_samp] .= N_sim[N_samp] .+ 1.0
    end

    thres = sort(vec(N_sim), rev=true)[L] # Lth maximum value
    N_null[N_sim .< thres] .= 0.0
    N_null[N_sim .> thres] .= 1.0

    i = sample(findall(N_sim .== thres), Int64(L-sum(N_null)), replace=false)
    N_null[i] .= 1.0 
    
    # convert to unipartite network 
    N_null = convert(Matrix{Bool}, N_null)
    N_null = UnipartiteNetwork(N_null)
    return N_null
end