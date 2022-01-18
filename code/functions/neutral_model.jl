"""
neutral_model(A::Vector{Int64}, chain)
    n: Vector of abundances
    chain: Monte-Carlo chain of the flexible links model fitted to empirical data
Returns a probabilistic unipartite network with probabilities given by the product of relative abundances
"""
function neutral_model(n::Vector, chain)
    rel_n = n ./ sum(n) # vector of relative abundances
    na_mat = rel_n * rel_n' # neutral abundance matrix (proportional to adjacency matrix)
    
    # predict the number of links using the flexible links model (median predicted value)
    S = length(n) # species richness
    L = convert(Int64, round(median(predict_links(S, 1000, chain)))) # predicted number of links
    
    # draw L interactions using the neutral abundance matrix as weights
    # repeat 100 times and keep the L higher values
    N = zeros(Float64, S, S)

    for i in 1:100
        N_samp = sample(eachindex(na_mat), Weights(vec(na_mat)), L, replace=false)
        N[N_samp] .= N[N_samp] .+ 1.0
    end

    thres = sort(vec(N), rev=true)[L] # Lth maximum value
    N[N .< thres] .= 0
    N[N .>= thres] .= 1
    
    # convert to unipartite network 
    N = convert(Matrix{Bool}, N)
    N = UnipartiteNetwork(N)
    return N
end

    
