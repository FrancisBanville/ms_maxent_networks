"""
neutral_model(A::Vector{Int64}, chain)
    n: Vector of abundances
    chain: Monte-Carlo chain of the flexible links model fitted to empirical data
Returns a probabilistic unipartite network with probabilities given by the product of relative abundances
"""
function neutral_model(n::Vector, chain)
    rel_n = n ./ sum(n) # vector of relative abundances
    N = rel_n * rel_n' # neutral abundance matrix (proportional to adjacency matrix)
    
    # predict the connectance using the flexible links model (median predicted value)
    S = length(n) # species richness
    L = convert(Int64, round(median(predict_links(S, 1000, chain)))) # predicted number of links
    C = L / S^2 # connectance
    
    # threshold the neutral abundance matrix using the predicted connectance
    N[diagind(N)] .= 0 # set diagonal elements to zero
    q = quantile(vec(N), 1-C) # find threshold
    N[N .< q] .= 0 
    N[N .>= q] .= 1
    
    # convert to unipartite network 
    N = convert(Matrix{Bool}, N)
    N = UnipartiteNetwork(N)
    return N
end


