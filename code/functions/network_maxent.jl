"""
svd_entropy(A::Matrix{Bool}) 
    A: Adjacency matrix of a unipartite directed simple network
Returns the SVD-entropy of the adjacency matrix of a deterministic network using a Singular Value Decomposition
"""
function svd_entropy(A::Matrix{Bool})
    # compute the SVD-entropy 
    F = svd(A)
    Λ = F.S[1:rank(A)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

"""
network_maxent(N::UnipartiteNetwork, nsteps::Int64, T::Float64)
    N: Unipartite directed simple network
    nsteps: Number of steps
Returns the adjacency matrix of maximum SVD-entropy constrained by the joint degree sequence of N using a simulating annealing algorithm 
"""
function network_maxent(N::UnipartiteNetwork, nsteps::Int64)
    A = convert(Matrix, N.edges) # convert UnipartiteNetwork to dense matrix format
    rmg = matrixrandomizer(A) # matrix generator object

    # initial random matrix that has the same joint degree sequence as N
    A0 = convert(Matrix{Bool}, rand(rmg))
    
    # initial temperature of the simulating annealing algorithm
    T = exp(-10/nsteps)  

    # simulating annealing algorithm
    for i in 1:nsteps
        # propose a new constrained random matrix and compute the difference in SVD-entropy 
        A1 = convert(Matrix{Bool}, rand(rmg))
        delta_entropy = svd_entropy(A1) - svd_entropy(A0)
        # accept if the difference is positive or with a probability p if it's negative
        if delta_entropy > 0 
            A0 = A1 
        else 
            p = exp(delta_entropy/T)
            P = rand(Uniform(0, 1))
            if P < p
                A0 = A1
            end
        end

        # update the temperature
        T = T * exp(-10/nsteps)  
    end
    max_entropy = svd_entropy(A0)
    Amax = (A = A0, entropy = max_entropy)
    return Amax
end

