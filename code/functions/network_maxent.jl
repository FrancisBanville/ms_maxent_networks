"""
svd_entropy(N::UnipartiteNetwork)
    N: Unipartite directed simple network
Returns the (truncated) SVD-entropy of the adjacency matrix of a deterministic network 
"""
function svd_entropy(N::UnipartiteNetwork)
    A = Array(N.edges)
    F = svd(A) # execute the singular value decomposition 
    R = rank(A) # find the rank of the matrix
    F.S[(R+1):end] .= 0.0 # truncate the SVD at the rank of the matrix
    λ = (F.S ./ sum(F.S))[1:R] # singular values
    return -sum(λ .* log.(λ)) * 1 / log(R) # SVD-entropy
end

"""
swap_degreeN(N::UnipartiteNetwork)
    N: Unipartite directed simple network
Returns a network with two swapped interactions (or the same network if couldn't find two interactions to be swapped)
"""
function swap_degreeN(N::UnipartiteNetwork)
    A = Array(N.edges)
    iy = findall(A .== 1) # find all interactions in the network
    i1, i2 = StatsBase.sample(iy, 2, replace=false) # choose two randomly
    n1, n2 = CartesianIndex(i1[1], i2[2]), CartesianIndex(i2[1], i1[2]) # the predator of the first eat the prey of the second interaction (and reversely)

    count = 0 
    # propose other interactions to swap if at least one is already in the network 
    while (n2 ∈ iy)|(n1 ∈ iy)
        i1, i2 = StatsBase.sample(iy, 2, replace=false)
        n1, n2 = CartesianIndex(i1[1], i2[2]), CartesianIndex(i2[1], i1[2])
        count = count + 1
        if count >= 1000 # stop searching if couldn't find two new interactions after many tries
            break
            return UnipartiteNetwork(A)
        end
    end
    # swap interactions
    A[i1] = false
    A[i2] = false
    A[n1] = true
    A[n2] = true
    return UnipartiteNetwork(A)
end

"""
network_maxent(N::UnipartiteNetwork, nsteps::Int64)
    N: Unipartite directed simple network
    nsteps: Number of steps
Returns the adjacency matrix of maximum SVD-entropy constrained by the joint degree sequence of N using a simulating annealing algorithm 
"""
function network_maxent(N::UnipartiteNetwork, nsteps::Int64)
    A0 = copy(N)

    # initial vector for SVD-entropies
    entropies = zeros(Float64, nsteps)
    best = svd_entropy(A0)
    entropies[1] = best

    # initial temperature of the simulating annealing algorithm
    T = 0.2

    # simulating annealing algorithm
    for i in 2:nsteps
        # propose a new constrained random matrix and compute the difference in SVD-entropy 
        A1 = swap_degreeN(A0)
        candidate = svd_entropy(A1)
        delta = candidate - best
        # accept if the difference is positive or with a probability p if it's negative
        if delta > 0 
            A0 = A1 
            best = candidate
        else 
            p = exp(delta/T)
            P = rand(Uniform(0, 1))
            if P < p
                A0 = A1
                best = candidate
            end
        end
        entropies[i] = best
        # update the temperature
        T = T * 0.99
    end
    Amax = (A = A0, entropy = maximum(entropies), entropies = entropies)
    return Amax
end
