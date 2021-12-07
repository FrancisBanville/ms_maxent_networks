# Define custom functions

"""
predict_links(S::T, n::T, chain) where {T <: Int64}
    S: number of species
    n: number of simulations
    chain: Monte-Carlo chain of the flexible links model fitted to empirical data
Returns n values of predicted numbers of links
"""
function predict_links(S::T, n::T, chain) where {T <: Int64}
    predicted_links = zeros(Int64, n)
    p = get_params(chain[200:end,:,:])
    for i in 1:n
        # select parameter values randomly
        p_i = rand(1:length(p.μ))
        μ, ϕ = p.μ[p_i], p.ϕ[p_i]
        # predict number of links 
        predicted_links[i] = rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))) + (S-1)
    end
    return predicted_links
  end


"""
p_k(S::T, L::T, k::T) where {T <: Int64}
    S: number of species
    L: number of links
    k: degree (number of species interacting with a given species)
Returns the probability that a species has a degree k (derived from maximum entropy - analytical solution)
"""
function p_k(S::T, L::T, k::T) where {T <: Int64}
    # mean degree constraint
    kavg = 2*L/S 
    # probability that a species has a degree of k
    c = 1/(kavg - 1)
    r = (kavg - 1)/kavg
    p_k = c*(r^k)
    return p_k
end


"""
dd_maxent_prob(S::T, L::T) where {T <: Int64}
    S: number of species
    L: number of links
Returns the degree distribution of maximum entropy (probability that a species has a degree k)
"""
function dd_maxent_prob(S::T, L::T) where {T <: Int64}
    dd_maxent_prob = [p_k(S, L, k) for k in 1:S]
    return dd_maxent_prob
end
  

"""
dd_maxent(S::T, L::T) where {T <: Int64}
    S: number of species
    L: number of links
Returns the expected degree for all species (ordered from lowest to highest)
"""
function dd_maxent(S::T, L::T) where {T <: Int64}
    # degree distribution of maximum entropy
    dd_maxent_proba = dd_maxent_prob(S, L)
    # cumulative number of species expected to have a degree of k or lower
    dd_maxent_cum = cumsum(dd_maxent_proba) .* S
    # number of species expected to have a degree of k
    dd_maxent = zeros(Int64, S)
    K = 1:S
    for k in K
      dd_maxent[k] = convert(Int64, round(dd_maxent_cum[k] - sum(dd_maxent)))
    end
    # expected degree for all species (ordered from lowest to highest)
    dd_maxent = vcat(fill.(K, dd_maxent)...)
    return dd_maxent
end


"""
p_kin_kout(S::T, L::T, kin::T, kout::T) where {T <: Int64}
    S: number of species
    L: number of links
    kin: in-degree (number of predators of a species)
    Kout: out-degree (number of preys of a species)
Returns the probability that a species has kin predators and kout preys (derived from maximum entropy - analytical solution)
"""
function p_kin_kout(S::T, L::T, kin::T, kout::T) where {T <: Int64}
    # a species cannot have a total degree k = 0 or higher than the number of species
    if (kin + kout == 0 || kin + kout > S)
        p_kin_kout = 0
    else
        # mean degree constraint
        kavg = 2L/S 
        # probability that a species has an in-degree of kin and an out-degree of kout
        k = kin + kout
        c = 1/(kavg - 1)
        r = (kavg - 1)/kavg
        p_kin_kout = binomial(BigInt(k), BigInt(kin))*c*(r/2)^k
    end
    return p_kin_kout
end
  

"""
jdd_maxent_prob(S::T, L::T) where {T <: Int64}
    S: number of species
    L: number of links
Returns the joint degree distribution of maximum entropy (probability that a species has an in-degree of kin and an out-degree of kout)
"""
function jdd_maxent_prob(S::T, L::T) where {T <: Int64}
    jdd_maxent_prob = zeros(Float64, S+1, S+1)
    # probability that a species has kin predators and kout preys for all values of kin and kout
    for kin in 0:S, kout in 0:S
      jdd_maxent_prob[kin+1, kout+1] = p_kin_kout(S, L, kin, kout)
    end
    return jdd_maxent_prob
end


"""
indd_maxent(S::T, L::T) where {T <: Int64}
    S: number of species
    L: number of links
Returns the expected in-degree for all species (ordered from lowest to highest)
"""
function indd_maxent(S::T, L::T) where {T <: Int64}
    # joint degree distribution of maximum entropy
    jdd_maxent_proba = jdd_maxent_prob(S, L)
    # cumulative number of species expected to have an in-degree kin or lower
    indd_maxent_cum = cumsum(sum.(eachcol(jdd_maxent_proba))) .* S
    # number of species expected to have an in-degree kin for all possible values of kin
    indd_maxent = zeros(Int64, S+1)
    Kin = 0:S
    for kin in Kin
      indd_maxent[kin+1] = convert(Int64, round(indd_maxent_cum[kin+1] - sum(indd_maxent)))
    end
    # expected in-degree for all species (ordered from lowest to highest)
    indd_maxent = vcat(fill.(Kin, indd_maxent)...)
    return indd_maxent
end


"""
matrix_rowcolsum(indd::T, outdd::T) where {T <: Vector{Int64}}
    indd: in-degree distribution 
    outdd: out-degree distributon 
Returns an adjacency matrix constrained by the in-degree and out-degree distributions
"""
function matrix_rowcolsum(indd::T, outdd::T) where {T <: Vector{Int64}}
    n = length(outdd)
    A = zeros(Bool, n, n)
    # row and column sums
    row_remain = copy(outdd)
    col_remain = copy(indd)
    # assign interactions iteratively until no more can be assigned
    while sum(row_remain) > 0
        # identify a species pair that could interact
        row_min = minimum(row_remain[row_remain .> 0])
        col_max = maximum(col_remain[col_remain .> 0])
        i = first(findall(row_remain .== row_min))
        j = first(findall(col_remain .== col_max))
        # assign the interaction
        A[i, j] = 1
        # adjust row and column sums
        row_remain[i] = row_remain[i] - 1
        col_remain[j] = col_remain[j] - 1
    end
    return A
end


"""
svd_entropy(A::T) where {T <: Matrix{Bool}}
    A: adjacency matrix
Returns the entropy of the adjacency matrix of a deterministic network using the singular values from a Singualr Value Decomposition
"""
function svd_entropy(A::T) where {T <: Matrix{Bool}}
    F = svd(A)
    Λ = F.S[1:rank(A)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end


"""
adjacency_matrix_maxent(A::T) where {T <: Matrix{Bool}}
    A: adjacency matrix
    n: number of permutations
Returns the adjacency matrix of maximum entropy constrained by row and column sums
"""
function adjacency_matrix_maxent(A::T, n::Int64) where {T <: Matrix{Bool}}
    # matrix generator object
    rmg = matrixrandomizer(A)
    # adjacency matrix of maximum entropy identified so far
    Amax = A
    svd_entropy_max = svd_entropy(A)
    for i in 1:n
        # generate a new matrix with the same row and column sums
        Anew = convert(Matrix{Bool}, rand(rmg))
        # SVD-entropy of that new matrix
        svd_entropy_new = svd_entropy(Anew)
        # keep that matrix if it has a higer entropy than all matrices checked so far
        if svd_entropy_new > svd_entropy_max
            Amax = Anew
            svd_entropy_max = svd_entropy_new
        end
    end
    return Amax
end


"""
predict_adjacency_matrix(S::T, L::T) where {T <: Int64}
    S: number of species
    L: number of links
    n: number of permutations
Returns the adjacency matrix of maximum entropy 
"""
function predict_adjacency_matrix(S::T, L::T, n::T) where {T <: Int64}
    # degree distribution of maximum entropy
    dd = dd_maxent(S, L)
    # in-degree distribution of maximum entropy
    indd = indd_maxent(S, L)
    # out-degree distribution
    outdd = dd .- indd
    # making sure the in-degree and out-degree distributions have the same sum
    m = length(dd)
    outdd[m] = last(outdd) + (sum(indd) - sum(outdd))
    # generate one matrix constrained by the in-degree and out-degree distributions
    A = matrix_rowcolsum(indd, outdd)
    # select the adjacency matrix with maximum SVD-entropy after n permutations
    A = adjacency_matrix_maxent(A, n)
    A = transpose(A)
    return A
end