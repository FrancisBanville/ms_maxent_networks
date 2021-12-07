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