include(joinpath("joint_degree_dist_maxent.jl")) # the degree distribution is obtained from the joint degree distribution of maximum entropy

"""
degree_dist_maxent(S::Int64, L::Int64) 
    S: number of species
    L: number of links
Returns the degree distribution of maximum entropy given S and L
"""
function degree_dist_maxent(S::Int64, L::Int64) 
    joint_degree_dist = joint_degree_dist_maxent(S, L) # get the joint degree distribution of maximum entropy given S and L
    
    degree_dist = zeros(S+1) # create vector for the degree distribution (S+1 entries because degrees go from 0 to S)

    # compute probability for all possible degrees (degree k = in-degree + out-degree)
    # there is many combinations of in and out-degrees for a given degree 
    # these combinations are on the antidiagonals of k by k submatrices
    for k in 0:S
        degree_dist[k+1] = sum(joint_degree_dist[k-i+1, i+1] for i in 0:k)
    end
    return degree_dist                                 
end

