"""
joint_degree_seq(N::UnipartiteNetwork) 
    N: Unipartite directed simple network
Returns the number of preys and predators for all species in N (i.e. the joint degree sequence)
"""
function joint_degree_seq(N::UnipartiteNetwork)
    # get the joint degree sequence 
    kin = vec(sum(N.edges, dims=1)) # number of predators for all species (column sums)
    kout = vec(sum(N.edges, dims=2)) # number of preys for all species (row sums)
    k = (kin = kin, kout = kout) # combine in and out-degree sequences
    return k
end


