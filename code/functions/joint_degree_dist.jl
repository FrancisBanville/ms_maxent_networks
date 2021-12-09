include(joinpath("joint_degree_seq.jl")) # the joint degree distribution is retrieved from the joint degree sequence

"""
joint_degree_dist(N::UnipartiteNetwork) 
    N: Unipartite directed simple network
Returns the proportion of species having a given number of preys and predators
"""
function joint_degree_dist(N::UnipartiteNetwork)
    S = richness(N) # number of species in N
    k = joint_degree_seq(N) # joint degree sequence of N

    # count the number of species having kin predators and kout preys
    k_count = zeros(Int64, S+1, S+1)
    for kin in 0:S
        for kout in 0:S
            k_count[kin+1, kout+1] = sum((k.kin .== kin) .& (k.kout .== kout))
        end
    end
    
    k_prop = k_count ./ S # proportion of species having kin predators and kout preys
    return k_prop
end
