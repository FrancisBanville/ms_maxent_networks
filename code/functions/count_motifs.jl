"""
count_motifs(N::UnipartiteNetwork) 
    N: unipartite network
Returns the proportion that each motif was found in the unipartite network
"""
function count_motifs(N::UnipartiteNetwork) 
    # list and count motifs
    motifs = unipartitemotifs() 
    motifs_count = [length(find_motif(N, motifs[i])) for i in 1:13]
    # returns the proportion of each motif for all networks
    return motifs_count ./ sum(motifs_count)
end