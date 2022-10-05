"""
MxSim(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the average of all species’ largest similarity index 
"""
function MxSim(N::UnipartiteNetwork)
      AJS_N = AJS(N) # Additive Jaccard similarity between all species pairs
      length_AJS_N = length(AJS_N)
      if length_AJS_N == 0
            return missing
      else 
            max_AJS = []
            try
                  # find the maximum similarity for every species
                  for i in 1:richness(N) 
                    spi = []
                        for j in 1:length_AJS_N
                              if in(AJS_N[j][1])(species(N)[i])
                               push!(spi, AJS_N[j][2])
                              end
                        end
                        if length(spi) > 0
                         push!(max_AJS, maximum(spi)[1])
                        end
                  end
            catch
                  return missing
            end
            return mean(max_AJS) # return average similarity index 
      end
end


