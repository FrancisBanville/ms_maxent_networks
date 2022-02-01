"""
simulate_degrees(JDD::Matrix{Float64})
    JDD: joint degree distribution 
Returns a simulated vector of in and out degrees using the joint degree distribution as weight
"""
function simulate_degrees(JDD::Matrix{Float64})
      S = size(JDD, 1) - 1 # number of species 
      deg = findall(JDD .>= 0) # get cartesian indices 
      deg_samp = sample(deg, Weights(vec(JDD)), S, replace=true) # select species degrees randomly 
      # get in and out degrees
      kin = zeros(Int64, S)
      kout = zeros(Int64, S)
      for i in 1:S
            kout[i] = deg_samp[i][1] - 1 # we substract one because of degree 0
            kin[i] = deg_samp[i][2] - 1
      end
      return (kin = kin, kout = kout)
end