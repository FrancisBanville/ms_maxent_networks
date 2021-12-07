"""
predict_links(S::Int64, n::Int64, chain) 
    S: number of species
    n: number of links that will be predicted
    chain: Monte-Carlo chain of the flexible links model fitted to empirical data
Returns n values of predicted numbers of links
"""
function predict_links(S::Int64, n::Int64, chain) 
    predicted_links = zeros(Int64, n) # create vector for the predicted links
    p = get_params(chain[200:end,:,:]) # get fitted parameters of the flexible links model 
    # simulate n number of links
    for i in 1:n
        # select parameter values randomly
        p_i = rand(1:length(p.μ)) 
        μ, ϕ = p.μ[p_i], p.ϕ[p_i]
        # predict number of links 
        predicted_links[i] = rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))) + (S-1)
    end
    return predicted_links
  end