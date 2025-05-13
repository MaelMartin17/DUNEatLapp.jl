
"""
function greet_DUNEatLapp()
Function to test if the package works
There is no entry and the output is Hello DUNEatLapp!
"""
function greet_DUNEatLapp()
    return "Hello DUNEatLapp!"
end

"""
	moving_window_filter(w_in::Vector, window::Int)
Function to apply a moving average window filter
It accepts a waveform and the size of the window to perform the average
It returns a new waveform
"""
function moving_window_filter(w_in::Vector{<:Real}, window::Int)
    w_out = zeros(length(w_in))
    half_window = window ÷ 2
    for i = 1 : 1 : half_window
        w_out[i] = sum(w_in[i:i+half_window]) / (half_window + 1)
    end
    for i = half_window + 1 : 1 : length(w_in)-half_window
        w_out[i] = sum(w_in[i-half_window:i+half_window]) / (window + 1)
    end
    for i = length(w_in)- half_window + 1 : 1 : length(w_in)
        w_out[i] = sum(w_in[i-half_window:i]) / (half_window + 1)
    end
    return w_out
end

"""
function apply_E_resolution(True_E_data::Vector, E_resolution::Int)
Function to apply an energy resolution on data. This resolution comes from the MicroBoone experiment and it is by default at 10% for 1 MeV.
To adjust the resolution, you have to take in account the 10 % and add or subtract to obtain the the desired value.
It accepts a vector for the True_E_data and an integer for the E_resolution. It returns a vector. 
"""
function apply_E_resolution(True_E_data::Vector,E_resolution::Int)
    distorted_E = zeros(length(True_E_data))
    for event = 1 : 1 : length(True_E_data)
        μ = True_E_data[event]
        σ = μ*sqrt((3.1/μ)*(3.1/μ) + (6.4/sqrt(μ))*(6.4/sqrt(μ)) + 7.3*7.3 )/100 +  E_resolution*μ/100
        d = Normal(μ,σ)        
        distorted_E[event] = rand(d)
    end
    return distorted_E
end

"""
function apply_std_E_resolution(True_E_data::Vector,E_resolution::Real)
function to apply a sigma/E = resolution / sqrt(E) energy resolution
It assumes that the True_E_data is in MeV, otherwise the results will be incorrect
"""
function apply_std_E_resolution(True_E_data::Vector,E_resolution::Real)
    distorted_E = zeros(length(True_E_data))
    for event = 1 : 1 : length(True_E_data)
        μ = True_E_data[event]
        σ = E_resolution * sqrt(μ)
        d = Normal(μ,σ)        
        distorted_E[event] = rand(d)
    end
    return distorted_E
end

"""
	 get_sampling(h::Histogram,n_samples::Int,res::Float64)
Funtion to sample from a histogram and apply an energy resolution
It accepts an histogram, the number of required samples n_samples and the resolution you want to apply res
"""
function get_sampling(h::Histogram,n_samples::Int,res::Float64)
    E_distribution = zeros(n_samples)
    items = collect(get_bin_centers(h))
    weights = h.weights .+ 1e-9
    for i = 1 : 1 : n_samples
        iSample = sample(items, Weights(weights), 1)[1]
        E_distribution[i] =  iSample 
    end
    return apply_std_E_resolution(E_distribution,res)
end

"""
    get_sampling(bin_centers::Vector{Float64},bin_weights::Vector{Float64},n_samples::Int,res::Float64)
Funtion to sample from a distribution and apply an energy resolution
It accepts a distribution x,y, where x are the values and y are the weights, the number of required samples n_samples and the resolution you want to apply res
"""
function get_sampling(bin_centers::Vector{Float64},bin_weights::Vector{Float64},n_samples::Int,res::Float64)
    E_distribution = zeros(n_samples)
    items = bin_centers
    weights = bin_weights .+ 1e-9
    for i = 1 : 1 : n_samples
        sampled_bin = sample(items, Weights(weights), 1)[1]
        E_distribution[i] =  sampled_bin 
    end
    return apply_std_E_resolution(E_distribution,res)
end

"""
 	get_bin_centers(h::Histogram)
Function to get the bin centers of a given histogram
It accepts an histogram and return the bin centers
"""
function get_bin_centers(h::Histogram)
    bin_edges = h.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2
    return bin_centers
end
