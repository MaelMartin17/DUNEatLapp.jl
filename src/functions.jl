#_______________________________________________________________________________________________________________________
"""
function greet_DUNEatLapp()
Function to test if the package works
There is no entry and the output is Hello DUNEatLapp!
"""
function greet_DUNEatLapp()
    return "Hello DUNEatLapp!"
end
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
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
#_______________________________________________________________________________________________________________________
"""
function Compute_Sun_initialVector(Event::Matrix)
function that takes for input a output.txt file from Marley and it computes the vector of the initial neutrino.
It returns a list of vectors.
"""
function Compute_Sun_initialVector(Event::Matrix)
    Liste_iniV = []
    for i in 1:length(Event[:,1])
        if Event[i,1] == 2
            Vector_nu = [Event[i+1,3], Event[i+1,4], Event[i+1,5]]
            push!(Liste_iniV, Vector_nu)
        end
    end
    return  Liste_iniV
end
#_______________________________________________________________________________________________________________________
"""
    rescale_neutron(h_neutron_side, side; n_runs::Int64=30,evts_per_run::Int64=500000)
Function to rescale the neutron histogram based on the specified side and simulation parameters.

# Arguments
- `h_neutron_side`: type::Histogram of the energy spectrum of the neutrons coming from the specified side.
- `side`: string indicating the side of the detector ("Bottom", "LongWall", or "ShortWall").
- `n_runs`: number of simulation runs.
- `evts_per_run`: number of events per run.

# Returns
- `counts`: Rescaled counts for the neutron histogram.
"""
function rescale_neutron(h_neutron_side, side; n_runs::Int64=30,evts_per_run::Int64=500000)

    ##############################
    # DUNE parameters for neutron flux at the FD Cavern
    # change this to adapt the scenario
    # We expect to have a value between 1 (best scenario) and 9 (worst scenario)
    neutron_flux = 2.2e-6#n/cm2/s 
    #############################
    # Surface
    DUNE_area_long = 12e2*64e2
    DUNE_area_short = 12e2*12e2
    # Time 
    number_of_years = 10*u"yr"
    days_per_year = 365
    hours_per_day = 24
    minutes_per_hour = 60 
    seconds_per_minute = 60
    seconds_in_period = uconvert(u"s",number_of_years)

    total_rate_Bottom = 2 * neutron_flux * DUNE_area_long * ustrip(seconds_in_period)
    total_rate_long_wall = 2 * neutron_flux * DUNE_area_long * ustrip(seconds_in_period)
    total_rate_short_wall = 2 * neutron_flux * DUNE_area_short * ustrip(seconds_in_period)

    if side == "Bottom"
        #Bottom
        rate_per_second_Bottom   = total_rate_Bottom/uconvert(u"s",10*u"yr")
        time_sim_neutrons_Bottom = n_runs*evts_per_run/rate_per_second_Bottom
        scaling_factor_Bottom    = uconvert(u"s",10*u"yr")/time_sim_neutrons_Bottom
        counts = rand.(Poisson.(h_neutron_side.weights .* scaling_factor_Bottom));
    elseif side == "LongWall"
        #Long Wall
        rate_per_second_longW   = total_rate_long_wall/uconvert(u"s",10*u"yr")
        time_sim_neutrons_longW = n_runs*evts_per_run/rate_per_second_longW
        scaling_factor_longW    = uconvert(u"s",10*u"yr")/time_sim_neutrons_longW
        counts = rand.(Poisson.(h_neutron_side.weights .* scaling_factor_longW));
    elseif side == "ShortWall"
        #Short Wall
        rate_per_second_shortW   = total_rate_short_wall/uconvert(u"s",10*u"yr")
        time_sim_neutrons_shortW = n_runs*evts_per_run/rate_per_second_shortW
        scaling_factor_shortW    = uconvert(u"s",10*u"yr")/time_sim_neutrons_shortW
        counts = rand.(Poisson.(h_neutron_side.weights .* scaling_factor_shortW));
    else
        error("Invalid side. Choose from 'Bottom', 'LongWall', or 'ShortWall'.")
    end
    return counts
end
#_______________________________________________________________________________________________________________________
