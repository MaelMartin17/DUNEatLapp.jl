using DataFrames, Random, Distributions

function greet_DUNEatLapp()
    return "Hello DUNEatLapp!"
end

function get_evts_index(df::DataFrame)
    event_number = df[1,:evt]
    evt_counter = 1
    evts_info = []
    start_evt_index = 1
    for i = 1 : 1 : length(df[!,:evt])
        if  event_number == df[i,:evt]
            continue
        else
            end_evt_index = i - 1
            push!(evts_info,[event_number start_evt_index end_evt_index])
            event_number = df[i,:evt]
            evt_counter += 1
            start_evt_index = i
        end
    end
    #evt_counter += 1
    start_evt_index = end_evt_index + 1
    end_evt_index = length(df[!,:evt])
    push!(evts_info,[event_number start_evt_index end_evt_index])
    return vcat(evts_info...)
end

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
