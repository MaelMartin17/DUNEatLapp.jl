using DataFrames, Random, Distributions, CSV

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

function get_n_primaries(my_file::String)
    df_primary = CSV.read(my_file, DataFrame,comment="#",header=["evt","pdg","E","x","y","z"])
    n_neutrons_in_file = df_primary[end,1]
    return n_neutrons_in_file
end

#=function get_rate_neutron_captures_Ar(my_file::String,name_primary::String)
    df_neutrons = CSV.read(my_file, DataFrame,comment="#",header=["evt","proc","Z","A","pdg","E","x","y","z","t"])
    n_neutrons = get_n_primaries(name_primary)
    n_capture_Ar = 0
    n_capture_Ar_fidu = 0

    n_evt = []
    n_evt_fidu = []
    fidu = 100
    for i = 1 : 1 : length(df_neutrons[:,1])
         if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18
            n_capture_Ar += 1
        end
        if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18 && abs(df_neutrons[i,:x]) < (3100 - fidu) && abs(df_neutrons[i,:y]) < (755 - fidu) && abs(df_neutrons[i,:z]) < (700 - fidu)
            n_capture_Ar_fidu += 1
        end
    end
    return n_capture_Ar/n_neutrons, n_capture_Ar_fidu/n_neutrons
end=#
