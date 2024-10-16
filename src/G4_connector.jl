"""
function get_evts_index(df::DataFrame)
Function to get the index fo start and end of each event in a ulalap.csv file
It accepts a dataframe and returns a matrix with column corresponding to the: number_of_event  index_start_evt  index_ends_evt
"""
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

"""
function get_n_primaries(my_file::String)
Function to get the number of primaries generated
"""
function get_n_primaries(my_file::String)
    df_primary = CSV.read(my_file, DataFrame,comment="#",header=["evt","pdg","E","x","y","z"])
    n_primaries_in_file = df_primary[end,1]+1
    return n_primaries_in_file
end
