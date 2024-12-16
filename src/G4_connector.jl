"""
function get_evts_index(df::DataFrame)
Function to get the index of start and end of each event in a ulalap.csv file
It accepts a dataframe and returns a matrix with column corresponding to the: number_of_event  index_start_evt  index_ends_evt
"""
function get_evts_index(df::DataFrame)
    # Initialize variables
    evts_info = Vector{Vector{Int}}()
    start_idx = 1
    event_number = df[1, :evt]
    evt_counter = 1

    # Iterate through rows to identify event boundaries
    for i in 2:size(df, 1)
        if df[i, :evt] != event_number
            # Record the event's index range
            push!(evts_info, [evt_counter, start_idx, i - 1])
            # Update for the new event
            evt_counter += 1
            start_idx = i
            event_number = df[i, :evt]
        end
    end

    # Handle the last event
    push!(evts_info, [evt_counter, start_idx, size(df, 1)])

    return vcat(evts_info...)
end

"""
function get_n_primaries(my_file::String)
Function to get the number of primaries generated
It accepts the name as string of the primary file produced by ulalap and returns the number of primaries that have been used for the simulations
"""
function get_n_primaries(my_file::String)
    df_primary = CSV.read(my_file, DataFrame,comment="#",select=[1],header=["evt","pdg","E","x","y","z"])
    n_primaries_in_file = df_primary[end,:evt]+1
    return n_primaries_in_file
end

"""
function get_primary_vertex(my_file::String)
function to get the primary info from a csv file into a data frame
It accepts the name of the file and returns a data frame
"""
function get_primary_vertex(my_file::String)
    df_primary = CSV.read(my_file, DataFrame,comment="#",header=["evt","pdg","E","x","y","z"])
    return df_primary
end
