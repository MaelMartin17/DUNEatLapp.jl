"""
function get_evts_index(df::DataFrame)
Function to get the index of start and end of each event in a ulalap.csv file
It accepts a dataframe and returns a matrix with column corresponding to the: number_of_event  index_start_evt  index_ends_evt
"""
function get_evts_index(df::DataFrame)
    #Initialize variables
    event_number = df[1,:evt]
    evts_info = Vector{Vector{Int}}()
    start_evt_index = 1
    # Iterate through rows to identify event boundaries
    for i = 1 : 1 : length(df[!,:evt])
        if  event_number == df[i,:evt]
            continue
        else
            end_evt_index = i - 1
    	    # Record the event's index range
            push!(evts_info,[event_number, start_evt_index, end_evt_index])
    	    # Update for the new event
            event_number = df[i,:evt]
            start_evt_index = i
        end
    end
    # Handle the last event
    start_evt_index = end_evt_index + 1
    end_evt_index = length(df[!,:evt])
    push!(evts_info,[event_number, start_evt_index, end_evt_index])
    return transpose(hcat(evts_info...))
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


"""
function get_hits_in_active_LAr_TPC(df_evt_all_hits::DataFrame,fidu::Float64 = 0.)
Function to filter the hits of a given event and select only hits in active region of the active volume
It uses the fact that the CRPs/cathode represents a surface of 6000x1300 cm2
For fiducialization purposes a given fiducialization can also be used (in cm)
"""
function get_hits_in_active_LAr_TPC(df_evt_all_hits::DataFrame, fidu::Float64 = 0.0)
    active_LAr_x = 3000.0
    active_LAr_y = 727.0
    active_LAr_z = 625.0

    active_x = active_LAr_x - fidu
    active_y = active_LAr_y - fidu
    active_z = active_LAr_z - fidu

    df_active = df_evt_all_hits[ (abs.(df_evt_all_hits.x) .< active_x) .&& 
                                 (abs.(df_evt_all_hits.y) .< active_y) .&& 
                                 (abs.(df_evt_all_hits.z) .< active_z), :]

    return df_active
end

"""
function get_hits_in_inactive_LAr(df_evt_all_hits::DataFrame)
Function to the hits in the inactive LAr between the field cage and the primary membrane of the cryostat
"""
function get_hits_in_inactive_LAr(df_evt_all_hits::DataFrame)
    active_LAr_x = 3000.0
    active_LAr_y = 727.0
    active_LAr_z = 625.0

    df_inactive = df_evt_all_hits[ (abs.(df_evt_all_hits.z) .> active_LAr_z) .|| 
                                   (abs.(df_evt_all_hits.x) .> active_LAr_x) .|| 
                                   (abs.(df_evt_all_hits.y) .> active_LAr_y), :]
    return df_inactive
end

"""
function get_hits_in_active_LAr_TPC_Lateral_fidu(df_evt_all_hits::DataFrame,fidu_x::Float64, fidu_z::Float64)
Function to filter the hits of a given event and select only hits in active region of the active volume
It uses the fact that the CRPs/cathode represents a surface of 6000x1300 cm2
For fiducialization purposes a given fiducialization can also be used (in cm)
"""
function get_hits_in_active_LAr_TPC_Lateral_fidu(df_evt_all_hits::DataFrame, fidu_x::Float64, fidu_z::Float64)
    active_LAr_x = 3000.0
    active_LAr_y = 727.0
    active_LAr_z = 625.0

    active_x = active_LAr_x - fidu_x
    active_y = active_LAr_y 
    active_z = active_LAr_z - fidu_z

    df_active = df_evt_all_hits[ (abs.(df_evt_all_hits.x) .< active_x) .&& 
                                 (abs.(df_evt_all_hits.y) .< active_y) .&& 
                                 (abs.(df_evt_all_hits.z) .< active_z), :]

    return df_active
end
