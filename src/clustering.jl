"""
function get_clusters_energy_of_evt(data_Ar::DataFrame,radius::Float64)
function to make the information of energy of each cluster of one event
It receives the dataframe from G4 of a single event and a radius in cm and returns the energy of each cluster
"""
function get_clusters_energy_of_evt(data_Ar::DataFrame,radius::Float64)
    clusters_info_E = Float64[]
    #if Edep is too small only one cluster
    if length(data_Ar[:,2]) <= 3
        Edep = sum(data_Ar[:,:E])
        Ecluster = Edep
        push!(clusters_info_E,Ecluster)
    else
        #apply clustering algorithm
        clustering = dbscan(Matrix(permutedims(data_Ar[:,2:4])), radius, min_neighbors = 1, min_cluster_size = 1)
        for a in clustering.clusters
            Ecluster = 0.
            for index_c in a.core_indices
                Ecluster +=data_Ar[index_c,:E]
            end
            push!(clusters_info_E,Ecluster)
        end
    end
    return clusters_info_E
end

"""
function process_neutron_file(df::DataFrame,radius::Float64,n_Ar_info::Vector{Int64},t_n_Ar_info::Vector{<:Real})
function to process a neutron file. It needs the ulalap file, the radius to use in the cluster, and the info about the capture
it will return the energy of clusters as well as the maximum energy of each cluster
"""
function process_clustering_neutron_file(df::DataFrame,radius::Float64,n_Ar_info::Vector{Int64},t_n_Ar_info::Vector{<:Real})
    clusters_info_E = Vector{<:AbstractFloat}[]
    clusters_info_max_E = AbstractFloat[]
    clusters_info_E_pre_post = Vector{<:AbstractFloat}[]
    clusters_info_max_E_pre_post = AbstractFloat[]
    clusters_info_E_other = Vector{<:AbstractFloat}[]
    clusters_info_max_E_other = AbstractFloat[]
    Index_evts = get_evts_index(df)
    #loop over all the events in the file
    for i in 1:1:length(Index_evts[:,1])
        i_evt = Index_evts[i,1] + 1 #since the index in geant4 starts at 0
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar_all = df[first:last,:]
        #When the capture takes place in Argon
        if n_Ar_info[i_evt] == 1
            t_capture_Ar = t_n_Ar_info[i_evt]
            #n-Ar capture gammas
            data_Ar = data_Ar_all[ (t_capture_Ar - 100) .< data_Ar_all[:,:t] .< (t_capture_Ar + 100),:]
            cl_all = get_clusters_energy_of_evt(data_Ar,radius)
            push!(clusters_info_E,cl_all)
            push!(clusters_info_max_E,maximum(cl_all))
            #had ellastic excitation gammas before n capture
            data_pre = data_Ar_all[data_Ar_all[:,:t] .< (t_capture_Ar - 100),:]
            cl_all = get_clusters_energy_of_evt(data_pre,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
            #n-Ar radioactivation gammas
            data_post = data_Ar_all[data_Ar_all[:,:t] .> (t_capture_Ar + 100),:]
            cl_all = get_clusters_energy_of_evt(data_post,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
        #When the capture take place elsewhere but with deposited energy in Argon
        else
            cl_all = get_clusters_energy_of_evt(data_Ar_all,radius)
            push!(clusters_info_E_other,cl_all)
            push!(clusters_info_max_E_other,maximum(cl_all))
        end
    end
    return vcat(clusters_info_E...), clusters_info_max_E, vcat(clusters_info_E_pre_post...), clusters_info_max_E_pre_post, vcat(clusters_info_E_other...), clusters_info_max_E_other
end
