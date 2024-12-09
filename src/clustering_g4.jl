"""
function get_clusters_energy_of_evt(data_Ar::DataFrame,radius::Float64)
function to make the information of energy of each cluster of one event
It receives the dataframe from G4 of a single event and a radius in cm and returns the energy of each cluster
"""
function get_clusters_energy_of_evt(data_Ar::DataFrame, radius::Float64)
    # Collect energy contributions for each cluster
    cluster_energies = Float64[]

    # If there are few data points, consider all as a single cluster
    if size(data_Ar, 1) <= 3
        total_energy = sum(data_Ar.E)
        push!(cluster_energies, total_energy)
    else
        # Apply DBSCAN clustering on the spatial coordinates (columns 2 to 4)
        spatial_coords = Matrix(permutedims(df_evt[:, [:x, :y, :z]]))
        clustering = dbscan(spatial_coords, radius, min_neighbors=1, min_cluster_size=1)

        # Loop over clusters
        for iCluster in clustering.clusters
            #Loop over all hits in cluster to compute total cluster energy
            jIndex = iCluster.core_indices
            Ecluster = sum(data_Ar[jIndex,:E])
            push!(cluster_energies,Ecluster)
        end
    end

    return cluster_energies
end


"""
function get_clusters_vertex_and_energy_of_evt(data_Ar::DataFrame, radius::Float64)
Function to get the energy and average position of each cluster
It receives the event data frame and the radius for the clsutering
"""
function get_clusters_vertex_and_energy_of_evt(data_Ar::DataFrame, radius::Float64)
    # Collect energy and coordinate information for each cluster
    cluster_energies = AbstractFloat[]
    cluster_centers = Vector{<:AbstractFloat}[]  # To store the average coordinates (x, y, z) of each cluster

    # If there are few data points, consider all as a single cluster
    if size(data_Ar, 1) <= 3
        total_energy = sum(data_Ar.E)
        avg_coords = [mean(data_Ar.x), mean(data_Ar.y), mean(data_Ar.z)]
        push!(cluster_energies, total_energy)
        push!(cluster_centers, avg_coords)
    else
        # Apply DBSCAN clustering on the spatial coordinates (:x, :y, :z)
        spatial_coords = Matrix(permutedims(data_Ar[:, [:x, :y, :z]]))
        clustering = dbscan(spatial_coords, radius, min_neighbors=1, min_cluster_size=1)

        # Loop over clusters
        for iCluster in clustering.clusters
            #Loop over all hits in cluster to compute total cluster energy
            jIndex = iCluster.core_indices
            cluster_data = data_Ar[jIndex, :]
            Ecluster = sum(cluster_data.E)
            avg_coords = [mean(cluster_data.x), mean(cluster_data.y), mean(cluster_data.z)]
            push!(cluster_energies,Ecluster)
            push!(cluster_centers,avg_coords)
        end
    end

    return cluster_energies, cluster_centers
end


"""
function process_clustering_neutron_file(df::DataFrame,radius::Float64,n_Ar_info::Vector{Int64},t_n_Ar_info::Vector{<:Real})
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

"""
function full_process_clustering_neutron_file(df::DataFrame,radius::Float64,n_Ar_info::Vector{Int64},t_n_Ar_info::Vector{<:Real})
Function to process a neutron run and return cluster information about energy and position
"""
function full_process_clustering_neutron_file(df::DataFrame,radius::Float64,n_Ar_info::Vector{Int64},t_n_Ar_info::Vector{<:Real})
    clusters_info_E = Vector{<:AbstractFloat}[]
    clusters_info_max_E = AbstractFloat[]
    clusters_info_E_bar = []
    clusters_info_E_pre_post = Vector{<:AbstractFloat}[]
    clusters_info_E_pre_post_bar = []
    clusters_info_max_E_pre_post = AbstractFloat[]
    clusters_info_E_other = Vector{<:AbstractFloat}[]
    clusters_info_E_other_bar = []
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
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_Ar,radius)
            push!(clusters_info_E,cl_all)
            push!(clusters_info_max_E,maximum(cl_all))
            push!(clusters_info_E_bar,pos_r)
            #had ellastic excitation gammas before n capture
            data_pre = data_Ar_all[data_Ar_all[:,:t] .< (t_capture_Ar - 100),:]
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_pre,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
            push!(clusters_info_E_pre_post_bar,pos_r)
            #n-Ar radioactivation gammas
            data_post = data_Ar_all[data_Ar_all[:,:t] .> (t_capture_Ar + 100),:]
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_post,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
            push!(clusters_info_E_pre_post_bar,pos_r)
        #When the capture take place elsewhere but with deposited energy in Argon
        else
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_Ar_all,radius)
            push!(clusters_info_E_other,cl_all)
            push!(clusters_info_max_E_other,maximum(cl_all))
            push!(clusters_info_E_other_bar,pos_r)
        end
    end
    return vcat(clusters_info_E...), clusters_info_max_E, clusters_info_E_bar, vcat(clusters_info_E_pre_post...), clusters_info_max_E_pre_post, clusters_info_E_pre_post_bar, vcat(clusters_info_E_other...), clusters_info_max_E_other, clusters_info_E_other_bar
end
