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
        spatial_coords = Matrix(permutedims(data_Ar[:, [:x, :y, :z]]))
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
#_______________________________________________________________________________________________________________________

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
#_______________________________________________________________________________________________________________________

"""
	get_clusters_vertex_and_energy_of_evt(data_Ar::SubDataFrame{DataFrame, DataFrames.Index, UnitRange{Int64}}, radius::Float64)
Same as before to ensure compatability with @view usage
"""
function get_clusters_vertex_and_energy_of_evt(data_Ar::SubDataFrame{DataFrame, DataFrames.Index, UnitRange{Int64}}, radius::Float64)
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
#_______________________________________________________________________________________________________________________

"""
function cluster_energy_Max(df::DataFrame,radius::Float64)
function to get the cluster with the highest energy.
It accepts a DataFrame for df and a Float for radius (in centimeters). It returns a DataFrame with the number of the event and the energy of the cluster.
"""
function cluster_energy_Max(df::DataFrame,radius::Float64)
    df_Info = DataFrame(evt = Int32[], E_max = Float32[])
    Index_evts = get_evts_index(df)
    for i in 1:1:length(Index_evts[:,1])
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar = df[first:last,:]
        if length(data_Ar[:,2]) > 3
            clustering = dbscan(Matrix(permutedims(data_Ar[:,2:4])), radius, min_neighbors = 1, min_cluster_size = 1)
            E_c = 0
            for a in clustering.clusters
                Ep = 0.
                for index_c in a.core_indices
                    Ep +=data_Ar[index_c,:E]
                end
                if Ep > E_c                
                    E_c = Ep
                end        
            end
            push!(df_Info,[data_Ar[1,:evt], E_c])
        else
            push!(df_Info,[data_Ar[1,:evt], sum(data_Ar[:,:E])])
        end
    end
    return df_Info
end
#_______________________________________________________________________________________________________________________

"""
function Condition_Cluster_Max_Emin(df_Ula::DataFrame,radius::Float64,Emin::Float64,limit::Float64,option::Bool)
function to select the clusters with the highest-energy and with no event between them and the limit variable. All clusters under Emin are removed.
It accepts a DataFrame for df_Ula, a Float for radius and limit (in centimeters), a float for Emin (in keV) and a Boolen for option. 
If option is True, it returns a DataFrame with the number of the event and the energy of the cluster. If option is false, it returns the number of events rejected and the number of initial events.
"""
function Condition_Cluster_Max_Emin(df_Ula::DataFrame,radius::Float64,Emin::Float64,limit::Float64,option::Bool)
    df_Info = DataFrame(evt = Int32[], E_max = Float32[])
    Index_evts = get_evts_index(df_Ula)
    nbr_evt_rejected = 0
    for i in 1:1:length(Index_evts[:,1])
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar = df_Ula[first:last, :]
        if length(data_Ar[:,2]) > 3
            clustering = dbscan(Matrix(permutedims(data_Ar[:,2:4])), radius, min_neighbors = 1, min_cluster_size = 1)
            data = []
            for a in clustering.clusters
                Ep = 0.
                x_moy = 0
                y_moy = 0
                z_moy = 0
                for index_c in a.core_indices
                    Ep    += data_Ar[index_c,:E]
                    x_moy += data_Ar[index_c,:x]
                    y_moy += data_Ar[index_c,:y]
                    z_moy += data_Ar[index_c,:z]
                end
                if Ep > Emin
	            push!(data,[Ep x_moy/length(a.core_indices) y_moy/length(a.core_indices) z_moy/length(a.core_indices)])
	    	end
            end
            if length(data[:,1])>0
	        data = vcat(data...)
        	data_bis = data[(data[:,1] .!= maximum(data[:,1])) .& (data[:,1] .> Emin), :]
        	condition = true
                for i in 1:1:length(data_bis[:,1])
                    dist = sqrt((data[argmax(data[:,1]),2]-data_bis[i,2])^2+(data[argmax(data[:,1]),3]-data_bis[i,3])^2+(data[argmax(data[:,1]),4]-data_bis[i,4])^2)
                    if dist < limit
                        condition = false 
                        nbr_evt_rejected += 1
                        break
                    else 
                        continue
                    end
                end
            else 
                nbr_evt_rejected += 1
                condition = false
            end
            if option == true && condition == true 
                push!(df_Info,[data_Ar[1,:evt], maximum(data[:,1])])
            end
        end
    end
    if option == true 
        return df_Info
    else
        return nbr_evt_rejected, length(Index_evts[:,1])
    end
end

#_______________________________________________________________________________________________________________________

"""
function nbr_cluster(df::DataFrame,radius::Float64)
function to get the number of cluster of an ULALAP event.
It accepts a DataFrame for df_Ula, a Float for radius. 
It returns a DataFrame with the event index and the number of cluster.
"""
function nbr_cluster(df::DataFrame,radius::Float64)
    df_Info = DataFrame(evt = Int32[], nbr_cl = Int32[])
    Index_evts = get_evts_index(df)
    for i in 1:1:length(Index_evts[:,1])
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar = df[first:last,:]
        if length(data_Ar[:,2]) > 3
            clustering = dbscan(Matrix(permutedims(data_Ar[:,2:4])), radius, min_neighbors = 1, min_cluster_size = 1)
            push!(df_Info,[data_Ar[1,:evt], length(clustering.clusters)])
        end
    end
    return df_Info
end
#_______________________________________________________________________________________________________________________

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
        data_Ar_all = @view df[first:last,:] #view to avoid copy
        #When the capture takes place in Argon
        if n_Ar_info[i_evt] == 1
            t_capture_Ar = t_n_Ar_info[i_evt]
            #n-Ar capture gammas
	    # Define the time window
            lower_bound = t_capture_Ar - 100
            upper_bound = t_capture_Ar + 100
	    #Select time of neutron capture to avoid mixing with post capture activity from radioactivation
	    data_Ar = data_Ar_all[(data_Ar_all[:, :t] .> lower_bound) .& (data_Ar_all[:, :t] .< upper_bound), :]
            cl_all = get_clusters_energy_of_evt(data_Ar,radius)
            push!(clusters_info_E,cl_all)
            push!(clusters_info_max_E,maximum(cl_all))
            #had ellastic excitation gammas before n capture
            data_pre = data_Ar_all[data_Ar_all[:,:t] .< lower_bound,:]
            cl_all = get_clusters_energy_of_evt(data_pre,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
            #n-Ar radioactivation gammas
            data_post = data_Ar_all[data_Ar_all[:,:t] .> upper_bound,:]
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
#_______________________________________________________________________________________________________________________

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
        data_Ar_all = @view df[first:last,:] #view to avoid copying
        #When the capture takes place in Argon
        if n_Ar_info[i_evt] == 1
            t_capture_Ar = t_n_Ar_info[i_evt]
            #n-Ar capture gammas
	    # Define the time window
            lower_bound = t_capture_Ar - 100
            upper_bound = t_capture_Ar + 100
	    # Select rows within the time window of neutron capture
            data_Ar = data_Ar_all[(data_Ar_all[:, :t] .> lower_bound) .& (data_Ar_all[:, :t] .< upper_bound), :]
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_Ar,radius)
            push!(clusters_info_E,cl_all)
            push!(clusters_info_max_E,maximum(cl_all))
            push!(clusters_info_E_bar,pos_r)
            #had ellastic excitation gammas before n capture
            data_pre = data_Ar_all[data_Ar_all[:,:t] .< lower_bound,:]
            cl_all, pos_r = get_clusters_vertex_and_energy_of_evt(data_pre,radius)
            push!(clusters_info_E_pre_post,cl_all)
            push!(clusters_info_max_E_pre_post,maximum(cl_all))
            push!(clusters_info_E_pre_post_bar,pos_r)
            #n-Ar radioactivation gammas
            data_post = data_Ar_all[data_Ar_all[:,:t] .> upper_bound,:]
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
    return clusters_info_E, clusters_info_max_E, clusters_info_E_bar, clusters_info_E_pre_post, clusters_info_max_E_pre_post, clusters_info_E_pre_post_bar, clusters_info_E_other, clusters_info_max_E_other, clusters_info_E_other_bar
end
#_______________________________________________________________________________________________________________________
