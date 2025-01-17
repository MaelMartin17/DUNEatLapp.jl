
"""
function greet_DUNEatLapp()
Function to test if the package works
There is no entry and the output is Hello DUNEatLapp!
"""
function greet_DUNEatLapp()
    return "Hello DUNEatLapp!"
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
        end
    end
    return df_Info
end

"""
function Condition_Cluster_Max(df_Ula::DataFrame,radius::Float64,limit::Float64,option::Bool)
function to select the clusters with the highest-energy and with no event between them and the limit variable.
It accepts a DataFrame for df_Ula, a Float for radius and limit (in centimeters) and a Boolen for option. 
If option is True, it returns a DataFrame with the number of the event and the energy of the cluster. If option is false, it returns the number of events rejected and the number of initial events.
"""
function Condition_Cluster_Max(df_Ula::DataFrame,radius::Float64,limit::Float64,option::Bool)
    #New_df_Ula = df_Ula[df_Ula[:,5] .> Emin, :]
    df_Info = DataFrame(evt = Int32[], E_max = Float32[])
    Index_evts = get_evts_index(df_Ula)
    nbr_evt_rejected = 0
    for i in 1:1:length(Index_evts[:,1])
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar = df_Ula[first:last,:]
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
                push!(data,[Ep x_moy/length(a.core_indices) y_moy/length(a.core_indices) z_moy/length(a.core_indices)])
            end
            data = vcat(data...)
            data_bis = data[data[:,1] .!= maximum(data[:,1]) , :]
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
                push!(data,[Ep x_moy/length(a.core_indices) y_moy/length(a.core_indices) z_moy/length(a.core_indices)])
            end
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

"""
function energy_all_clusters(df::DataFrame,radius::Float64)
function to get the energy of all clusters of the ULALAP events.
It accepts a DataFrame for df_Ula, a Float for radius. 
It returns a DataFrame with the event index and the energy of the clusters.
"""
function energy_all_clusters(df::DataFrame,radius::Float64)
    df_Info = DataFrame(evt = Int32[], E_cl = Float32[])
    Index_evts = get_evts_index(df)
    for i in 1:1:length(Index_evts[:,1])
        first = Index_evts[i,2]
        last  = Index_evts[i,3]
        data_Ar = df[first:last,:]
        if length(data_Ar[:,2]) > 3
            clustering = dbscan(Matrix(permutedims(data_Ar[:,2:4])), radius, min_neighbors = 1, min_cluster_size = 1)
            for a in clustering.clusters
                Ep = 0.
                for index_c in a.core_indices
                    Ep += data_Ar[index_c,:E]
                end
                push!(df_Info,[data_Ar[1,:evt],Ep])
            end
        end
    end
    return df_Info
end
