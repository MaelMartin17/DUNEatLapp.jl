using DataFrames, Random, Distributions, CSV, Clustering


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


function get_rate_neutron_captures_Ar(my_file::String,name_primary::String,fidu::Real=100.)
    df_neutrons = CSV.read(my_file, DataFrame,comment="#",header=["evt","proc","Z","A","pdg","E","x","y","z","t"])
    n_neutrons = get_n_primaries(name_primary)
    n_capture_Ar = 0
    n_capture_Ar_fidu = 0
    FD_x_size = 3100.
    FD_y_size = 755.
    FD_z_size = 700.

    for i = 1 : 1 : length(df_neutrons[!,1])
         if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18
            n_capture_Ar += 1
        end
        if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18 && abs(df_neutrons[i,:x]) < (FD_x_size - fidu) && abs(df_neutrons[i,:y]) < (FD_y_size - fidu) && abs(df_neutrons[i,:z]) < (FD_z_size - fidu)
            n_capture_Ar_fidu += 1
        end
    end
    return n_capture_Ar/n_neutrons, n_capture_Ar_fidu/n_neutrons
end


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
