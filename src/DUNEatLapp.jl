__precompile__(true)

module DUNEatLapp

using DelimitedFiles
using DataFrames
using Random
using Distributions
using CSV
using Clustering

export greet_DUNEatLapp
export get_evts_index
export apply_E_resolution
export apply_std_E_resolution 
export get_n_primaries
export get_rate_neutron_captures_Ar
export get_info_neutron_captures_Ar
export get_primary_vertex
export cluster_energy_Max
export Condition_Cluster_Max
export get_clusters_energy_of_evt
export get_clusters_vertex_and_energy_of_evt
export process_clustering_neutron_file
export get_capture_position
export is_n_capture_on_Ar
export Single_Hits_lardon
export nbr_cluster
include("functions.jl")
include("clustering_g4.jl")
include("G4_connector.jl")
include("neutron_capture_info.jl")
include("analyse_lardon.jl")

end
