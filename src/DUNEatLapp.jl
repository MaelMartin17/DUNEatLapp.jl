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
export cluster_energy_Max
export Condition_Cluster_Max
export get_capture_position
include("functions.jl")
include("G4_connector.jl")
include("neutron_capture_info.jl")

end
