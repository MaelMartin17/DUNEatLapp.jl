using DUNEatLapp
using Test
using DataFrames, Random, Distributions, CSV, Clustering, ProgressBars

@testset "DUNEatLapp.jl" begin
	
	df = DataFrame(evt = 1:10, b = rand(10))

	@test DUNEatLapp.greet_DUNEatLapp() == "Hello DUNEatLapp!"
	@test DUNEatLapp.greet_DUNEatLapp() != "Hello world!"
	@test DUNEatLapp.get_evts_index(df) != 2
end
