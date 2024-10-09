using DUNEatLapp
using Test

@testset "DUNEatLapp.jl" begin
	@test DUNEatLapp.greet_DUNEatLapp() == "Hello DUNEatLapp!"
	@test DUNEatLapp.greet_DUNEatLapp() != "Hello world!"
end
