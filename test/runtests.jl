using BVPInterface
using Test, LinearAlgebra

@testset "BVPInterface.jl" begin
    include("bvpm2_test.jl")
    include("bvpsol_test.jl")
    include("coldae_test.jl")
    include("colmod_test.jl")
    include("colnew_test.jl")
end