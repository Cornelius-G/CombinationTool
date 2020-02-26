import Compat.Test
using CombinationTool


Test.@testset "Package CombinationTool" begin
    include("covariance.jl")
    include("combination.jl")
end
