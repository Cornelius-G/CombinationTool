import Compat.Test

push!(LOAD_PATH, "1_BAT2.0")
using CombinationTool


Test.@testset "Package CombinationTool" begin
    include("covariance.jl")
    include("combination.jl")
end
