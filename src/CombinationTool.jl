__precompile__(true)

module CombinationTool

using BAT, StatsBase, IntervalSets
using LinearAlgebra


include("datahandling.jl")
include("combination.jl")
include("callBAT.jl")
include("utils.jl")

end # module
