#using CombinationTool
using BAT
using LinearAlgebra
using ValueShapes
using IntervalSets
using Plots

push!(LOAD_PATH, "1_BAT2.0")
using CombinationTool

include("inputs.jl")

m = createmodel(observables, measurements, correlations)
#CombinationTool.printmodel(m)

likelihood = CombinationDensity(m)

prior = NamedTupleDist(
    p1 = [-10.0..10.0],
    p2 = [-20.0..20.0]
)

posterior = PosteriorDensity(likelihood, prior)

nchains = 4
nsamples = 10^4
algorithm = MetropolisHastings()

samples, stats = bat_sample(posterior, (nsamples, nchains), algorithm)
