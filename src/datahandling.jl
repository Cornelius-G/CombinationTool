export Model
export Parameter
export Measurement
export Observable
export Uncertainties
export Correlation
export createmodel

using FunctionWrappers
import FunctionWrappers: FunctionWrapper

const Uncertainties=Dict

struct RealRealFunc
    f::FunctionWrapper{Real,Tuple{Real}}
end
(cb::RealRealFunc)(v) = cb.f(v)


struct Model
    observable_names::Vector{String}
    observable_functions::Vector{RealRealFunc}

    measurement_names::Vector{String}
    measurement_observables::Vector{Int}
    measurement_values::Vector{Real}
    active_measurements::Vector{Int} # indices of the active measurements

    uncertainty_names::Vector{String}
    uncertainties::Array{Real, 2}

    correlationmatrices::Vector{Array{Real,2}}
end


struct Observable
    name::String
    func::Function
    # TODO min::Real
    # TODO max::Real
end


struct Measurement
    name::String
    obs::String
    value::Float64
    uncertainties::Dict{String, Float64}# TODO: Float -> Real?
    activity::Bool
end


# constructor with default value active=true
function Measurement(
    name::String,
    obs::String,
    value::Float64,
    uncertainties::Dict{String, Float64};
    active=true
)
    Measurement(name, obs, value, uncertainties, active)
end


struct Correlation
    name::String
    correlationmatrix::Array{Float64}
    active::Bool
end

# constructor with default value active=true
function Correlation(
    name::String,
    correlationmatrix::Array{Float64};
    active=true
)
    Correlation(name, correlationmatrix, active)
end



function createmodel(
    observables::AbstractArray{Observable},
    measurements::AbstractArray{Measurement},
    correlations::AbstractArray{Correlation}
)

    observable_names, observable_functions = createobservables(observables)

    measurement_names, measurement_observables, measurement_values, active_measurements = createmeasurements(measurements, observable_names)

    uncertainty_names, uncertainties, correlationmatrices = createuncertainties(measurements, correlations, active_measurements)

    Model(observable_names,
            observable_functions,
            measurement_names,
            measurement_observables,
            measurement_values,
            active_measurements,
            uncertainty_names,
            uncertainties,
            correlationmatrices)
end



function createobservables(observables::AbstractArray{Observable})
    nobs = length(observables)

    names = Vector{String}(undef, nobs)
    funcs = Vector{RealRealFunc}(undef, nobs)

    for i in 1:nobs
        names[i] = observables[i].name
        funcs[i] = RealRealFunc(observables[i].func)
    end

    duplicates = findfirstduplicate(names)
    if(duplicates[1])
        throw(ArgumentError("Observable with the name \"" * duplicates[2] * "\" already exists."))
    end

    names, funcs
end


function getobservableindex(name::String, observable_names::Vector{String})
    index = findfirst(x -> x==name, observable_names)
    if index == 0
        throw(ArgumentError("Observable with the name \"" * name * "\" does not exists."))
    end
    index
end


function createmeasurements(
    measurements::AbstractArray{Measurement},
    observable_names::Vector{String}
)

    nmeas = length(measurements)

    nactives=0
    for i in 1:nmeas
        if(measurements[i].activity)
            nactives += 1
        end
    end

    names = Vector{String}(undef,  nactives)
    measuredobs = Vector{Int}(undef, nactives)
    values = Vector{Real}(undef, nactives)
    actives = Vector{Int}(undef, nactives)

    a=0
    for i in 1:nmeas
        if(measurements[i].activity)
            a += 1
            actives[a] = i
            names[a] = measurements[i].name #TODO: check for name duplicates (only when same observable)
            values[a] = measurements[i].value
            measuredobs[a] = getobservableindex(measurements[i].obs, observable_names)
        end
    end

    names, measuredobs, values, actives
end


function createuncertainties(
    measurements::AbstractArray{Measurement},
    correlations::AbstractArray{Correlation},
    active_meas::AbstractVector{Int}
)
    n_active_meas = length(active_meas)
    active_unc = Vector{Correlation}()

    for i in 1:length(correlations)
        if(correlations[i].active)
            push!(active_unc, correlations[i]) #TODO: is push! ok?
        end
    end

    n_active_unc = length(active_unc)

    uncertainties = Array{Real}(undef, n_active_meas, n_active_unc)
    uncertainty_names = Array{String}(undef, n_active_unc)
    corr_array =  Array{Array{Real}}(undef, n_active_unc)

    # fill array of uncertainties
    for i in 1:n_active_meas
        for j in 1:n_active_unc
            uncertainties[i,j] = measurements[active_meas[i]].uncertainties[active_unc[j].name]
        end
    end

    # fill uncertainty names and correlationmatrices
    for i in 1:n_active_unc
        uncertainty_names[i] = active_unc[i].name
        corr_array[i] = active_unc[i].correlationmatrix[active_meas, active_meas]
    end

    uncertainty_names, uncertainties, corr_array
end



# TODO: remove
function printmodel(m)
    println("\nObservables:")
    println(m.observable_names)
    println(m.observable_functions)

    println("\nMeasurements:")
    println(m.measurement_names)
    println(m.measurement_observables)
    println(m.measurement_values)

    println("\nUncertainties:")
    println(m.uncertainty_names)
    println(m.uncertainties)
    println(m.correlationmatrices)

end
#=======================================================================#
