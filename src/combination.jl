function calculate_covariancematrix(
    uncertainties::AbstractArray{Real, 2},
    correlationmatrices::Vector{Array{Real,2}}
)

    nmeas= size(uncertainties, 1)
    nunc = size(uncertainties, 2)
    covariance = zeros(nmeas, nmeas)

    for i in 1:nunc
        for j in 1:nmeas
            for k in j:nmeas
                if(j == k)
                    covariance[j, k] += uncertainties[j, i]^2
                else
                    covariance[j, k] += correlationmatrices[i][j, k] * uncertainties[j, i] * uncertainties[k, i]
                    covariance[k, j]  = covariance[j, k]
                end
            end
        end
    end

    return covariance
end



function combinemeasurements(
    m::Model,
    invcov::Matrix{<:Real},
    sqrtdetcov::Real,
    parameters
)

    result::Float64=0.0

    nmeas=length(m.measurement_values)

    for i in 1:nmeas
        r1 = m.measurement_values[i] - m.observable_functions[m.measurement_observables[i]].f.obj.x(parameters)[1]
        for j in (i+1):nmeas
            r2 = m.measurement_values[j] - m.observable_functions[m.measurement_observables[j]].f.obj.x(parameters)[1]
            result += r1 * invcov[i,j] * r2
        end
        result += 0.5 * r1 * invcov[i,i] * r1
    end
    
    final_result::Float64 = -result - log((2*π)^(0.5*nmeas) * sqrtdetcov) 
    return  final_result
end
