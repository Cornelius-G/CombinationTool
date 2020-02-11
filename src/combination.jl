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

function eval_block_array(params,BA::Array{ObservableBlock})
    vcat(eval_block.(Ref(params),BA)...)
end

function eval_block(params,B::ObservableBlock)
    if(B.ObsandCoeffs.Coeff == [[] for i in 1:length(B.ObsandCoeffs.Coeff)])
        vcat(B.f.(Ref(params))...)
    else
        vcat(B.f.(Ref(params),B.ObsandCoeffs.Coeff)...)
    end
end


function combinemeasurements(
    m::Model,
    invcov::Matrix{<:Real},
    sqrtdetcov::Real,
    parameters
)

    result::Float64=0.0

    nmeas=length(m.measurement_values)
    
    r1 = eval_block_array(parameters,m.observable_blocks) - m.measurement_values
    
    final_result = -0.5 * transpose(r1)* invcov * r1
    
    return  final_result
end
