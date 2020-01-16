export CombinationDensity


struct CombinationDensity <: AbstractDensity
    model::Model
    invcov::Matrix{Real}
    sqrtdetcov::Real
end

function CombinationDensity(m::Model)
    covariance = calculate_covariancematrix(m.uncertainties, m.correlationmatrices)
    sqrtdetcov = sqrt(abs(det(covariance)))

    CombinationDensity(m, inv(covariance), sqrtdetcov)
end


function BAT.density_logval(
    dens::CombinationDensity,
    params
)

    combinemeasurements(
        dens.model,
        dens.invcov,
        dens.sqrtdetcov,
        params
    ) 
end 