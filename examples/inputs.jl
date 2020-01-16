#============= Observables =============================================#
function obs2(params)
    2*params.p1.^2 + 4*params.p2
end

observables = [
    Observable("Obs1", params -> params.p1.^4-params.p2.^2),
    Observable("Obs2", obs2)
]


#============= Measurements ============================================#
measurements = [
    Measurement("Meas1","Obs1", 42.0, Uncertainties("stat"=>1.1,
                                                    "syst"=>1.2)),

    Measurement("Meas2","Obs2", 22.0, Uncertainties("stat"=>2.1,
                                                    "syst"=>2.2))
]


#============= Correlations ============================================#
correlations = [
    Correlation("stat", [1.0 0.0
                         0.0 1.0], active=false),

    Correlation("syst", [1.0 0.5
                         0.5 1.0])
]
#=======================================================================#
