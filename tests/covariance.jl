using Compat.Test
using CombinationTool
import CombinationTool.calculate_covariancematrix

@testset "CovarianceMatrix" begin

    @testset "CovarianceMatrix: 2 Measurements, 2 Uncertainties" begin
        observables = [
            Observable("Obs1", params -> params[1].^4-params[2].^2)
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                             "syst1" => 5.0)),

            Measurement("Meas2", "Obs1", 39.0, Uncertainties("stat" => 2.5,
                                                             "syst1" => 4.0))
        ]

        correlations = [
            Correlation("stat", [1.0 0.0;
                                 0.0 1.0]),

            Correlation("syst1", [1.0 0.5;
                                  0.5 1.0])
        ]

        m = createmodel(observables, measurements, correlations)

        @test calculate_covariancematrix(m.uncertainties, m.correlationmatrices)==[29.0 10.0; 10.0 22.25]
    end;


    @testset "CovarianceMatrix: 3 Measurements, 3 Uncertainties" begin
        observables = [
             Observable("Obs1", params -> params[1].^4-params[2].^2)
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                             "syst1" => 5.0,
                                                             "syst2" => 3.0)),

            Measurement("Meas2", "Obs1", 39.0, Uncertainties("stat" => 2.5,
                                                             "syst1" => 4.0,
                                                             "syst2" => 1.2)),

            Measurement("Meas3", "Obs1", 24.0, Uncertainties("stat" => 1.2,
                                                             "syst1" => 2.3,
                                                             "syst2" => 0.75)),
        ]

        correlations = [
            Correlation("stat", [1.0 0.0 0.0;
                                 0.0 1.0 0.0;
                                 0.0 0.0 1.0]),

            Correlation("syst1", [1.0 0.5 0.0;
                                  0.5 1.0 0.0;
                                  0.0 0.0 1.0]),

             Correlation("syst2", [1.0 -0.3  0.2;
                                  -0.3  1.0 -0.5;
                                   0.2 -0.5  1.0])
        ]

        m = createmodel(observables, measurements, correlations)
        @test calculate_covariancematrix(m.uncertainties, m.correlationmatrices)≈[38.0 8.92 0.45; 8.92 23.69 -0.45; 0.45 -0.45 7.2925]

    end

    @testset "CovarianceMatrix: 3 Measurements (1 deactivated), 3 Uncertainties" begin
        observables = [
             Observable("Obs1", params -> params[1].^4-params[2].^2)
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                             "syst1" => 5.0,
                                                             "syst2" => 3.0)),

            Measurement("Meas2", "Obs1", 39.0, Uncertainties("stat" => 2.5,
                                                             "syst1" => 4.0,
                                                             "syst2" => 1.2), false),

            Measurement("Meas3", "Obs1", 24.0, Uncertainties("stat" => 1.2,
                                                             "syst1" => 2.3,
                                                             "syst2" => 0.75)),
        ]

        correlations = [
            Correlation("stat", [1.0 0.0 0.0
                                 0.0 1.0 0.0
                                 0.0 0.0 1.0]),

            Correlation("syst1", [1.0 0.5 0.0;
                                  0.5 1.0 0.0;
                                  0.0 0.0 1.0]),

             Correlation("syst2", [1.0 -0.3  0.2;
                                  -0.3  1.0 -0.5;
                                   0.2 -0.5  1.0])
        ]

        m = createmodel(observables, measurements, correlations)
        @test calculate_covariancematrix(m.uncertainties, m.correlationmatrices)≈[38.0 0.45; 0.45 7.2925]

    end

end;# CovarianceMatrix
