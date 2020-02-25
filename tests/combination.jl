using Compat.Test
using CombinationTool
using LinearAlgebra
import CombinationTool.calculate_covariancematrix
import CombinationTool.combinemeasurements


@testset "CombineMeasurements" begin

    @testset "Combination: 2 Obs, 2 Meas" begin
        observables = [
            Observable("Obs1", params -> params[1].^4-params[2].^2),
            Observable("Obs2", params -> 2*params[1].^2 + 4*params[2])
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                             "syst1" => 5.0)),

            Measurement("Meas2", "Obs2", 39.0, Uncertainties("stat" => 2.5,
                                                             "syst1" => 4.0))
        ]

        correlations = [
            Correlation("stat", [1.0 0.0;
                                 0.0 1.0]),

            Correlation("syst1", [1.0 0.5;
                                  0.5 1.0])
        ]

        m = createmodel(observables, measurements, correlations)

        covariance=calculate_covariancematrix(m.uncertainties, m.correlationmatrices)
        sqrtdetcov=sqrt(abs(det(covariance)))

        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [3.694, 1.15]), digits=8) == -373.40688529
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [1, -2]), digits=8) == -71.61944062
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [-3.2, 5.02]), digits=8) == -13.80088051
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [-6, 0.01]), digits=8) == -30800.19990017
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [28.3, -21.2]), digits=8) ≈ -8.36405819e9
    end;


    @testset "Combination: 2 Obs, 3 Meas" begin
        observables = [
            Observable("Obs1", params -> params[1].^4-params[2].^2),
            Observable("Obs2", params -> 2*params[1].^2 + 4*params[2])
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                                     "syst1" => 5.0,
                                                                     "syst2" => 3.0)),

            Measurement("Meas2", "Obs2", 39.0, Uncertainties("stat" => 2.5,
                                                                     "syst1" => 4.0,
                                                                     "syst2" => 1.2)),

            Measurement("Meas3", "Obs2", 24.0, Uncertainties("stat" => 1.2,
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

        covariance=calculate_covariancematrix(m.uncertainties, m.correlationmatrices)

        sqrtdetcov=sqrt(abs(det(covariance)))

        # comparison with EFTfitter v1.03
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [3.694, 1.15]), digits=4) == -264.3656
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [2.611,1.952]), digits=4) == -7.2546
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [3.277, 1.68]), digits=4) == -61.3434
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [3.139, 1.008]), digits=4) == -39.4537

    end;

    @testset "Combination: 2 Obs, 3 Meas (1 deactivated)" begin
        observables = [
            Observable("Obs1", params -> params[1].^4-params[2].^2),
            Observable("Obs2", params -> 2*params[1].^2 + 4*params[2])
        ]

        measurements = [
            Measurement("Meas1", "Obs1", 53.0, Uncertainties("stat" => 2.0,
                                                             "syst1" => 5.0)),

            Measurement("Meas2", "Obs2", 39.0, Uncertainties("stat" => 2.5,
                                                             "syst1" => 4.0)),

            Measurement("Meas3", "Obs2", 139.0, Uncertainties("stat" => 9.5,
                                                             "syst1" => 6.0), active=false)
        ]

        correlations = [
            Correlation("stat", [1.0 0.0 0.5;
                                 0.0 1.0 0.0;
                                 0.5 0.0 1.0]),

            Correlation("syst1", [1.0 0.5 0.0;
                                  0.5 1.0 0.3;
                                  0.0 0.3 1.0])
        ]

        m = createmodel(observables, measurements, correlations)

        covariance=calculate_covariancematrix(m.uncertainties, m.correlationmatrices)
        sqrtdetcov=sqrt(abs(det(covariance)))

        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [3.694, 1.15]), digits=8) == -373.40688529
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [1, -2]), digits=8) == -71.61944062
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [-3.2, 5.02]), digits=8) == -13.80088051
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [-6, 0.01]), digits=8) == -30800.19990017
        @test round(combinemeasurements(m, inv(covariance), sqrtdetcov, [28.3, -21.2]), digits=8) ≈ -8.36405819e9
    end;

end;# CombineMEasurements
