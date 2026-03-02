using DGGRIDRunner
using Test

@testset "DggridRunner.jl" begin
    @testset "DGGRIDParams" begin
        include("test_dggrid_params.jl")
    end
    @testset "Point Output" begin
        include("test_point_output.jl")
    end
    @testset "DGGRIDRunner Lib" begin
    #   # NOTE: This test is might need to be disabled due to a dynamic library loading issue with DGGRID7_jll.
    #   # It might be possible to run locally using the test/run_with_env.sh script.
        include("test_dggrid_runner_lib.jl")
    end
end
