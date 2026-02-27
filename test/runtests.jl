using DggridRunner
using Test

@testset "DggridRunner.jl" begin
    @testset "DggridParams" begin
        include("test_dggrid_params.jl")
    end
    @testset "Point Output" begin
        include("test_point_output.jl")
    end
    # @testset "Runner Lib" begin
    #     # NOTE: This test is currently disabled due to a dynamic library loading issue with DGGRID7_jll.
    #     # It can be run locally using the test/run_with_env.sh script.
    #     # include("test_dggrid_runner_lib.jl")
    # end
end
