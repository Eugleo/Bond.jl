using Test

include("../src/Mass.jl")

@testset "Mass tests" begin
    @testset "Atomic and molecular weights" begin
        @test Mass.mass("H") ≈ 1.00782503207
        @test Mass.mass("H2") ≈ 2.01565006414
        @test Mass.mass("H2O") ≈ 18.0105646837
        @test Mass.mass("NH3") ≈ 17.02654910101
        @test Mass.mass("C6H12O6") ≈ 180.0633881022
    end

    @testset "Protein mass" begin
        @test Mass.mass(Mass.Types.Protein("EUGLEO")) ≈ 834.3026405049
        @test Mass.mass(Mass.Types.Protein("HITCHHIKERSGUIDE")) ≈ 1924.81168589048
    end
end