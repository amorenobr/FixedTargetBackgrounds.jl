using FixedTargetBackgrounds
using Test

@testset "FixedTargetBackgrounds.jl" begin
    @testset "Boost - SHiP 400 GeV" begin
        b = Boost(400.0)
        @test b.β   ≈ 0.9977  atol = 1e-4
        @test b.γ   ≈ 14.62   atol = 1e-2
        @test b.y   ≈ 3.375   atol = 1e-2
        @test b.ecm ≈ 27.43   atol = 1e-2
        @test ship_boost.e_beam == 400.0
    end

    @testset "lab_eta - CMFrame boosts, LabFrame is identity" begin
        # A particle already at CM rapidity 0 (pz=0) lands near Y_boost in the lab.
        η_cm = lab_eta(ship_frame, 0.0, 1.0, 1.0)
        @test η_cm > 3.0

        # LabFrame applies no boost: same (pz,e,pT) gives the raw lab η.
        lab = LabFrame(400.0)
        pz, e, pT = 5.0, 6.0, 1.0
        @test lab_eta(lab, pz, e, pT) ≈ atanh(pz / sqrt(pT^2 + pz^2))
        # And it must differ from the boosted CM results for the same inputs.
        @test lab_eta(lab, pz, e, pT) != lab_eta(ship_frame, pz, e, pT)
    end

    @testset "in_window" begin
        @test in_window(3.0)
        @test in_window(1.0)                    # inclusive lower edge
        @test in_window(5.0)                    # inclusive upper edge
        @test !in_window(0.5)
        @test !in_window(5.5)
        @test in_window(2.0, (1.5, 2.5))        # custom window
    end

end
