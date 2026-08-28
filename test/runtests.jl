using FixedTargetBackgrounds
using Test
using PYTHIA8

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

    @testset "nucleus_pdg_id" begin
        @test nucleus_pdg_id(42, 96)  == 1000420960     # Mo-96 (the SHiP target)
        @test nucleus_pdg_id(82, 208) == 1000822080     # Pb-208
        @test nucleus_pdg_id(74, 184) == 1000741840     # W-184
    end

    if get(ENV, "FTB_RUN_PYTHIA", "false") == "true"
        @testset "Pythia - charm p+p (200 events)" begin
            pythia = PYTHIA8.Pythia()
            configure_beams!(pythia; idA = 2212, idB = 2212)
            configure_charm!(pythia)
            @test PYTHIA8.init(pythia)
            cfg, res = measure_charm(pythia, ship_frame; n_events = 200)
            @test res.n_gen > 0
            @test res.n_D > 0
            @test 1.0 < res.D_per_event < 3.0   # canonical ~1.78
        end
    end

end
