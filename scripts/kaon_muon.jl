# =============================================================================================
# scripts/kaon_muon.jl - kaon→μ background for p+Mo via superposition
#
# Reproduces the reference ../SHiP_Pythia/analysis_kaon_muon.jl using the
# FixedTargetBackgrounds package.
# Run: julia --project=. scripts/kaon_muon.jl
# =============================================================================================
using FixedTargetBackgrounds
using PYTHIA8

const n_events = 50_000         # min-bias -> cheap per event; match the reference stats

function run_kaon(idB::Int)
    p = PYTHIA8.Pythia()
    configure_beams!(p; idA = 2212, idB = idB)
    configure_kaon!(p)
    PYTHIA8.init(p)
    _, r = measure_kaon_muons(p, ship_frame; n_events = n_events)
    return r
end

pp = run_kaon(2212)     # p+p
pn = run_kaon(2112)     # p+n

# Superpose per-collision yields to Mo-96 (42 p : 54 n). For kaons, p+p and p+n
# genuinely differ (isospin: the neutron's valence content shifts K±/K⁰ rates),
# so the 42:54 weighting actually matters here, unlike charm
mu_per_evt     = superpose(pp.mu_per_event, pn.mu_per_event, mo96)
acc_per_evt    = superpose(pp.n_acc / pp.n_gen, pn.n_acc / pn.n_gen, mo96)
plus_per_evt   = superpose(pp.n_muplus / pp.n_gen, pn.n_muplus / pn.n_gen, mo96)
minus_per_evt  = superpose(pp.n_muminus / pp.n_gen, pn.n_muminus / pn.n_gen, mo96)

pct(x) = round(100x; digits = 1)

println("\n--- Kaon→μ p+Mo (42:54 p+p/p+n superposition), N=$n_events each ---")
println("μ/collision (per inelastic collision) : ", round(mu_per_evt; digits = 3))
println("  in acceptance (η_lab 1-5)           : ", round(acc_per_evt; digits = 3),
        "  (", pct(acc_per_evt / mu_per_evt), "%)")
println("  μ⁺/μ⁻ ratio                         : ", round(plus_per_evt / minus_per_evt; digits = 2))
println("\nReference (analysis_kaon_muon.jl): ~0.59 μ/coll, ~0.49 in acc (~82%), μ⁺/μ⁻ ≈ 1.2")
