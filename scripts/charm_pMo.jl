# =============================================================================================
# scripts/charm_pMo.jl - charm-μ background for p+Mo via superposition.
#
# Reproduces the D-meson level numbers of the reference
# ../SHiP_Pythia/analysis_charm_pMo.jl using the FixedTargetBackgrounds package.
# Run: julia --project=. scripts/charm_pMo.jl
# =============================================================================================
using FixedTargetBackgrounds
using PYTHIA8

const n_events = 20_000

function run_charm(idB::Int)
    p = PYTHIA8.Pythia()
    configure_beams!(p; idA = 2212, idB = idB)
    configure_charm!(p)
    PYTHIA8.init(p)
    _, r = measure_charm(p, ship_frame; n_events = n_events)
    return r
end

pp = run_charm(2212)    # p+p
pn = run_charm(2112)    # p+n

# Per-event yields, superposed to Mo-96 (42 p : 54 n)
# Ratios (acceptance, species) are combined by superposing numerator and denominator separately,
# never average the ratios directly
D_per_evt     = superpose(pp.D_per_event, pn.D_per_event, mo96)
acc_per_evt   = superpose(pp.n_acc / pp.n_gen, pn.n_acc / pn.n_gen, mo96)
D0_per_evt   = superpose(pp.n_D0 / pp.n_gen, pn.n_D0 / pn.n_gen, mo96)
Dp_per_evt   = superpose(pp.n_Dp / pp.n_gen, pn.n_Dp / pn.n_gen, mo96)
Ds_per_evt   = superpose(pp.n_Ds / pp.n_gen, pn.n_Ds / pn.n_gen, mo96)

pct(x) = round(100x; digits = 1)

println("\n--- Charm p+Mo (42:54 p+p/p+n superposition), N=$n_events each ---")
println("D/event (per charm collision) : ", round(D_per_evt; digits = 3))
println("  in acceptance (η_lab 1-5)   : ", pct(acc_per_evt / D_per_evt), "%")
println("  species D⁰/D⁺/Ds            : ",
        pct(D0_per_evt / D_per_evt), "/", pct(Dp_per_evt / D_per_evt), "/",
        pct(Ds_per_evt / D_per_evt), " %")
println("\nReference (analysis_charm_pMo.jl): 1.78 D/evt, ~80.7% D acc, ~59/31/9%")
