# =============================================================================================
# scripts/superposition_vs_angantyr.jl
#
# Reproduces the p+p/p+n superposition vs TRUE p+Mo (Angantyr). Charm should be nulcear
# invariant (~1.0x); soft kaon yield should rise ~1.8x (wounded nucleons); (N_ch) 11.6 → 21.4
# Run: julia --project=. scripts/superposition_vs_angantyr.jl
# =============================================================================================
using FixedTargetBackgrounds
using PYTHIA8

const n_sup = 20_000    # superposition (p+p, p+n): cheap
const n_ang = 5_000     # Angantyr p+Mo: much slower, fewer events
const mo_id = nucleus_pdg_id(42, 96)

# --- Builders (angantyr=true → register Mo-96 + nuclear σ-fit knobs) ---
function build(idB; angantyr, process!::Function)
    p = PYTHIA8.Pythia()
    angantyr && register_mo96!(p)
    configure_beams!(p; idA = 2212, idB = idB)
    process!(p)
    angantyr && configure_angantyr!(p)
    PYTHIA8.init(p)
    return p
end

run_charm(idB; angantyr, n) = measure_charm(build(idB; angantyr, process! = configure_charm!), ship_frame; n_events = n)[2]
run_kaon(idB; angantyr, n) = measure_kaon_muons(build(idB; angantyr, process! = configure_kaon!), ship_frame; n_events = n)[2]
run_mult(idB; angantyr, n) = measure_multiplicity(build(idB; angantyr, process! = configure_kaon!); n_events = n)[2]

# --- Charm ---
cpp = run_charm(2212; angantyr = false, n = n_sup)
cpn = run_charm(2112; angantyr = false, n = n_sup)
cmo = run_charm(mo_id; angantyr = true, n = n_ang)
charm_sup = superpose(cpp.D_per_event, cpn.D_per_event, mo96)

# --- Kaon ---
kpp = run_kaon(2212; angantyr = false, n = n_sup)
kpn = run_kaon(2112; angantyr = false, n = n_sup)
kmo = run_kaon(mo_id; angantyr = true, n = n_ang)
kaon_sup = superpose(kpp.mu_per_event, kpn.mu_per_event, mo96)
kaon_sup_acc = superpose(kpp.n_acc / kpp.n_gen, kpn.n_acc / kpn.n_gen, mo96)

# --- Multiplicity ---
mpp = run_mult(2212; angantyr = false, n = n_sup)
mmo = run_mult(mo_id; angantyr = true, n = n_ang)

rd(x) = round(x; digits = 3)
rr(x) = round(x; digits = 2)

println("\n===== Superposition (42:54 p+p/p+n) vs true Angantyr p+Mo =====\n")

println("CHARM (D/event)")
println("  superposition : ", rd(charm_sup))
println("  Angantyr p+Mo : ", rd(cmo.D_per_event))
println("  ratio ang/sup : ", rr(cmo.D_per_event / charm_sup), "   (expect ~1.0, hard process, nuclear invariant)\n")

println("KAON (μ/collision)")
println("  superposition : ", rd(kaon_sup),             "  in-acc  ", rd(kaon_sup_acc))
println("  Angantyr p+Mo : ", rd(kmo.mu_per_event),     "  in-acc  ", rd(kmo.n_acc / kmo.n_gen))
println("  ratio ang/sup : ", rr(kmo.mu_per_event / kaon_sup), "   (expect ~1.8, wounded nucleon enhancement)\n")

println("(N_ch)")
println("  p+p (single N-N) : ", rd(mpp.mean_nch))
println("  Angantyr p+Mo    : ", rd(mmo.mean_nch))
println("  ratio            : ", rr(mmo.mean_nch / mpp.mean_nch), "\n")

println("Reference: charm ~1.0x, kaon ~1.8x, (N_ch) 11.6→21.4")
println("Note: the A^(1/3)≈4.6 heuristic overestimates; measured enhancement is ~1.8x")
