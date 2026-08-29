# scripts/spectra_plots.jl - Plots the spectra via the Makie extension.
# Run: julia --project=scripts scripts/spectra_plots.jl
using FixedTargetBackgrounds
using PYTHIA8
using CairoMakie        # Loading a backend activates the plotting extension

const n = 20_000

function run_charm(idB)
    p = PYTHIA8.Pythia(); configure_beams!(p; idA = 2212, idB = idB); configure_charm!(p)
    PYTHIA8.init(p);
    return measure_charm(p, ship_frame; n_events = n)[2]
end

function run_kaon(idB)
    p = PYTHIA8.Pythia(); configure_beams!(p; idA = 2212, idB = idB); configure_kaon!(p)
    PYTHIA8.init(p);
    return measure_kaon_muons(p, ship_frame; n_events = n)[2]
end

cpp = run_charm(2212)
kpp = run_kaon(2212)

save("charm_eta.png", plot_eta(cpp; labels = ["charm p+p"]))
save("charm_pT.png", plot_pT(cpp; labels = ["charm p+p"], xmax = 5))
save("kaon_eta.png", plot_eta(kpp; labels = ["K→μ p+p"]))
save("kaon_pT.png", plot_pT(kpp; labels = ["K→μ p+p"], xmax = 2))

println("Saved charm_eta.png, charm_pT.png, kaon_eta.png, kaon_pT.png" )
