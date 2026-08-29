module FixedTargetBackgrounds

using PYTHIA8

include("kinematics.jl")
include("pythia_setup.jl")
include("measurements.jl")
include("superposition.jl")

export m_proton, ship_eta_window, ship_boost, ship_frame
export Boost, GenerationFrame, CMFrame, LabFrame
export lab_eta, in_window

export nucleus_pdg_id, register_nucleus!, register_mo96!
export configure_beams!, configure_angantyr!, configure_charm!, configure_kaon!
export default_ecm

export measure_charm, measure_kaon_muons, measure_multiplicity
export d_meson_ids, muon_id

export NuclearComposition, w_pp, w_pn, superpose, mo96

# Plotting entry points. Methods live in ext/FixedTargetBackgroundsMakieExt.jl, which loads automatically
# once a Makie backend (CairoMakie/GLMakie) is imported
"""
    plot_eta(results...; labels, edges, per_event=true, window=ship_eta_window) -> Figure

Overlay lab-frame pseudorapidity spectra (dN/dη per event) of one or more `measure_*`result NamedTuples.
Requires a Makie backend to be loaded.
"""
function plot_eta end

"""
    plot_pT(results...; labels, edges) -> Figure

Overlay transverse-momentum spectra (unit-area density, log-y). Requires a Makie backend.
"""
function plot_pT end

"""
    plot_multiplicity(results...; labels, edges) -> Figure

Overlay charged-multiplicity distributions from `measure_multiplicity`. Requires a backend.
"""
function plot_multiplicity end

export plot_eta, plot_pT, plot_multiplicity

end
