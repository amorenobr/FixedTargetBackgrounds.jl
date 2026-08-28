module FixedTargetBackgrounds

using PYTHIA8

include("kinematics.jl")
include("pythia_setup.jl")
include("measurements.jl")

export m_proton, ship_eta_window, ship_boost, ship_frame
export Boost, GenerationFrame, CMFrame, LabFrame
export lab_eta, in_window

export nucleus_pdg_id, register_nucleus!, register_mo96!
export configure_beams!, configure_angantyr!, configure_charm!, configure_kaon!
export default_ecm

export measure_charm, measure_kaon_muons, measure_multiplicity
export d_meson_ids, muon_id

end
