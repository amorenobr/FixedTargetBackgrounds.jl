module FixedTargetBackgrounds

using PYTHIA8

include("kinematics.jl")
include("pythia_setup.jl")

export m_proton, ship_eta_window, ship_boost, ship_frame
export Boost, GenerationFrame, CMFrame, LabFrame
export lab_eta, in_window

export nucleus_pdg_id, register_nucleus!, register_mo96!
export configure_beams!, configure_angantyr!, configure_Charm!, configure_kaon!
export default_ecm

end
