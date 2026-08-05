# =============================================================================================
# kinematics.jl - fixed-target boost + generation frame
#
# Two orthogonal axes:
#       * Boost:                - the physical setup (beam energy, beam/target mass) → β, γ,
#                                 rapidity Y_boost, √s. General across any foxed-target energy.
#       * GenerationFrame       - How the event was generated:
#                                 CMFrame → generated in the CM frame, must be boosted to lab
#                                 (the SHiP workflow: Beams:eCM, then boost).
#                                 LabFrame → generated directly in the lab/fixed-target frame
#                                 (Pythia frameType 2/3). No boost applied.
# Downstream code calls lab_eta(frame, ...) polymorphically, so switching an experiment's
# generation frame needs no change in the mearument loops.
# =============================================================================================

const m_proton = 0.938                  # GeV
const ship_eta_window = (1.0, 5.0)      # forward window proxy for the SHiP aperture

"""
    Boost(e_beam; m_beam=m_proton, m_target=m_proton)

Fixed-target kinematics for a beam of energy `e_beam` (GeV) on a target at rest.
Stores β, γ, the boost rapidity `y`, and the CM energy `ecm`(√s).
"""
struct Boost
    e_beam::Float64
    m_beam::Float64
    m_target::Float64
    β::Float64
    γ::Float64
    y::Float64
    ecm::Float64
end

function Boost(e_beam::Real; m_beam::Real = m_proton, m_target::Real = m_proton)
    ecm = sqrt(m_beam^2 + m_target^2 + 2 * m_target * e_beam)
    γ = (e_beam + m_target) / ecm
    β = e_beam / (e_beam + m_target)
    y = atanh(β)
    return Boost(e_beam, m_beam, m_target, β, γ, y, ecm)
end

"""
    Generation Frame

Abstract type for the frame in which events are generated. Concrete subtypes: [`CMFrame`](@ref)
(boost to lab) and [`LabFrame`](@ref) (identity).
"""
abstract type GenerationFrame
end

"""
    CMFrame(boost::Boost)
    CMFrame(e_beam; kwargs...)

Events generated in the CM frame (`Beams:eCM`); `lab_eta` boosts them to the lab. This is the SHiP
workflow. Carries the [`Boost`](@ref) so √s and Y_boost are on hand.
"""
struct CMFrame <: GenerationFrame
    boost::Boost
end
CMFrame(e_beam::Real; kwargs...) = CMFrame(Boost(e_beam;kwargs...))

"""
    LabFrame/boost::Boost)
    LabFrame(e_beam; kwargs...)

Events generated directly in the lab/fixed-target frame; `lab_eta` is the identity (no boost).
Carries the [`Boost`](@ref) too, so √s / Y_boost remain available for reporting even though they are
not apply to η.
"""
struct LabFrame <: GenerationFrame
    boost::Boost
end
LabFrame(e_beam::Real; kwargs...) = LabFrame(Boost(e_beam; kwargs...))

"""
    lab_eta(frame, pz, e, pT)

Pseudorapidity in the lab frame for a particle with longitudinal momentum `pz`, energy `e`, and transverse
momentum `pT`. For [`CMFrame`](@ref) the momentum is boosted to the lab first; for [`LabFrame`](@ref) it
is used as is.
"""
function lab_eta(f::CMFrame, pz::Real, e::Real, pT::Real)
    pzl = f.boost.γ * (pz + f.boost.β * e)
    return atanh(pzl / sqrt(pT^2 + pzl^2))
end
lab_eta(::LabFrame, pz::Real, e::Real, pT::Real) = atanh(pz / sqrt(pT^2 + pz^2))

"""
    in window(η, window=ship_eta_window)

Whether pseudorapidity `η` falls in the (inclusive) acceptance `window`. The SHiP default is a
forward-window proxy, not the true trakcer apperture.
"""
in_window(η::Real, window::Tuple{<:Real,<:Real} = ship_eta_window) = window[1] <= η <= window[2]

# --- SHiP canonical defaults ---
const ship_boost = Boost(400.0)
const ship_frame = CMFrame(ship_boost)
