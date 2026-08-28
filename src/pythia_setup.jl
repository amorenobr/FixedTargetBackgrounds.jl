# =============================================================================================
# pythia_setup.jl - Pythia8 beam and process configutaion
#
# Thin, composable wrappers over readString
#       ⋅ register_nucleus!(Z, A)       - the generalized Mo-96 fix (Mo, W, any nculeus)
#       ⋅ configure_beams!              - p+p, p+n, or p+A in the CM frame
#       ⋅ configure_charm!              - charm-enriched HardQCD, D mesons frozen
#       ⋅ configure_kaon!               - softQCD inelastic, K decays on, π/n frozen
#       ⋅ configure_angantyr!           - nuclear cross-section fit knobs for p+A
# Each return the mutated `pythia`so calls can be chained
# =============================================================================================


const default_ecm = 27.4                        # √s (GeV) for SHiP 400 GeV p on fix target
const atomic_mass_unit = 0.9314941              # GeV - u; nuclear m0 ≈ A⋅u (ignores binding E)

"""
    nucleus_pdg_id(Z, A) -> Int

PDG nuclear code `10-LZZZAAAI`(L=0 strange, I=0 isomer). Mo-96 → `1000420960`.
"""
nucleus_pdg_id(Z::Integer, A::Integer) = 1_000_000_000 + Z * 10_000 + A * 10

"""
    register_nucleus!(pythia, Z, A; m0=A*atomic_mass_unit, name="nuc\$(z)_\$(A)") -> Int

Register a nucleus so Angantyr accepts it as a beam. Pythia's `ParticleData.xml` predefines only ~16
nuclei (Mo-96 is not one), so any other target must be added before `init`. REturns the PDG id. The `m0`
default `A⋅u` ignores nuclear binding energy (a sub-permille efect on the beam kinematics); pass `m0`
for an exact mass.
"""
function register_nucleus!(pythia, Z::Integer, A::Integer; m0::Real = A * atomic_mass_unit,
        name::AbstractString = "nuc$(Z)_$(A)")
    id = nucleus_pdg_id(Z, A)
    # fields: id:new = name antiName spinType(=1) chargeType(=3Z) colType(=0) m0
    PYTHIA8.readString(pythia, "$id:new = $name $(name)bar 1 $(3 * Z) 0 $m0")
    return id
end

"""
    register_mo96!(pythia) -> Int

SHiP's molybdenum target with the canonical mass (89.3346 GeV) that reproduces the reference Angantyr
numbers exactly.
"""
register_mo96!(pythia) = register_nucleus!(pythia, 42, 96; m0 = 89.3346, name = "96Mo")

"""
    configure_beams!(pythia; idA=2212, idB=2212, ecm=default_ecm, quiet=true) -> pythia

CM-frame beams (`frameType = 1`). `idB` may be a nucleon (2212 p / 2112 n) or a registered nucleus id.
For a nucleus, register it first and add [`configure_angantyr!`](@ref).
"""
function configure_beams!(pythia; idA::Integer = 2212, idB::Integer = 2212, ecm::Real = default_ecm,
        quiet::Bool = true)
    PYTHIA8.readString(pythia, "Beams:frameType = 1")
    PYTHIA8.readString(pythia, "Beams:idA = $idA")
    PYTHIA8.readString(pythia, "Beams:idB = $idB")
    PYTHIA8.readString(pythia, "Beams:eCM = $ecm")
    quiet && PYTHIA8.readString(pythia, "Print:quiet = on")
    return pythia
end

"""
    configure_angantyr!(pythia; sigfit_par="10.03,15.60,7.72") -> pythia

Angantyr nuclear cross-section knobs for p+A: skip the per-init σ refit (`SigFitNGen = 0`) and resue a
converged parameter set (`SigFitDefPar`).
"""
function configure_angantyr!(pythia; sigfit_par::AbstractString = "10.03,15.60,7.72")
    PYTHIA8.readString(pythia, "HeavyIon:SigFitNGen = 0")
    PYTHIA8.readString(pythia, "HeavyIon:SigFitDefPar = $sigfit_par")
    return pythia
end

"""
    configure_charm!(pythia) -> pythia

Charm-enriched HardQCD (`gg/qqbar → ccbar`). Freezes D°, D⁺, Ds (`mayDecay = off`) so they are final
state and countable. The step that, if missed, yields "0 D".
"""
function configure_charm!(pythia)
    PYTHIA8.readString(pythia, "HardQCD:gg2ccbar = on")
    PYTHIA8.readString(pythia, "HardQCD:qqbar2ccbar = on")
    for id in("421", "411", "431")              # D°, D⁺, Ds
        PYTHIA8.readString(pythia, "$id:mayDecay = off")
    end
    return pythia
end

"""
    configure_kaon!(pythia) -> pythia

Min-bias SoftQCD inelastic with kaon decays enabled (K± and K°_L are stable by default in Pythia) and
π±/neutrons frozen, to isolate muons from kaon decays.
"""
function configure_kaon!(pythia)
    PYTHIA8.readString(pythia, "SoftQCD:inelastic = on")
    PYTHIA8.readString(pythia, "321:mayDecay = on")             # K± → μν
    PYTHIA8.readString(pythia, "130:mayDecay = on")             # K°_L → πμν
    PYTHIA8.readString(pythia, "211:mayDecay = off")            # freeze π±
    PYTHIA8.readString(pythia, "2112:mayDecay = off")           # freeze neutrons
    return pythia
end
