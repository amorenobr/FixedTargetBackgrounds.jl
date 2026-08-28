# =============================================================================================
# superposition.jl - p+A as a weighted p+p / p+n combination
#
# In the superposition approximation, the beam proton strikes one target nucleon, a proton with
# probability Z/A, a neutron with probability N/A. So a per-collision yield for p+A is the Z/A:
# N/A weighted average of the p+p and p+n yields. This is the approximation the true Angantyr
# p+Mo replaces, it is exact for hard charm, and ~1.8x low for soft kaon production.
# =============================================================================================

"""
    NuclearComposition(Z, N)
    NuclearComposition(; Z, A)

Proton (`Z`) and neutron (`N`) content of a target nucleus. Construct by mass number with
`NuclearComposition(Z=42, A=96)` (N = A - Z).
"""
struct NuclearComposition
    Z::Int
    N::Int
end
NuclearComposition(; Z::Integer, A::Integer) = NuclearComposition(Z, A - Z)

"proton fraction Z/A"
w_pp(c::NuclearComposition) = c.Z / (c.Z + c.N)
"neutron fraction N/A"
w_pn(c::NuclearComposition) = c.N / (c.Z + c.N)

"""
    superpose(x_pp, x_pn, comp) -> Float64

Z/A: N/A weighted average of a p+p and p+n per-collision quantity `x`. Generic in `x`, apply
it to any per-event yield (D/event, μ/event, a species sub-yield, an accepted yield) and
combine ratios by superposing numerator and denominator separately.
"""
superpose(x_pp::Real, x_pn::Real, comp::NuclearComposition) = w_pp(comp) * x_pp + w_pn(comp) * x_pn

# SHiP's molybdenum target: Mo-96 (Z=42, N=54) -> weights 42:54
const mo96 = NuclearComposition(Z = 42, A = 96)
