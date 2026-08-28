# =============================================================================================
# measurements.jl - generator-level event loop counters
#
# Each measurement takes a CONFIGURED + init()'d Pythia object (built with the helpers in
# pythia_setup.jl) and a GenerationFrame, runs n_events, and returns (config, results) - 
# self-describing NamedTuples. Scalar rates are the canonical outputs (regression-tested);
# the return η_lab / pT vectors let callers bin spectra without pulling a plotting dependency
# into the core.
# =============================================================================================

const d_meson_ids = (421, 411, 431)                     # D⁰, D⁺, Ds
const muon_id = 13

"""
    measure_charm(pythia, frame; n_events, window=ship_eta_window) -> (config, results)

Count final state D mesons (charm proxy) per generated event. D mesons must be frozen
(`configure_charm!`does this) or they decay and none are counted.

`results`: `n_gen`, `n_D`, `n_acc`, species counts `n_D0`/`n_Dp`/`n_Ds`, derived
`D_per_event`/`accept_frac`, and spectra vectors `η_lab`/`p_T`.
"""
function measure_charm(pythia, frame::GenerationFrame; n_events::Integer, window = ship_eta_window)
    n_gen = 0; n_D = 0; n_acc = 0; n_D0 = 0; n_Dp = 0; n_Ds = 0
    η_lab = Float64[]; pT = Float64[]

    for _ in 1:n_events
        PYTHIA8.next(pythia) || continue
        n_gen += 1
        event = PYTHIA8.event(pythia)
        for j in 1:PYTHIA8.size(event)
            p = event[j]
            PYTHIA8.isFinal(p) || continue
            pid = abs(PYTHIA8.id(p))
            pid in d_meson_ids || continue

            pT_val = PYTHIA8.pT(p)
            η = lab_eta(frame, PYTHIA8.pz(p), PYTHIA8.e(p), pT_val)
            push!(pT, pT_val); push!(η_lab, η)
            n_D += 1
            in_window(η, window) && (n_acc += 1)
            pid == 421 ? (n_D0 += 1) : pid == 411 ? (n_Dp += 1) : (n_Ds += 1)
        end
    end

    config = (; measurement = :charm, n_events, window, frame)
    results = (; n_gen, n_D, n_acc, n_D0, n_Dp, n_Ds,
               D_per_event = n_D / max(n_gen, 1),
               accept_frac = n_acc / max(n_D, 1),
               η_lab, pT)
    return (config, results)
end

"""
    measure_kaon_muons(pythia, frame; n_events, window=ship_eta_window) -> (config, results)

Count final state muons (from kaon decays) per generated collision. Use with `configure_kaon!`
(K decays on, π/n frozen)

`results`: `n_gen`, `n_mu`, `n_acc`, charge split `n_muplus`/`n_muminus`, derived
`mu_per_event`/`accept_frac`, and spectra vectors `η_lab`/`p_T`/`charge`.
"""
function measure_kaon_muons(pythia, frame::GenerationFrame; n_events::Integer, window = ship_eta_window)
    n_gen = 0; n_mu = 0; n_acc = 0; n_muplus = 0; n_muminus
    η_lab = Float64[]; pT = Float64[]; charge = Int[]
    
    for _ in 1:n_events
        PYTHIA8.next(pythia) || continue
        n_gen += 1
        event = PYTHIA8.event(pythia)
        for j in 1:PYTHIA8.size(event)
            p = event[j]
            PYTHIA8.isFinal(p) || continue
            id = PYTHIA8.id(p)
            abs(id) == muon_id || continue

            pT_val = PYTHIA8.pT(p)
            η = lab_eta(frame, PYTHIA8.pz(p), PYTHIA8.e(p), pT_val)
            push!(pT, pT_val); push!(η_lab, η)
            n_mu += 1
            in_window(η, window) && (n_acc += 1)
            # id = -13 is μ⁺ (anti-muon), +13 is μ⁻
            id == -13 ? (n_muplus += 1; push!(charge, +1)) : (n_muminus += 1; push!(charge, -1))
        end
    end

    config = (; measurement = :kaon_muons, n_events, window, frame)
    results = (; n_gen, n_mu, n_acc, n_muplus, n_muminus,
               mu_per_event = n_mu / max(n_gen, 1),
               accept_frac = n_acc / max(n_mu, 1),
               η_lab, pT, charge)
    return (config, results)
end

"""
    measure_multiplicity(pythia; n_events) -> (config, results)

Mean charged final state multiplicity (N_ch) per event (p+p ≈ 12, Angantyr p+Mo ≈ 21). `results`: `n_gen`,
`mean_ch`, and the per-event vector `nch`.
"""
function measure_multiplicity(pythia, n_events::Integer)
    n_gen = 0; total_ch = 0; nch = Int[]
   
    for _ in 1:n_events
        PYTHIA8.next(pythia) || continue
        n_gen += 1
        event = PYTHIA8.event(pythia)
        c = 0
        for j in 1:PYTHIA8.size(event)
            p = event[j]
            (PYTHIA8.isFinal(p) && PYTHIA8.isCharged(p)) || continue
            c += 1
        end
        push!(nch, c); total_h += c
    end

    config = (; measurement = :multiplicity, n_events)
    results = (; n_gen, mean_nch = total_ch / max(n_gen, 1), nch)
    return (config, results)
end
