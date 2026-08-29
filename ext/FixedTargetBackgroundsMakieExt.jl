module FixedTargetBackgroundsMakieExt

using FixedTargetBackgrounds
using FixedTargetBackgrounds: ship_eta_window
using Makie

# Minimal histogram over a range of bin edges
function _hist(x, edges)
    n = length(edges) -1
    centers = [(edges[i] + edges[i + 1]) / 2 for i in 1:n]
    counts = zeros(Float64, n)
    for v in x
        i = searchsortedlast(edges, v)
        1 <= i <= n && (counts[i] += 1)
    end
    return centers, counts
end

_labels(n, given) = given === nothing ? ["series $i" for i in 1:n] : collect(given)


function FixedTargetBackgrounds.plot_eta(results::NamedTuple...; labels = nothing, edges = range(-2, 10; length = 121),
        per_event::Bool = true, window = ship_eta_window)
    labs = _labels(length(results), labels)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "η (lab frame)",
              ylabel = per_event ? "dN/dη per event" : "count",
              title = "Lab-frame Pseudorapidity")
    bw = step(edges)
    for (r, lab) in zip(results, labs)
        c, n = _hist(r.η_lab, edges)
        y = per_event ? n ./ max(r.n_gen, 1) ./ bw : n
        lines!(ax, c, y; linewidth = 2, label = lab)
    end
    vlines!(ax, [window[1], window[2]]; color = :gray, linestyle = :dash, label = "acceptance")
    axislegend(ax)
    return fig
end

function FixedTargetBackgrounds.plot_pT(results::NamedTuple...; labels = nothing, edges = range(0, 8; length = 81), xmax = nothing)
    labs = _labels(length(results), labels)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "pT [GeV]",
              ylabel = "probability density", yscale = log10,
              title = "Transverse Momentum")
    bw = step(edges)
    for (r, lab) in zip(results, labs)
        c, n = _hist(r.pT, edges)
        area = sum(n) * bw
        y = n ./ max(area, eps())
        m = y .> 0              # log-y: drop empty bins
        lines!(ax, c[m], y[m]; linewidth = 2, label = lab)
    end
    xmax !== nothing && xlims!(ax, 0, xmax)
    axislegend(ax)
    return fig
end

function FixedTargetBackgrounds.plot_multiplicity(results::NamedTuple...; labels = nothing, edges = range(0, 60; length = 61))
    labs = _labels(length(results), labels)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "N_ch",
              ylabel = "fraction of events",
              title = "Charge Multiplicity")
    for (r, lab) in zip(results, labs)
        c, n = _hist(r.nch, edges)
        y = n ./ max(sum(n), 1)
        stairs!(ax, c, y; linewidth = 2, label = lab, step = :center)
    end
    axislegend(ax)
    return fig
end

end # module
