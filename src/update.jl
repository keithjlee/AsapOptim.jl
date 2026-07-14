"""
    updatemodel(p::OptParams, x) -> Model

Materialize a design vector back into the reference `Asap.Model`: write the
evaluated node positions and element sections into the model, re-solve, and
return it. The returned model is the ordinary mutable definition object —
ready for post-processing, force recovery, visualization, or export.

(The model's topology is untouched, so the frozen analysis cache is reused;
this is a numeric re-assembly plus one factorization.)
"""
function updatemodel(p::OptParams, x::AbstractVector)
    X, _, sections, _, ends = _design_state(x, p)
    for (i, node) in enumerate(p.model.nodes)
        node.position = SVector{3,Float64}(X[1, i], X[2, i], X[3, i])
    end
    for (i, el) in enumerate(p.model.elements)
        p.amask[i] && (el.section = sections[i])
        if ends !== nothing && el isa FrameElement &&
           (p.jslot1[i] != 0 || p.jslot2[i] != 0)
            el.ends = ends[i]
        end
    end
    # changed end conditions live in the frozen analysis cache (and can alter
    # DOF activity), so joint-variable designs need a re-process
    solve!(p.model; reprocess = ends !== nothing)
    return p.model
end
