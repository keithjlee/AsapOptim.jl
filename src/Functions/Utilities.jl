"""
    replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Replace the values of `values[indices]` with the values in `newvalues` in a differentiable way. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})
    
    # v2 = copy(values)
    # v2[indices] .= newvalues

    values[indices] .= newvalues

    return values
end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement.

df/dnewvalues = df/dvalues ⋅ dvalues / dnewvalues = v̄ ⋅ dvalues / dnewvalues

Where v̄ = [nvalues × 1], dvalues/dnewvalues = [nvalues × nnewvalues] so:

df/dnewvalues = (dvalues/dnewvalues)ᵀv̄ = [nnewvalues × 1]

Is simply the values of v̄ at the indices of the new values.
"""
function ChainRulesCore.rrule(::typeof(replacevalues), values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

    v = replacevalues(values, indices, newvalues)

    function replacevalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, replacevalues_pullback 
end


"""
    addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

Add the values of `increments` to the current values in `values` at `indices`. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    # v2 = copy(values)
    # v2[indices] .+= increments

    values[indices] .+= increments

    return values
end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement
"""
function ChainRulesCore.rrule(::typeof(addvalues), values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    v = addvalues(values, indices, increments)

    function addvalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, addvalues_pullback 
end