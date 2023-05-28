function cleartrace!(params::TrussOptParams)
    empty!(params.losstrace)
    empty!(params.valtrace)
end