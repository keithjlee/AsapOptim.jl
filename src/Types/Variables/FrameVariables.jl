const valid_properties = [:A, :Ix, :Iy, :J]

"""
    SectionVariable <: IndependentVariable

A variable tied to a geometric section property of an element. 

```
SectionVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
SectionVariable(element::Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
SectionVariable(element::Element, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
```
"""
mutable struct SectionVariable <: IndependentVariable
    i::Int64
    property::Symbol
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Float64

    function SectionVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
        @assert in(property, valid_properties) "Valid property symbols: :A, :Ix, :Iy, :J"

        new(elementindex, property, value, lowerbound, upperbound, 0)
    end

    function SectionVariable(element::Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
        @assert in(property, valid_properties) "Valid property symbols: :A, :Ix, :Iy, :J"

        new(element.elementID, property, value, lowerbound, upperbound, 0)
    end

    function SectionVariable(element::Element, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
        @assert in(property, valid_properties) "Valid property symbols: :A, :Ix, :Iy, :J"

        value = getproperty(element.section, property)
        @assert lowerbound ≤ value ≤ upperbound "Element property value not inside defined bounds"

        new(element.elementID, property, value, lowerbound, upperbound, 0)
    end
end