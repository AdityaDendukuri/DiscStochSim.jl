"""
    RectLatticeBoundaryCondition(x::CartesianIndex, bounds::Tuple{<:Real,<:Real})

Return `true` if every component of `x` lies in `[bounds[1], bounds[2]]`.
"""
function RectLatticeBoundaryCondition(x::CartesianIndex, bounds::Tuple{<:Real,<:Real})
    lo, hi = bounds
    all(c -> lo <= c <= hi, Tuple(x))
end

"""
    RectLatticeBoundaryCondition(x::CartesianIndex, bound::Real)

Return `true` if every component of `x` lies in `[0, bound]`.
"""
function RectLatticeBoundaryCondition(x::CartesianIndex, bound::Real)
    all(c -> 0 <= c <= bound, Tuple(x))
end
