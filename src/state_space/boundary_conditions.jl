
# boundary condition function
function RectLatticeBoundaryCondition(X::CartesianIndex, bound::Tuple)
  all_pos = all(Tuple(X) .> 0)
  lower = all(Tuple(X) .>= bound[1])
  upper = all(Tuple(X) .<= bound[2])
  lower && upper && all_pos
end

function RectLatticeBoundaryCondition(X::CartesianIndex, bound::Int)
  lower = all(Tuple(X) .> 0)
  upper = all(Tuple(X) .< bound)
  lower && upper
end

export RectLatticeBoundaryCondition
