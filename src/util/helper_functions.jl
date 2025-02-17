"""
    CombineSets(X::Set{Set{Element}}) where {Element}

TBW
"""
function CombineSets(X)
  reduce(union!, X)
end

"""
    FindElement(X::Element, states::Set{Element}) where {Element}

TBW
"""
function FindElement(X::Element, states::Union{Set{Element},Vector{Element}}) where {Element}
  findfirst(x -> x == X, states |> collect)
end



"""
    FindElement(X::Element, states::Set{Element}) where {Element}

TBW
"""
function FindElement(X::Set{Element}, states::Set{Element}) where {Element}
  map(X |> collect) do x
    FindElement(x, states)
  end
end
export FindElement

"""
    FilterEmptyValues(X::Vector{Union{Nothing,T}}) where {T}

TBW
"""
function FilterEmptyValues(X::Vector{Union{Nothing,T}}) where {T}
  filter(x -> !isnothing(x), X) |> Vector
end

"""
    FilterEmptyValues(X::Vector{T}) where {T}

TBW
"""
function FilterEmptyValues(X::Vector{T}) where {T}
  filter(x -> !isnothing(x), X) |> Vector{T}
end

"""
    Make1D(X::Vector{Vector{ElementType}}) where {ElementType}

TBW
"""

function Make1D(X::Vector{Vector})
  vcat(X...)
end
"""
    Make1D(X::Vector{Vector{ElementType}}) where {ElementType}

TBW
"""
function Make1D(X::Vector{Vector{ElementType}}) where {ElementType}
  vcat(X...)
end

"""
    PrintMatSorted(A, X)

TBW
"""
function PrintMatSorted(A, X)
  idxs = sortperm([x[1] for x ∈ collect(X)])
  display(A[idxs, idxs])
end
export PrintMatSorted

"""
    PurgeSpacePerCent(p::Vector{T}, d::Number) where {T}

TBW
"""
function FindLowestValuesPercent(p::Vector{T}, d::Number) where {T}
  x_decimal = d / 100
  num_elements = ceil(Int, length(p) * x_decimal)
  sortperm(p)[1:num_elements]
end

"""
    PurgeSpacePerCent(p::Vector{T}, d::Number) where {T}

TBW
"""
function FindLowestValuesPercent(p::Vector{T}, d::Number) where {T}
  x_decimal = d / 100
  num_elements = ceil(Int, length(p) * x_decimal)
  sortperm(p)[1:num_elements]
end



function findLowestValuesPercent_naive(p::Vector{T}, d::Real) where {T<:Real}
  target_fraction = d / 100
  sorted_indices = sortperm(p)
  running_sum = zero(T)
  selected_indices = Int[]
  for idx in sorted_indices
    if running_sum + p[idx] > target_fraction
      break
    end
    push!(selected_indices, idx)
    running_sum += p[idx]
  end
  return selected_indices
end


