# Currently no refinement on min/max... will add later

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `max(::MC, ::MC)`. Note that this function is not
invertible. If `b` dominates `c`, the  `b = b ∩ a`. If `c` dominates `b`,  `c = c ∩ a`.
No update otherwise. TODO: See if constructing an reverse from the identity
`max(x,y) = 0.5*(x + y + abs(x - y))` yields a tigher relaxation in this case.
"""
function max_rev(a::MC{N,T}, b::MC{N,T}, c::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a, a)
    if hi(b) <= lo(c)
        c = c ∩ a
    elseif lo(b) >= hi(c)
        b = b ∩ a
    end
    return a, b, c
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `min`. Note that this function is not
invertible. If `b` dominates `c`, the  `b = b ∩ a`. If `c` dominates `b`,  `c = c ∩ a`.
No update otherwise.  TODO: See if constructing an reverse from the identity
`min(x, y) = 0.5*(x + y - abs(x - y))` yields a tigher relaxation in this case.
"""
function min_rev(a::MC{N,T}, b::MC{N,T}, c::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a, a)
    if hi(b) <= lo(c)
        b = b ∩ a
    elseif lo(b) >= hi(c)
        c = c ∩ a
    end
    a, b, c
end
