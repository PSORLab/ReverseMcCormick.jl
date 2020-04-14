# Currently no refinement on min/max... will add later

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `max`.
"""
function max_rev(a::MC{N,T}, b::MC{N,T}, c::MC{N,T}) where {N, T<:RelaxTag}
    aIntv, bIntv, cIntv = IntervalContractors.max_rev(a.intv, b.Intv, c.Intv)
    A = a ∩ MC{N,T}(aIntv)
    B = b ∩ MC{N,T}(bIntv)
    C = c ∩ MC{N,T}(cIntv)
    A,B,C
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `min`.
"""
function min_rev(a::MC{N,T}, b::MC{N,T}, c::MC{N,T}) where {N, T<:RelaxTag}
    aIntv, bIntv, cIntv = IntervalContractors.min_rev(a.intv, b.Intv, c.Intv)
    A = a ∩ MC{N,T}(aIntv)
    B = b ∩ MC{N,T}(bIntv)
    C = c ∩ MC{N,T}(cIntv)
    A,B,C
end
