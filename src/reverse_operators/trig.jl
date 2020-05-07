"""
$(FUNCTIONNAME)

Reverse McCormick operator for `sin`.
"""
function sin_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if !isempty(a)
        aintv, bintv = sin_rev(a.Intv, b.Intv)
        a = a ∩ MC{N,T}(aintv)
        b = b ∩ MC{N,T}(bintv)
    else
        b = isempty(b)
    end
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `cos`.
"""
function cos_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if !isempty(a)
        aintv, bintv = cos_rev(a.Intv, b.Intv)
        a = a ∩ MC{N,T}(aintv)
        b = b ∩ MC{N,T}(bintv)
    else
        b = isempty(b)
    end
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `tan`.
"""
function tan_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if !isempty(a)
        aintv, bintv = tan_rev(a.Intv, b.Intv)
        a = a ∩ MC{N,T}(aintv)
        b = b ∩ MC{N,T}(bintv)
    else
        b = isempty(b)
    end
    a, b
end


# trivial definitions
for f in (:sec_rev, :csc_rev, :cot_rev, :secd_rev, :cscd_rev, :cotd_rev)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            x = empty(x)
        end
        y, x
    end
end
