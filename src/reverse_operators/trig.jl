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
        return a, a
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
        return a, a
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
        return a, a
    end
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `sind`.
"""
function sind_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = sin_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `cosd`.
"""
function cosd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = cos_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `tand`.
"""
function tand_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = tan_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

# trivial definitions
for f in (:sec_rev, :csc_rev, :cot_rev, :secd_rev, :cscd_rev, :cotd_rev)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            return y, y
        end
        y, x
    end
end
