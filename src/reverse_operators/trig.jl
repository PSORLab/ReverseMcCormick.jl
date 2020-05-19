"""
$(SIGNATURES)

Reverse McCormick operator for `sin`.
"""
function sin_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    aintv, bintv = sin_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `cos`.
"""
function cos_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    aintv, bintv = cos_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `tan`.
"""
function tan_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    aintv, bintv = tan_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `sind`.
"""
function sind_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = sin_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `cosd`.
"""
function cosd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = cos_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `tand`.
"""
function tand_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = tan_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `sec`.
"""
function sec_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ asec(a)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `csc`.
"""
function csc_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ acsc(a)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `csc`.
"""
function cot_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ acot(a)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `secd`.
"""
function secd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = sec_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `cscd`.
"""
function cscd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = csc_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `cscd`.
"""
function cotd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = cot_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end
