const HALF_PI = Interval{Float64}(pi)/2

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asin`.
"""
function asin_rev(y::MC, x::MC)
    if isempty(y)
        return y, y
    else
        y = y ∩ Interval{Float64}(-HALF_PI.lo, HALF_PI.hi)
        x = sin(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acos`.
"""
function acos_rev(y::MC, x::MC)
    if isempty(y)
        return y, y
    else
        y = y ∩ Interval{Float64}(0.0, HALF_PI.hi)
        x = x ∩ cos(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atan`.
"""
function atan_rev(y::MC, x::MC)
    if isempty(y)
        return y, y
    else
        y = y ∩ Interval{Float64}(-HALF_PI.lo, HALF_PI.hi)
        x = x ∩ tan(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asind`.
"""
function asind_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = asin_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acosd`.
"""
function acosd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = acos_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atand`.
"""
function atand_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = atan_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

# trivial definitions
for f in (:asec_rev, :acsc_rev, :acot_rev)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            return y, y
        end
        y, x
    end
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asecd`.
"""
function asecd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = asec_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acscd`.
"""
function acscd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = acsc_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acotd`.
"""
function acotd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    if isempty(a)
        return a, a
    end
    atemp, btemp = acot_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(bintv)
    a, b
end
