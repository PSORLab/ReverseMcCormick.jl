const HALF_PI = Interval{Float64}(pi)/2

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asin`.
"""
function asin_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(-HALF_PI.lo, HALF_PI.hi)
    x = sin(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acos`.
"""
function acos_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, HALF_PI.hi)
    x = x ∩ cos(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atan`.
"""
function atan_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(-HALF_PI.lo, HALF_PI.hi)
    x = x ∩ tan(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asind`.
"""
function asind_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = asin_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acosd`.
"""
function acosd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = acos_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atand`.
"""
function atand_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = atan_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asec`.
"""
function asec_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ sec(a)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acsc`.
"""
function acsc_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ csc(a)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acot`.
"""
function acot_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    a = a ∩ Interval(0.0, pi)
    b = b ∩ cot(a)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asecd`.
"""
function asecd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = asec_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acscd`.
"""
function acscd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = acsc_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acotd`.
"""
function acotd_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    isempty(a) && (return a, a)
    atemp, btemp = acot_rev(a, deg2rad(b))
    btemp = b ∩ rad2deg(btemp)
    a = a ∩ atemp
    b = b ∩ MC{N,T}(btemp)
    a, b
end
