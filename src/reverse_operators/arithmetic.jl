# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# McCormick.jl
# A McCormick operator library in Julia
# See https://github.com/PSORLab/ReverseMcCormick.jl
#############################################################################
# src/reverse_operators/arithmetic.jl
# Contains definitions of reverse +, -, /, *, promotions, conversion, one,
# zero.
#############################################################################

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `real(b)`.
"""
function real_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ a
    a, b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `one(b)`.
"""
function one_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    1.0 ∉ a    && (return a, empty(b))
    a, b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `zero(b)`.
"""
function zero_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    0.0 ∉ a    && (return a, empty(b))
    a, b
end

const RAD2DEG = Interval{Float64}(pi)/Interval{Float64}(180.0)
"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `rad2deg(b)`. That is, it is the reverse
contractor of `a = (180/pi)*b`
"""
function rad2deg_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = RAD2DEG*a ∩ b
    a, b
end


const DEG2RAD = Interval{Float64}(180.0)/Interval{Float64}(pi)

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `deg2rad(b)`. That is, it is the reverse
contractor of `a = (pi/180)*b`
"""
function deg2rad_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = DEG2RAD*a ∩ b
    a, b
end


"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `b` +`c`
"""
function plus_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    b_new = b ∩ (a - c)
    c_new = c ∩ (a - b)
    a, b_new, c_new
end
function plus_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if c ∈ a - b
        b = b ∩ (a - c)
        return a, b, MC{N,T}(c)
    end
    empty(a), empty(a), empty(a)
end
function plus_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    a, b, c = plus_rev(a, b, c)
    a, c, b
end
function plus_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ a
    a, b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `b`- `c`
"""
function minus_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    b_new = b ∩ (a + c)
    c_new = c ∩ (b - a)
    a, b_new, c_new
end
function minus_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    b_new = b ∩ (a + c)
    a, b_new, MC{N,T}(c)
end

function minus_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    c_new = c ∩ (b - a)
    a, MC{N,T}(b), c_new
end


"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `-b``
"""
function minus_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b_new = b ∩ (-a)
    return a, b_new
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `b`*`c`
"""
function mult_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    c_new = (0.0 ∉ b) ? (c ∩ (a/b)) : c
    b_new = (0.0 ∉ c) ? (b ∩ (a/c)) : b
    a, b_new, c_new
end

function mult_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if !iszero(c)
        b = b ∩ (a*inv(c))
    elseif 0.0 ∉ a
        return empty(a), empty(a), empty(a)
    end
    a, b, MC{N,T}(c)
end

function mult_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if !iszero(c)
        b = b ∩ (a*inv(c))
    elseif 0.0 ∉ a
        return empty(a), empty(a), empty(a)
    end
    a, MC{N,T}(c), b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `b`/`c`
"""
function div_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    b = b ∩ (a * c)
    if ~isempty(b) && 0.0 ∉ a
        c = c ∩ (b / a)
    end
    a, b, c
end

function div_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    b = b ∩ (a * c)
    a, b, c
end

function div_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    b ∉ (a * c) && (return a, a, a)
    if 0.0 ∉ a
        c = c ∩ (b / a)
    end
    a, b, c
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `inv(b)`
"""
function inv_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    0.0 ∉ a && (b = b ∩ inv(a))
    a, b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `b`^`c`
"""
function power_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    (isempty(b) || isempty(c)) && (return a, a, a)
    0.0 ∉ c && (b = b ∩ a^inv(c))
    if 0.0 ∉ a
        if 0.0 < lo(b) < Inf
            blog = log(b)
            if 0.0 ∉ blog
                c = c ∩ (log(a) / blog)
            end
        end
    end
    a, b, c
end
function power_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if !isempty(c)
        if (0.0 ∉ a) && (b > zero(C))
            blog = log(b)
            if 0.0 ∉ blog
                c = c ∩ (log(a)/blog)
            end
        end
    else
        a = empty(a)
    end
    a, MC{N,T}(b), c
end

function int_power_rev(a::MC{N,T}, b::MC{N,T}, n::Int, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if n == 2
        root = sqrt(a)
        if lo(b) > 0.0
            b = b ∩ root
        else
            b = b ∩ -root
        end
    elseif iseven(n)
        root = a^(1/n)
        if lo(b) > 0.0
            b = b ∩ root
        else
            b = b ∩ -root
        end
    elseif isodd(n)
        if lo(b) > 0.0
            b = b ∩ (a ∩ Interval(0, Inf))^(1/n)
        else
            b = b ∩ -(((-a) ∩ Interval(0, Inf))^(1/n))
        end
    end
    a, b, MC{N,T}(c)
end

function power_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    (isempty(b) || 0.0 > lo(b)) && (return a, a, a)
    if isone(-c)
        a, b = inv_rev(a,b)
        return a, b, MC{N,T}(c)
    elseif !iszero(c)
        if isinteger(c)
            return int_power_rev(a, b, Int(c), c)
        elseif lo(b) > 0.0
            b = b ∩ a^(one(C)/c)
        end
    elseif iszero(c) && 1.0 ∉ a
        return empty(a), empty(a), empty(a)
    end
    a, b, MC{N,T}(c)
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `sqrt(b)`. That is
`b = b ∩ a^2`
"""
function sqrt_rev(a::MC, b::MC)
    isempty(a) && (return a, a)
    b = b ∩ (a^2)
    a, b
end
sqr_rev(f, x) = power_rev(f, x, 2)

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `abs(b)`
"""
function abs_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y, x
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `step(b)`
"""
function step_rev(a::MC, b::MC)
    a_lo = a.Intv.lo
    a_hi = a.Intv.hi
    if isempty(a) || 0.0 > a_hi || 1.0 < a_lo ||
        (0.0 < a_lo && a_hi < 1.0)
        b = empty(b)
        return a, b
    elseif 0.0 ∈ a && 1.0 ∉ a
        b = b ∩ Interval(-Inf, 0.0)
    elseif -1.0 ∉ a && 1.0 ∈ a
        b = b ∩ Interval(0.0, Inf)
    end
    a, b
end

"""
$(SIGNATURES)

Creates reverse McCormick contractor for `a` = `sign(b)`.
"""
function sign_rev(a::MC, b::MC)
    a_lo = a.Intv.lo
    a_hi = a.Intv.hi
    if isempty(a) || -1.0 > a_hi || 1.0 < a_lo || (-1.0 < a_lo && a_hi < 1.0)
        b = empty(b)
        return a, b
    elseif -1.0 ∈ a && 1.0 ∉ a
        b = b ∩ Interval(-Inf, 0.0)
    elseif -1.0 ∉ a && 1.0 ∈ a
        b = b ∩ Interval(0.0, Inf)
    end
    a, b
end
