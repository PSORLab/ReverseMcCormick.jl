"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `real(b)`.
"""
function real_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    else
        b = b ∩ a
    end
    a, b
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `one(b)`.
"""
function one_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    elseif 1.0 ∉ a
        b = empty(b)
    end
    a, b
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `zero(b)`.
"""
function zero_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    elseif 0.0 ∉ a
        b = empty(b)
    end
    a, b
end

const RAD2DEG = Interval{Float64}(pi)/Interval{Float64}(180.0)
"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `rad2deg(b)`. That is, it is the reverse
contractor of `a = (180/pi)*b`
"""
function rad2deg_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    b = RAD2DEG*a ∩ b
    a, b
end


const DEG2RAD = Interval{Float64}(180.0)/Interval{Float64}(pi)
"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `deg2rad(b)`. That is, it is the reverse
contractor of `a = (pi/180)*b`
"""
function deg2rad_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    b = DEG2RAD*a ∩ b
    a, b
end


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b` +`c`
"""
function plus_rev(a::MC, b::MC, c::MC)
    if isempty(a)
        return a, a, a
    end
    bold = b
    b = b ∩ (a - c)
    c = c ∩ (a - bold)
    a, b, c
end
function plus_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if c ∈ a - b
        b = b ∩ (a - c)
        return a, b, MC{N,T}(c)
    end
    empty(MC{N,T}), empty(MC{N,T}), c
end
function plus_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    a, b, c = plus_rev(a, b, c)
    a, c, b
end
function plus_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    b = b ∩ a
    a, b
end


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`- `c`
"""
function minus_rev(a::MC, b::MC, c::MC)
    if isempty(a)
        return a, a, a
    end
    bold = b
    b = b ∩ (a + c)
    c = c ∩ (bold - a)
    a, b, c
end
function minus_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    b = b ∩ (a + c)
    a, b, MC{N,T}(c)
end

function minus_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    c = c ∩ (b - a)
    a, MC{N,T}(b), c
end


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `-b``
"""
function minus_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    b = b ∩ -a
    return a, b
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`*`c`
"""
function mul_rev(a::MC, b::MC, c::MC)
    if isempty(a)
        return a, a, a
    end
    bflag = 0.0 ∉ b
    if bflag
        temp1 = a / b
        ((0.0 ∉ a) || bflag) && (c = c ∩ temp1)
    end
    cflag = 0.0 ∉ c
    if cflag
        temp2 = a / c
        ((0.0 ∉ a) || cflag) && (b = b ∩ temp2)
    end
    a, b, c
end
function mul_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    mul_rev(a, b, MC{N,T}(c))
end
function mul_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    mul_rev(a, MC{N,T}(c), b)
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`/`c`
"""
function div_rev(a::MC, b::MC, c::MC)
    if isempty(a)
        return a, a, a
    end
    b = b ∩ (a * c)
    if ~isempty(b) && 0.0 ∉ a
        c = c ∩ (b / a)
    end
    a, b, c
end
function div_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    b = b ∩ (a * c)
    a, b, c
end
function div_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    if b ∉ (a * c)
        return a, a, a
    end
    if 0.0 ∉ a
        c = c ∩ (b / a)
    end
    a, b, c
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `inv(b)`
"""
function inv_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    0.0 ∉ a && (b = b ∩ inv(a))
    a, b
end
inv_rev(a,b) = inv_rev(promote(a,b)...)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`^`c`
"""
function power_rev(a::MC, b::MC, c::MC)
    if isempty(a)
        return a, a, a
    end
    if !isempty(b) && !isempty(c)
        0.0 ∉ c && (b = b ∩ (a^(inv(c))))
        if 0.0 ∉ a
            if (0.0 < lo(b) < Inf)
                blog = log(b)
                if 0.0 ∉ blog
                    c = c ∩ (log(a) / blog)
                end
            end
        end
    else
        a = empty(a)
    end
    a, b, c
end
function power_rev(a::MC{N,T}, b::C, c::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    if !isempty(c)
        if (0.0 ∉ a) && (b > zero(C))
            blog = log(b)
            c = c ∩ (log(a)/blog)
        end
    else
        a = empty(a)
    end
    a, MC{N,T}(b), c
end
function power_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    if isempty(a)
        return a, a, a
    end
    if !isempty(b) && 0.0 < lo(b) < Inf
            b = b ∩ exp(-c*log(b))
    else
        a = empty(a)
    end
    a, b, MC{N,T}(c)
end

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `sqrt(b)`. That is
`b = b ∩ a^2`
"""
function sqrt_rev(a::MC, b::MC)
    if isempty(a)
        return a, a
    end
    b = b ∩ a^2
    a, b
end
sqr_rev(f, x)  = power_rev(f, x, 2)


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `abs(b)`
"""
function abs_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    end
    y, x
end

"""
$(FUNCTIONNAME)

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
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `sign(b)`.
"""
function sign_rev(a::MC, b::MC)
    a_lo = a.Intv.lo
    a_hi = a.Intv.hi
    if isempty(a) || -1.0 > a_hi || 1.0 < a_lo ||
        (-1.0 < a_lo && a_hi < 1.0)
        b = empty(b)
        return a, b
    elseif -1.0 ∈ a && 1.0 ∉ a
        b = b ∩ Interval(-Inf, 0.0)
    elseif -1.0 ∉ a && 1.0 ∈ a
        b = b ∩ Interval(0.0, Inf)
    end
    a, b
end
