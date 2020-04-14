"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b` +`c`
"""
function plus_rev(a::MC, b::MC, c::MC)  # a = b + c
    bold = b
    b = b ∩ (a - c)
    c = c ∩ (a - bold)
    a, b, c
end
plus_rev(a,b,c) = plus_rev(promote(a,b,c)...)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`- `c`
"""
function minus_rev(a::MC, b::MC, c::MC)  # a = b - c
    bold = b
    b = b ∩ (a + c)
    c = c ∩ (bold - a)
    a, b, c
end
function minus_rev(a::MC, b::MC)  # a = -b
    b = b ∩ -a
    return a, b
end
minus_rev(a,b,c) = minus_rev(promote(a,b,c)...)
minus_rev(a,b) = minus_rev(promote(a,b)...)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`*`c`
"""
function mul_rev(a::MC, b::MC, c::MC)  # a = b * c
    bflag = (0.0 ∉ b.Intv)
    if bflag
        temp1 = a / b
        ((0.0 ∉ a.Intv) || bflag) && (c = c ∩ temp1)
    end
    cflag = (0.0 ∉ c.Intv)
    if cflag
        temp2 = a / c
        ((0.0 ∉ a.Intv) || cflag) && (b = b ∩ temp2)
    end
    a,b,c
end
mul_rev(a::MC{N,T},b::MC{N,T},c::Float64) where {N, T<:RelaxTag} = mul_rev(a,b,MC{N,T}(c))
mul_rev(a::MC{N,T},b::Float64,c::MC{N,T}) where {N, T<:RelaxTag} = mul_rev(a,MC{N,T}(b),c)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`/`c`
"""
function div_rev(a::MC, b::MC, c::MC)  # a = b / c
    b = b ∩ (a * c)
    if (~isempty(b) && ~in(0.0, a.Intv))
        c = c ∩ (b / a)
    end
    a,b,c
end
div_rev(a,b,c) = div_rev(promote(a,b,c)...)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `inv(b)`
"""
function inv_rev(a::MC, b::MC)  # a = inv(b)
    ~in(0.0, a.Intv) && (b = b ∩ inv(a))
    a,b
end
inv_rev(a,b) = inv_rev(promote(a,b)...)

"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `b`^`c`
"""
function power_rev(a::MC, b::MC, c::MC)  # a = b^c
    if (~isempty(b) && ~isempty(c))
        ~in(0.0, c.Intv) && (b = b ∩ (a^(inv(c))))
        if ~in(0.0, a.Intv)
            if (0.0 < b.Intv.lo < Inf)
                blog = log(b)
                if ~in(0.0, blog.Intv)
                    c = c ∩ (log(a) / blog)
                end
            end
        end
    else
        a = empty(a)
    end
    a,b,c
end
power_rev(a,b,c) = power_rev(promote(a,b,c)...)


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `sqrt(b)`
"""
function sqrt_rev(a::MC, b::MC)  # a = sqrt(b)
    A, B = sqrt_rev(a,b)
    b = B ∩ (a^2)
    a,b
end
sqr_rev(f, x)  = power_rev(f,x,2)


"""
$(FUNCTIONNAME)

Creates reverse McCormick contractor for `a` = `abs(b)`
"""
function abs_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    AIntv, BIntv = abs_rev(a, b)
    A = a ∩ MC{N,T}(aIntv)
    B = b ∩ MC{N,T}(bIntv)
    A,B
end
