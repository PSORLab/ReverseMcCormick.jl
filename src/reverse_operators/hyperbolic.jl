"""
$(FUNCTIONNAME)

Reverse McCormick operator for `sinh`.
"""
function sinh_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    else
        x = x ∩ asinh(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `cosh`.
"""
function cosh_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    else
        y = y ∩ Interval{Float64}(1.0, Inf)
        x = x ∩ acosh(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `tanh`.
"""
function tanh_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    else
        y = y ∩ Interval{Float64}(-1.0, 1.0)
        x = x ∩ atanh(y)
    end
    y, x
end

# trivial definitions
for f in (:sech_rev, :csch_rev, :coth_rev)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            x = empty(x)
        end
        y, x
    end
end
