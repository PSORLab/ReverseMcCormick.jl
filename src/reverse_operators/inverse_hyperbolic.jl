"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asinh`.
"""
function asinh_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    else
        x = x ∩ sinh(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acosh`.
"""
function acosh_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
    else
        y = y ∩ Interval{Float64}(0.0, Inf)
        x = x ∩ cosh(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atanh`.
"""
function atanh_rev(y::MC, x::MC)
    x = x ∩ tanh(y)
    y, x
end

# trivial definitions
for f in (:asech_rev, :acsch_rev, :acoth_rev)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            return y, empty(x)
        end
        y, x
    end
end
