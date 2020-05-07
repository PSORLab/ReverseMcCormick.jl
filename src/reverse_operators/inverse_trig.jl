const HALF_PI = Interval{Float64}(pi)/2

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asin`.
"""
function asin_rev(y::MC, x::MC)
    if isempty(y)
        return y, empty(x)
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
        return y, empty(x)
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
        return y, empty(x)
    else
        y = y ∩ Interval{Float64}(-HALF_PI.lo, HALF_PI.hi)
        x = x ∩ tan(y)
    end
    y, x
end

# trivial definitions
for f in (:asec_rev, :acsc_rev, :acot_rev, :asecd, :acscd, :acotd)
    @eval function ($f)(y::MC, x::MC)
        if isempty(y)
            return y, empty(x)
        end
        y, x
    end
end
