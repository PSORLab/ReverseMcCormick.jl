half_pi(::T) = Interval{T}(pi)/2

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asin`.
"""
function asin_rev(y::MC, x::MC)
    h = lo(half_pi(Float64))
    y = y ∩ Interval{Float64}(-h, h)
    if ~isempty(y)
        x = sin(y)
    end
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acos`.
"""
function acos_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(0.0, hi(half_pi(Float64)))
    if ~isempty(y)
        x = x ∩ cos(y)
    end
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atan`.
"""
function atan_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(-lo(half_pi(Float64)), hi(half_pi(Float64)))
    if ~isempty(y)
        x = x ∩ tan(y)
    end
    y,x
end

# trivial definitions
for f in (:asec_rev, :acsc_rev, :acot_rev, :asecd, :acscd, :acotd)
    @eval function ($f)(y::MC, x::MC)
        isempty(y) && x = empty(x)
        y, x
    end
end
