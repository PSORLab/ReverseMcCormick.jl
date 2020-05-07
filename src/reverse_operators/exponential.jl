"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp`.
"""
function exp_rev(y::MC, x::MC)
    if lo(y) > 0.0
        y = y ∩ Interval{Float64}(0.0, Inf)
        x = x ∩ log(y)
    else
        x = empty(x)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp2`.
"""
function exp2_rev(y::MC, x::MC)
    if lo(y) > 0.0
        y = y ∩ Interval{Float64}(0.0, Inf)
        x = x ∩ log2(y)
    else
        x = empty(x)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp10`.
"""
function exp10_rev(y::MC, x::MC)
    if lo(y) > 0.0
        y = y ∩ Interval{Float64}(0.0, Inf)
        x = x ∩ log10(y)
    else
        x = empty(x)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `expm1`.
"""
function expm1_rev(y::MC, x::MC)
    if lo(y) > -1.0
        y = y ∩ Interval{Float64}(-1.0, Inf)
        x = x ∩ log1p(y)
    else
        x = empty(x)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log`.
"""
function log_rev(y::MC, x::MC)
    if isempty(y)
        x = empty(x)
    else
        x = x ∩ exp(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log2`.
"""
function log2_rev(y::MC, x::MC)
    if isempty(y)
        x = empty(x)
    else
        x = x ∩ exp2(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log10`.
"""
function log10_rev(y::MC, x::MC)
    if isempty(y)
        x = empty(x)
    else
        x = x ∩ exp10(y)
    end
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log1p`.
"""
function log1p_rev(y::MC, x::MC)
    if isempty(y)
        x = empty(x)
    else
        x = x ∩ expm1(y)
    end
    y, x
end
