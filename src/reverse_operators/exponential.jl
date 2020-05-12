"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp`.
"""
function exp_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp2`.
"""
function exp2_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log2(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `exp10`.
"""
function exp10_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log10(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `expm1`.
"""
function expm1_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(-1.0, Inf)
    x = x ∩ log1p(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log`.
"""
function log_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log2`.
"""
function log2_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp2(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log10`.
"""
function log10_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp10(y)
    y, x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `log1p`.
"""
function log1p_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ expm1(y)
    y, x
end
