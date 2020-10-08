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
# src/reverse_operators/exponential.jl
# Contains definitions of reverse McCormick operators of exp, exp2, exp10, expm1
# log, log2, log10, log1p
#############################################################################

"""
$(SIGNATURES)

Reverse McCormick operator for `exp`.
"""
function exp_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `exp2`.
"""
function exp2_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log2(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `exp10`.
"""
function exp10_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ log10(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `expm1`.
"""
function expm1_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(-1.0, Inf)
    x = x ∩ log1p(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `log`.
"""
function log_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `log2`.
"""
function log2_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp2(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `log10`.
"""
function log10_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ exp10(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `log1p`.
"""
function log1p_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ expm1(y)
    y, x
end
