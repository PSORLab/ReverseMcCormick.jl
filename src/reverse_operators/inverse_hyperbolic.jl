# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# ReverseMcCormick.jl
# A Reverse McCormick operator library in Julia
# See https://github.com/PSORLab/ReverseMcCormick.jl
#############################################################################
# src/reverse_operators/inverse_hyperbolic.jl
# Contains definitions of reverse McCormick operators of asinh, acosh, atanh,
# asech, acsch, acoth
#############################################################################

"""
$(SIGNATURES)

Reverse McCormick operator for `asinh`.
"""
function asinh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ sinh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `acosh`.
"""
function acosh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ cosh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `atanh`.
"""
function atanh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ tanh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `asech`.
"""
function asech_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, Inf)
    x = x ∩ sech(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `acsch`.
"""
function acsch_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ csch(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `acoth`.
"""
function acoth_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ coth(y)
    y, x
end
