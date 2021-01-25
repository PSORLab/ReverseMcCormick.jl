# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# ReverseMcCormick.jl
# A Reverse McCormick operator library in Julia
# See https://github.com/PSORLab/ReverseMcCormick.jl
#############################################################################
# src/reverse_operators/hyperbolic.jl
# Contains definitions of reverse McCormick operators of sinh, cosh, tanh,
# sech, csch, coth
#############################################################################


"""
$(SIGNATURES)

Reverse McCormick operator for `sinh`.
"""
function sinh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ asinh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `cosh`.
"""
function cosh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(1.0, Inf)
    x = x ∩ acosh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `tanh`.
"""
function tanh_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(-1.0, 1.0)
    x = x ∩ atanh(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `sech`.
"""
function sech_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    y = y ∩ Interval{Float64}(0.0, 1.0)
    x = x ∩ asech(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `csch`.
"""
function csch_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ acsch(y)
    y, x
end

"""
$(SIGNATURES)

Reverse McCormick operator for `coth`.
"""
function coth_rev(y::MC, x::MC)
    isempty(y) && (return y, y)
    x = x ∩ coth(y)
    y, x
end
