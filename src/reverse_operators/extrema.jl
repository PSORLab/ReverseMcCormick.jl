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
# src/reverse_operators/extrema.jl
# Contains definitions of reverse min and max.
#############################################################################

"""
$(SIGNATURES)

Reverse McCormick operator for `max(::MC, ::MC)`. Note that this function is not
invertible. If `b` dominates `c`, the  `b = b ∩ a`. If `c` dominates `b`,  `c = c ∩ a`.
No update otherwise. TODO: See if constructing an reverse from the identity
`max(x,y) = 0.5*(x + y + abs(x - y))` yields a tigher relaxation in this case.
"""
function max_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    if hi(b) <= lo(c)
        c = c ∩ a
    elseif lo(b) >= hi(c)
        b = b ∩ a
    end
    return a, b, c
end

function max_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if c < lo(a)
        b = b ∩ a
    end
    a, b, MC{N,T}(c)
end

function max_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if c < lo(a)
        b = b ∩ a
    end
    a, MC{N,T}(c), b
end

"""
$(SIGNATURES)

Reverse McCormick operator for `min`. Note that this function is not
invertible. If `b` dominates `c`, the  `b = b ∩ a`. If `c` dominates `b`,  `c = c ∩ a`.
No update otherwise.  TODO: See if constructing an reverse from the identity
`min(x, y) = 0.5*(x + y - abs(x - y))` yields a tigher relaxation in this case.
"""
function min_rev(a::MC, b::MC, c::MC)
    isempty(a) && (return a, a, a)
    if hi(b) <= lo(c)
        b = b ∩ a
    elseif lo(b) >= hi(c)
        c = c ∩ a
    end
    a, b, c
end

function min_rev(a::MC{N,T}, b::MC{N,T}, c::C) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if c > hi(a)
        b = b ∩ a
    end
    a, b, MC{N,T}(c)
end

function min_rev(a::MC{N,T}, c::C, b::MC{N,T}) where {N, T<:RelaxTag, C<:NumberNotRelax}
    isempty(a) && (return a, a, a)
    if c > hi(a)
        b = b ∩ a
    end
    a, MC{N,T}(c), b
end
