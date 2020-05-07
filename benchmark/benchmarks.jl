# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# McCormick.jl
# A McCormick operator library in Julia
# See https://github.com/PSORLab/McCormick.jl
#############################################################################
# benchmark/benchmarks.jl
# Runs benchmarks on all supported McCormick operators.
# Allocations are expected for the following operations:
# - MC(Float64, Interval, 1): Nonstatic definition of static type.
# - exp2, exp10, tanh, atanh, acosh: Underlying interval operation may allocate.
#############################################################################

using BenchmarkTools, IntervalArithmetic

# don't really want to measure interval performance, so we're using the
# fastest interval rounding mode
setrounding(Interval, :none)
using McCormick, ReverseMcCormick

const SUITE = BenchmarkGroup()

for T in (NS, Diff, MV)
    begin
    end
end
