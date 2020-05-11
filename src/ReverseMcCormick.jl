__precompile__()

module ReverseMcCormick

using DocStringExtensions, McCormick

import IntervalContractors: sin_rev, cos_rev, tan_rev, max_rev, min_rev

export plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sec_rev, csc_rev, cot_rev, asec_rev, acsc_rev, acot_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       sech_rev, csch_rev, coth_rev, asech_rev, acsch_rev, acoth_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev, real_rev, one_rev, zero_rev,
       deg2rad_rev, rad2deg_rev, step_rev, sign_rev

const NumberNotRelax = Union{Int16, Int32, Int64, Float16, Float32, Float64, Interval{Float64}}

include("reverse_operators/reverse.jl")

end
