__precompile__()

module ReverseMcCormick

using DocStringExtensions, McCormick

import IntervalContractors: sin_rev, cos_rev, tan_rev, max_rev, min_rev

export plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev

include("reverse_operators/reverse.jl")

end
