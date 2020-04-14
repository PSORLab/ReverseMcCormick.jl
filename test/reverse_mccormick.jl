# REVERSE OPERATORS

function MC_1_is_equal(y, x, tol)
    bool1 = isapprox(y.cc,x.cc,atol=tol)
    bool2 = isapprox(y.cv,x.cv,atol=tol)
    bool3 = isapprox(y.cv_grad[1], x.cv_grad[1], atol=tol)
    bool4 = isapprox(y.cc_grad[1], x.cc_grad[1], atol=tol)
    bool5 = isapprox(y.Intv.lo, x.Intv.lo, atol=tol)
    bool6 = isapprox(y.Intv.hi, x.Intv.hi, atol=tol)
    return (bool1 && bool2 && bool3 && bool4 && bool5 && bool6)
end

@testset "Reverse Multiplication" begin

    # THE BINARY OPERATOR
    a = MC{1,NS}(1.0, Interval{Float64}(0.4,3.0), 1)
    b = MC{1,NS}(Interval{Float64}(-10.0,-1.0))
    c = MC{1,NS}(2.0, Interval{Float64}(1.1,4.5), 1)

    aout1, bout1, cout1 = mul_rev(a,b,c)

    #@test bout1.Intv.lo == -10.0
    #@test bout1.Intv.hi == -1.0
    #@test cout1.Intv.lo == 1.1
    #@test cout1.Intv.hi == 4.5

    bout1 = MC{1,NS}(1.0, Interval{Float64}(0.4,3.0), 1)
    cout1 = MC{1,NS}(Interval{Float64}(-10.0,-1.0))
    aout1 = bout1*cout1

    aout1_a, bout1_a, cout1_a = mul_rev(aout1, bout1, cout1)

    MC_1_is_equal(aout1_a, aout1, 0.00001)
    MC_1_is_equal(bout1_a, bout1, 0.00001)
    MC_1_is_equal(cout1_a, cout1, 0.00001)

    bout2 = MC{1,NS}(1.0, Interval{Float64}(0.4,3.0), 1)
    cout2 = MC{1,NS}(Interval{Float64}(-10.0,-1.0))
    aout2 = 0.3*bout1*cout1+1.0

    aout2_a, bout2_a, cout2_a = mul_rev(aout2, bout2, cout2)

    MC_1_is_equal(aout2_a, aout2, 0.00001)
    MC_1_is_equal(bout2_a, bout2, 0.00001)
    @test cout2_a.Intv.lo == -10.0
    @test cout2_a.Intv.hi == -1.0
#    @test cout2_a.cv == -6.32
    @test cout2_a.cc == -1.0
    #@test cout2_a.cv_grad[1] == -8.38
    @test cout2_a.cc_grad[1] == 0.0

    # WITH FLOAT

end

#=
@testset "Reverse Addition" begin
end
=#

#=
@testset "Reverse Division" begin
end
=#

#=
@testset "Reverse Subtraction" begin
end
=#

@testset "Reverse Exponential" begin
    a = MC{1,NS}(1.0,Interval{Float64}(0.4,3.0),1)
    expa = exp(a)*1.1
    y,x = exp_rev(expa,a)

    @test MC_1_is_equal(expa, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test isapprox(x.Intv.lo, 0.49531, atol=1E-5)
    @test x.Intv.hi == 3.0
    @test x.cc_grad[1] == 1.0
    @test x.cv_grad[1] == 1.0

    exp2a = exp2(a)*1.1
    y,x = exp2_rev(exp2a,a)

    @test MC_1_is_equal(exp2a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test isapprox(x.Intv.lo, 0.537503, atol=1E-5)
    @test x.Intv.hi == 3.0
    @test x.cc_grad[1] == 1.0
    @test x.cv_grad[1] == 1.0

    exp10a = exp10(a)*1.1
    y,x = exp10_rev(exp10a,a)

    @test MC_1_is_equal(exp10a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test isapprox(x.Intv.lo, 0.441392, atol=1E-5)
    @test x.Intv.hi == 3.0
    @test x.cc_grad[1] == 1.0
    @test x.cv_grad[1] == 1.0

    expm1a = expm1(a)*1.1
    y,x = expm1_rev(expm1a,a)

    @test MC_1_is_equal(expm1a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test isapprox(x.Intv.lo, 0.432436, atol=1E-5)
    @test x.Intv.hi == 3.0
    @test x.cc_grad[1] == 1.0
    @test x.cv_grad[1] == 1.0
end

@testset "Reverse Logarithm" begin

   a = MC{1,NS}(9.5,9.6,Interval{Float64}(8.0,11.0))
   b = MC{1,NS}(9.5,9.6,Interval{Float64}(9.4,9.8))
   log1a = log(b)
   y1,x1 = log_rev(log1a,a)

   @test isapprox(9.5, x1.cv, atol=1E-5)
   @test isapprox(9.6, x1.cc, atol=1E-5)
   @test isapprox(9.399999999999997, x1.Intv.lo, atol=1E-4)
   @test isapprox(9.800000000000002, x1.Intv.hi, atol=1E-4)
   @test isapprox(0.0, x1.cv_grad[1], atol=1E-4)
   @test isapprox(0.0, x1.cc_grad[1], atol=1E-4)

   @test isapprox(2.2511278633761, y1.cv, atol=1E-5)
   @test isapprox(2.2617630, y1.cc, atol=1E-5)
   @test isapprox(2.240709689275958, y1.Intv.lo, atol=1E-4)
   @test isapprox(2.2823823856765264, y1.Intv.hi, atol=1E-4)
   @test isapprox(0.0, y1.cv_grad[1], atol=1E-4)
   @test isapprox(0.0, y1.cc_grad[1], atol=1E-4)

   a = MC{1,NS}(9.5,9.6,Interval{Float64}(8.0,11.0))
   b = MC{1,NS}(9.5,9.6,Interval{Float64}(9.4,9.8))
   log2a = log2(b)
   y2,x2 = log2_rev(log2a,a)

   @test isapprox(3.247691004899668, y2.cv, atol=1E-5)
   @test isapprox(3.2630344, y2.cc, atol=1E-5)
   @test isapprox(3.232660756790275, y2.Intv.lo, atol=1E-4)
   @test isapprox(3.292781749227846, y2.Intv.hi, atol=1E-4)
   @test isapprox(0.0, y2.cv_grad[1], atol=1E-4)
   @test isapprox(0.0, y2.cc_grad[1], atol=1E-4)

   a = MC{1,NS}(9.5,9.6,Interval{Float64}(8.0,11.0))
   b = MC{1,NS}(9.5,9.6,Interval{Float64}(9.4,9.8))
   log10a = log10(b)
   y3,x3 = log10_rev(log10a,a)

   @test isapprox(0.9776524091228977, y3.cv, atol=1E-5)
   @test isapprox(0.9822712, y3.cc, atol=1E-5)
   @test isapprox(0.9731278535996987, y3.Intv.lo, atol=1E-4)
   @test isapprox(0.991226075692495, y3.Intv.hi, atol=1E-4)
   @test isapprox(0.0, y3.cv_grad[1], atol=1E-4)
   @test isapprox(0.0, y3.cc_grad[1], atol=1E-4)

   a = MC{1,NS}(9.5,9.6,Interval{Float64}(8.0,11.0))
   b = MC{1,NS}(9.5,9.6,Interval{Float64}(9.4,9.8))
   log1pa = log1p(b)
   y4,x4 = log1p_rev(log1pa,a)

   @test isapprox(2.3512408881430384, y4.cv, atol=1E-5)
   @test isapprox(2.3608540, y4.cc, atol=1E-5)
   @test isapprox(2.3418058061473266, y4.Intv.lo, atol=1E-4)
   @test isapprox(2.3795461341301745, y4.Intv.hi, atol=1E-4)
   @test isapprox(0.0, y4.cv_grad[1], atol=1E-4)
   @test isapprox(0.0, y4.cc_grad[1], atol=1E-4)
end

@testset "Reverse Trignometric" begin
end

@testset "Reverse Hyperbolic" begin

   a = MC{1,NS}(2.5,2.6,Interval{Float64}(2.0,3.0))
   b = MC{1,NS}(2.4,2.6,Interval{Float64}(2.2,2.8))
   sinh1 = sinh(b)
   y1,x1 = sinh_rev(sinh1,a)
   @test isapprox(y1.cv, 5.466229213676094, atol=1E-7)
   @test isapprox(y1.cc, 6.946980626335908, atol=1E-7)
   @test isapprox(y1.Intv.lo, 4.457105170535894, atol=1E-7)
   @test isapprox(y1.Intv.hi, 8.191918354235915, atol=1E-7)
   @test isapprox(y1.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y1.cc_grad[1], 0.0, atol=1E-7)

   coshb = cosh(b)
   y2,x2 = cosh_rev(coshb,a)
   @test isapprox(y2.cv, 5.556947166965507, atol=1E-7)
   @test isapprox(y2.cc, 7.024455054206832, atol=1E-7)
   @test isapprox(y2.Intv.lo, 4.567908328898228, atol=1E-7)
   @test isapprox(y2.Intv.hi, 8.252728416861133, atol=1E-7)
   @test isapprox(y2.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y2.cc_grad[1], 0.0, atol=1E-7)

   tanh1 = tanh(b)
   y3,x3 = tanh_rev(tanh1,a)
   @test isapprox(y3.cv, 0.9813725934213438, atol=1E-7)
   @test isapprox(y3.cc, 0.9890274022010992, atol=1E-7)
   @test isapprox(y3.Intv.lo, 0.9757431300314515, atol=1E-7)
   @test isapprox(y3.Intv.hi, 0.9926315202011281, atol=1E-7)
   @test isapprox(y3.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y3.cc_grad[1], 0.0, atol=1E-7)

   @test isapprox(x1.cv, 2.5, atol=1E-7)
   @test isapprox(x1.cc, 2.6, atol=1E-7)
   @test isapprox(x1.Intv.lo, 2.19999, atol=1E-5)
   @test isapprox(x1.Intv.hi, 2.80001, atol=1E-5)

   a = MC{1,NS}(2.5,2.6,Interval{Float64}(1.3,4.5))
   b = MC{1,NS}(2.4,2.6,Interval{Float64}(2.3,3.4))

   asinhb = asinh(b)
   y1,x1 = asinh_rev(asinhb,a)
   @test isapprox(y1.cv, 1.6036967920467529, atol=1E-7)
   @test isapprox(y1.cc, 1.6837431439977444, atol=1E-7)
   @test isapprox(y1.Intv.lo, 1.570278543484978, atol=1E-7)
   @test isapprox(y1.Intv.hi, 1.9378792776645006, atol=1E-7)
   @test isapprox(y1.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y1.cc_grad[1], 0.0, atol=1E-7)

   acoshb = acosh(b)
   y2,x2 = acosh_rev(acoshb,a)
   @test isapprox(y2.cv, 1.5131824386442316, atol=1E-7)
   @test isapprox(y2.cc, 1.6094379124341003, atol=1E-7)
   @test isapprox(y2.Intv.lo, 1.47504478124142, atol=1E-7)
   @test isapprox(y2.Intv.hi, 1.894559012672298, atol=1E-7)
   @test isapprox(y2.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y2.cc_grad[1], 0.0, atol=1E-7)

   a = MC{1,NS}(0.31,0.34,Interval{Float64}(-0.2,0.5))
   b = MC{1,NS}(0.31,0.34,Interval{Float64}(0.3,0.4))

   atanhb = atanh(b)
   y3,x3 = atanh_rev(atanhb,a)
   @test isapprox(y3.cv, 0.3205454093019461, atol=1E-7)
   @test isapprox(y3.cc, 0.35517133459930783, atol=1E-7)
   @test isapprox(y3.Intv.lo, 0.3095196042031117, atol=1E-7)
   @test isapprox(y3.Intv.hi, 0.42364893019360184, atol=1E-7)
   @test isapprox(y3.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(y3.cc_grad[1], 0.0, atol=1E-7)

   @test isapprox(x3.cv, 0.31, atol=1E-7)
   @test isapprox(x3.cc, 0.34, atol=1E-7)
   @test isapprox(x3.Intv.lo, 0.29999999999999993, atol=1E-7)
   @test isapprox(x3.Intv.hi, 0.4000000000000001, atol=1E-7)
   @test isapprox(x3.cv_grad[1], 0.0, atol=1E-7)
   @test isapprox(x3.cc_grad[1], 0.0, atol=1E-7)

end
