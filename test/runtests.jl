using CalcDist
using Test

# tolerance level
TOL8 = 1.0e-8
TOL6 = 1.0e-6
TOL4 = 1.0e-4
TOL = TOL8

@testset "CalcDist.jl" begin
   @testset "Binomial" begin
      # binompdf
      @test_throws DomainError binompdf(10, -0.1, 4)
      @test_throws DomainError binompdf(-5, 0.2, 3)
         # single
      @test only(binompdf(10, 0.5, 5)) == 0.24609375
      @test only(binompdf(20, 0.95, 8)) ≈ 2.0403e-11 atol=TOL
         # vector
      @test binompdf(5, 0.6, [3,4,5]) ≈ [0.3456, 0.2592, 0.07776] atol=TOL
         # missing x
      @test binompdf(5, 0.5) == [0.03125, 0.15625, 0.3125, 0.3125, 0.15625, 0.03125]

      # binomcdf
      @test_throws DomainError binomcdf(10, -0.1, 4)
      @test_throws DomainError binomcdf(-5, 0.2, 3)
         # single
      @test only(binomcdf(5, 0.6, 3)) ≈ 0.66304 atol=TOL
         # vector
      @test binomcdf(5, 0.6, [3,4,5]) ≈ [0.66304, 0.92224, 1.0] atol=TOL
         # missing x
      @test binomcdf(4, 0.25) ≈ [0.31640625, 0.73828125, 0.94921875, 0.99609375, 1.0] atol=TOL
   end

   @testset "Poisson" begin
      # poissonpdf
      @test_throws DomainError poissonpdf(-2, 5)
         # single
      @test only(poissonpdf(6, 10)) ≈ 0.0413030934 atol=TOL
         # vector
      @test poissonpdf(6, [8,9,10]) ≈ [0.10325773, 0.06883849, 0.04130309] atol=TOL

      # poissoncdf
      @test_throws DomainError poissoncdf(-4, 7)
         # single
      @test only(poissoncdf(0.1875, 2)) ≈ 0.99904486 atol=TOL
         # vector
      @test poissoncdf(0.126, [0,1,2,3]) ≈ [0.88161485, 0.99269832, 0.99969658, 0.99999050] atol=TOL
   end

   @testset "Geometric" begin
      # geometpdf
      @test_throws DomainError geometpdf(-0.5, 6)
      @test_throws DomainError geometpdf(1.2, 6)
         # single
      @test only(geometpdf(0.4, 6)) ≈ 0.031104 atol=TOL
         # vector
      @test geometpdf(0.25, [2,3,4]) == [0.1875, 0.140625, 0.10546875]

      # geometcdf
      @test_throws DomainError geometcdf(-0.4, [4,5,6])
      @test_throws DomainError geometcdf(1.4, [4,5,6])
         # single
      @test only(geometcdf(0.5, 2)) == 0.75
         # vector
      @test geometcdf(0.5, [1,2,3]) == [0.5, 0.75, 0.875]
   end

   @testset "Counting" begin
      # nCr
      @test_throws DomainError nCr(-9, 4)
      @test_throws DomainError nCr(9, [4,-5,6])
         # single
      @test only(nCr(9,4)) == 126
         # vector
      @test nCr(9, [5,4,3]) == [126, 126, 84]
      @test nCr([8,7,6], 3) == [56, 35, 20]

      # nPr
      @test_throws DomainError nCr(-8, 4)
      @test_throws DomainError nCr(8, [2,-4,5])
         # single
      @test only(nPr(7,3)) == 210
         # vector
      @test nPr(7, [1,2,3]) == [7, 42, 210]
      @test nPr([5,6,7], 2) == [20, 30, 42]
   end
end
