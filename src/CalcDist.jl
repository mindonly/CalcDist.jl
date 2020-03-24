module CalcDist

import StatsFuns: normpdf, normcdf, norminvcdf,
                  tdistpdf, tdistcdf, tdistinvcdf,
                  chisqpdf, chisqcdf, chisqinvcdf,
                  fdistpdf, fdistcdf, fdistinvcdf

export normalpdf, normalcdf, invNorm,
       tpdf, tcdf, invT,
       X2pdf, X2cdf, invX2,
       Fpdf, Fcdf, invF,
       binompdf, binomcdf,
       poissonpdf, poissoncdf,
       geometpdf, geometcdf,
       nCr, nPr

##
# Continuous
##

"""
   normalpdf(x, μ = 0, σ = 1) \n
      x: value
      μ: mean
      σ: standard deviation
   probability density function (Normal)
"""
function normalpdf(x, μ = 0, σ = 1)
   return normpdf(μ, σ, x)
end

"""
   normalcdf(lower, upper, μ = 0, σ = 1) \n
      lower: lower boundary
      upper: upper boundary
          μ: mean
          σ: standard deviation
   cumulative distribution function (Normal)
"""
function normalcdf(lower, upper, μ = 0, σ = 1)
   return normcdf(μ, σ, upper) - normcdf(μ, σ, lower)
end

"""
   invNorm(area, μ = 0, σ = 1) \n
      area: probability area
         μ: mean
         σ: standard deviation
   inverse cumulative distribution function (Normal)
"""
function invNorm(area, μ = 0, σ = 1)
   return norminvcdf(μ, σ, area)
end

"""
   tpdf(x, ν) \n
      x: value
      ν: degrees of freedom
   probability density function (Student's t)
"""
function tpdf(x, ν)
   return tdistpdf(ν, x)
end

"""
   tcdf(lower, upper, ν) \n
      lower: lower boundary
      upper: upper boundary
          ν: degrees of freedom
   cumulative distribution function (Student's t)
"""
function tcdf(lower, upper, ν)
   return tdistcdf(ν, upper) - tdistcdf(ν, lower)
end

"""
   invT(area, ν) \n
      area: probability area
         ν: degrees of freedom
   inverse cumulative distribution function (Student's t)
"""
function invT(area, ν)
   return tdistinvcdf(ν, area)
end

"""
   X2pdf(x, ν) \n
      x: value
      ν: degrees of freedom
   probability density function (Chi-squared)
"""
function X2pdf(x, ν)
   return chisqpdf(ν, x)
end

"""
   X2cdf(lower, upper, ν) \n
      lower: lower boundary
      upper: upper boundary
          ν: degrees of freedom
   cumulative distribution function (Chi-squared)
"""
function X2cdf(lower, upper, ν)
   return chisqcdf(ν, upper) - chisqcdf(ν, lower)
end

"""
   invX2(area, ν) \n
      area: probability area
         ν: degrees of freedom
   inverse cumulative distribution function (Chi-squared)
"""
function invX2(area, ν)
   return chisqinvcdf(ν, area)
end

"""
   Fpdf(x, nν, dν) \n
       x: value
      nν: numerator degrees of freedom
      dν: denominator degrees of freedom
   probability density function (F)
"""
function Fpdf(x, nν, dν)
   return fdistpdf(nν, dν, x)
end

"""
   Fcdf(lower, upper, nν, dν) \n
      lower: lower boundary
      upper: upper boundary
         nν: numerator degrees of freedom
         dν: denominator degrees of freedom
   cumulative distribution function (F)
"""
function Fcdf(lower, upper, nν, dν)
   return fdistcdf(nν, dν, upper) - fdistcdf(nν, dν, lower)
end

"""
   invF(area, nν, dν) \n
      area: probability area
        nν: numerator degrees of freedom
        dν: denominator degrees of freedom
   inverse cumulative distribution function (F)
"""
function invF(area, nν, dν)
   fdistinvcdf(nν, dν, area)
end

##
# Discrete
##

"""
   binompdf(n, p, x = nothing) \n
      n: trials
      p: probability of success
      x: value
   probability mass function (Binomial)
"""
function binompdf(n, p, x = nothing)
   if n ≤ 0
      throw(DomainError(n, "\t[n > 0] must be true."))
   end
   if p < 0 || p > 1
      throw(DomainError(p, "\t[0 ≤ p ≤ 1] must be true."))
   end

   if isnothing(x)
      x = range(0, stop=n)
   end

   probs = Vector{Real}()

   for i in x
      push!(probs, binomial(n, i) * p^i * (1-p)^(n-i))
   end

   return probs
end

"""
   binomcdf(n, p, x = nothing) \n
      n: trials
      p: probability of success
      x: value
   cumulative distribution function (Binomial)
"""
function binomcdf(n, p, x = nothing)
   if n ≤ 0
      throw(DomainError(n, "\t[n > 0] must be true."))
   end
   if p < 0 || p > 1
      throw(DomainError(p, "\t[0 ≤ p ≤ 1] must be true."))
   end

   if isnothing(x)
      x = range(0, stop=n)
   end

   probs = Vector{Real}()

   for i in x
      push!(probs, sum([ binomial(n, i) * p^i * (1-p)^(n-i) for i in 0:i ]))
   end

   return probs
end

"""
   poissonpdf(λ, x) \n
      λ: mean, variance
      x: value
   probability mass function (Poisson)
"""
function poissonpdf(λ, x)
   if λ ≤ 0
      throw(DomainError(λ, "\t[λ > 0] must be true."))
   end
   
   probs = Vector{Real}()

   for i in x
      push!(probs, ℯ^(-λ) * (λ^i / factorial(i)))
   end
   
   return probs
end

"""
   poissoncdf(λ, x) \n
      λ: mean, variance
      x: value
   cumulative distribution function (Poisson)
"""
function poissoncdf(λ, x)
   if λ ≤ 0
      throw(DomainError(λ, "\t[λ > 0] must be true."))
   end
   
   probs = Vector{Real}()

   for i in x
      push!(probs, ℯ^(-λ) * sum([ λ^i / factorial(i) for i in 0:i ]))
   end

   return probs
end

"""
   geometpdf(p, x) \n
      p: probability of success
      x: trial on which 1st success occurs
   probability density function (Geometric)
"""
function geometpdf(p, x)
   if p < 0 || p > 1
      throw(DomainError(p, "\t[0 ≤ p ≤ 1] must be true."))
   end

   probs = Vector{Real}()

   for i in x
      push!(probs, p * (1-p)^(i-1))
   end

   return probs
end

"""
   geometcdf(p, x) \n
      p: probability of success
      x: trial on which 1st success occurs
   cumulative distribution function (Geometric)
"""
function geometcdf(p, x)
   if p < 0 || p > 1
      throw(DomainError(p, "\t[0 ≤ p ≤ 1] must be true."))
   end

   probs = Vector{Real}()

   for i in x
      push!(probs, 1 - (1-p)^i)
   end

   return probs
end

##
# Counting
##

"""
   nCr(n, r) \n
      return number of combinations of n items taken r at a time. n and r must be
      non-negative. both n and r can be vectors.
"""
function nCr(n, r)
   if n < 0
      throw(DomainError(n, "\t[n ≥ 0] must be true."))
   end
   
   combs = Vector{Real}()

   for i in r
      if i < 0
         throw(DomainError(i, "\t[r ≥ 0] must be true."))
      end
      push!(combs, factorial(n) / ( factorial(i) * factorial(n-i) ))
   end

   return combs
end

function nCr(n::Vector{Int}, r)
   if r < 0
      throw(DomainError(r, "\t[r ≥ 0] must be true."))
   end
   
   combs = Vector{Real}()

   for i in n
      if i < 0
         throw(DomainError(i, "\t[n ≥ 0] must be true."))
      end
      push!(combs, factorial(i) / ( factorial(r) * factorial(i-r) ))
   end

   return combs
end

"""
   nPr(n, r)
      return number of permutations of n items taken r at a time. n and r must be
      non-negative. both n and r can be vectors.
"""
function nPr(n, r)
   if n < 0
      throw(DomainError(n, "\t[n ≥ 0] must be true."))
   end
   
   perms = Vector{Real}()

   for i in r
      if i < 0
         throw(DomainError(i, "\t[r ≥ 0] must be true."))
      end
      push!(perms, factorial(n) / factorial(n-i))
   end

   return perms
end

function nPr(n::Vector{Int}, r)
   if r < 0
      throw(DomainError(r, "\t[r ≥ 0] must be true."))
   end
   
   perms = Vector{Real}()

   for i in n
      if i < 0
         throw(DomainError(i, "\t[n ≥ 0] must be true."))
      end
      push!(perms, factorial(i) / factorial(i-r)) 
   end

   return perms
end

end # module
