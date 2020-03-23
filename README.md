# CalcDist.jl

A Julia package implementing probability distribution functions and other bits according to the
interface found on some commercially-available graphing calculators.

### Distributions
```
* normalpdf(x, μ, σ)
* normalcdf(lower, upper, μ, σ)
* invNorm(area, μ, σ)
* tpdf(x, ν)
* tcdf(lower, upper, ν)
* X2pdf(x, ν)
* X2cdf(lower, upper, ν)
* Fpdf(x, nν, dν)
* Fcdf(lower, upper, nν, dν)
* binompdf(n, p, x)
* binomcdf(n, p, x)
* poissonpdf(λ, x)
* poissoncdf(λ, x)
* geometpdf(p, x)
* geometcdf(p, x)
```

### Counting
```
* nCr(n, r)
* nPr(n, r)
```
