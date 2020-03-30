abstract type Tax end

net(x, T::Tax)=x-tax(x, T)

# * Gouveia and Strauss (1994)
struct GouveiaStrauss <: Tax
    maxmargtax :: Float64
    coefexp    :: Float64
    coefcons   :: Float64
end

function GouveiaStrauss(;maxmargtax=0.275,
                        coefexp=5.,
                        coefcons=0.3)
    GouveiaStrauss(maxmargtax, coefexp, coefcons)
end


tax(y, gs::GouveiaStrauss)=_gouveiastrauss(y, gs.maxmargtax, gs.coefexp, gs.coefcons)



@doc doc"
    _gouveiastrauss(y, φ₀, φ₁, φ₂)

Inspired by Gouveia and Strauss (1994)

Parametrizes income taxes with a differentiable function characterized by three
parameters (the phis). If
```math
  \tau_y(y) = \phi_0 \left[ y - \left( y^{-\phi_1} + \phi_2  \right)^{-\frac{1}{\phi_1}} \right]
```

The $\phi_0$ coefficient is the limit $\lim \tau_y^{\prime}(y)$ when
$y\to\infty$. The other coefficients govern the degree of progressiveness of the
function. (For example, when $\phi_2 \to 0$, no one is taxed.)
"
_gouveiastrauss(y, φ₀, φ₁, φ₂)=φ₀*(y-(y^(-φ₁)+φ₂)^(-1/φ₁))


# * Simple linear schedule

struct AdValorem <: Tax
    t::Float64
end

tax(y, lav::AdValorem)=y*lav.t



# * Null tax


struct NullTax <: Tax end
tax(y, ::NullTax)=zero(y)
