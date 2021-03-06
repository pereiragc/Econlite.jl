abstract type AbstractAggregator end


#=
INSTRUCTIONS FOR EXTENDING ABSTRACT CES AGGREGATOR

Subtypes should have:
- A type `shares` whose length determines the number of aggregated entities (e.g. goods)
- If you  want to use `aggr`, `gradient_aggr`, `partial_aggr` defined for AbstractCESAggregator:
+ `_aggr(x, t::SubType, NoSkip)`  should evaluate the aggregator at (iterable) point x.
                                  `aggr` will default to this specification.
+ `_aggr(x, t::SubType, leaveout::Int64)` should evaluate the aggregator leaving out element `j`

+ `_partial_aggr(x, j, Z, t::SubType)` should evaluate the aggregator's j-th partial derivative at,
                                       point x.
                                       The argument Z is the evaluated aggregator. Even if it's not
                                       not necessary, leave it there anyhow.
Otherwise, you have to define those methods yourself.

Moreover, users might expect the following:
+ `invpartial_aggr(x,j,k,t::SubType)` should solve the equation ∂ⱼF(ξ, x₋ⱼ) = k for ξ
=#
abstract type AbstractCESAggregator <: AbstractAggregator end
ngoods(aa::AbstractCESAggregator)=length(aa.shares)
share(j, aa::AbstractCESAggregator)=aa.shares[j]



aggr(x, aa::AbstractCESAggregator)=_aggr(x, aa, NoSkip)
# aggr(x, aa::AbstractCESAggregator, leaveout)=_aggr(x, aa, leaveout)

" Partial derivative of aggregator wrt coordinate `j`"
partial_aggr(x,j,aa::AbstractCESAggregator)=_partial_aggr(x,j,aggr(x, aa),aa)

function gradient_aggr!(gr, x, aa::AbstractCESAggregator)
    ev = aggr(x, aa)
    @inbounds for g in 1:length(aa)
        gr[g] = _partial_aggr(x, g, ev, aa)
    end
end

"Aggregator gradient (allocates) "
function gradient_aggr(x,aa::AbstractCESAggregator)
    gr = fill(0., ngoods(aa))
    gradient_aggr!(gr, x, aa)
    return gr
end

"Aggregator gradient (returns tuple)"
function gradient_aggr_t(x,aa::AbstractCESAggregator)
    ev = aggr(x, aa)
    map(g -> _partial_aggr(x, g, ev, aa), 1:ngoods(aa))
end


struct NoSkip end
# Inlining breaks the code. Seems to make execution (a bit) slower. Perhaps
# revert to regular checking?
@noinline cond_include(x, g, neutral, ::Type{NoSkip})=x
@noinline cond_include(x, g, neutral, n::Int64)=ifelse(g == n, neutral, x)
# @noinline cond_include(x, g, neutral, nn::NTuple{N, Int64}) where {N} = ifelse(g in nn, neutral, x)









# ARMINGTON AGGREGATOR (longest section) ---------------------------------------
struct Armington{N} <: AbstractCESAggregator
    es::Float64 # elasticity of substitution
    shares::NTuple{N, Float64}
    function Armington(es, shares)
        es == 1 && error("Please use Cobb-Douglas aggregator")
        any(shares .< 0.) && error("Negative shares not allowed")
        new{length(shares)}(es, Tuple(shares))
    end
end

function _aggr(x, aa::Armington, leaveout)
    @unpack shares, es = aa
    pow = 1 - 1/es

    r = 0.
    @inbounds for g in 1:ngoods(aa)
        # Include only if not left out
        r += cond_include(shares[g]*x[g]^pow, g, 0., leaveout)
    end
    return r^(1/pow)
end


doc"""
    invpartial_aggr(x,j,k,a::Armington)
Solve the equation
```math
\partial_j A(x_j, x_{-j}) = k
```
for xⱼ.

Note that entry x[j] is redundant.
"""
function invpartial_aggr(x,j,k,aa::Armington)
    pow = 1 - 1/aa.es
    vv = aggr(x, aa, j)

    t1 = (k/share(j, aa))^(aa.es-1)-share(j, aa)
    t1 ^= -1/pow

    return vv*t1
end

" Specialize inverse partial aggregate for 2 goods "
function invpartial_aggr(x, j, k, aa::Armington{2})
    pow = aa.es/(aa.es - 1)
    pow2 =  aa.es - 1

    w_mine = share(j, aa)
    w_other = share(3-j, aa)
    x_other = x[3-j]

    # Just check the formula... :(
    r = w_other/((k/w_mine)^pow2 - w_mine)
    r ^= pow
    r *= x_other
end


doc"""

Take advantage of the fact that each partial is

```math
\partial_j \mathcal A = \omega_j x_j^{-\frac{1}{\xi}}
\left[ \sum \omega_i x_i^{1-\frac{1}{\xi}}
      \right]^{\frac{1/\xi}{1 - 1/\xi}}
```

The term in square brackets corresponds to ``\mathcal A(x)^{\frac{1}{\xi}}``. In
`_partial_aggr`, we name this `aggr_val_e` (e stands for "along with exponent").

The extra argument for the function value is useful when computing the gradient
(need only compute it once).

"""
function _partial_aggr(x, j, aggr_val_e, aa::Armington)
    pow = 1/aa.es
    r = x[j]^(-pow) * share(j, aa)
    r *= aggr_val_e^pow
end














# COBB-DOUGLAS AGGREGATOR ------------------------------------------------------

# This is a special case of Armington when shares add up to one, and the
# elasticity of substitution parameter equals 1.

struct CobbDouglas{N} <: AbstractCESAggregator
    shares::NTuple{N, Float64}
    function CobbDouglas(shares)
        !all(shares .>= 0) && error("Please specify only positive shares")
        # !(sum(shares) == 1.) && error("Shares must add up to one")
        new{length(shares)}(shares)
    end
end

CobbDouglas(a::Float64)=CobbDouglas((a, 1-a))



function _aggr(x, cbd::CobbDouglas, leaveout)
    r = one(eltype(x))
    for i in 1:ngoods(cbd)
        r *= cond_include(x[i]^cbd.shares[i], i, 1., leaveout)
    end
    return r
end

function invaggr(x,j,y,cbd::CobbDouglas)
    r = y
    r /=  _aggr(x, cbd, j) # leave out entry j
    r ^= 1/share(j, cbd)
end

function invpartial_aggr(x,j,k,cbd::CobbDouglas)
    (k/shares(j, cbd)/_aggr(x, cbd, j))^(1/(shares(j, cbd)-1))
end

_partial_aggr(x, j, aggr_val, cbd::CobbDouglas)=shares(j, cbd)*aggr_val/x[j]

# Specialize 2 goods
aggr(x, cbd::CobbDouglas{2})=x[1]^share(1, cbd)*x[2]^share(2, cbd)

function invaggr(x, j, y, cbd::CobbDouglas{2})
    j_other = 3 - j
    r = y/x[j_other]^share(j_other, cbd)
    r ^= 1/share(j, cbd)
end

function partial_aggr(x, j, cbd::CobbDouglas{2})
    j_other = 3 - j
    share_mine = share(j, cbd)
    share_other = share(j_other, cbd)
    x_other = x[j_other]

    share_mine * x[j]^(share_mine-one(share_mine)) * x_other^share_other
end

function invpartial_aggr(x,j,k,cbd::CobbDouglas{2})
    aux = share(j, cbd) * x[3-j]^share(3-j, cbd)
    (k / aux)^(1/(share(j,cbd) - 1))
end







# LINEAR AGGREGATOR   ----------------------------------------------------------
# This when the elasticity of substitution is infinite.

struct AggLinear{N} <: AbstractCESAggregator
    shares::NTuple{N, Float64}
end

# aggr(x, la::AggLinear)=sum(la.shares.*x)
function aggr(x, la::AggLinear)
    for g in 1:ngoods(la)
        la.shares.*x
    end
end

invpartial_aggr(x, j, k, la::AggLinear)=error("Can't invert constant function")
partial_aggr(x,j,la::AggLinear)=la.shares[j]
gradient_aggr(x, la::AggLinear)=collect(la.shares)
