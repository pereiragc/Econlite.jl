abstract type AbstractAggregator end


abstract type AbstractCESAggregator <: AbstractAggregator end

is_crs(::AbstractCESAggregator)=true
ngoods(aa::AbstractCESAggregator)=length(aa.shares)

struct Armington{N} <: AbstractCESAggregator
    es::Float64 # elasticity of substitution
    shares::NTuple{N, Float64}
    function Armington(es, shares)
        es == 1 && error("Please use Cobb-Douglas aggregator")
        any(shares .< 0.) && error("Negative shares not allowed")
        new{length(shares)}(es, Tuple(shares))
    end
end


struct NoSkip end
aggr(x, aa::Armington)=_aggr(x, aa.shares, 1-1/aa.es, NoSkip)
aggr(x, aa::Armington, leaveout)=_aggr(x, aa.shares, 1-1/aa.es, leaveout)

function _aggr(x, shares, pow, leaveout)
    r = 0.
    @inbounds for g in 1:length(shares)
        # Include only if not left out
        r += cond_include(shares[g]*x[g]^pow , g, leaveout)
    end
    return r^(1/pow)
end

@noinline cond_include(x, g, ::Type{NoSkip})=x
@noinline cond_include(x, g, n::Int64)=ifelse(g == n, 0., x)
@noinline cond_include(x, g, nn::NTuple{N, Int64}) where {N} =ifelse(g in nn, 0., x)


" Partial derivative of aggregator wrt coordinate `j`"
partial_aggr(x,j,aa::Armington)=_partial_aggr(x,j,aggr(x, aa),aa)


"Aggregator gradient"
function gradient_aggr(x,aa::Armington)
    gr = fill(0., size(x))
    _gradient_aggr!(gr, x, aa)
    return gr
end


" Mutating version of aggregator gradient (`gr` and `x` should be same dimension) "
function gradient_aggr!(gr,x,aa::Armington)
    aggr_val_e = aggr(x, aa)
    for j in eachindex(gr)
        gr[j]=_partial_aggr(x,j,aggr_val_e,aa)
    end
    nothing
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

    t1 = (k/aa.shares[j])^(aa.es-1)-aa.shares[j]
    t1 ^= -1/pow

    return vv*t1
end

" Specialize inverse partial aggregate for 2 goods "
function invpartial_aggr(x, j, k, aa::Armington{2})
    pow = aa.es/(aa.es - 1)
    pow2 =  aa.es - 1

    w_mine = aa.shares[j]
    w_other = aa.shares[3-j]
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
    r = x[j]^(-pow) * aa.shares[j]
    r *= aggr_val_e^pow
end
