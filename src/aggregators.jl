abstract type AbstractAggregator end


abstract type AbstractCESAggregator <: AbstractAggregator end
ngoods(aa::AbstractCESAggregator)=length(aa.shares)
share(j, aa::AbstractCESAggregator)=aa.shares[j]


struct NoSkip end

aggr(x, aa::AbstractCESAggregator)=_aggr(x, aa, NoSkip)
aggr(x, aa::AbstractCESAggregator, leaveout)=_aggr(x, aa, leaveout)

# Seems to make execution (a bit) slower. Perhaps revert to regular checking?
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

" Partial derivative of aggregator wrt coordinate `j`"
partial_aggr(x,j,aa::Armington)=_partial_aggr(x,j,aggr(x, aa),aa)


"Aggregator gradient"
function gradient_aggr(x,aa::Armington)
    gr = fill(0., ngoods(aa))
    gradient_aggr!(gr, x, aa)
    return gr
end


" Mutating version of aggregator gradient (`gr` and `x` should be same dimension) "
function gradient_aggr!(gr,x,aa::Armington)
    aggr_val_e = aggr(x, aa)
    @inbounds for j in eachindex(gr)
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














# COBB-DOUGLAS AGGREGATOR ------------------------------------------------------
# This is a special case of Armington when shares add up to one.

struct CobbDouglas{N} <: AbstractCESAggregator
    shares::NTuple{N, Float64}
    function CobbDouglas(shares)
        !all(shares .>= 0) && error("Please specify only positive shares")
        # !(sum(shares) == 1.) && error("Shares must add up to one")
        new{length(shares)}(shares)
    end
end

CobbDouglas(a::Float64)=CobbDouglas((a, 1-a))



function _aggr(x, cobbdoug::CobbDouglas, leaveout)
    r = one(eltype(x))
    for i in 1:ngoods(cobbdoug)
        r *= cond_include(x[i]^cobbdoug.shares[i], i, 1., leaveout)
    end
    return r
end


partial_aggr(x, j, cobbdoug::CobbDouglas)=_partial_aggr(x,j,aggr(x,cobbdoug),cobbdoug)

function gradient_aggr(x,cobbdoug::CobbDouglas)
    gr = similar(x)
    gradient_aggr!(gr, x, cobbdoug)
    return gr
end

function invpartial_aggr(x,j,k,cobbdoug::CobbDouglas)
    (k/cobbdoug.shares[j]/aggr(x, cobbdoug, leaveout=j))^(1/(cobbdoug.shares[j]-1))
end

_partial_aggr(x, j, aggr_val, cobbdoug::CobbDouglas)=cobbdoug.shares[j]*aggr_val/x[j]
function gradient_aggr!(gr,x,cobbdoug::CobbDouglas)
    aggr_val = aggr(x, cobbdoug)
    for i in eachindex(x)
        gr[i]=_partial_aggr(x, i, aggr_val, cobbdoug)
    end
    nothing
end

# Specialize 2 goods
aggr(x, cobbdoug::CobbDouglas{2})=x[1]^cobbdoug.shares[1]*x[2]^cobbdoug.shares[2]
function partial_aggr(x, j, cobbdoug::CobbDouglas{2})
    jother = 3 - j
    share_mine = cobbdoug.shares[j]
    share_other = cobbdoug.shares[jother]
    x_other = x[jother]

    share_mine * x[j]^(share_mine-one(share_mine)) * x_other^share_other
end

function invpartial_aggr(x,j,k,cobbdoug::CobbDouglas{2})
    aux = share(j, cobbdoug) * x[3-j]^share(3-j, cobbdoug)
    (k / aux)^(1/(share(j,cobbdoug) - 1))
end