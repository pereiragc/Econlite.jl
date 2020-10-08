abstract type AbstractSingleValuedUtility end
abstract type AbstractIsoElastic <: AbstractSingleValuedUtility end

dutil(c, ps::AbstractIsoElastic)=c^(-rra(ps))
invdutil(m, ps::AbstractIsoElastic)=m^(-1/rra(ps))

"
   IsoElastic

IsoElastic specification for a utility with one good.

The parametrization is either
```math
U(c) = \frac{c^{1-\alpha}-1}{1-\alpha}
```

or

```math
U(c) = \frac{c^{1-\alpha}}{1-\alpha}
```

(depending on whether IsoElasticVar2 is passed to `util`.)

The field `elast` corresponds to α.
"
struct IsoElastic <: AbstractIsoElastic
    elast::Float64
    function IsoElastic(elast)
        elast == 1 && error("Please use log utility specification")
        new(elast)
    end
end
struct LogUtility <: AbstractIsoElastic end
function IsoElastic(;a::Float64)
    a == 1 && return LogUtility()
    return(IsoElastic(a))
end


# Methods for iso-elastic types
struct IsoElasticVar2 end

util(c, ps::IsoElastic, ::Type{IsoElasticVar2})=c^(1-ps.elast)/(1-ps.elast)
util(c, ps::IsoElastic)=util(c,ps,IsoElasticVar2)-1/(1-ps.elast)
invutil(u, ps::IsoElastic)=(u*(1-ps.elast)+1)^(1/(1-ps.elast))

util(c, ::LogUtility)=log(c)
invutil(u, ::LogUtility)=exp(u)

rra(ps::IsoElastic)=ps.elast
rra(::LogUtility)=1.

" Constant absolute risk aversion specification "
struct CARA <: AbstractSingleValuedUtility
    a::Float64
end

util(c, ps::CARA)=(1 - exp(-ps.a * c))/ps.a
invutil(u, ps::CARA)=(-1/ps.a)*log(1 - ps.a*u)
dutil(c, ps::CARA)=exp(-ps.a * c)
invdutil(m, ps::CARA)=-log(m)/ps.a



" Risk neutral specification "
struct RiskNeutral <: AbstractIsoElastic end

util(c, ps::RiskNeutral)=c
invutil(u, ps::RiskNeutral)=u
dutil(c, ps::RiskNeutral)=one(c)
invdutil(c, ps::RiskNeutral)=error("Tried to invert constant marginal utility from risk neutral agent")
