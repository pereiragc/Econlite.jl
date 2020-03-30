#
# GENERALIZED ARGMAX FUNCTION
# Implements the EVT-I business...


errmsg = "Please specify a log sum behavior. Usage is
`generalized_argmax!(y,x,sig, SumNormalize)` or
`generalized_argmax!(y,x,sig,SumRaw)`"

generalized_argmax!(y,x,sig)=error(errmsg)

function generalized_argmax!(y, x, sig, logsumspec)
    # `logsumspec` controls whether we subtract by greatest element of x or not
    # when computing stuff

    expsum = vexpsum(x, sig, logsumspec)  # Note: might be a number or a tuple
    # (depending on `logsumspec`)

    @inbounds for i in eachindex(x)
        y[i] = donorm(x[i], expsum, sig, logsumspec)
    end

    return genlog(expsum, sig, logsumspec)
end


struct SumNormalize end
struct SumRaw end

genlog(ss,sig,::Type{SumNormalize})=ss[2]+sig*log(ss[1])
genlog(ss,sig,::Type{SumRaw})=sig*log(ss)


vexpsum(x, sig, ::Type{SumNormalize})=vexpsum_norm(x,sig)
vexpsum(x, sig, ::Type{SumRaw})=vexpsum_raw(x,sig)

function vexpsum_norm(x::AbstractArray, sig)
    s = zero(eltype(x))
    maxel = maximum(x)
    @inbounds for i in eachindex(x)
        xnorm = x[i]-maxel
        s += exp(xnorm/sig)
    end
    return s, maxel
end

function vexpsum_raw(x::AbstractArray, sig)
    s = zero(eltype(x))
    @inbounds for i in eachindex(x)
        s += exp(x[i]/sig)
    end
    return s
end

donorm(xx, ss, sig, ::Type{SumNormalize})=exp((xx - ss[2])/sig)/ss[1]
donorm(xx, ss, sig, ::Type{SumRaw})=exp(xx/sig)/ss
