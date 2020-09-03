"""
A `SolutionObject` is typically a struct with multiple Arrays being iterated.

Example:

```
struct VFIObject{T} <: SolutionObject
     vfun::Array{Float64, T}
     pol::Array{Int64, T}
end
```

Importantly, you have to define an `unpack` method for the subtypes of
`SolutionObject`. That should return a Tuple with the objects that are used to
determine convergence.

In the example above, if one is doing conventional VFI, they would do

```
unpack(vfiobj::VFIObject)=(vfiobj.vfun,)
```

"""
abstract type SolutionObject end

errmsg = "Please define an `unpack` method for "
unpack(x::SolutionObject)=error(errmsg * "$(typeof(x))")


" Interface for governing whether and how many of the iterations are displayed. "
abstract type IterationDisplay end


" Display whenever current iteration is a multiple of `mod` "
struct DisplayMod <: IterationDisplay
    mod :: Int64
end

display_iter(i, dm::DisplayMod) =
    (mod(i, dm.mod) == 0) && println("Iteration $i")


struct NoDisplay <: IterationDisplay end

display_iter(i, ::NoDisplay)=nothing

@with_kw struct IterationSpec{T <: IterationDisplay}
    maxiter::Int64 = 3000
    tol::Float64 = 1e-8
    show_convergence::Bool = true
    iterdisp::T = NoDisplay()
    lastiter_showdist::Bool = false
end


" Conservative specification for value function iteration "
const _vfi_conserv = IterationSpec(maxiter=3000,
                             tol = 1e-9,
                             iterdisp = DisplayMod(100),
                             lastiter_showdist=false)



function basic_message(niter, loop_cont, iterspec)
    if iterspec.show_convergence
        if loop_cont
            println("Stopped due to maximum iterations (N = $(iterspec.maxiter))")
        else
            println("Converged after $niter iterations")
        end
    end
end

"""
    contract!(container, xiter, other_args, F!, vfi::IterationSpec)

Simple interface for contraction mappings when all iteration objects
relevant to a given problem are stored in a container object, but such
container object might also include other fields.

A common use case of this is a regular VFI procedure. `contract!` allows
you to store both policy and value functions in the same object (`container`)
but to only iterate on the value function. In regular VFI, the policy function
is obtained as a "residual" of the VFI step, and it might be efficient to save it.

`container` will typically be a 'container type' that stores objects in a value
function iteration procedure. For example, container might contain a (array
representation of) a value function and a policy function.

`xiter` is typically a subcollection of the objects in `container` that will be
actually iterated on.

`F!` should be the function that takes iteration objects and executes the VFI
step.

**IMPORTANT: Has side-effects on both `container` and `xiter`**

"""
function contract!(container, xiter, other_args, F!, vfi::IterationSpec=_vfi_conserv)
    niter = 0
    loop_cont = true

    while loop_cont && (niter <= vfi.maxiter)
        niter += 1

        # Execute iteration
        F!(container, xiter, other_args...)

        # Evaluate convergence
        loop_cont = ! within_tolerance(container,
                                       xiter, vfi.tol)

        # Diagnostics
        display_iter(niter, vfi.iterdisp)

        # Replace `xiter` with container's updated version
        transfer!(xiter, container)
    end

    if vfi.lastiter_showdist
        niter += 1
        F!(container, xiter, other_args...)
        dd = maxdist_unroll(xiter, container)
        println(" Maximum distance: $dd ")
    end

    basic_message(niter, loop_cont, vfi)
end

unpack(x::Tuple)=x

transfer!(target::SolutionObject,
          source::SolutionObject)=transfer!(unpack(target), unpack(source))


@unroll function transfer!(tup_target::Tuple, tup_source::Tuple)
    i = 0
    @unroll for mu in tup_target
        i+=1
        copyto!(mu, tup_source[i])
    end
end

transfer!(xiter::AbstractArray, container)=copyto!(xiter, container)

within_tolerance(a,b,tol)=error("`within_tolerance` not defined for $(typeof(a))")
within_tolerance(a::T,b::T,tol) where {T <: SolutionObject} =
    within_tolerance(unpack(a), unpack(b), tol)

#
# Implements efficient tolerance check for l^∞
function within_tolerance(a::AbstractArray,b,tol)
    for i in eachindex(b)
        if !within_tolerance(a[i],b[i],tol)
            return false
        end
    end
    return true
end


@inline within_tolerance(a::Float64, b, tol)=abs(b-a) <= tol
@inline within_tolerance(a::Int64, b::Int64, tol)=(a == b)

@unroll function within_tolerance_unroll(tup1::Tuple, tup2::Tuple, tol)
    i = 0

    @unroll for X in tup1
        i += 1
        Xpr = tup2[i]

        if !within_tolerance(X, Xpr, tol)
            return false
        end
    end

    return true
end


within_tolerance(tup1::Tuple, tup2::Tuple, tol)=within_tolerance_unroll(tup1, tup2, tol)


maxdist_unroll(x1::T, x2::T) where {T <: SolutionObject} =
    maxdist_unroll(unpack(x1), unpack(x2))

@unroll function maxdist_unroll(tup1, tup2)
    i = 0
    m = -Inf
    @unroll for X in tup1
        i += 1
        s = maximum(abs.(X .- tup2[i]))
        m = max(s, m)
    end

    m
end
