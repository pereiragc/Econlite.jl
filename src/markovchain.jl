abstract type AbstractMarkovChain{T} <: AbstractArray{T, 1} end
supp(mc::AbstractMarkovChain)=mc.support
Base.getindex(mc::AbstractMarkovChain, i)=supp(mc)[i]
Base.length(mc::AbstractMarkovChain)=length(supp(mc))
Base.IndexStyle(::Type{<:AbstractMarkovChain}) = IndexLinear()
transition(MC::AbstractMarkovChain)=MC.transition


" Unconditional mean of a MarkovChain "
Statistics.mean(mc::AbstractMarkovChain)=dot(markov_invariant(mc), supp(mc))



"""

Markov Chain type with "arbitrary" support (`T`) and pre-stored cumulative
transition matrix. Assumes that **rows** add up to one.

"""
struct MarkovChain{T} <: AbstractMarkovChain{T}
    support::Vector{T}
    transition::Matrix{Float64}
    transition_cumu::Matrix{Float64}
end

transition(MC::MarkovChain, i, j)=transition(MC)[i,j] # Make rows add up to one

" Simples constructor for any iterable `v`. Note: no sanity check for dimensions."
MarkovChain(v, M::Matrix)=MarkovChain([v[i] for i in 1:length(v)], M,
                                      cumsum(M, dims=2))



# Note: will also work with `Aggregate` (& anything for which `transition` is implemented)
function expect_markov(v::AbstractArray, current_state, MC)
    e = 0.

    for next_state in 1:length(MC)
        e += v[next_state] * transition(MC, current_state, next_state)
    end

    e
end

" TODO: unrolled version for arbitrary tuple of arrays "
function expect_markov(v1::AbstractArray, v2::AbstractArray, current_state, MC)
    e1 = e2 = 0.

    for next_state in 1:length(MC)
        e1 += v1[next_state] * transition(MC, current_state, next_state)
        e2 += v2[next_state] * transition(MC, current_state, next_state)
    end

    e1, e2
end

markov_invariant(mc::MarkovChain, tol=1e-12)=markov_invariant(mc.transition,tol)

" Returns eigenvector associated with A's last eigenvalue "
markov_invariant(A, tol)=begin
    ed = eigen(transpose(A))

    !within_tolerance(ed.values[end], 1, tol) &&
        error("No invariant distribution found")

    evec = ed.vectors[:, end]

    makeprob!(evec)

    return evec
end

makeprob!(u)=(u ./= sum(u))




"""

Simulate a Markov Chain  for `len` periods, starting with index `i0`

The `rand_tmp` vector should contain `len` draws from a uniform distribution.

Overwrites vector `path_prealloc`.

"""
function markov_genpath!(path_prealloc, mc::MarkovChain, i0, len, rand_tmp)
    for t in 1:len
        i0 = draw_next(mc, i0, rand_tmp[t])
        path_prealloc[t] = i0
    end
end


function markov_genpath(mc::MarkovChain, i0, len, seed=1234)
    seed!(seed)
    rand_tmp=rand(len)
    path_prealloc = fill(0, len)
    markov_genpath!(path_prealloc, mc, i0, len, rand_tmp)
    return path_prealloc
end



" Guaranteed to not be `nothing`; use with EXTREME caution "
function findfirstg(bvec)
    idx = 1

    for i in 1:length(bvec)
        if bvec[i]
            idx = i
            break
        end
    end

    return idx
end

" Guaranteed to not be `nothing`; use with EXTREME caution "
function findfirstpositive(nvec)
    j = 0
    while (nvec[j+1] < 0) && (j+1 < length(nvec))
        j += 1
    end
    j+1
end


draw_next(mc::MarkovChain, i0, r)=begin
    @views findfirstpositive(mc.transition_cumu[i0, :] .- r)
end

draw_next(mc::MarkovChain, i0)=draw_next(mc, i0, rand())
