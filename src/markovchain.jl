"""
     IndexedMarkovChain{T}

Representation of a Markov Chain type with "arbitrary" support (`T`). Assumes
that **rows** add up to one.

"""
abstract type IndexedMarkovChain{T} <: AbstractArray{T, 1} end
Base.getindex(mc::IndexedMarkovChain, i)=supp(mc)[i]
Base.length(mc::IndexedMarkovChain)=length(supp(mc))
Base.size(mc::IndexedMarkovChain)=(length(mc), )
Base.IndexStyle(::Type{<:IndexedMarkovChain}) = IndexLinear()
supp(mc::IndexedMarkovChain, i)=error("Method `supp` must be implemented for $(typeof(MC)) ")
transition(mc::IndexedMarkovChain, i, j)=transition(mc)[i,j]


" Unconditional mean of a MarkovChain "
Statistics.mean(mc::IndexedMarkovChain)=dot(markov_invariant(mc), supp(mc))


"Markov Chain with pre-stored cumulative transition matrix"
struct MarkovChain{T} <: IndexedMarkovChain{T}
    support::Vector{T}
    transition::Matrix{Float64}
    transition_cumu::Matrix{Float64}
end
supp(mc::MarkovChain)=mc.support
transition(MC::MarkovChain)=MC.transition # Make rows add up to one
transition_cumu(MC::MarkovChain)=MC.transition_cumu


" Simples constructor for any iterable `v`. Note: no sanity check for dimensions."
MarkovChain(v, M::Matrix)=MarkovChain([v[i] for i in 1:length(v)], M,
                                      cumsum(M, dims=2))


" Same as `MarkovChain`, but pre-store invariant distribution "
struct MarkovChainExtended{T} <: IndexedMarkovChain{T}
    mc::MarkovChain{T}
    invariant::Array{Float64, 1}
end

MarkovChainExtended(mc::MarkovChain)=
    MarkovChainExtended(mc, markov_invariant(mc))
MarkovChainExtended(v, M::Matrix, tol=1e-12)=
    MarkovChainExtended(MarkovChain(v, M))


supp(mc::MarkovChainExtended)=supp(mc.mc)
transition(mc::MarkovChainExtended)=transition(mc.mc)
transition_cumu(mc::MarkovChainExtended)=transition_cumu(mc.mc)


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

markov_invariant(mc::MarkovChainExtended)=mc.invariant
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
function markov_genpath!(path_prealloc, mc, i0, len, rand_tmp)
    for t in 1:len
        i0 = draw_next(mc, i0, rand_tmp[t])
        path_prealloc[t] = i0
    end
end

function markov_genpath_seeded(mc, i0, len, seed=1234)
    rand_tmp=rand(len)
    path_prealloc = fill(0, len)
    markov_genpath!(path_prealloc, mc, i0, len, rand_tmp)
    return path_prealloc
end


function markov_genpath(mc, i0, len)
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


draw_next(mc::IndexedMarkovChain, i0, r)=begin
    @views findfirstpositive(transition_cumu(mc)[i0, :] .- r)
end

draw_next(mc::IndexedMarkovChain, i0)=draw_next(mc, i0, rand())



indeptensor(mc0::IndexedMarkovChain, mc1::IndexedMarkovChain)=MarkovChainExtended(vec(collect(Base.product(mc0, mc1))),
                                                                                         Base.kron(transition(mc1), transition(mc0)))
