module Econlite

using Parameters
using Markdown
using Unrolled
using LinearAlgebra
using Statistics
import Random.seed!



# Single good preferences
export IsoElastic, LogUtility, CARA, RiskNeutral
export util, invutil, dutil, invdutil
include("utility.jl")


# Aggregators
export aggr, invaggr, partial_aggr, invpartial_aggr, gradient_aggr, gradient_aggr!
export Armington, CobbDouglas, AggLinear
include("aggregators.jl")


# Taxes
export net, gross, tax
export GouveiaStrauss, AdValorem, NullTax
include("taxes.jl")

# Maximum
export generalized_argmax!
export SumNormalize, SumRaw
include("maximum.jl")


# Generate grid
export gridpoints
include("gridpoints.jl")


# Iteration
export SolutionObject, IterationSpec, DisplayMod, NoDisplay
export within_tolerance, contract!
include("iteration.jl")


# Markov chains
export IndexedMarkovChain, MarkovChain, MarkovChainExtended
export expect_markov, transition, markov_invariant, draw_next,
    markov_genpath, markov_genpath_seeded, markov_genpath!, mean,
    indeptensor, makemean1!
include("markovchain.jl")

export rouwenhorst, rouwenhorst_mc
include("rouwenhorst.jl")


end
