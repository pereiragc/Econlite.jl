module Econlite

using Parameters
using Markdown



# Single good preferences
export IsoElastic, LogUtility, CARA
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

end # module
