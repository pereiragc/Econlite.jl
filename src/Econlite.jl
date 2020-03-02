module Econlite

using Parameters
using Markdown



# Single good preferences
export IsoElastic, LogUtility, CARA
export util, invutil, dutil, invdutil
include("utility.jl")


# Aggregators
export aggr, partial_aggr, invpartial_aggr, gradient_aggr, gradient_aggr!
export Armington, CobbDouglas, AggLinear
include("aggregators.jl")


# Taxes
export net, tax
export GouveiaStrauss, AdValorem
include("taxes.jl")

end # module
