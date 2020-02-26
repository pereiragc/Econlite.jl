module Econlite

using Parameters
using Markdown



# Single good preferences
export IsoElastic, LogUtility, CARA
export util, invutil, dutil, invdutil
include("utility.jl")


export aggr, partial_aggr, invpartial_aggr, Armington
include("aggregators.jl")

end # module
