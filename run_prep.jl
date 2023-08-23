using Dates
using Statistics
using NLsolve
using Printf
using Distributed
using Plots

# load some handy auxiliary functions
include("loadfunctions.jl")

# load main functions
include("algo.jl")
include("hh.jl")
include("firm.jl")