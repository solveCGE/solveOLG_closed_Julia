#------------------------------------------------------------------------------#
# "solveOLG_closed_Julia"
# 
# Solves an AK-OLG-Model, closed economy, with income effects in Julia
# Philip Schuster, August, 2023
#
# run.jl: runs everything (first run_prep.jl then run_sim.jl)
#------------------------------------------------------------------------------#

include("run_prep.jl")
include("run_sim.jl") # run this line again to measure time without precompilation

