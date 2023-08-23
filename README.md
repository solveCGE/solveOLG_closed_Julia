# solveOLG closed economy in Julia
Solves a simple AK-OLG-model for closed economy in Julia

## About
Shows how to solve a simple deterministic overlapping-generations model (OLG) of Auerbauch-Kotlikoff type, solving for the transition path between two steady-states. The code used global variables. Cohorts can be computed in parallel.

A model description can be found here: <https://github.com/solveCGE/solveOLG_doc>.

## How to run
Parameters can be set in `calib.jl`. Policy shocks are defined in `run_sim.jl` (or just uncomment some of the predefined exemplary shocks). The model is then solved by just running `run.jl`:

```bash
julia -t4 run.jl
```
from the terminal or in the Julia REPL (`julia -t4`)

```julia
include("run.jl")
```
The number of threads has to be specified when calling Julia using the `-t` flag (i.e. the code above uses 4 threads). Running `run_sim.jl` twice allows taking time without precompilation.

## Author
Philip Schuster
