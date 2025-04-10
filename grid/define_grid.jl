cd(@__DIR__)
import Pkg; Pkg.activate(".")
using NCDatasets

# Several useful functions
include("grids.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    domain = ARGS[1]
    grid_name = ARGS[2]

    println("domain = $domain")
    println("grid_name = $grid_name")

    # Define the grid
    define_my_grid(domain,grid_name)
end
