cd(@__DIR__)
import Pkg; Pkg.activate(".")
using NCDatasets

# Several useful functions
include("grids.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    domain = ARGS[1]
    grid_name = ARGS[2]

    println("domain, grid_name = $domain, $grid_name")

    # Define the grid
    grid = generate_grid(domain,grid_name)

    # Write to output file
    filename = joinpath("../maps","grid_$(grid_name).nc")
    grid_write_nc(grid,filename)

end
