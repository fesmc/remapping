cd(@__DIR__)
import Pkg; Pkg.activate(".")
using NCDatasets

# Several useful functions
include("grids.jl")

function define_my_grid(domain,grid_name,griddes)

    #grid_info = read_cdo_grid_file("../maps/grid_GRL-PAL-32KM.txt","GRL-PAL-32KM","Greenland")
    grid_info = read_cdo_grid_file(griddes,grid_name,domain)

    # Create a new file
    filename = "maps/grid_$(grid_info["grid_name"]).nc"
    xc, yc = define_grid_nc(grid_info, filename)

    lat2D, lon2D = projected_to_latlon(grid_info, xc, yc)
    area = cell_areas(lat2D, lon2D)

    write_2d_variable(filename,"lat2D",lat2D)
    write_2d_variable(filename,"lon2D",lon2D)
    write_2d_variable(filename,"area",area)

    println("File written: $filename")
end

if abspath(PROGRAM_FILE) == @__FILE__
    domain = ARGS[1]
    grid_name = ARGS[2]
    griddes = "maps/grid_$grid_name.txt"

    println("domain = $domain")
    println("grid_name = $grid_name")
    println("griddes = $griddes")

    # Define the grid
    define_my_grid(domain,grid_name,griddes)
end
