cd(@__DIR__)
import Pkg; Pkg.activate(".")

using Shapefile
using CairoMakie
using GeoInterface
using GeometryBasics
using GeometryOps

using NCDatasets

# Several useful functions
include("../grid/grids.jl")

# Load shapefile (change path to your own shapefile)
#table = Shapefile.Table("marineregions/World_EEZ_v12_20231025/eez_boundaries_v12.shp")
table = Shapefile.Table("marineregions/World_EEZ_v12_20231025/eez_v12.shp")

# Find regions of interest
kk = []
for (k,reg) in enumerate(table.GEONAME)
    if occursin("Greenland",reg)
        push!(kk,k)
    end
end


# Extract geometries
geoms = GeoInterface.geometry.(table)

begin
    # Start the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect())

    #for geom in geoms
    #    lines!(ax,geom, color=:black)
    #end

    k = kk[1]
    geom = geoms[k]
    lines!(ax,geom, color=:red)

    pt = Point2f(-45.0, 70.0)
    scatter!(ax,pt,color=:green)
end
fig


# Test point (longitude, latitude)
k = kk[1]
geom = geoms[k]
pt = Point2f(-45.0, 70.0)
GeometryOps.contains(geom, pt)



# Assuming `geom` is your Polygon object
coords = GeoInterface.coordinates(geom)  # Extract coordinates

# If it's a multipolygon or polygon with holes, you might have a nested structure
# For a multipolygon, each polygon can be isolated like this:

if eltype(coords) <: Tuple  # Simple polygon with one outer boundary
    rings = [coords]  # Just one region
else
    rings = coords  # Multipolygon: List of polygons
end

# Now, `rings` contains individual polygon parts (rings)

# Example: Plot the first ring of the first polygon
ring = rings[1][1]  # Outer boundary of the first polygon (or the first multipolygon part)
xs = first.(ring)
ys = last.(ring)

# Plot the isolated region
fig = Figure()
ax = Axis(fig[1, 1]; aspect=DataAspect())
lines!(ax, xs, ys, color=:black)
fig


pt = Point2f(-45.0, 70.0)
polygon = GeometryOps.Polygon(Point2f.(ring))
GeometryOps.contains(polygon, pt)

### Saving a region ###

# Define the grid
grid = generate_grid("Greenland","GRL-PAL-16KM")

nx, ny = grid["nx"], grid["ny"]
mask = fill(0.0,nx,ny)

polygon = GeometryOps.Polygon(Point2f.(ring))

for i in 1:nx, j in 1:ny
    pt = Point2f(grid["lon2D"][i,j],grid["lat2D"][i,j])
    if GeometryOps.contains(polygon, pt)
        mask[i,j] = 1.0
    end
end

# Write a new file with the mask
filename_out = "../out/$(grid["grid_name"])_REGIONS.nc"
grid_write_nc(grid,filename_out)
write_2d_variable(filename_out,"mask",mask)

begin
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect())
    lines!(ax, polygon, color=:black)
    kk = findall(mask .> 0.0)
    scatter!(ax,grid["lon2D"][kk],grid["lat2D"][kk])
    fig
end

begin
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect())
    heatmap!(ax,grid["xc"],grid["yc"],mask)
    fig
end

