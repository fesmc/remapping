# remapping
Repository to share scripts and methods to remap original datasets onto regional ice-sheet grids.

## First steps

To remap data to a new grid, it is necessary to define the grid. If this has not been done externally, then the following steps will typically need to be performed before proceeding with processing and remapping a particular dataset.

### Step 1: define the grid description file

Currently remapping is performed via `cdo`. Therefore we need a grid description file that can be read by `cdo`. This file contains projection information, resolution between grid points and axis locations. This file should be stored in the `maps` folder. See examples there for how to define it. Importantly, the grid description filename should follow this approach `grid_GRIDNAME.txt`, where `GRIDNAME` is the name of the current grid of interest. This is typically named with an acronym for the region (like `GRL` for Greenland) and the resolution (like `16KM`).

The grid description file can also be produced via `cdo` directly, if a NetCDF file is available with this grid defined. In that case, run the command `cdo griddes FILE.nc > maps/grid_GRIDNAME.txt`.

### Step 2: generate a grid file

Once a grid description file has been made, we can create an actual NetCDF file with the grid information defined. This will include the projection information, Cartesian axes (`xc`, `yc`), lat/lon information (`lat2D`, `lon2D`) and the grid cell area (`area`).

This can be created by going into the `grid` subfolder and running the script `julia define_grid.jl DOMAIN GRIDNAME` with the variables modified appropriately for your current grid. So to get a grid file for the Greenland domain on the grid GRL-32KM run:

```bash
julia grid/define_grid.jl Greenland GRL-32KM
```

This will produce a file in the output folder named "out/GRL-32KM_grid.nc", as well as a file with the standard grid name in the maps folder "maps/grid_GRL-32KM.nc". The latter is convenient for remapping via cdo commands later, so that the target nc file is always together with the grid description file.

### Step 3: generate grid remapping files

Once a grid description file is available and a NetCDF grid file, now we can use `cdo` to calculate conservative remapping weights between our target grid and a source grid.

```bash
TO DO
```
