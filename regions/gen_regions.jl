cd(@__DIR__)
import Pkg; Pkg.activate(".")
using NCDatasets

function read_cdo_grid_file(filepath::String)
    grid_dict = Dict{String, Any}()

    for line in eachline(filepath)
        # Remove leading/trailing whitespace
        line = strip(line)

        # Skip empty lines or comment lines (if any)
        isempty(line) && continue
        startswith(line, "#") && continue

        # Split by '=' and strip whitespace
        if occursin("=", line)
            key, value = strip.(split(line, "=", limit=2))

            # Try to parse value as number (Int or Float), otherwise leave as String
            if occursin(r"^-?\d+\.?\d*$", value)
                parsed_value = parse(Float64, value)
                # Convert to Int if it’s actually an integer
                if isinteger(parsed_value)
                    parsed_value = Int(parsed_value)
                end
                grid_dict[key] = parsed_value
            else
                grid_dict[key] = value
            end
        end
    end

    return grid_dict
end

function write_regions(grid_info::Dict{String,Any}, filename::String)
    # Extract basic grid info
    xsize = grid_info["xsize"]
    ysize = grid_info["ysize"]
    xfirst = grid_info["xfirst"]
    xinc = grid_info["xinc"]
    yfirst = grid_info["yfirst"]
    yinc = grid_info["yinc"]
    xname = grid_info["xname"]
    yname = grid_info["yname"]
    xunits = get(grid_info, "xunits", "")
    yunits = get(grid_info, "yunits", "")

    # Construct coordinate vectors
    xc = xfirst:xinc:(xfirst + xinc * (xsize - 1))
    yc = yfirst:yinc:(yfirst + yinc * (ysize - 1))

    # Create NetCDF file
    ds = NCDataset(filename, "c")

    # Define dimensions
    defDim(ds, xname, xsize)
    defDim(ds, yname, ysize)

    # Define coordinate variables
    var_xc = defVar(ds, xname, Float64, (xname,))
    var_yc = defVar(ds, yname, Float64, (yname,))
    var_xc.attrib["units"] = xunits
    var_yc.attrib["units"] = yunits
    var_xc[:] = xc
    var_yc[:] = yc

    # Define the projection variable 'crs'
    crs = defVar(ds, "crs", Int32, ())
    crs.attrib["grid_mapping_name"] = grid_info["grid_mapping_name"]
    crs.attrib["straight_vertical_longitude_from_pole"] = grid_info["straight_vertical_longitude_from_pole"]
    crs.attrib["latitude_of_projection_origin"] = grid_info["latitude_of_projection_origin"]
    crs.attrib["standard_parallel"] = grid_info["standard_parallel"]
    crs.attrib["false_easting"] = grid_info["false_easting"]
    crs.attrib["false_northing"] = grid_info["false_northing"]
    crs.attrib["semi_major_axis"] = grid_info["semi_major_axis"]
    crs.attrib["inverse_flattening"] = grid_info["inverse_flattening"]

    close(ds)
end

function write_2d_variable(
    filename::String,
    varname::String,
    data::AbstractMatrix;
    xdim::String = "xc",
    ydim::String = "yc",
)

    
    

    # Open NetCDF file in append mode
    ds = Dataset(filename, "a")

    # Define the variable (infer data type from matrix)
    dtype = eltype(data)
    var = defVar(ds, varname, dtype, (xdim, ydim))

    # Write the data
    var[:, :] = data

    # Add optional attributes
    # for (k, v) in attribs
    #     var.attrib[k] = v
    # end

    var.attrib["grid_mapping"] = "crs"
    
    close(ds)
end

function gen_mask(xc,yc)
    mask = fill(0.0,(size(xc,1),size(yc,1)))
    return mask
end

grid_info = read_cdo_grid_file("../maps/grid_GRL-PAL-32KM.txt")

# Create a new file
filename = "empty_mask.nc"
write_regions(grid_info, filename)

nc = NCDataset(filename)
xc = nc["xc"][:]
yc = nc["yc"][:]
close(nc)

mask = gen_mask(xc,yc);

write_2d_variable(filename,"mask",mask)