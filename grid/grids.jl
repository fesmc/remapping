cd(@__DIR__)
using NCDatasets
using Proj
using GeographicLib
using LinearAlgebra

function read_cdo_grid_file(filepath::String,grid_name::String,domain::String)
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

    # Define projection string based on grid_dict
    grid_dict["proj_str"] = "+proj=stere " *
               "+lat_0=$(grid_dict["latitude_of_projection_origin"]) " *
               "+lat_ts=$(grid_dict["standard_parallel"]) " *
               "+lon_0=$(grid_dict["straight_vertical_longitude_from_pole"]) " *
               "+k=1 +x_0=$(grid_dict["false_easting"] * 1000) +y_0=$(grid_dict["false_northing"] * 1000) " * # Convert to meters
               "+a=$(grid_dict["semi_major_axis"]) +rf=$(grid_dict["inverse_flattening"]) +units=m +no_defs" # Units set to meters

    # Add domain and grid_name too
    grid_dict["domain"] = domain
    grid_dict["grid_name"] = grid_name

    return grid_dict
end

function define_grid_nc(grid_info::Dict{String,Any}, filename::String)
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

    return collect(xc), collect(yc)
end

function projected_to_latlon(grid_info::Dict, xc::Vector{<:Real}, yc::Vector{<:Real})
    
    # Define the source and destination CRS
    source_crs = grid_info["proj_str"]
    dest_crs = "EPSG:4326" # Latitude/longitude

    # Create the transformation object
    transformation = Proj.Transformation(source_crs, dest_crs; always_xy=true)

    # Create meshgrid of x/y (convert km to m)
    X = repeat(xc' .* 1000, length(yc), 1)  # (ysize, xsize) in meters
    Y = repeat(yc  .* 1000, 1, length(xc))  # (ysize, xsize) in meters

    # Apply the transformation
    lonlat = transformation.(X,Y)

    # Extract the longitude and latitude
    lon2D = [coord[1] for coord in lonlat]
    lat2D = [coord[2] for coord in lonlat]

    # Reshape to 2D arrays (nx, ny)
    lon2D = reshape(lon2D, length(yc), length(xc))
    lat2D = reshape(lat2D, length(yc), length(xc))

    lon2D = Matrix(lon2D')
    lat2D = Matrix(lat2D')

    return lat2D, lon2D
end

function geodesic_polygon_area(p)

    lons = first.(p);
    lats = last.(p);

    pg = Polygon(lons,lats)

    n, perimeter, area = properties(pg)

    return area
end

function extend2D(var)
    nx, ny = size(var)
    vare = fill(NaN,nx+2,ny+2)
    vare[2:nx+1,2:ny+1] .= var

    vare[1,2:ny+1]    = var[1,:] .- (var[2,:].-var[1,:])
    vare[nx+2,2:ny+1] = var[nx,:] .+ (var[nx,:].-var[nx-1,:])
    vare[:,1]    = vare[:,2] .- (vare[:,3].-vare[:,2])
    vare[:,ny+2] = vare[:,ny+1] .+ (vare[:,ny+1].-vare[:,ny])
    
    # Fill the corners
    vare[1, 1]       = var[1, 1] .- (var[2, 1] .- var[1, 1]) .- (var[1, 2] .- var[1, 1])
    vare[1, ny+2]    = var[1, ny] .- (var[2, ny] .- var[1, ny]) .+ (var[1, ny] .- var[1, ny-1])
    vare[nx+2, 1]    = var[nx, 1] .+ (var[nx, 1] .- var[nx-1, 1]) .- (var[nx, 2] .- var[nx, 1])
    vare[nx+2, ny+2] = var[nx, ny] .+ (var[nx, ny] .- var[nx-1, ny]) .+ (var[nx, ny] .- var[nx, ny-1])

    return vare
end

function cell_areas(lat2D::Matrix{Float64}, lon2D::Matrix{Float64})
    nx, ny = size(lat2D)
    area = fill(0.0, nx, ny)

    lon2De = extend2D(lon2D)
    lat2De = extend2D(lat2D)
    
    for j in 2:ny+1, i in 2:nx+1
        
        im1 = i-1
        ip1 = i+1
        jm1 = j-1
        jp1 = j+1

        # Bottom-left
        lon1 = 0.25*(lon2De[i,j]+lon2De[im1,j]+lon2De[im1,jm1]+lon2De[i,jm1])
        lat1 = 0.25*(lat2De[i,j]+lat2De[im1,j]+lat2De[im1,jm1]+lat2De[i,jm1])

        # Bottom-right
        lon2 = 0.25*(lon2De[i,j]+lon2De[ip1,j]+lon2De[ip1,jm1]+lon2De[i,jm1])
        lat2 = 0.25*(lat2De[i,j]+lat2De[ip1,j]+lat2De[ip1,jm1]+lat2De[i,jm1])

        # Top-right
        lon3 = 0.25*(lon2De[i,j]+lon2De[ip1,j]+lon2De[ip1,jp1]+lon2De[i,jp1])
        lat3 = 0.25*(lat2De[i,j]+lat2De[ip1,j]+lat2De[ip1,jp1]+lat2De[i,jp1])

        # Top-left
        lon4 = 0.25*(lon2De[i,j]+lon2De[im1,j]+lon2De[im1,jp1]+lon2De[i,jp1])
        lat4 = 0.25*(lat2De[i,j]+lat2De[im1,j]+lat2De[im1,jp1]+lat2De[i,jp1])

        vertices = [
            (lon1,lat1),
            (lon2,lat2),
            (lon3,lat3),
            (lon4,lat4),
            (lon1,lat1)     # Close the polygon
        ]

        area[i-1, j-1] = geodesic_polygon_area(vertices)

    end

    # Recalculate area around the border using extended values
    area = extend2D(area[2:nx-1,2:ny-1])

    return area
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