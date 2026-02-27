module DGGRIDRunnerLib

include("dggrid_params.jl")
import .DGGRIDParams

# re-export DGGRIDParams module
export DGGRIDParams

# export functions
export prep_output_stats
export parse_output_stats
export run_output_stats
export run_dggrid_simple
export prep_generate_grid_whole_earth
export prep_generate_grid_coarse_cells
export prep_generate_grid_clip_region, prep_generate_grid_clip_cells
export grid_gen_convenience!

# z3_coarse_and_convenience_test()
# dryrun()

export prep_output_stats

import DGGRID7_jll
import ArchGDAL

const libs_paths = DGGRID7_jll.LIBPATH_list
const dggrid_exec = DGGRID7_jll.get_dggrid_path()

# maybe update to make flexible between platforms, Apple/macosx, Linux, Windows
if Sys.isapple()
    ENV["DYLD_LIBRARY_PATH"] = join(libs_paths, ":")
elseif Sys.islinux()
    ENV["LD_LIBRARY_PATH"] = join(libs_paths, ":")
elseif Sys.iswindows()
    ENV["PATH"] = join(libs_paths, ";") * ";" * ENV["PATH"]
end

# test version
function dryrun()
    lineinfo = readchomp(`$dggrid_exec -h`)
    println("DGGRID Executable found: $lineinfo")
end

# ---------------------------------------------
# output_stats operation functions
# ---------------------------------------------

function prep_output_stats(dggs_type::String, resolution::Int)
    # grid stats table
    # dggrid_operation OUTPUT_STATS
    # dggs_type ISEA7H
    # dggs_res_spec 20
    # precision 7 (default)
    # verbosity 0 (default)
    params = DGGRIDParams.DGGRIDMetafile()
    DGGRIDParams.add_parameter!(params, "dggrid_operation", "OUTPUT_STATS")
    DGGRIDParams.add_parameter!(params, "dggs_type", dggs_type)
    # DGGRIDParams.add_parameter!(params, "dggs_res_specify_type", "SPECIFIED")
    DGGRIDParams.add_parameter!(params, "dggs_res_spec", resolution)
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        for err in is_ok[2]
            println(" - $err")
        end
        # error("Invalid DGGRID parameters for OUTPUT_STATS")
    else
        println("DGGRID parameters validated successfully.")
    end
    return (true, params)
end

function parse_output_stats(lines::Vector{<:AbstractString})
    # Skip header line (Julia is 1-based array index)
    data_lines = lines[4:end]

    # headers are resolution, num_cells, area_km2, cls_km

    # Dictionary to store results
    # Dict{Int, NamedTuple{(:num_cells, :area_km2, :area_m2, :cls_km, :cls_m), Tuple{Int, Float64, Float64, Float64, Float64}}}
    stats = Dict{Int, NamedTuple}()
    
    for line in data_lines
        # Split by whitespace and filter out empty strings
        parts = filter(!isempty, split(line))
        
        if length(parts) >= 4
            resolution = parse(Int, parts[1])
            num_cells = parse(Int64, replace(parts[2], "," => ""))
            area_km2 = parse(Float64, replace(parts[3], "," => ""))
            cls_km = parse(Float64, replace(parts[4], "," => ""))
            
            # Convert to meters with full precision
            area_m2 = area_km2 * 1_000_000.0  # km² to m²
            cls_m = cls_km * 1_000.0           # km to m
            
            # Store as NamedTuple for easy access
            stats[resolution] = (
                num_cells = num_cells,
                area_km2 = area_km2,
                area_m2 = area_m2,
                cls_km = cls_km,
                cls_m = cls_m
            )
        end
    end
    return stats
end

function run_output_stats(params::DGGRIDParams.DGGRIDMetafile; temp_prefix::String = "")
    # PipeBuffer is a type of IOBuffer
    # create temporary metafile
    metafile = if temp_prefix == ""
            tempname()
        else
            tempname(temp_prefix)
    end

    DGGRIDParams.write_metafile(params, metafile)
    println("Created temporary metafile at: $metafile")

    shell_io = PipeBuffer()
    mycommand = `$dggrid_exec $metafile`
    println("Running DGGRID with command: $mycommand")

    run(mycommand, devnull, shell_io, stderr)
    stats_io_out = readlines(shell_io)
    needed_lines = String[]
    skip_out = true
    for line in stats_io_out
        if occursin("Earth Radius", line)
            skip_out = false
            # push!(needed_lines, line)
        end
        if !skip_out
            push!(needed_lines, line)
            println(line)
        end
    end
    stats = parse_output_stats(needed_lines)
    println("Parsed stats: $stats")
    return stats
end

# optionally provide prefix for temporary files as paramter
function run_dggrid_simple(params::DGGRIDParams.DGGRIDMetafile; temp_prefix::String = "")
    # create temporary metafile
    metafile = if temp_prefix == ""
            tempname()
        else
            tempname(temp_prefix)
    end

    DGGRIDParams.write_metafile(params, metafile)
    println("Created temporary metafile at: $metafile")
    
    # run dggrid with metafile
    cmd = `$dggrid_exec $metafile`
    println("Running DGGRID with command: $cmd")
    run(cmd)
    
    # optionally delete metafile
    rm(metafile; force=true)
    println("Deleted temporary metafile.")
end

# ---------------------------------------------
# GENERATE_GRID functions
# 3 typical use cases:
# - whole grid generation given a dggs_type and resolution
# - coarse grid generation given a dggs_type, resolution, and a list of coarser cells and their resolution
# - clipped region grid generation given a dggs_type, resolution, and a spatial file defining the clip region
# whole and coarse cells don't require creating input files
# clip_region_files requires function to write a spatial file as input
# all functions require defining a path to and then reading in the geodata and cleaning up in- and output files
# ---------------------------------------------

"""
    prep_generate_grid_whole_earth(dggs_type::String,
                                   resolution::Int; 
                                   output_format::String="GPKG",
                                   output_dir::String="",
                                   address_type::String="SEQNUM",
                                   densification::Int=0,
                                   point_output::Bool=false)

Prepare parameters for generating a complete global grid at specified resolution.

# Arguments
- `dggs_type::String`: DGG preset type (e.g., "ISEA7H", "IGEO7")
- `resolution::Int`: Grid resolution level (0-20 for ISEA7H/IGEO7)
- `output_format::String`: GDAL format ("GPKG", "FlatGeobuf", "Arrow", "GeoJSON", "ESRI Shapefile")
- `output_dir::String`: Output directory (empty string uses temp directory)
- `address_type::String`: Cell addressing scheme ("SEQNUM", "Q2DI", "HIERNDX")
- `densification::Int`: Number of additional points per cell edge (default 0)
- `point_output::Bool`: If true, output cell center points instead of cell boundaries (default false)

# Returns
- `(success::Bool, params::DGGRIDMetafile, output_path::String)`: Validation status, metafile parameters, and output file path

# Example
```julia
# Generate cell boundaries
success, params, output_path = prep_generate_grid_whole_earth("ISEA7H", 3, output_format="GPKG")
if success
    run_dggrid_simple(params)
    # Read output from output_path
end

# Generate cell center points
success, params, output_path = prep_generate_grid_whole_earth("ISEA7H", 3, point_output=true)
if success
    run_dggrid_simple(params)
    # Read point output from output_path
end
```
"""
function prep_generate_grid_whole_earth(dggs_type::String,
                                        resolution::Int; 
                                        output_format::String="GPKG",
                                        output_dir::String="",
                                        address_type::String="SEQNUM",
                                        densification::Int=0,
                                        point_output::Bool=false)
    
    # Create metafile with basic parameters
    params = DGGRIDParams.DGGRIDMetafile()
    DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
    DGGRIDParams.add_parameter!(params, "dggs_type", dggs_type)
    DGGRIDParams.add_parameter!(params, "dggs_res_spec", resolution)
    
    # Set clip type for whole earth
    DGGRIDParams.add_parameter!(params, "clip_subset_type", "WHOLE_EARTH")

    # Generate output file path
    output_prefix = if output_dir == ""
        tempname()
    else
        joinpath(output_dir, "grid_$(dggs_type)_res$(resolution)")
    end

    # Determine output file extension
    ext_map = Dict(
        "GPKG" => ".gpkg",
        "FlatGeobuf" => ".fgb",
        "Arrow" => ".arrow",
        "GeoJSON" => ".geojson",
        "ESRI Shapefile" => ".shp"
    )
    output_path = output_prefix * get(ext_map, output_format, ".gpkg")
    
    # Configure output format based on point_output flag
    if point_output
        # Configure point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "point_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "point_output_file_name", output_path)
        
        # Disable cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "NONE")
    else
        # Configure cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "cell_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "cell_output_file_name", output_path)
        
        # Disable point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "NONE")
    end
    
    # Set address type for output labels
    if address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "output_hier_ndx_form", "INT64")
    elseif address_type == "SEQNUM"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "GLOBAL_SEQUENCE")
    elseif address_type == "Q2DI"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "Q2DI")
    end
    
    # Set densification
    if densification > 0
        DGGRIDParams.add_parameter!(params, "densification", densification)
    end
    
    # Validate parameters
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        println("ERROR: Invalid DGGRID parameters for GENERATE_GRID (whole earth):")
        for err in is_ok[2]
            println(" - $err")
        end
        return (false, params, "")
    # else
    #     println("DGGRID parameters validated successfully for whole earth grid generation.")
    end
    
    return (true, params, output_path)
end

"""
    prep_generate_grid_coarse_cells(dggs_type::String,
                                    resolution::Int, 
                                    coarse_res::Int,
                                    coarse_cells::Vector{<:AbstractString};
                                    output_format::String="GPKG",
                                    output_dir::String="",
                                    address_type::String="SEQNUM",
                                    densification::Int=0,
                                    clip_densification::Int=1,
                                    point_output::Bool=false)

Prepare parameters for generating grid cells within specified coarser resolution cells.

# Arguments
- `dggs_type::String`: DGG preset type (e.g., "ISEA7H", "IGEO7")
- `resolution::Int`: Target grid resolution level
- `coarse_res::Int`: Resolution of coarse clipping cells (must be < resolution and > 0)
- `coarse_cells::Vector{<:AbstractString}`: Cell identifiers (sequence numbers or addresses) of coarse cells to use as clipping regions
- `output_format::String`: GDAL format ("GPKG", "FlatGeobuf", "Arrow", "GeoJSON", "ESRI Shapefile")
- `output_dir::String`: Output directory (empty string uses temp directory)
- `address_type::String`: Cell addressing scheme ("SEQNUM", "Q2DI", "HIERNDX")
- `densification::Int`: Number of additional points per cell edge (default 0)
- `clip_densification::Int`: Densification for coarse cell boundaries (default 1)
- `point_output::Bool`: If true, output cell center points instead of cell boundaries (default false)

# Returns
- `(success::Bool, params::DGGRIDMetafile, output_path::String)`: Validation status, metafile parameters, and output file path

# Example
```julia
# Generate resolution 5 cell boundaries within coarse resolution 2 cells 1, 2, and 3
success, params, output_path = prep_generate_grid_coarse_cells("ISEA7H", 5, 2, string.([1, 2, 3]))
if success
    run_dggrid_simple(params)
    # Read output from output_path
end

# Generate resolution 5 cell center points
success, params, output_path = prep_generate_grid_coarse_cells("ISEA7H", 5, 2, string.([1, 2, 3]), point_output=true)
if success
    run_dggrid_simple(params)
    # Read point output from output_path
end
```
"""
function prep_generate_grid_coarse_cells(dggs_type::String, 
                                         resolution::Int, 
                                         coarse_res::Int,
                                         coarse_cells::Vector{<:AbstractString};
                                         output_format::String="GPKG",
                                         output_dir::String="",
                                         address_type::String="SEQNUM",
                                         densification::Int=0,
                                         clip_densification::Int=1,
                                         point_output::Bool=false)
    
    # Validate coarse resolution
    if coarse_res <= 0 || coarse_res >= resolution
        println("ERROR: coarse_res must be > 0 and < resolution ($resolution)")
        return (false, DGGRIDParams.DGGRIDMetafile(), "")
    end
    
    if isempty(coarse_cells)
        println("ERROR: coarse_cells cannot be empty")
        return (false, DGGRIDParams.DGGRIDMetafile(), "")
    end
    
    # Create metafile with basic parameters
    params = DGGRIDParams.DGGRIDMetafile()
    DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
    DGGRIDParams.add_parameter!(params, "dggs_type", dggs_type)
    DGGRIDParams.add_parameter!(params, "dggs_res_spec", resolution)
    
    # Set clip type for coarse cells
    DGGRIDParams.add_parameter!(params, "clip_subset_type", "COARSE_CELLS")
    DGGRIDParams.add_parameter!(params, "clip_cell_res", coarse_res)
    
    # Convert seqnums to space-delimited string
    # we can have both clip_cell_seqnums when address_type is SEQNUM and clip_cell_addresses when address_type is HIERNDX or Q2DI
    # in the mean time the data vector we pass in can alsways be called coarse_cells
    # seqnums are integers but addresses are strings
    addresses_str = join(coarse_cells, " ")
    if address_type == "SEQNUM"
        DGGRIDParams.add_parameter!(params, "clip_cell_seqnums", addresses_str)
    elseif address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "clip_cell_addresses", addresses_str)
        DGGRIDParams.add_parameter!(params, "input_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "INT64")
    elseif address_type == "Q2DI"
        # we don't know if that is defined behaviour
    end
    DGGRIDParams.add_parameter!(params, "clip_cell_densification", clip_densification)
    
    # Generate output file path
    output_prefix = if output_dir == ""
        tempname()
    else
        joinpath(output_dir, "grid_$(dggs_type)_res$(resolution)_coarse$(coarse_res)")
    end
    
    # Determine output file extension
    ext_map = Dict(
        "GPKG" => ".gpkg",
        "FlatGeobuf" => ".fgb",
        "Arrow" => ".arrow",
        "GeoJSON" => ".geojson",
        "ESRI Shapefile" => ".shp"
    )
    output_path = output_prefix * get(ext_map, output_format, ".gpkg")
    
    # Configure output format based on point_output flag
    if point_output
        # Configure point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "point_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "point_output_file_name", output_path)
        
        # Disable cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "NONE")
    else
        # Configure cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "cell_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "cell_output_file_name", output_path)
        
        # Disable point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "NONE")
    end
    
    # Set address type for output labels
    if address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "output_hier_ndx_form", "INT64")
    elseif address_type == "SEQNUM"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "GLOBAL_SEQUENCE")
    elseif address_type == "Q2DI"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "Q2DI")
    end
    
    # Set densification
    if densification > 0
        DGGRIDParams.add_parameter!(params, "densification", densification)
    end
    
    # Validate parameters
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        println("ERROR: Invalid DGGRID parameters for GENERATE_GRID (coarse cells):")
        for err in is_ok[2]
            println(" - $err")
        end
        return (false, params, "")
    # else
    #     println("DGGRID parameters validated successfully for coarse cells grid generation.")
    end
    
    return (true, params, output_path)
end

"""
    prep_generate_grid_clip_region(dggs_type::String,
                                   resolution::Int, 
                                   clip_file::String;
                                   output_format::String="GPKG",
                                   output_dir::String="",
                                   address_type::String="SEQNUM",
                                   densification::Int=0,
                                   geodetic_densify::Float64=0.0,
                                   use_holes::Bool=false,
                                   point_output::Bool=false)

Prepare parameters for generating grid cells clipped to a spatial region defined by a file.

# Arguments
- `dggs_type::String`: DGG preset type (e.g., "ISEA7H", "IGEO7")
- `resolution::Int`: Grid resolution level
- `clip_file::String`: Path to clipping polygon file (Shapefile, GeoPackage, GeoJSON, etc.)
- `output_format::String`: GDAL format ("GPKG", "FlatGeobuf", "Arrow", "GeoJSON", "ESRI Shapefile")
- `output_dir::String`: Output directory (empty string uses temp directory)
- `address_type::String`: Cell addressing scheme ("SEQNUM", "Q2DI", "HIERNDX")
- `densification::Int`: Number of additional points per cell edge (default 0)
- `geodetic_densify::Float64`: Max arc length in degrees for clipping polygon segments (0.0 = no densification)
- `use_holes::Bool`: Handle holes in clipping polygons (default false)
- `point_output::Bool`: If true, output cell center points instead of cell boundaries (default false)

# Returns
- `(success::Bool, params::DGGRIDMetafile, output_path::String)`: Validation status, metafile parameters, and output file path

# Example
```julia
# Generate cell boundaries clipped to a region defined in a shapefile
success, params, output_path = prep_generate_grid_clip_region(
    "ISEA7H", 5, "/path/to/region.shp",
    geodetic_densify=1.0
)
if success
    run_dggrid_simple(params)
    # Read output from output_path
end

# Generate cell center points clipped to a region
success, params, output_path = prep_generate_grid_clip_region(
    "ISEA7H", 5, "/path/to/region.shp",
    point_output=true
)
if success
    run_dggrid_simple(params)
    # Read point output from output_path
end
```
"""
function prep_generate_grid_clip_region(dggs_type::String, 
                                        resolution::Int, 
                                        clip_file::String;
                                        output_format::String="GPKG",
                                        output_dir::String="",
                                        address_type::String="SEQNUM",
                                        densification::Int=0,
                                        geodetic_densify::Float64=0.0,
                                        use_holes::Bool=false,
                                        point_output::Bool=false)
    
    # Validate clip file exists
    if !isfile(clip_file)
        println("ERROR: Clip file does not exist: $clip_file")
        return (false, DGGRIDParams.DGGRIDMetafile(), "")
    end
    
    # Create metafile with basic parameters
    params = DGGRIDParams.DGGRIDMetafile()
    DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
    DGGRIDParams.add_parameter!(params, "dggs_type", dggs_type)
    DGGRIDParams.add_parameter!(params, "dggs_res_spec", resolution)
    
    # Set clip type for spatial file
    DGGRIDParams.add_parameter!(params, "clip_subset_type", "GDAL")
    DGGRIDParams.add_parameter!(params, "clip_region_files", clip_file)
    
    # Set geodetic densification if specified
    if geodetic_densify > 0.0
        DGGRIDParams.add_parameter!(params, "geodetic_densify", geodetic_densify)
    end
    
    # Set hole handling
    if use_holes
        DGGRIDParams.add_parameter!(params, "clip_using_holes", true)
    end
    
    # Generate output file path
    output_prefix = if output_dir == ""
        tempname()
    else
        clip_basename = splitext(basename(clip_file))[1]
        joinpath(output_dir, "grid_$(dggs_type)_res$(resolution)_$(clip_basename)")
    end

    # Determine output file extension
    ext_map = Dict(
        "GPKG" => ".gpkg",
        "FlatGeobuf" => ".fgb",
        "Arrow" => ".arrow",
        "GeoJSON" => ".geojson",
        "ESRI Shapefile" => ".shp"
    )
    output_path = output_prefix * get(ext_map, output_format, ".gpkg")
    
    # Configure output format based on point_output flag
    if point_output
        # Configure point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "point_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "point_output_file_name", output_path)
        
        # Disable cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "NONE")
    else
        # Configure cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "cell_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "cell_output_file_name", output_path)
        
        # Disable point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "NONE")
    end
    
    # Set address type for output labels
    if address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "output_hier_ndx_form", "INT64")
    elseif address_type == "SEQNUM"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "GLOBAL_SEQUENCE")
    elseif address_type == "Q2DI"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "Q2DI")
    end
    
    # Set densification
    if densification > 0
        DGGRIDParams.add_parameter!(params, "densification", densification)
    end
    
    # Validate parameters
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        println("ERROR: Invalid DGGRID parameters for GENERATE_GRID (clip region):")
        for err in is_ok[2]
            println(" - $err")
        end
        return (false, params, "")
    # else
    #     println("DGGRID parameters validated successfully for clipped region grid generation.")
    end
    
    return (true, params, output_path)
end


"""
    prep_generate_grid_clip_region(dggs_type::String,
                                   resolution::Int, 
                                   clip_file::String;
                                   output_format::String="GPKG",
                                   output_dir::String="",
                                   address_type::String="SEQNUM",
                                   densification::Int=0,
                                   geodetic_densify::Float64=0.0,
                                   use_holes::Bool=false,
                                   point_output::Bool=false)

Prepare parameters for generating grid cells clipped to a spatial region defined by a file.

# Arguments
- `dggs_type::String`: DGG preset type (e.g., "ISEA7H", "IGEO7")
- `resolution::Int`: Grid resolution level
- `clip_file::String`: Path to clipping text file containing a list of cell ids within the same resolution
- `output_format::String`: GDAL format ("GPKG", "FlatGeobuf", "Arrow", "GeoJSON", "ESRI Shapefile")
- `output_dir::String`: Output directory (empty string uses temp directory)
- `output_address_type::String`: Cell addressing scheme for output ("SEQNUM", "Q2DI", "HIERNDX")
- `input_address_type::String`: Cell addressing scheme for input ("SEQNUM", "Q2DI", "HIERNDX")
- `densification::Int`: Number of additional points per cell edge (default 0)
- `geodetic_densify::Float64`: Max arc length in degrees for clipping polygon segments (0.0 = no densification)
- `point_output::Bool`: If true, output cell center points instead of cell boundaries (default false)

# Returns
- `(success::Bool, params::DGGRIDMetafile, output_path::String)`: Validation status, metafile parameters, and output file path

# Example
```julia
# Generate cells  clipped to a region defined by a list of cell IDs in a text file
success, params, output_path = prep_generate_grid_clip_cells(
    "ISEA7H", 5, "/path/to/cellids.txt",
    geodetic_densify=1.0
)
if success
    run_dggrid_simple(params)
    # Read output from output_path
end

# Generate cell center points clipped to a region
success, params, output_path = prep_generate_grid_clip_cells(
    "ISEA7H", 5, "/path/to/cellids.txt",
    point_output=true
)
if success
    run_dggrid_simple(params)
    # Read point output from output_path
end
```
"""
function prep_generate_grid_clip_cells( dggs_type::String,
                                        resolution::Int, 
                                        clip_file::String;
                                        output_format::String="GPKG",
                                        output_dir::String="",
                                        output_address_type::String="SEQNUM",
                                        input_address_type::String="SEQNUM",
                                        densification::Int=0,
                                        geodetic_densify::Float64=0.0,
                                        point_output::Bool=false)
    
    # Validate clip file exists
    if !isfile(clip_file)
        println("ERROR: Clip file does not exist: $clip_file")
        return (false, DGGRIDParams.DGGRIDMetafile(), "")
    end
    
    # Create metafile with basic parameters
    params = DGGRIDParams.DGGRIDMetafile()
    DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
    DGGRIDParams.add_parameter!(params, "dggs_type", dggs_type)
    DGGRIDParams.add_parameter!(params, "dggs_res_spec", resolution)
    
    # Set clip type for spatial file
    DGGRIDParams.add_parameter!(params, "clip_subset_type", "INPUT_ADDRESS_TYPE")
    DGGRIDParams.add_parameter!(params, "clip_region_files", clip_file)
    DGGRIDParams.add_parameter!(params, "input_address_type", input_address_type)
    
    if input_address_type == "SEQNUM"
        # pass
        
    elseif input_address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "input_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "INT64")
    elseif input_address_type == "Q2DI"
        # we don't know if that is defined behaviour
    end

    # Set geodetic densification if specified
    if geodetic_densify > 0.0
        DGGRIDParams.add_parameter!(params, "geodetic_densify", geodetic_densify)
    end
    
    # Generate output file path
    output_prefix = if output_dir == ""
        tempname()
    else
        clip_basename = splitext(basename(clip_file))[1]
        joinpath(output_dir, "grid_$(dggs_type)_res$(resolution)_$(clip_basename)")
    end

    # Determine output file extension
    ext_map = Dict(
        "GPKG" => ".gpkg",
        "FlatGeobuf" => ".fgb",
        "Arrow" => ".arrow",
        "GeoJSON" => ".geojson",
        "ESRI Shapefile" => ".shp"
    )
    output_path = output_prefix * get(ext_map, output_format, ".gpkg")
    
    # Configure output format based on point_output flag
    if point_output
        # Configure point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "point_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "point_output_file_name", output_path)
        
        # Disable cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "NONE")
    else
        # Configure cell output
        DGGRIDParams.add_parameter!(params, "cell_output_type", "GDAL")
        DGGRIDParams.add_parameter!(params, "cell_output_gdal_format", output_format)
        DGGRIDParams.add_parameter!(params, "cell_output_file_name", output_path)
        
        # Disable point output
        DGGRIDParams.add_parameter!(params, "point_output_type", "NONE")
    end
    
    # Set address type for output labels
    if output_address_type == "HIERNDX"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "HIERNDX")
        if (dggs_type in ["IGEO7", "ISEA7H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z7")
        elseif (dggs_type in ["ISEA3H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "Z3")
        elseif (dggs_type in ["ISEA4H"])
            DGGRIDParams.add_parameter!(params, "output_hier_ndx_system", "ZORDER")
        end
        DGGRIDParams.add_parameter!(params, "output_hier_ndx_form", "INT64")
    elseif output_address_type == "SEQNUM"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "GLOBAL_SEQUENCE")
    elseif address_type == "Q2DI"
        DGGRIDParams.add_parameter!(params, "output_cell_label_type", "OUTPUT_ADDRESS_TYPE")
        DGGRIDParams.add_parameter!(params, "output_address_type", "Q2DI")
    end
    
    # Set densification
    if densification > 0
        DGGRIDParams.add_parameter!(params, "densification", densification)
    end
    
    # Validate parameters
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        println("ERROR: Invalid DGGRID parameters for GENERATE_GRID (clip region):")
        for err in is_ok[2]
            println(" - $err")
        end
        return (false, params, "")
    # else
    #     println("DGGRID parameters validated successfully for clipped region grid generation.")
    end
    
    return (true, params, output_path)
end


function grid_gen_convenience!(params::DGGRIDParams.DGGRIDMetafile; longitude_wrap_mode="UNWRAP_EAST", dggs_vert0_lon=11.2)
    # 
    # if dggs_type IGEO7 or ISEA7H and HIERNDX address type, switch to digit_string
    # input_hier_ndx_form digit_string
    # if dggs_type IGEO7 or ISEA7H or ISEA3H or ISEA4H and HIERNDX address type, set authalic prep params
    # dggs_vert0_lon 11.2
    # dggs_vert0_lat 58.2825
    # generally we also like to set
    # longitude_wrap_mode UNWRAP_EAST
    # if using shapefile we need to add
    # shapefile_id_field_length 22
    dggs_type = DGGRIDParams.get_parameter(params, "dggs_type")
    address_type = DGGRIDParams.get_parameter(params, "output_address_type")
    if (dggs_type in ["IGEO7", "ISEA7H", "ISEA3H", "ISEA4H"]) && (address_type == "HIERNDX")
        # DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "digit_string")
        DGGRIDParams.add_parameter!(params, "output_hier_ndx_form", "digit_string")

        operation = DGGRIDParams.get_parameter(params, "dggrid_operation")
        clip_type = DGGRIDParams.get_parameter(params, "clip_subset_type")
        if (operation == "GENERATE_GRID") && ( clip_type == "COARSE_CELLS")
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "digit_string")
        elseif (operation == "GENERATE_GRID") && ( clip_type == "INPUT_ADDRESS_TYPE")
            if DGGRIDParams.get_parameter(params, "input_address_type") == "HIERNDX"
                DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "digit_string")
            end
        end
        if (operation == "TRANSFORM_ADDRESSES")
            DGGRIDParams.add_parameter!(params, "input_hier_ndx_form", "digit_string")
        end
    end

    DGGRIDParams.add_parameter!(params, "dggs_vert0_lon", dggs_vert0_lon)
    DGGRIDParams.add_parameter!(params, "dggs_vert0_lat", 58.2825)
    DGGRIDParams.add_parameter!(params, "longitude_wrap_mode", longitude_wrap_mode)

    cell_output_type = DGGRIDParams.get_parameter(params, "cell_output_type")
    if cell_output_type == "SHAPEFILE"
        DGGRIDParams.add_parameter!(params, "shapefile_id_field_length", 22)
    end

    # Validate parameters
    is_ok = DGGRIDParams.validate_metafile(params)
    if !is_ok[1]
        println("ERROR: Invalid DGGRID parameters for grid_gen_convenience update!:")
        for err in is_ok[2]
            println(" - $err")
        end
        return false
    # else
    #     println("DGGRID parameters validated successfully for grid_gen_convenience update!")
    end
    return true
end

# ---------------------------------------------
# dev test section
# ---------------------------------------------

function z3_coarse_and_convenience_test()
    dggs_type = "ISEA3H"
    resolution = 4
    coarse_res = 2
    coarse_cells = ["0330"]
    output_format = "FlatGeobuf" # GPKG
    output_dir = ""
    address_type = "HIERNDX"
    densification = 0
    point_output = false

    success, params, output_path = prep_generate_grid_coarse_cells(
        dggs_type, resolution, coarse_res, coarse_cells;
        output_format=output_format,
        output_dir=output_dir,
        address_type=address_type,
        densification=densification,
        point_output=point_output
    )

    validated = grid_gen_convenience!(params)

    # just print out results
    println("Generated parameters... success: $success, validated: $validated, output_path: $output_path")
    for (key, value) in params.params
        println(" - $key : $value")
    end
    
end

end # module

# ---------------------------------------------
# if name is main equivalent in Julia
# ---------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    import .DGGRIDRunnerLib
    
    DGGRIDRunnerLib.dryrun()

    DGGRIDRunnerLib.z3_coarse_and_convenience_test()
end

# dggs_type IGEO7
# dggs_vert0_lon 11.2
# dggs_vert0_lat 58.2825
# input_hier_ndx_form digit_string
# longitude_wrap_mode UNWRAP_EAST
# when using shapefile we need
# shapefile_id_field_length 22
