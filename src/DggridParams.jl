"""
DGGRID Metafile Parameter Library

This module provides a comprehensive parameter library for DGGRID metafiles.
It defines all parameters with their types, allowed values, defaults, and validation.

Based on DGGRID documentation and focused on:
- Operations: GENERATE_GRID, TRANSFORM_POINTS, OUTPUT_STATS
- Address forms: SEQNUM, Q2DI, HIERNDX (Z7, Z7_STRING, INT64, DIGIT_STRING)
- DGG presets: ISEA7H, IGEO7
"""
module DGGRIDParams

export DGGRIDMetafile, Parameter, ParameterType
export create_metafile, add_parameter!, get_parameter, validate_metafile, remove_parameter!
export write_metafile, read_metafile, get_parameter_info, list_parameters

# Parameter types
# STRINGLIST is a STRING with space-delimited values
@enum ParameterType begin
    BOOLEAN
    INTEGER
    DOUBLE
    STRING
    STRINGLIST
    CHOICE
end

const max_resolution = Dict(
    "CUSTOM" => 35,
    "ISEA7H" => 20,
    "IGEO7"  => 20,
    "ISEA3H" => 34,
    "ISEA4H" => 29,
    "ISEA43H" => 22,
    "ISEA4T" => 29,
    "ISEA4D" => 29
)

"""
    Parameter

Represents a single DGGRID metafile parameter with its metadata.
"""
struct Parameter
    name::String
    type::ParameterType
    default::Any
    allowed_values::Union{Nothing, Vector{String}}  # For CHOICE types
    description::String
    used_when::String  # Conditions when this parameter is used
end

"""
    DGGRIDMetafile

Container for DGGRID metafile parameters.
"""
mutable struct DGGRIDMetafile
    params::Dict{String, Any}
    
    DGGRIDMetafile() = new(Dict{String, Any}())
end

"""
    add_parameter!(metafile::DGGRIDMetafile, name::String, value::Any)

Add or update a parameter in the metafile.
"""
function add_parameter!(metafile::DGGRIDMetafile, name::String, value::Any)
    metafile.params[name] = value
    return metafile
end

"""
    remove_parameter!(metafile::DGGRIDMetafile, name::String)

Remove a parameter from the metafile.
"""
function remove_parameter!(metafile::DGGRIDMetafile, name::String)
    delete!(metafile.params, name)
    return metafile
end

"""
    get_parameter(metafile::DGGRIDMetafile, name::String, default=nothing)

Get a parameter value from the metafile.
"""
function get_parameter(metafile::DGGRIDMetafile, name::String, default=nothing)
    return get(metafile.params, name, default)
end

"""
    create_metafile(; kwargs...)

Create a new metafile with the given parameters.
"""
function create_metafile(; kwargs...)
    metafile = DGGRIDMetafile()
    for (key, value) in kwargs
        add_parameter!(metafile, String(key), value)
    end
    return metafile
end

# ============================================================================
# PARAMETER DEFINITIONS
# ============================================================================

"""
All DGGRID parameters with their metadata.
Focused on GENERATE_GRID, TRANSFORM_POINTS, OUTPUT_STATS operations.
"""
const PARAMETER_DEFINITIONS = Dict{String, Parameter}(
    
    # ========================================================================
    # GENERAL PARAMETERS
    # ========================================================================
    
    "dggrid_operation" => Parameter(
        "dggrid_operation",
        CHOICE,
        "GENERATE_GRID",
        ["GENERATE_GRID", "TRANSFORM_POINTS", "OUTPUT_STATS"],
        "Specifies the operation to be performed",
        "always"
    ),
    
    "precision" => Parameter(
        "precision",
        INTEGER,
        7,
        nothing,
        "Number of digits to right of decimal point for floating point output",
        "always"
    ),
    
    "verbosity" => Parameter(
        "verbosity",
        INTEGER,
        0,
        nothing,
        "Amount of debugging output (0-3)",
        "always"
    ),
    
    # ========================================================================
    # DGG SPECIFICATION PARAMETERS
    # ========================================================================
    
    "dggs_type" => Parameter(
        "dggs_type",
        CHOICE,
        "ISEA7H",
        ["CUSTOM", "ISEA7H", "IGEO7", "ISEA3H", "ISEA4H"],
        "Preset DGG type",
        "always"
    ),
    
    "dggs_topology" => Parameter(
        "dggs_topology",
        CHOICE,
        "HEXAGON",
        ["HEXAGON"],
        "Cell shape (only HEXAGON supported in this wrapper)",
        "dggs_type is CUSTOM"
    ),
    
    "dggs_proj" => Parameter(
        "dggs_proj",
        CHOICE,
        "ISEA",
        ["ISEA"],
        "Projection (only ISEA supported in this wrapper)",
        "dggs_type is CUSTOM"
    ),
    
    "dggs_aperture" => Parameter(
        "dggs_aperture",
        INTEGER,
        7,
        nothing,
        "DGGS aperture (3, 4, or 7)",
        "dggs_aperture_type is PURE"
    ),
    
    "dggs_aperture_type" => Parameter(
        "dggs_aperture_type",
        CHOICE,
        "PURE",
        ["PURE"],
        "Aperture sequence type (only PURE supported in this wrapper)",
        "dggs_topology is HEXAGON"
    ),
    
    "dggs_res_spec" => Parameter(
        "dggs_res_spec",
        INTEGER,
        9,
        nothing,
        "DGG resolution (0-35)",
        "dggs_res_specify_type is SPECIFIED"
    ),
    
    "dggs_res_specify_type" => Parameter(
        "dggs_res_specify_type",
        CHOICE,
        "SPECIFIED",
        ["SPECIFIED", "CELL_AREA", "INTERCELL_DISTANCE"],
        "How DGG resolution is specified",
        "always"
    ),
    
    "dggs_res_specify_area" => Parameter(
        "dggs_res_specify_area",
        DOUBLE,
        100.0,
        nothing,
        "Desired cell area in square kilometers",
        "dggs_res_specify_type is CELL_AREA"
    ),
    
    "dggs_res_specify_intercell_distance" => Parameter(
        "dggs_res_specify_intercell_distance",
        DOUBLE,
        100.0,
        nothing,
        "Desired intercell distance in kilometers",
        "dggs_res_specify_type is INTERCELL_DISTANCE"
    ),
    
    "dggs_res_specify_rnd_down" => Parameter(
        "dggs_res_specify_rnd_down",
        BOOLEAN,
        true,
        nothing,
        "Round down (true) or up (false) to nearest resolution",
        "dggs_res_specify_type is CELL_AREA or INTERCELL_DISTANCE"
    ),
    
    # ========================================================================
    # DGG ORIENTATION PARAMETERS
    # ========================================================================
    
    "dggs_orient_specify_type" => Parameter(
        "dggs_orient_specify_type",
        CHOICE,
        "SPECIFIED",
        ["SPECIFIED", "RANDOM", "REGION_CENTER"],
        "How DGG orientation is specified",
        "always"
    ),
    
    "dggs_vert0_lon" => Parameter(
        "dggs_vert0_lon",
        DOUBLE,
        11.25,
        nothing,
        "Longitude of icosahedron vertex 0 in degrees",
        "dggs_orient_specify_type is SPECIFIED"
    ),
    
    "dggs_vert0_lat" => Parameter(
        "dggs_vert0_lat",
        DOUBLE,
        58.28252559,
        nothing,
        "Latitude of icosahedron vertex 0 in degrees",
        "dggs_orient_specify_type is SPECIFIED"
    ),
    
    "dggs_vert0_azimuth" => Parameter(
        "dggs_vert0_azimuth",
        DOUBLE,
        0.0,
        nothing,
        "Azimuth from vertex 0 to vertex 1 in degrees",
        "dggs_orient_specify_type is SPECIFIED"
    ),
    
    "region_center_lon" => Parameter(
        "region_center_lon",
        DOUBLE,
        0.0,
        nothing,
        "Longitude of study region center in degrees",
        "dggs_orient_specify_type is REGION_CENTER"
    ),
    
    "region_center_lat" => Parameter(
        "region_center_lat",
        DOUBLE,
        0.0,
        nothing,
        "Latitude of study region center in degrees",
        "dggs_orient_specify_type is REGION_CENTER"
    ),
    
    # ========================================================================
    # EARTH RADIUS PARAMETERS
    # ========================================================================
    
    "proj_datum" => Parameter(
        "proj_datum",
        CHOICE,
        "WGS84_AUTHALIC_SPHERE",
        ["WGS84_AUTHALIC_SPHERE", "WGS84_MEAN_SPHERE", "CUSTOM_SPHERE"],
        "Earth radius datum",
        "always"
    ),
    
    "proj_datum_radius" => Parameter(
        "proj_datum_radius",
        DOUBLE,
        6371.007180918475,
        nothing,
        "Custom earth radius in kilometers",
        "proj_datum is CUSTOM_SPHERE"
    ),
    
    # ========================================================================
    # CELL OUTPUT PARAMETERS
    # ========================================================================
    
    "cell_output_type" => Parameter(
        "cell_output_type",
        CHOICE,
        "GDAL",
        ["NONE", "SHAPEFILE", "KML", "GEOJSON", "GDAL", "GDAL_COLLECTION"],
        "Cell boundary output file format",
        "dggrid_operation is GENERATE_GRID or TRANSFORM_POINTS"
    ),
    
    "cell_output_file_name" => Parameter(
        "cell_output_file_name",
        STRING,
        "cells",
        nothing,
        "Cell boundary output file name prefix",
        "cell_output_type is not NONE"
    ),
    
    "cell_output_gdal_format" => Parameter(
        "cell_output_gdal_format",
        STRING,
        "GPKG",
        nothing,
        "GDAL output format (GPKG, FlatGeobuf, Arrow, GeoJSON, ESRI Shapefile, KML)",
        "cell_output_type is GDAL or GDAL_COLLECTION"
    ),
    
    "point_output_type" => Parameter(
        "point_output_type",
        CHOICE,
        "NONE",
        ["NONE", "SHAPEFILE", "KML", "GEOJSON", "GDAL", "GDAL_COLLECTION", "TEXT"],
        "Cell point output file format",
        "dggrid_operation is GENERATE_GRID or TRANSFORM_POINTS"
    ),
    
    "point_output_file_name" => Parameter(
        "point_output_file_name",
        STRING,
        "centers",
        nothing,
        "Cell point output file name prefix",
        "point_output_type is not NONE"
    ),
    
    "point_output_gdal_format" => Parameter(
        "point_output_gdal_format",
        STRING,
        "GPKG",
        nothing,
        "GDAL output format for points",
        "point_output_type is GDAL or GDAL_COLLECTION"
    ),
    
    "output_cell_label_type" => Parameter(
        "output_cell_label_type",
        CHOICE,
        "GLOBAL_SEQUENCE",
        ["GLOBAL_SEQUENCE", "ENUMERATION", "OUTPUT_ADDRESS_TYPE"],
        "Output form for cell indexes",
        "dggrid_operation is GENERATE_GRID or TRANSFORM_POINTS"
    ),
    
    "densification" => Parameter(
        "densification",
        INTEGER,
        0,
        nothing,
        "Number of additional points per edge (0 = no densification)",
        "cell_output_type is not NONE"
    ),
    
    "longitude_wrap_mode" => Parameter(
        "longitude_wrap_mode",
        CHOICE,
        "WRAP",
        ["WRAP", "UNWRAP_EAST", "UNWRAP_WEST"],
        "How to handle longitude for cells crossing anti-meridian",
        "always"
    ),
    
    "unwrap_points" => Parameter(
        "unwrap_points",
        BOOLEAN,
        true,
        nothing,
        "Output point longitudes using longitude_wrap_mode",
        "always"
    ),
    
    "max_cells_per_output_file" => Parameter(
        "max_cells_per_output_file",
        INTEGER,
        0,
        nothing,
        "Maximum cells per output file (0 = no limit)",
        "cell_output_type is not NONE"
    ),
    
    "shapefile_id_field_length" => Parameter(
        "shapefile_id_field_length",
        INTEGER,
        11,
        nothing,
        "Number of digits in Shapefile cell index strings",
        "cell_output_type or point_output_type is SHAPEFILE"
    ),
    
    # ========================================================================
    # ADDRESS TYPE PARAMETERS
    # ========================================================================
    
    "input_address_type" => Parameter(
        "input_address_type",
        CHOICE,
        "GEO",
        ["GEO", "SEQNUM", "Q2DI", "HIERNDX"],
        "Cell address form in input files",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "output_address_type" => Parameter(
        "output_address_type",
        CHOICE,
        "SEQNUM",
        ["GEO", "SEQNUM", "Q2DI", "HIERNDX"],
        "Cell address form in output",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "input_hier_ndx_system" => Parameter(
        "input_hier_ndx_system",
        CHOICE,
        "Z7",
        ["Z7", "Z3", "ZORDER"],
        "Hierarchical indexing system for input (Z7 for aperture 7)",
        "input_address_type is HIERNDX"
    ),
    
    "output_hier_ndx_system" => Parameter(
        "output_hier_ndx_system",
        CHOICE,
        "Z7",
        ["Z7", "Z3", "ZORDER"],
        "Hierarchical indexing system for output (Z7 for aperture 7)",
        "output_address_type is HIERNDX"
    ),
    
    "input_hier_ndx_form" => Parameter(
        "input_hier_ndx_form",
        CHOICE,
        "INT64",
        ["INT64", "DIGIT_STRING"],
        "Index representation for input (INT64 or DIGIT_STRING)",
        "input_address_type is HIERNDX"
    ),
    
    "output_hier_ndx_form" => Parameter(
        "output_hier_ndx_form",
        CHOICE,
        "INT64",
        ["INT64", "DIGIT_STRING"],
        "Index representation for output (INT64 or DIGIT_STRING)",
        "output_address_type is HIERNDX"
    ),
    
    # ========================================================================
    # GRID GENERATION PARAMETERS
    # ========================================================================
    
    "clip_subset_type" => Parameter(
        "clip_subset_type",
        CHOICE,
        "WHOLE_EARTH",
        ["WHOLE_EARTH", "SHAPEFILE", "GDAL", "COARSE_CELLS", "SEQNUMS", "INPUT_ADDRESS_TYPE"],
        "How to determine portion of DGG to generate",
        "dggrid_operation is GENERATE_GRID"
    ),
    
    "clip_region_files" => Parameter(
        "clip_region_files",
        STRING,
        "",
        nothing,
        "Space-delimited list of clipping polygon files",
        "clip_subset_type is SHAPEFILE or GDAL or SEQNUMS or INPUT_ADDRESS_TYPE"
    ),
    
    "clip_cell_res" => Parameter(
        "clip_cell_res",
        INTEGER,
        1,
        nothing,
        "Resolution of coarse clipping cells",
        "clip_subset_type is COARSE_CELLS"
    ),
    
    "clip_cell_seqnums" => Parameter(
        "clip_cell_seqnums",
        STRING,
        "",
        nothing,
        "Space-delimited sequence numbers of coarse clipping cells",
        "clip_subset_type is COARSE_CELLS"
    ),

    "clip_cell_addresses" => Parameter(
        "clip_cell_addresses",
        STRING,
        "",
        nothing,
        "Space-delimited address labels of coarse clipping cells",
        "clip_subset_type is COARSE_CELLS"
    ),
    
    "clip_cell_densification" => Parameter(
        "clip_cell_densification",
        INTEGER,
        1,
        nothing,
        "Points per edge for clipping cell densification",
        "clip_subset_type is COARSE_CELLS"
    ),
    
    "geodetic_densify" => Parameter(
        "geodetic_densify",
        DOUBLE,
        0.0,
        nothing,
        "Maximum arc length in degrees for clipping polygon segments (0 = no densification)",
        "clip_subset_type is SHAPEFILE or GDAL"
    ),
    
    "clip_using_holes" => Parameter(
        "clip_using_holes",
        BOOLEAN,
        false,
        nothing,
        "Handle holes in clipping polygons",
        "clip_subset_type is GDAL"
    ),
    
    "clipper_scale_factor" => Parameter(
        "clipper_scale_factor",
        INTEGER,
        1000000,
        nothing,
        "Coarseness of polygon intersection grid",
        "clip_subset_type is SHAPEFILE or GDAL"
    ),
    
    "update_frequency" => Parameter(
        "update_frequency",
        INTEGER,
        100000,
        nothing,
        "Number of cells tested before status update",
        "dggrid_operation is GENERATE_GRID"
    ),
    
    "output_first_seqnum" => Parameter(
        "output_first_seqnum",
        INTEGER,
        0,
        nothing,
        "First cell sequence number to generate",
        "clip_subset_type is WHOLE_EARTH"
    ),
    
    "output_last_seqnum" => Parameter(
        "output_last_seqnum",
        INTEGER,
        typemax(Int),
        nothing,
        "Last cell sequence number to generate",
        "clip_subset_type is WHOLE_EARTH"
    ),
    
    # bin_coverage, not supported in 8.42 for GENERATE_GRID
    "bin_coverage" => Parameter(
        "bin_coverage",
        STRING,
        "GLOBAL",
        ["GLOBAL", "PARTIAL"],
        "Allows DGGRID to determine how to trade-off speed vs. memory usage | not supported in 8.42 for GENERATE_GRID?",
        "dggrid_operation is GENERATE_GRID"
    ),

    # ========================================================================
    # TRANSFORM POINTS PARAMETERS
    # ========================================================================
    
    "input_file_name" => Parameter(
        "input_file_name",
        STRING,
        "input.txt",
        nothing,
        "Input file name for point transformation",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "input_delimiter" => Parameter(
        "input_delimiter",
        STRING,
        " ",
        nothing,
        "Delimiter character for input files",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "output_file_name" => Parameter(
        "output_file_name",
        STRING,
        "output.txt",
        nothing,
        "Output file name for transformed points",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "output_delimiter" => Parameter(
        "output_delimiter",
        STRING,
        " ",
        nothing,
        "Delimiter character for output files",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    "output_file_type" => Parameter(
        "output_file_type",
        CHOICE,
        "TEXT",
        ["NONE", "TEXT"],
        "Text output file type",
        "dggrid_operation is TRANSFORM_POINTS"
    ),
    
    # ========================================================================
    # NEIGHBOR AND CHILDREN OUTPUT
    # ========================================================================
    
    "neighbor_output_type" => Parameter(
        "neighbor_output_type",
        CHOICE,
        "NONE",
        ["NONE", "TEXT", "GDAL_COLLECTION"],
        "Output cell neighbors",
        "dggs_topology is HEXAGON"
    ),
    
    "neighbor_output_file_name" => Parameter(
        "neighbor_output_file_name",
        STRING,
        "nbr",
        nothing,
        "Neighbors output file name",
        "neighbor_output_type is TEXT"
    ),
    
    "children_output_type" => Parameter(
        "children_output_type",
        CHOICE,
        "NONE",
        ["NONE", "TEXT", "GDAL_COLLECTION"],
        "Output cell spatial children",
        "always"
    ),
    
    "children_output_file_name" => Parameter(
        "children_output_file_name",
        STRING,
        "chd",
        nothing,
        "Children output file name",
        "children_output_type is TEXT"
    ),
    
    # ========================================================================
    # COLLECTION OUTPUT
    # ========================================================================
    
    "collection_output_file_name" => Parameter(
        "collection_output_file_name",
        STRING,
        "collection",
        nothing,
        "Collection output file name prefix",
        "any output type is GDAL_COLLECTION"
    ),
    
    "collection_output_gdal_format" => Parameter(
        "collection_output_gdal_format",
        STRING,
        "GPKG",
        nothing,
        "GDAL format for collection output",
        "any output type is GDAL_COLLECTION"
    ),
)

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

"""
    validate_parameter(name::String, value::Any)

Validate a parameter value against its definition.
Returns (is_valid::Bool, error_message::String)
"""
function validate_parameter(name::String, value::Any)
    if !haskey(PARAMETER_DEFINITIONS, name)
        return (false, "Unknown parameter: $name")
    end
    
    param_def = PARAMETER_DEFINITIONS[name]
    
    # Type checking
    if param_def.type == BOOLEAN
        if !(value isa Bool)
            return (false, "$name must be a Boolean (true/false)")
        end
    elseif param_def.type == INTEGER
        if !(value isa Integer)
            return (false, "$name must be an Integer")
        end
    elseif param_def.type == DOUBLE
        if !(value isa Real)
            return (false, "$name must be a Real number")
        end
    elseif param_def.type == STRING || param_def.type == STRINGLIST
        if !(value isa AbstractString)
            return (false, "$name must be a String")
        end
    elseif param_def.type == CHOICE
        value_str = uppercase(string(value))
        if param_def.allowed_values !== nothing && !(value_str in param_def.allowed_values)
            return (false, "$name must be one of: $(join(param_def.allowed_values, ", "))")
        end
    end
    
    return (true, "")
end

"""
    validate_metafile(metafile::DGGRIDMetafile)

Validate all parameters in a metafile.
Returns (is_valid::Bool, errors::Vector{String})
"""
function validate_metafile(metafile::DGGRIDMetafile)
    errors = String[]
    
    for (name, value) in metafile.params
        is_valid, error_msg = validate_parameter(name, value)
        if !is_valid
            push!(errors, error_msg)
        end
    end

    # Additional validation for resolution based on dggs_type
    if haskey(metafile.params, "dggs_type") && haskey(metafile.params, "dggs_res_spec")
        dggs_type = metafile.params["dggs_type"]
        dggs_res_spec = metafile.params["dggs_res_spec"]
        if haskey(max_resolution, dggs_type)
            max_res = max_resolution[dggs_type]
            if dggs_res_spec > max_res
                push!(errors, "Resolution $dggs_res_spec exceeds maximum allowed for $dggs_type ($max_res)")
            end
        end
    end
    
    return (isempty(errors), errors)
end

# ============================================================================
# METAFILE I/O FUNCTIONS
# ============================================================================

"""
    write_metafile(metafile::DGGRIDMetafile, filepath::String)

Write a metafile to disk in DGGRID format.
"""
function write_metafile(metafile::DGGRIDMetafile, filepath::String)
    open(filepath, "w") do io
        println(io, "# DGGRID Metafile")
        println(io, "# Generated by DGGRIDParams (DggridRunner.jl)")
        println(io, "")
        
        # Write parameters in a logical order
        param_order = [
            "dggrid_operation",
            "dggs_type",
            "dggs_topology",
            "dggs_proj",
            "dggs_aperture_type",
            "dggs_aperture",
            "dggs_res_specify_type",
            "dggs_res_spec",
            "dggs_res_specify_area",
            "dggs_res_specify_intercell_distance",
            "dggs_res_specify_rnd_down",
            "dggs_orient_specify_type",
            "dggs_vert0_lon",
            "dggs_vert0_lat",
            "dggs_vert0_azimuth",
            "region_center_lon",
            "region_center_lat",
            "proj_datum",
            "proj_datum_radius",
            "input_address_type",
            "input_hier_ndx_system",
            "input_hier_ndx_form",
            "output_address_type",
            "output_hier_ndx_system",
            "output_hier_ndx_form",
            "output_cell_label_type",
            "cell_output_type",
            "cell_output_file_name",
            "cell_output_gdal_format",
            "point_output_type",
            "point_output_file_name",
            "point_output_gdal_format",
            "clip_subset_type",
            "clip_region_files",
            "clip_cell_res",
            "clip_cell_seqnums",
            "clip_cell_densification",
            "geodetic_densify",
            "clip_using_holes",
            "clipper_scale_factor",
            "densification",
            "longitude_wrap_mode",
            "unwrap_points",
            "shapefile_id_field_length",
            "max_cells_per_output_file",
            "update_frequency",
            "output_first_seqnum",
            "output_last_seqnum",
            "input_file_name",
            "input_delimiter",
            "output_file_name",
            "output_delimiter",
            "output_file_type",
            "neighbor_output_type",
            "neighbor_output_file_name",
            "children_output_type",
            "children_output_file_name",
            "collection_output_file_name",
            "collection_output_gdal_format",
            "precision",
            "verbosity",
        ]
        
        # Write ordered parameters
        for param_name in param_order
            if haskey(metafile.params, param_name)
                value = metafile.params[param_name]
                write_parameter(io, param_name, value)
            end
        end
        
        # Write any remaining parameters not in the order list
        for (param_name, value) in metafile.params
            if !(param_name in param_order)
                write_parameter(io, param_name, value)
            end
        end
    end
    
    return filepath
end

"""
    write_parameter(io::IO, name::String, value::Any)

Write a single parameter to an IO stream in DGGRID format.
"""
function write_parameter(io::IO, name::String, value::Any)
    if value isa Bool
        # DGGRID uses TRUE/FALSE
        value_str = value ? "TRUE" : "FALSE"
    elseif value isa AbstractString
        value_str = value
    else
        value_str = string(value)
    end
    
    println(io, "$name $value_str")
end

"""
    get_parameter_info(name::String)

Get information about a parameter.
"""
function get_parameter_info(name::String)
    if haskey(PARAMETER_DEFINITIONS, name)
        return PARAMETER_DEFINITIONS[name]
    else
        return nothing
    end
end

"""
    list_parameters(; operation::Union{Nothing,String}=nothing)

List all available parameters, optionally filtered by operation.
"""
function list_parameters(; operation::Union{Nothing,String}=nothing)
    params = collect(keys(PARAMETER_DEFINITIONS))
    sort!(params)
    
    if operation !== nothing
        # Filter by operation (simplified - could be more sophisticated)
        filter!(p -> contains(PARAMETER_DEFINITIONS[p].used_when, operation) || 
                     PARAMETER_DEFINITIONS[p].used_when == "always", params)
    end
    
    return params
end

# read metafile function that takes a filepath and returns a DGGRIDMetafile
function read_metafile(filepath::String)
    metafile = DGGRIDMetafile()
    open(filepath, "r") do io
        for line in eachline(io)
            # Skip empty lines and comments
            if isempty(strip(line)) || startswith(strip(line), "#")
                continue
            end
            # split the line at the first whitespace
            # first part is the parameter name, and the rest second part is the value
            parts = split(line, r"\s+", limit=2)
            if length(parts) == 2
                name = parts[1]
                value = parts[2]
                # Convert value to appropriate type if needed
                if haskey(PARAMETER_DEFINITIONS, name)
                    param_info = PARAMETER_DEFINITIONS[name]
                    if param_info.type == BOOLEAN
                        value = lowercase(value) == "true"
                    elseif param_info.type == INTEGER
                        value = parse(Int, value)
                    elseif param_info.type == DOUBLE
                        value = parse(Float64, value)
                    end
                end
                metafile.params[name] = value
            end
        end
    end
    return metafile
end

end # module DGGRIDParams
