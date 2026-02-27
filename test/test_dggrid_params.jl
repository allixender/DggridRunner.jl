"""
Test suite for dggrid_params.jl

Tests the DGGRID parameter library including:
- Parameter definitions and metadata
- Metafile creation and manipulation
- Parameter validation
- Metafile I/O operations
"""

using Test

# Include the module to test
# include("DggridParams.jl")
using DggridRunner
using .DggridParams

@testset "DggridParams Tests" begin
    
    # ========================================================================
    # Test Parameter Definitions
    # ========================================================================
    
    @testset "Parameter Definitions" begin
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "dggrid_operation")
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "dggs_type")
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "dggs_res_spec")
        
        # Test operation parameter
        op_param = DggridParams.PARAMETER_DEFINITIONS["dggrid_operation"]
        @test op_param.type == DggridParams.CHOICE
        @test "GENERATE_GRID" in op_param.allowed_values
        @test "TRANSFORM_POINTS" in op_param.allowed_values
        @test "OUTPUT_STATS" in op_param.allowed_values
        
        # Test DGG type parameter
        dgg_param = DggridParams.PARAMETER_DEFINITIONS["dggs_type"]
        @test dgg_param.type == DggridParams.CHOICE
        @test "ISEA7H" in dgg_param.allowed_values
        @test "IGEO7" in dgg_param.allowed_values
        
        # Test address type parameters
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "input_address_type")
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "output_address_type")
        addr_param = DggridParams.PARAMETER_DEFINITIONS["output_address_type"]
        @test "SEQNUM" in addr_param.allowed_values
        @test "Q2DI" in addr_param.allowed_values
        @test "HIERNDX" in addr_param.allowed_values
        
        # Test hierarchical index parameters
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "output_hier_ndx_system")
        @test haskey(DggridParams.PARAMETER_DEFINITIONS, "output_hier_ndx_form")
        hier_form = DggridParams.PARAMETER_DEFINITIONS["output_hier_ndx_form"]
        @test "INT64" in hier_form.allowed_values
        @test "DIGIT_STRING" in hier_form.allowed_values
    end
    
    # ========================================================================
    # Test Metafile Creation
    # ========================================================================
    
    @testset "Metafile Creation" begin
        # Test empty metafile
        mf = DGGRIDMetafile()
        @test isempty(mf.params)
        
        # Test metafile with parameters
        mf2 = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H",
            dggs_res_spec=9
        )
        @test length(mf2.params) == 3
        @test mf2.params["dggrid_operation"] == "GENERATE_GRID"
        @test mf2.params["dggs_type"] == "ISEA7H"
        @test mf2.params["dggs_res_spec"] == 9
    end
    
    # ========================================================================
    # Test Parameter Manipulation
    # ========================================================================
    
    @testset "Parameter Manipulation" begin
        mf = DGGRIDMetafile()
        
        # Test adding parameters
        add_parameter!(mf, "dggrid_operation", "GENERATE_GRID")
        @test get_parameter(mf, "dggrid_operation") == "GENERATE_GRID"
        
        add_parameter!(mf, "dggs_res_spec", 10)
        @test get_parameter(mf, "dggs_res_spec") == 10
        
        # Test updating parameters
        add_parameter!(mf, "dggs_res_spec", 12)
        @test get_parameter(mf, "dggs_res_spec") == 12
        
        # Test getting non-existent parameter with default
        @test get_parameter(mf, "nonexistent", "default") == "default"
        @test get_parameter(mf, "nonexistent") === nothing
    end
    
    # ========================================================================
    # Test Parameter Validation
    # ========================================================================
    
    @testset "Parameter Validation" begin
        # Test valid parameters
        @test DggridParams.validate_parameter("dggrid_operation", "GENERATE_GRID")[1] == true
        @test DggridParams.validate_parameter("dggs_res_spec", 9)[1] == true
        @test DggridParams.validate_parameter("precision", 7)[1] == true
        @test DggridParams.validate_parameter("verbosity", 0)[1] == true
        @test DggridParams.validate_parameter("dggs_vert0_lon", 11.25)[1] == true
        @test DggridParams.validate_parameter("clip_using_holes", true)[1] == true
        
        # Test invalid parameter name
        @test DggridParams.validate_parameter("invalid_param", "value")[1] == false
        
        # Test invalid parameter types
        @test DggridParams.validate_parameter("dggs_res_spec", "not_an_int")[1] == false
        @test DggridParams.validate_parameter("precision", 7.5)[1] == false  # Should be integer
        @test DggridParams.validate_parameter("dggs_vert0_lon", "not_a_number")[1] == false
        @test DggridParams.validate_parameter("clip_using_holes", "not_a_bool")[1] == false
        
        # Test invalid choice values
        @test DggridParams.validate_parameter("dggrid_operation", "INVALID_OP")[1] == false
        @test DggridParams.validate_parameter("dggs_type", "INVALID_TYPE")[1] == false
        
        # Test valid choice values (case insensitive)
        @test DggridParams.validate_parameter("dggrid_operation", "generate_grid")[1] == true
        @test DggridParams.validate_parameter("dggs_type", "isea7h")[1] == true
    end
    
    @testset "Metafile Validation" begin
        # Test valid metafile
        mf = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H",
            dggs_res_spec=9,
            precision=7
        )
        is_valid, errors = validate_metafile(mf)
        @test is_valid == true
        @test isempty(errors)
        
        # Test invalid metafile
        mf_invalid = DGGRIDMetafile()
        add_parameter!(mf_invalid, "dggrid_operation", "INVALID_OP")
        add_parameter!(mf_invalid, "dggs_res_spec", "not_an_int")
        is_valid, errors = validate_metafile(mf_invalid)
        @test is_valid == false
        @test length(errors) == 2

        # Test max resolution validation
        mf_highres = create_metafile(
            dggrid_operation="OUTPUT_STATS",
            dggs_type="ISEA7H",
            dggs_res_spec=25  # Exceeds max for ISEA7H
        )
        is_valid, errors = validate_metafile(mf_highres)
        @test is_valid == false
        @test length(errors) == 1
    end
    
    # ========================================================================
    # Test Address Type Parameters
    # ========================================================================
    
    @testset "Address Type Parameters" begin
        # Test SEQNUM address type
        mf = create_metafile(
            output_address_type="SEQNUM"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test Q2DI address type
        mf = create_metafile(
            output_address_type="Q2DI"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test HIERNDX with Z7 and INT64
        mf = create_metafile(
            output_address_type="HIERNDX",
            output_hier_ndx_system="Z7",
            output_hier_ndx_form="INT64"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test HIERNDX with Z7 and DIGIT_STRING
        mf = create_metafile(
            output_address_type="HIERNDX",
            output_hier_ndx_system="Z7",
            output_hier_ndx_form="DIGIT_STRING"
        )
        @test validate_metafile(mf)[1] == true
    end
    
    # ========================================================================
    # Test DGG Presets
    # ========================================================================
    
    @testset "DGG Presets" begin
        # Test ISEA7H preset
        mf = create_metafile(
            dggs_type="ISEA7H",
            dggs_res_spec=9
        )
        @test validate_metafile(mf)[1] == true
        
        # Test IGEO7 preset
        mf = create_metafile(
            dggs_type="IGEO7",
            dggs_res_spec=9
        )
        @test validate_metafile(mf)[1] == true
    end
    
    # ========================================================================
    # Test Operations
    # ========================================================================
    
    @testset "Operations" begin
        # Test GENERATE_GRID operation
        mf = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H",
            dggs_res_spec=9,
            clip_subset_type="WHOLE_EARTH"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test TRANSFORM_POINTS operation
        mf = create_metafile(
            dggrid_operation="TRANSFORM_POINTS",
            dggs_type="ISEA7H",
            dggs_res_spec=9,
            input_address_type="GEO",
            output_address_type="SEQNUM",
            input_file_name="input.txt",
            output_file_name="output.txt"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test OUTPUT_STATS operation
        mf = create_metafile(
            dggrid_operation="OUTPUT_STATS",
            dggs_type="ISEA7H",
            dggs_res_spec=9
        )
        @test validate_metafile(mf)[1] == true
    end
    
    # ========================================================================
    # Test Output Formats
    # ========================================================================
    
    @testset "Output Formats" begin
        # Test GDAL output with GeoPackage
        mf = create_metafile(
            cell_output_type="GDAL",
            cell_output_gdal_format="GPKG",
            cell_output_file_name="cells"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test GDAL output with FlatGeobuf
        mf = create_metafile(
            cell_output_type="GDAL",
            cell_output_gdal_format="FlatGeobuf",
            cell_output_file_name="cells"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test GDAL output with Arrow
        mf = create_metafile(
            cell_output_type="GDAL",
            cell_output_gdal_format="Arrow",
            cell_output_file_name="cells"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test Shapefile output
        mf = create_metafile(
            cell_output_type="SHAPEFILE",
            cell_output_file_name="cells"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test KML output
        mf = create_metafile(
            cell_output_type="KML",
            cell_output_file_name="cells"
        )
        @test validate_metafile(mf)[1] == true
    end
    
    # ========================================================================
    # Test Clip Subset Types
    # ========================================================================
    
    @testset "Clip Subset Types" begin
        # Test WHOLE_EARTH
        mf = create_metafile(
            clip_subset_type="WHOLE_EARTH"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test COARSE_CELLS
        mf = create_metafile(
            clip_subset_type="COARSE_CELLS",
            clip_cell_res=3,
            clip_cell_seqnums="1 2 3 4"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test SHAPEFILE
        mf = create_metafile(
            clip_subset_type="SHAPEFILE",
            clip_region_files="region.shp"
        )
        @test validate_metafile(mf)[1] == true
        
        # Test GDAL
        mf = create_metafile(
            clip_subset_type="GDAL",
            clip_region_files="region.gpkg"
        )
        @test validate_metafile(mf)[1] == true
    end
    
    # ========================================================================
    # Test Metafile I/O
    # ========================================================================
    
    @testset "Metafile I/O" begin
        # Create a test metafile
        mf = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H",
            dggs_res_spec=9,
            clip_subset_type="WHOLE_EARTH",
            cell_output_type="GDAL",
            cell_output_gdal_format="GPKG",
            cell_output_file_name="test_cells",
            precision=7,
            verbosity=1
        )
        
        # Write to temporary file
        temp_file = tempname() * ".meta"
        write_metafile(mf, temp_file)
        
        # Check file exists
        @test isfile(temp_file)
        
        # Read and verify content
        content = read(temp_file, String)
        @test contains(content, "dggrid_operation GENERATE_GRID")
        @test contains(content, "dggs_type ISEA7H")
        @test contains(content, "dggs_res_spec 9")
        @test contains(content, "clip_subset_type WHOLE_EARTH")
        @test contains(content, "cell_output_type GDAL")
        @test contains(content, "cell_output_gdal_format GPKG")
        @test contains(content, "precision 7")
        @test contains(content, "verbosity 1")
        
        # Test boolean formatting
        mf_bool = create_metafile(
            clip_using_holes=true,
            unwrap_points=false
        )
        temp_file_bool = tempname() * ".meta"
        write_metafile(mf_bool, temp_file_bool)
        content_bool = read(temp_file_bool, String)
        @test contains(content_bool, "clip_using_holes TRUE")
        @test contains(content_bool, "unwrap_points FALSE")
        
        # Test reading metafile
        mf_read = read_metafile(temp_file)
        @test get_parameter(mf_read, "dggrid_operation") == "GENERATE_GRID"
        @test get_parameter(mf_read, "dggs_type") == "ISEA7H"
        @test get_parameter(mf_read, "dggs_res_spec") == 9
        @test get_parameter(mf_read, "clip_subset_type") == "WHOLE_EARTH"
        @test get_parameter(mf_read, "cell_output_type") == "GDAL"

        # Clean up
        rm(temp_file, force=true)
        rm(temp_file_bool, force=true)

    end
    
    @testset "Metafile read extern" begin
        filename = joinpath(@__DIR__, "metafile_example_global_grid")
        mf = read_metafile(filename)
        @test get_parameter(mf, "dggrid_operation") == "GENERATE_GRID"
        @test get_parameter(mf, "output_cell_label_type") == "OUTPUT_ADDRESS_TYPE"
        @test get_parameter(mf, "longitude_wrap_mode") == "UNWRAP_EAST"
    end

    # ========================================================================
    # Test Parameter Info Functions
    # ========================================================================
    
    @testset "Parameter Info Functions" begin
        # Test get_parameter_info
        info = DggridParams.get_parameter_info("dggrid_operation")
        @test info !== nothing
        @test info.name == "dggrid_operation"
        @test info.type == DggridParams.CHOICE
        
        # Test non-existent parameter
        @test DggridParams.get_parameter_info("nonexistent") === nothing
        
        # Test list_parameters
        all_params = DggridParams.list_parameters()
        @test length(all_params) > 0
        @test "dggrid_operation" in all_params
        @test "dggs_type" in all_params
        
        # Test filtered list
        gen_params = DggridParams.list_parameters(operation="GENERATE_GRID")
        @test length(gen_params) > 0
        @test "clip_subset_type" in gen_params
    end
    
    # ========================================================================
    # Test Complex Scenarios
    # ========================================================================
    
    @testset "Complex Scenarios" begin
        # Test complete GENERATE_GRID configuration
        mf = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H",
            dggs_res_spec=10,
            clip_subset_type="GDAL",
            clip_region_files="region.gpkg",
            geodetic_densify=1.0,
            cell_output_type="GDAL",
            cell_output_gdal_format="GPKG",
            cell_output_file_name="output_cells",
            point_output_type="GDAL",
            point_output_gdal_format="GPKG",
            point_output_file_name="output_points",
            output_address_type="HIERNDX",
            output_hier_ndx_system="Z7",
            output_hier_ndx_form="INT64",
            densification=2,
            precision=10,
            verbosity=2
        )
        is_valid, errors = validate_metafile(mf)
        @test is_valid == true
        @test isempty(errors)
        
        # Test complete TRANSFORM_POINTS configuration
        mf2 = create_metafile(
            dggrid_operation="TRANSFORM_POINTS",
            dggs_type="IGEO7",
            dggs_res_spec=12,
            input_address_type="GEO",
            output_address_type="HIERNDX",
            output_hier_ndx_system="Z7",
            output_hier_ndx_form="DIGIT_STRING",
            input_file_name="points.txt",
            output_file_name="transformed.txt",
            input_delimiter=",",
            output_delimiter=",",
            output_file_type="TEXT",
            precision=8
        )
        is_valid2, errors2 = validate_metafile(mf2)
        @test is_valid2 == true
        @test isempty(errors2)
    end
    
    # ========================================================================
    # Test Edge Cases
    # ========================================================================
    
    @testset "Edge Cases" begin
        # Test empty metafile validation
        mf_empty = DGGRIDMetafile()
        is_valid, errors = validate_metafile(mf_empty)
        @test is_valid == true  # Empty metafile is technically valid
        @test isempty(errors)
        
        # Test metafile with only invalid parameters
        mf_invalid = DGGRIDMetafile()
        add_parameter!(mf_invalid, "invalid_param1", "value1")
        add_parameter!(mf_invalid, "invalid_param2", "value2")
        # Note: validation only checks known parameters, unknown ones are ignored
        
        # Test numeric edge cases
        mf_numeric = create_metafile(
            dggs_res_spec=0,  # Minimum resolution
            precision=0,      # Minimum precision
            verbosity=3       # Maximum verbosity
        )
        @test validate_metafile(mf_numeric)[1] == true
        
        # Test large numbers
        mf_large = create_metafile(
            output_last_seqnum=999999999,
            clipper_scale_factor=10000000
        )
        @test validate_metafile(mf_large)[1] == true

        # test remove_parameter!
        mf_remove = create_metafile(
            dggrid_operation="GENERATE_GRID",
            dggs_type="ISEA7H"
        )
        @test haskey(mf_remove.params, "dggrid_operation")
        remove_parameter!(mf_remove, "dggrid_operation")
        @test !haskey(mf_remove.params, "dggrid_operation")
    end
end

println("\n✓ All tests completed!")
