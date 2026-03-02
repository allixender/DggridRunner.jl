"""
Test suite for dggrid_runner_lib.jl

Tests the DGGRID runner library including:
- operation preps and prep operations
- run operations and read in results
- Metafile I/O operations
"""

using Test

# Include the module to test
# include("DggridRunner.jl")
using DGGRIDRunner

@testset "DGGRIDRunner Tests" begin
    
    # ========================================================================
    # Test Basic functions
    # ========================================================================
    
    @testset "Basic Functions" begin
        @testset "output_stats function" begin
            dggs_type = "ISEA7H"
            resolution = 20
            is_ok = prep_output_stats(dggs_type, resolution)
            @test is_ok[1] == true
            @test typeof(is_ok[2]) == DGGRIDParams.DGGRIDMetafile

            shell_out = """Earth Radius: 6,371.0071809

Res                 # Cells        Area (km^2)      CLS (km)
  0                      12 51,006,562.1724089 8,199.5003701
  1                      72  7,286,651.7389156 3,053.2232428
  2                     492  1,040,950.2484165 1,151.6430095
  3                   3,432    148,707.1783452   435.1531492
  4                  24,012     21,243.8826207   164.4655799
  5                 168,072      3,034.8403744    62.1617764
  6               1,176,492        433.5486249    23.4949231
  7               8,235,432         61.9355178     8.8802451
  8              57,648,012          8.8479311     3.3564171
  9             403,536,072          1.2639902     1.2686064
 10           2,824,752,492          0.1805700     0.4794882
 11          19,773,267,432          0.0257957     0.1812295
 12         138,412,872,012          0.0036851     0.0684983
 13         968,890,104,072          0.0005264     0.0258899
 14       6,782,230,728,492          0.0000752     0.0097855
 15      47,475,615,099,432          0.0000107     0.0036986
 16     332,329,305,696,012          0.0000015     0.0013979
 17   2,326,305,139,872,072          0.0000002     0.0005284
 18  16,284,135,979,104,492          0.0000000     0.0001997
 19 113,988,951,853,731,432          0.0000000     0.0000755
 20 797,922,662,976,120,012          0.0000000     0.0000285"""

            parsed_stats = parse_output_stats(split(shell_out, '\n'))
            @test haskey(parsed_stats, resolution)

            @test haskey(parsed_stats[3], :num_cells) == true # should be symbol, not string
            @test parsed_stats[3].num_cells == 3432

            params = is_ok[2]
            stats_io = run_output_stats(params)
            @test haskey(stats_io[3], :num_cells) == true # should be symbol, not string
            @test stats_io[3].num_cells == 3432
            @test stats_io[20].num_cells == 797922662976120012

        end
    end
    
    # ========================================================================
    # Test GENERATE_GRID functions
    # ========================================================================
    
    @testset "GENERATE_GRID Functions" begin
        
        @testset "prep_generate_grid_whole_earth" begin
            # Test basic whole earth grid preparation
            dggs_type = "ISEA7H"
            resolution = 3
            
            success, params, output_path = prep_generate_grid_whole_earth(
                dggs_type, resolution,
                output_format="GPKG"
            )
            
            @test success == true
            @test typeof(params) == DGGRIDParams.DGGRIDMetafile
            @test endswith(output_path, ".gpkg")
            @test DGGRIDParams.get_parameter(params, "dggrid_operation") == "GENERATE_GRID"
            @test DGGRIDParams.get_parameter(params, "dggs_type") == dggs_type
            @test DGGRIDParams.get_parameter(params, "dggs_res_spec") == resolution
            @test DGGRIDParams.get_parameter(params, "clip_subset_type") == "WHOLE_EARTH"
            @test DGGRIDParams.get_parameter(params, "cell_output_type") == "GDAL"
            @test DGGRIDParams.get_parameter(params, "cell_output_gdal_format") == "GPKG"
            @test DGGRIDParams.get_parameter(params, "point_output_type") == "NONE"
            
            # Test with HIERNDX address type
            success2, params2, output_path2 = prep_generate_grid_whole_earth(
                "IGEO7", 2,
                address_type="HIERNDX",
                output_format="FlatGeobuf"
            )
            
            @test success2 == true
            @test endswith(output_path2, ".fgb")
            @test DGGRIDParams.get_parameter(params2, "output_cell_label_type") == "OUTPUT_ADDRESS_TYPE"
            @test DGGRIDParams.get_parameter(params2, "output_address_type") == "HIERNDX"
            @test DGGRIDParams.get_parameter(params2, "output_hier_ndx_system") == "Z7"
            
            # Test with densification
            success3, params3, output_path3 = prep_generate_grid_whole_earth(
                "ISEA7H", 2,
                densification=2
            )
            
            @test success3 == true
            @test DGGRIDParams.get_parameter(params3, "densification") == 2
        end
        
        @testset "prep_generate_grid_coarse_cells" begin
            # Test coarse cells grid preparation
            dggs_type = "ISEA7H"
            resolution = 5
            coarse_res = 2
            coarse_seqnums = string.([1, 2, 3])
            
            success, params, output_path = prep_generate_grid_coarse_cells(
                dggs_type, resolution, coarse_res, coarse_seqnums,
                output_format="GPKG"
            )
            
            @test success == true
            @test typeof(params) == DGGRIDParams.DGGRIDMetafile
            @test endswith(output_path, ".gpkg")
            @test DGGRIDParams.get_parameter(params, "dggrid_operation") == "GENERATE_GRID"
            @test DGGRIDParams.get_parameter(params, "clip_subset_type") == "COARSE_CELLS"
            @test DGGRIDParams.get_parameter(params, "clip_cell_res") == coarse_res
            @test DGGRIDParams.get_parameter(params, "clip_cell_seqnums") == "1 2 3"
            @test DGGRIDParams.get_parameter(params, "clip_cell_densification") == 1
            
            # Test validation: coarse_res must be < resolution
            success_fail, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 5, ["1"]
            )
            @test success_fail == false
            
            # Test validation: coarse_res must be > 0
            success_fail2, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 0, ["1"]
            )
            @test success_fail2 == false
            
            # Test validation: coarse_seqnums cannot be empty
            success_fail3, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 2, String[]
            )
            @test success_fail3 == false
            
            # Test with custom clip densification
            success2, params2, _ = prep_generate_grid_coarse_cells(
                "ISEA3H",
                6,
                3,
                ["01001"],
                address_type="HIERNDX",
                clip_densification=3
            )
            @test success2 == true
            @test DGGRIDParams.get_parameter(params2, "output_hier_ndx_system") == "Z3"
            @test DGGRIDParams.get_parameter(params2, "clip_cell_densification") == 3
            
            # Apply convenience function
            success2 = grid_gen_convenience!(params2)
            @test success2 == true
            
            # Check that convenience parameters were added
            @test DGGRIDParams.get_parameter(params2, "output_hier_ndx_form") == "digit_string"
            @test DGGRIDParams.get_parameter(params2, "input_hier_ndx_form") == "digit_string"
        end

        @testset "prep_generate_grid_coarse_cells with IGEO7" begin
            # Test coarse cells grid preparation
            dggs_type = "IGEO7"
            resolution = 5
            coarse_res = 2
            coarse_cell_labels = ["0012", "0013", "0014"]
            
            success, params, output_path = prep_generate_grid_coarse_cells(
                dggs_type,
                resolution,
                coarse_res,
                coarse_cell_labels,
                output_format="GPKG",
                address_type="HIERNDX"
            )
            
            @test success == true
            @test typeof(params) == DGGRIDParams.DGGRIDMetafile
            @test endswith(output_path, ".gpkg")
            @test DGGRIDParams.get_parameter(params, "dggrid_operation") == "GENERATE_GRID"
            @test DGGRIDParams.get_parameter(params, "clip_subset_type") == "COARSE_CELLS"
            @test DGGRIDParams.get_parameter(params, "clip_cell_res") == coarse_res
            @test DGGRIDParams.get_parameter(params, "clip_cell_seqnums") == nothing
            @test DGGRIDParams.get_parameter(params, "clip_cell_addresses") == "0012 0013 0014"
            @test DGGRIDParams.get_parameter(params, "clip_cell_densification") == 1
            
            # Test validation: coarse_res must be < resolution
            success_fail, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 5, ["1"]
            )
            @test success_fail == false
            
            # Test validation: coarse_res must be > 0
            success_fail2, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 0, ["1"]
            )
            @test success_fail2 == false
            
            # Test validation: coarse_seqnums cannot be empty
            success_fail3, _, _ = prep_generate_grid_coarse_cells(
                dggs_type, 5, 2, String[]
            )
            @test success_fail3 == false
            
            # Test with custom clip densification
            success2, params2, _ = prep_generate_grid_coarse_cells(
                "IGEO7", 6, 3, ["0110", "0120", "1310"],
                clip_densification=3
            )
            
            @test success2 == true
            @test DGGRIDParams.get_parameter(params2, "clip_cell_densification") == 3
        end
        
        @testset "prep_generate_grid_clip_region" begin
            # Create a temporary test file to simulate a clip region file
            test_clip_file = tempname() * ".geojson"
            write(test_clip_file, """
            {
              "type": "FeatureCollection",
              "features": [{
                "type": "Feature",
                "geometry": {
                  "type": "Polygon",
                  "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]]
                },
                "properties": {}
              }]
            }
            """)
            
            try
                dggs_type = "ISEA7H"
                resolution = 4
                
                success, params, output_path = prep_generate_grid_clip_region(
                    dggs_type, resolution, test_clip_file,
                    output_format="GPKG"
                )
                
                @test success == true
                @test typeof(params) == DGGRIDParams.DGGRIDMetafile
                @test endswith(output_path, ".gpkg")
                @test DGGRIDParams.get_parameter(params, "dggrid_operation") == "GENERATE_GRID"
                @test DGGRIDParams.get_parameter(params, "clip_subset_type") == "GDAL"
                @test DGGRIDParams.get_parameter(params, "clip_region_files") == test_clip_file
                
                # Test with geodetic densification
                success2, params2, _ = prep_generate_grid_clip_region(
                    "IGEO7", 5, test_clip_file,
                    geodetic_densify=1.0
                )
                
                @test success2 == true
                @test DGGRIDParams.get_parameter(params2, "geodetic_densify") == 1.0
                
                # Test with holes enabled
                success3, params3, _ = prep_generate_grid_clip_region(
                    "ISEA7H", 4, test_clip_file,
                    use_holes=true
                )
                
                @test success3 == true
                @test DGGRIDParams.get_parameter(params3, "clip_using_holes") == true
                
                success4, params4, _ = prep_generate_grid_clip_region(
                    "ISEA3H",
                    4, test_clip_file,
                    address_type="HIERNDX",
                    densification=2
                )
                @test success4 == true

                # convenience function should be applied
                success5 = grid_gen_convenience!(params4)
                @test success5 == true
                @test DGGRIDParams.get_parameter(params4, "output_hier_ndx_form") == "digit_string"
                @test DGGRIDParams.get_parameter(params4, "output_hier_ndx_system") == "Z3"
                @test DGGRIDParams.get_parameter(params4, "densification") == 2
                
            finally
                # Clean up test file
                rm(test_clip_file; force=true)
            end
            
            # Test validation: clip file must exist
            success_fail, _, _ = prep_generate_grid_clip_region(
                "ISEA7H", 4, "/nonexistent/file.shp"
            )
            @test success_fail == false

            # TODO: prep_generate_grid_clip_cells, write temp txt file with 3 rows 0012 0013 0014
            
        end
    end
    
    # ========================================================================
    # Test grid_gen_convenience! function
    # ========================================================================
    
    @testset "grid_gen_convenience! Function" begin
        
        @testset "HIERNDX address type with ISEA7H" begin
            # Create params with HIERNDX address type
            success, params, _ = prep_generate_grid_whole_earth(
                "ISEA7H", 3,
                address_type="HIERNDX"
            )
            
            @test success == true
            
            # Apply convenience function
            success = grid_gen_convenience!(params)
            @test success == true
            
            # Check that convenience parameters were added
            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_form") == "digit_string"
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lon") == 11.2
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lat") == 58.2825
            @test DGGRIDParams.get_parameter(params, "longitude_wrap_mode") == "UNWRAP_EAST"
        end
        
        @testset "HIERNDX address type with IGEO7" begin
            # Create params with HIERNDX address type
            success, params, _ = prep_generate_grid_whole_earth(
                "IGEO7", 4,
                address_type="HIERNDX"
            )
            
            @test success == true
            
            # Apply convenience function
            grid_gen_convenience!(params)
            
            # Check that convenience parameters were added
            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_form") == "digit_string"
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lon") == 11.2
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lat") == 58.2825
            @test DGGRIDParams.get_parameter(params, "longitude_wrap_mode") == "UNWRAP_EAST"
        end
        
        @testset "HIERNDX address type with ISEA3H" begin
            # Create params with HIERNDX address type
            success, params, _ = prep_generate_grid_whole_earth(
                "ISEA3H", 5,
                address_type="HIERNDX"
            )
            @test success == true

            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_system") == "Z3"
            
            # Apply convenience function
            success = grid_gen_convenience!(params)
            @test success == true
            
            # Check that convenience parameters were added
            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_form") == "digit_string"
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lon") == 11.2
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lat") == 58.2825
            @test DGGRIDParams.get_parameter(params, "longitude_wrap_mode") == "UNWRAP_EAST"
        end
        
        @testset "Non-HIERNDX address type (SEQNUM)" begin
            # Create params with SEQNUM address type
            success, params, _ = prep_generate_grid_whole_earth(
                "ISEA7H", 3,
                address_type="SEQNUM"
            )
            
            @test success == true
            
            # Apply convenience function
            grid_gen_convenience!(params)
            
            # Check that HIERNDX-specific parameters were NOT added
            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_form") === nothing
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lon") == 11.2
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lat") == 58.2825
            @test DGGRIDParams.get_parameter(params, "longitude_wrap_mode") == "UNWRAP_EAST"
        end
        
        @testset "Shapefile output type" begin
            # Create params with SHAPEFILE output type
            params = DGGRIDParams.DGGRIDMetafile()
            DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
            DGGRIDParams.add_parameter!(params, "dggs_type", "ISEA7H")
            DGGRIDParams.add_parameter!(params, "dggs_res_spec", 3)
            DGGRIDParams.add_parameter!(params, "cell_output_type", "SHAPEFILE")
            DGGRIDParams.add_parameter!(params, "output_address_type", "SEQNUM")
            
            # Apply convenience function
            grid_gen_convenience!(params)
            
            # Check that shapefile_id_field_length was added
            @test DGGRIDParams.get_parameter(params, "shapefile_id_field_length") == 22
        end
        
        @testset "Non-Shapefile output type (GDAL)" begin
            # Create params with GDAL output type
            success, params, _ = prep_generate_grid_whole_earth(
                "ISEA7H", 3,
                output_format="GPKG"
            )
            
            @test success == true
            
            # Apply convenience function
            grid_gen_convenience!(params)
            
            # Check that shapefile_id_field_length was NOT added
            @test DGGRIDParams.get_parameter(params, "shapefile_id_field_length") === nothing
        end
        
        @testset "Combined: HIERNDX + Shapefile" begin
            # Create params with both HIERNDX and SHAPEFILE
            params = DGGRIDParams.DGGRIDMetafile()
            DGGRIDParams.add_parameter!(params, "dggrid_operation", "GENERATE_GRID")
            DGGRIDParams.add_parameter!(params, "dggs_type", "IGEO7")
            DGGRIDParams.add_parameter!(params, "dggs_res_spec", 4)
            DGGRIDParams.add_parameter!(params, "output_address_type", "HIERNDX")
            DGGRIDParams.add_parameter!(params, "cell_output_type", "SHAPEFILE")
            
            # Apply convenience function
            grid_gen_convenience!(params)
            
            # Check that both sets of parameters were added
            @test DGGRIDParams.get_parameter(params, "output_hier_ndx_form") == "digit_string"
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lon") == 11.2
            @test DGGRIDParams.get_parameter(params, "dggs_vert0_lat") == 58.2825
            @test DGGRIDParams.get_parameter(params, "longitude_wrap_mode") == "UNWRAP_EAST"
            @test DGGRIDParams.get_parameter(params, "shapefile_id_field_length") == 22
        end
    end
end

