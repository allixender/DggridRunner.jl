"""
Test suite for point_output parameter in prep_generate_grid_* functions

Tests that the point_output parameter correctly switches between cell and point output
"""

using Test

# Include the module to test
# include("dggrid_runner_lib.jl")
using DggridRunner

@testset "Point Output Parameter Tests" begin
    
    # ========================================================================
    # Test prep_generate_grid_whole_earth with point_output
    # ========================================================================
    
    @testset "prep_generate_grid_whole_earth - cell output (point_output=false)" begin
        dggs_type = "ISEA7H"
        resolution = 3
        
        success, params, output_path = prep_generate_grid_whole_earth(
            dggs_type, resolution,
            output_format="GPKG",
            point_output=false
        )
        
        @test success == true
        @test typeof(params) == DggridParams.DGGRIDMetafile
        
        # Check that cell output is configured
        @test DggridParams.get_parameter(params, "cell_output_type") == "GDAL"
        @test DggridParams.get_parameter(params, "cell_output_gdal_format") == "GPKG"
        @test DggridParams.get_parameter(params, "cell_output_file_name") !== nothing
        
        # Check that point output is disabled
        @test DggridParams.get_parameter(params, "point_output_type") == "NONE"
        @test DggridParams.get_parameter(params, "point_output_gdal_format") === nothing
        @test DggridParams.get_parameter(params, "point_output_file_name") === nothing
    end
    
    @testset "prep_generate_grid_whole_earth - point output (point_output=true)" begin
        dggs_type = "ISEA7H"
        resolution = 3
        
        success, params, output_path = prep_generate_grid_whole_earth(
            dggs_type, resolution,
            output_format="GPKG",
            point_output=true
        )
        
        @test success == true
        @test typeof(params) == DggridParams.DGGRIDMetafile
        
        # Check that point output is configured
        @test DggridParams.get_parameter(params, "point_output_type") == "GDAL"
        @test DggridParams.get_parameter(params, "point_output_gdal_format") == "GPKG"
        @test DggridParams.get_parameter(params, "point_output_file_name") !== nothing
        
        # Check that cell output is disabled
        @test DggridParams.get_parameter(params, "cell_output_type") == "NONE"
        @test DggridParams.get_parameter(params, "cell_output_gdal_format") === nothing
        @test DggridParams.get_parameter(params, "cell_output_file_name") === nothing
    end
    
    @testset "prep_generate_grid_whole_earth - different formats with point_output" begin
        formats = ["GPKG", "FlatGeobuf", "GeoJSON"]
        
        for format in formats
            # Test with cell output
            success_cell, params_cell, _ = prep_generate_grid_whole_earth(
                "IGEO7", 2,
                output_format=format,
                point_output=false
            )
            @test success_cell == true
            @test DggridParams.get_parameter(params_cell, "cell_output_gdal_format") == format
            @test DggridParams.get_parameter(params_cell, "point_output_type") == "NONE"
            
            # Test with point output
            success_point, params_point, _ = prep_generate_grid_whole_earth(
                "IGEO7", 2,
                output_format=format,
                point_output=true
            )
            @test success_point == true
            @test DggridParams.get_parameter(params_point, "point_output_gdal_format") == format
            @test DggridParams.get_parameter(params_point, "cell_output_type") == "NONE"
        end
    end
    
    # ========================================================================
    # Test prep_generate_grid_coarse_cells with point_output
    # ========================================================================
    
    @testset "prep_generate_grid_coarse_cells - cell output (point_output=false)" begin
        dggs_type = "ISEA7H"
        resolution = 5
        coarse_res = 2
        coarse_cells = string.([1, 2, 3])
        
        success, params, output_path = prep_generate_grid_coarse_cells(
            dggs_type, resolution, coarse_res, coarse_cells,
            output_format="GPKG",
            point_output=false
        )
        
        @test success == true
        @test typeof(params) == DggridParams.DGGRIDMetafile
        
        # Check that cell output is configured
        @test DggridParams.get_parameter(params, "cell_output_type") == "GDAL"
        @test DggridParams.get_parameter(params, "cell_output_gdal_format") == "GPKG"
        @test DggridParams.get_parameter(params, "cell_output_file_name") !== nothing
        
        # Check that point output is disabled
        @test DggridParams.get_parameter(params, "point_output_type") == "NONE"
    end
    
    @testset "prep_generate_grid_coarse_cells - point output (point_output=true)" begin
        dggs_type = "ISEA7H"
        resolution = 5
        coarse_res = 2
        coarse_cells = string.([1, 2, 3])
        
        success, params, output_path = prep_generate_grid_coarse_cells(
            dggs_type, resolution, coarse_res, coarse_cells,
            output_format="FlatGeobuf",
            point_output=true
        )
        
        @test success == true
        @test typeof(params) == DggridParams.DGGRIDMetafile
        
        # Check that point output is configured
        @test DggridParams.get_parameter(params, "point_output_type") == "GDAL"
        @test DggridParams.get_parameter(params, "point_output_gdal_format") == "FlatGeobuf"
        @test DggridParams.get_parameter(params, "point_output_file_name") !== nothing
        
        # Check that cell output is disabled
        @test DggridParams.get_parameter(params, "cell_output_type") == "NONE"
    end
    
    @testset "prep_generate_grid_coarse_cells - HIERNDX with point_output" begin
        dggs_type = "IGEO7"
        resolution = 5
        coarse_res = 2
        coarse_cells = ["0012", "0013", "0014"]
        
        # Test with cell output
        success_cell, params_cell, _ = prep_generate_grid_coarse_cells(
            dggs_type, resolution, coarse_res, coarse_cells,
            address_type="HIERNDX",
            point_output=false
        )
        @test success_cell == true
        @test DggridParams.get_parameter(params_cell, "cell_output_type") == "GDAL"
        @test DggridParams.get_parameter(params_cell, "point_output_type") == "NONE"
        
        # Test with point output
        success_point, params_point, _ = prep_generate_grid_coarse_cells(
            dggs_type, resolution, coarse_res, coarse_cells,
            address_type="HIERNDX",
            point_output=true
        )
        @test success_point == true
        @test DggridParams.get_parameter(params_point, "point_output_type") == "GDAL"
        @test DggridParams.get_parameter(params_point, "cell_output_type") == "NONE"
    end
    
    # ========================================================================
    # Test prep_generate_grid_clip_region with point_output
    # ========================================================================
    
    @testset "prep_generate_grid_clip_region - cell and point output" begin
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
            
            # Test with cell output (default)
            success_cell, params_cell, output_path_cell = prep_generate_grid_clip_region(
                dggs_type, resolution, test_clip_file,
                output_format="GPKG",
                point_output=false
            )
            
            @test success_cell == true
            @test DggridParams.get_parameter(params_cell, "cell_output_type") == "GDAL"
            @test DggridParams.get_parameter(params_cell, "cell_output_gdal_format") == "GPKG"
            @test DggridParams.get_parameter(params_cell, "point_output_type") == "NONE"
            
            # Test with point output
            success_point, params_point, output_path_point = prep_generate_grid_clip_region(
                dggs_type, resolution, test_clip_file,
                output_format="GeoJSON",
                point_output=true
            )
            
            @test success_point == true
            @test DggridParams.get_parameter(params_point, "point_output_type") == "GDAL"
            @test DggridParams.get_parameter(params_point, "point_output_gdal_format") == "GeoJSON"
            @test DggridParams.get_parameter(params_point, "cell_output_type") == "NONE"
            
        finally
            # Clean up test file
            rm(test_clip_file; force=true)
        end
    end
    
    @testset "prep_generate_grid_clip_region - point output with geodetic_densify" begin
        # Create a temporary test file
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
            success, params, _ = prep_generate_grid_clip_region(
                "IGEO7", 5, test_clip_file,
                geodetic_densify=1.0,
                point_output=true
            )
            
            @test success == true
            @test DggridParams.get_parameter(params, "geodetic_densify") == 1.0
            @test DggridParams.get_parameter(params, "point_output_type") == "GDAL"
            @test DggridParams.get_parameter(params, "cell_output_type") == "NONE"
            
        finally
            rm(test_clip_file; force=true)
        end
    end
    
    # ========================================================================
    # Test backward compatibility (default behavior)
    # ========================================================================
    
    @testset "Backward compatibility - default is cell output" begin
        # When point_output is not specified, it should default to false (cell output)
        
        # Test whole_earth
        success1, params1, _ = prep_generate_grid_whole_earth("ISEA7H", 3)
        @test success1 == true
        @test DggridParams.get_parameter(params1, "cell_output_type") == "GDAL"
        @test DggridParams.get_parameter(params1, "point_output_type") == "NONE"
        
        # Test coarse_cells
        success2, params2, _ = prep_generate_grid_coarse_cells(
            "ISEA7H", 5, 2, string.([1, 2])
        )
        @test success2 == true
        @test DggridParams.get_parameter(params2, "cell_output_type") == "GDAL"
        @test DggridParams.get_parameter(params2, "point_output_type") == "NONE"
        
        # Test clip_region
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
            success3, params3, _ = prep_generate_grid_clip_region(
                "ISEA7H", 4, test_clip_file
            )
            @test success3 == true
            @test DggridParams.get_parameter(params3, "cell_output_type") == "GDAL"
            @test DggridParams.get_parameter(params3, "point_output_type") == "NONE"
        finally
            rm(test_clip_file; force=true)
        end
    end
end

println("\n✓ All point_output parameter tests completed successfully!")
