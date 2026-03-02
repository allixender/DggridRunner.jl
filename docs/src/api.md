# API Reference

## DGGRIDRunner

```@docs
prep_output_stats
parse_output_stats
run_output_stats
run_dggrid_simple
prep_generate_grid_whole_earth
prep_generate_grid_coarse_cells
prep_generate_grid_clip_region
prep_generate_grid_clip_cells
grid_gen_convenience!
```

## DGGRIDParams

```@docs
DGGRIDParams.DGGRIDParams
DGGRIDParams.DGGRIDMetafile
DGGRIDParams.Parameter
DGGRIDParams.create_metafile
DGGRIDParams.add_parameter!
DGGRIDParams.get_parameter
DGGRIDParams.remove_parameter!
DGGRIDParams.validate_metafile
DGGRIDParams.write_metafile
DGGRIDParams.read_metafile
DGGRIDParams.list_parameters
DGGRIDParams.get_parameter_info
```

## AuthalicConversion

```@docs
AuthalicConversion.AuthalicToWGS84
AuthalicConversion.WGS84ToAuthalic
AuthalicConversion.transform_and_unwrap
AuthalicConversion.transform_and_unwrap_point
AuthalicConversion.unwrap_polygon_lon!
AuthalicConversion.unwrap_polygon
```
