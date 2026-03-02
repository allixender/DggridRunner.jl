```@meta
EditURL = "../../../examples/coarse_cells_workflow.jl"
```

# Coarse-Cells Sub-Grid with HIERNDX Addressing

This example generates the fine-resolution children of a single IGEO7 parent
cell, using **HIERNDX** (hierarchical index) cell addressing.

HIERNDX encodes the full path from the coarsest cell down to the target cell
as a digit string, e.g. `"0001611"` for a resolution-5 IGEO7 cell.
`grid_gen_convenience!` applies the standard extra parameters required for
HIERNDX grids: digit-string index format, vertex-0 position, and longitude
unwrap mode.

## Extra packages

```julia
using Pkg
Pkg.add(["GeoDataFrames", "GeoFormatTypes"])
```

```julia
using DGGRIDRunner
import DGGRIDRunner.AuthalicConversion: AuthalicToWGS84, transform_and_unwrap
import GeoDataFrames as GDF
import GeoFormatTypes
using Base.Threads
```

## Step 1 — Prepare parameters

We request IGEO7 cells at resolution 10 that fall inside the resolution-5
parent cell `"0001611"`.  Setting `address_type = "HIERNDX"` makes DGGRID
label every output cell with its hierarchical index.

```julia
parent_cell_id = "0001611"

success, params, output_path = prep_generate_grid_coarse_cells(
    "IGEO7", 10, 5, [parent_cell_id];
    output_format = "GPKG",
    address_type  = "HIERNDX",
    point_output  = false,
)
```

## Step 2 — Apply convenience settings

`grid_gen_convenience!` adds the parameters that HIERNDX grids require but
that the individual `prep_*` functions do not set automatically:
* `output_hier_ndx_form  = "digit_string"`
* `dggs_vert0_lon / lat` (vertex-0 of the icosahedron)
* `longitude_wrap_mode   = "UNWRAP_EAST"`

```julia
success && grid_gen_convenience!(params)
```

## Step 3 — Run DGGRID

```julia
success && run_dggrid_simple(params)
```

## Step 4 — Read and convert geometry

Like all IGEO7 output, the raw coordinates are in the authalic CRS.
We apply the fused transform-and-unwrap pass and parallelise with `@threads`.

```julia
gdf = GDF.read(output_path)
println("Sub-grid cells generated: $(nrow(gdf))")

f = AuthalicToWGS84()
geometry_wgs84 = Vector{Any}(undef, nrow(gdf))

@threads for i in eachindex(gdf.geometry)
    geometry_wgs84[i] = transform_and_unwrap(f, gdf.geometry[i])
end

gdf.geometry = geometry_wgs84
```

## Step 5 — Write the WGS84 result

```julia
out_path = replace(output_path, ".gpkg" => "_wgs84.gpkg")
GDF.write(out_path, gdf; crs = GeoFormatTypes.EPSG(4326))
println("WGS84 sub-grid written to: $out_path")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

