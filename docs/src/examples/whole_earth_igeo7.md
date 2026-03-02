```@meta
EditURL = "../../../examples/whole_earth_igeo7.jl"
```

# Whole-Earth IGEO7 Grid: Generation and Authalic-to-WGS84 Conversion

This example generates a complete global IGEO7 grid at a given resolution and
converts the output geometry from the authalic CRS used internally by DGGRID
back to WGS84 (EPSG:4326) for use in standard GIS tooling.

**IGEO7** is an isometric equal-area DGGS.  Because DGGRID places cell vertices
on an authalic sphere, the raw output coordinates are *not* in WGS84 and must be
transformed before writing to a standard geospatial file.

## Extra packages

This example reads and writes spatial files with
[GeoDataFrames.jl](https://github.com/evetion/GeoDataFrames.jl) and
[GeoFormatTypes.jl](https://github.com/JuliaGeo/GeoFormatTypes.jl).
Install them once with:
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

## Step 1 — Prepare DGGRID parameters

`prep_generate_grid_whole_earth` validates the requested DGGS type and
resolution, builds a `DGGRIDMetafile`, and returns a temporary output path.
We use SEQNUM addressing here (the default), which assigns each cell a
sequential integer ID.

```julia
dggs_type  = "IGEO7"
resolution = 3

success, params, output_path = prep_generate_grid_whole_earth(
    dggs_type, resolution;
    output_format = "GPKG",
    address_type  = "SEQNUM",
)
```

## Step 2 — Run DGGRID

`run_dggrid_simple` writes the metafile to a temporary location, calls the
DGGRID binary bundled via `DGGRID7_jll`, and cleans up afterwards.
Output cells are written to `output_path`.

```julia
success && run_dggrid_simple(params)
```

## Step 3 — Read the raw output

The GeoPackage produced by DGGRID contains cell polygons whose coordinates
are in the authalic equatorial CRS, not in WGS84.

```julia
gdf = GDF.read(output_path)
println("Grid cells generated: $(nrow(gdf))")
```

## Step 4 — Convert authalic → WGS84

`AuthalicToWGS84` is a zero-allocation callable struct.
`transform_and_unwrap` applies the coordinate conversion *and* unwraps
longitudes across the ±180° dateline in a single pass per polygon — no
intermediate allocations.  We parallelise across cells with `@threads`.

```julia
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
println("WGS84 grid written to: $out_path")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

