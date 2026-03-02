```@meta
EditURL = "../../../examples/clip_cells_workflow.jl"
```

# Cell-Clipped Grid from a HIERNDX Cell-ID List

This example generates grid cells for a *specific set* of IGEO7 cells
identified by their HIERNDX addresses, using `prep_generate_grid_clip_cells`.

This is useful when you already know which cells you need — for instance when
you have computed relevant cell IDs from an upstream analysis and want to
materialise only those polygons rather than generating the entire grid.

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

## Step 1 — Write a cell-ID input file

`prep_generate_grid_clip_cells` expects the cell IDs as a plain-text file
with one ID per line.  Here we write a temporary file from a Julia vector.

```julia
cell_ids = ["0001615332", "0001615333", "0001615335"]

id_file = tempname() * ".txt"
write(id_file, join(cell_ids, "\n"))
```

## Step 2 — Prepare parameters

Both `input_address_type` and `output_address_type` are set to `"HIERNDX"` so
that DGGRID interprets the input IDs correctly and labels the output cells with
the same addressing scheme.

```julia
success, params, output_path = prep_generate_grid_clip_cells(
    "IGEO7", 8, id_file;
    output_format       = "GPKG",
    output_address_type = "HIERNDX",
    input_address_type  = "HIERNDX",
    point_output        = false,
)
```

## Step 3 — Apply convenience settings and run DGGRID

```julia
success && grid_gen_convenience!(params)
success && run_dggrid_simple(params)
```

The temporary cell-ID file is no longer needed once DGGRID has read it.

```julia
rm(id_file)
```

## Step 4 — Read and convert geometry

As with all IGEO7 output, coordinates are in the authalic CRS.
We apply the fused transform-and-unwrap pass in parallel.

```julia
gdf = GDF.read(output_path)
println("Clipped grid cells generated: $(nrow(gdf))")

f = AuthalicToWGS84()
geometry_wgs84 = Vector{Any}(undef, nrow(gdf))

@threads for i in eachindex(gdf.geometry)
    geometry_wgs84[i] = transform_and_unwrap(f, gdf.geometry[i])
end

gdf.geometry = geometry_wgs84
```

## Step 5 — Inspect results

For this small set of cells we simply print the converted geometries.
In a real workflow you would write to a file or pass `gdf` to a plotting function.

```julia
for i in eachindex(gdf.geometry)
    println("Cell $i ($(gdf[i, :name])): $(gdf.geometry[i])")
end
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

