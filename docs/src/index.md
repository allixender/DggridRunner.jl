# DGGRIDRunner.jl

A Julia wrapper for the [DGGRID](https://discreteglobalgrids.org/) CLI tool,
providing high-level functions to create and query Discrete Global Grid Systems (DGGS).

## Installation

```julia
using Pkg
Pkg.add("DGGRIDRunner")
```

## Quick Start

```julia
using DGGRIDRunner

# Generate an ISEA7H grid at resolution 3 (whole earth)
success, params, output_path = prep_generate_grid_whole_earth("ISEA7H", 3)
run_dggrid_simple(params)
# → output_path now points to a GeoPackage with the grid cells
```
