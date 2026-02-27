# DggridRunner

A Julia wrapper library to run highlevel functions of the DGGRID cli tool

[![Build Status](https://github.com/allixender/DggridRunner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/allixender/DggridRunner.jl/actions/workflows/CI.yml?query=branch%3Amain)

[DGGRID](https://www.discreteglobalgrids.org/software/) is a free software program for creating and manipulating Discrete Global Grids created and maintained by Kevin Sahr.

- [DGGRID Version 8.43 on GitHub](https://github.com/sahrk/DGGRID)
- [DGGRID User Manual](https://github.com/sahrk/DGGRID/blob/master/dggridManualV841.pdf)

[![Population Gridded](day-04-hexa.png)](https://twitter.com/allixender/status/1324055326111485959)


## Notes

Technically, you need the `dggrid` tool compiled available on the system. However, we have a [Julia Yggdrasil BinaryBuilder](https://github.com/JuliaPackaging/Yggdrasil/tree/master/D/DGGRID7) package available, `DGGRID7_jll` that you can install.

If you need that for other purposes, you could something like that:

```julia
import DGGRID7_jll

const libs_paths = DGGRID7_jll.LIBPATH_list
const dggrid_exec = DGGRID7_jll.get_dggrid_path()
```

We internally use that already, so ou are good to go