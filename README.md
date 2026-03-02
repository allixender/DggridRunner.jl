# DggridRunner

A Julia wrapper library to run highlevel functions of the DGGRID cli tool

[![Build Status](https://github.com/allixender/DggridRunner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/allixender/DggridRunner.jl/actions/workflows/CI.yml?query=branch%3Amain)

[DGGRID](https://www.discreteglobalgrids.org/software/) is a free software program for creating and manipulating Discrete Global Grids created and maintained by Kevin Sahr.

- [DGGRID Version 8.43 on GitHub](https://github.com/sahrk/DGGRID)
- [DGGRID User Manual](https://github.com/sahrk/DGGRID/blob/master/dggridManualV841.pdf)

[![Population Gridded](docs/src/day-04-hexa.png)](https://twitter.com/allixender/status/1324055326111485959)


## Inspiration

There is a very similar Python package with a longer history: [dggrid4py](https://github.com/allixender/dggrid4py). That package tries to abstract away the DGGRID parameters in order to give users an easier API. This [DggridRunner](https://github.com/allixender/DggridRunner.jl) Julia package goes down a different road and rather gives an easy access to better use the specific paramters for more fine-grained DGGRID usage.

## Notes

Technically, you need the `dggrid` tool compiled available on the system. However, we have a [Julia Yggdrasil BinaryBuilder](https://github.com/JuliaPackaging/Yggdrasil/tree/master/D/DGGRID7) package available, `DGGRID7_jll` that you can install.

If you need that for other purposes, you could something like that:

```julia
import DGGRID7_jll

DGGRID7_jll.dggrid() do dggrid_exec
    run(`$dggrid_exec -h`)
end
```

We internally use that already, so you are good to go