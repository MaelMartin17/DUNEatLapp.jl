# DUNEatLapp

This Julia package brings together the functions used at Lapp to study low-energy signals in the DUNE experiment.

## Installation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MaelMartin17.github.io/DUNEatLapp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MaelMartin17.github.io/DUNEatLapp.jl/dev/)
[![Build Status](https://github.com/MaelMartin17/DUNEatLapp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MaelMartin17/DUNEatLapp.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MaelMartin17/DUNEatLapp.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MaelMartin17/DUNEatLapp.jl)

```julia
Pkg.add("DUNEatLapp")
```
## Available functions

get_evts_index(df::DataFrame) : Provides indexes for each event in a Geant4 simulation 
