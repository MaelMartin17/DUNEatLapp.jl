# DUNEatLapp

This Julia package brings together the functions used at Lapp to study low-energy signals in the DUNE experiment.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MaelMartin17.github.io/DUNEatLapp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MaelMartin17.github.io/DUNEatLapp.jl/dev/)
[![Build Status](https://github.com/MaelMartin17/DUNEatLapp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MaelMartin17/DUNEatLapp.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MaelMartin17/DUNEatLapp.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MaelMartin17/DUNEatLapp.jl)

## Installation

```julia
Pkg.add("https://github.com/MaelMartin17/DUNEatLapp.jl")
```
## Description : 

clustering.jl => functions that cluster the deposits created after ULALAP simulation (https://github.com/lmanzanillas/ULALAP) and select some information about this clustering 

G4_connector.jl => functions that provide basic information about the ULALAP simulation

analyse_lardon.jl => functions that extract the data from a Lardon analysis

neutron_capture_info.jl => functions that give information about the neutron capture in the Liquid Argon, after ULALAP simulation

function.jl => Basic functions to apply as energy resolution, obtain the bin center, etc...



