# md_simple

[![Build Status](https://github.com/viniciuspiccoli/md_simple.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/viniciuspiccoli/md_simple.jl/actions/workflows/CI.yml?query=branch%3Amain)

#### in progress

to run the code:


```julia
import Pkg
Pkg.add("https://github.com/viniciuspiccoli/md_simple.jl")

using md_simple, Plots

# type of atom
const argon = LennardJonesParticle(
    3.4e-10,      # m // sigma
    1.65e-21,     # J // epsilon
    6.63e-26      # kg // mass
)

# simulation object (reduced units) - Number of particles, density, temp...
sim = Simulation{TwoD}(
    substance = argon,
    nx = 20, ny = 20, nz = 2,  # nz is ignored in 2D
    density = 0.6,
    temp = 2.28,
    rcut = 5.0,
    dt = 0.001
)

# run NVT simulation
results = run_md!(sim, 100.0, 100, "argon_2d.xyz")
```
#### plotting the RDF

```julia
using LaTeXStrings
r, g = compute_rdf("argon_2d.xyz")


``` 