# In a new file, e.g., "main.jl"
# Make sure you have added Plots and ProgressMeter to your Julia environment

#include("SimpleMD.jl")
using md_simple

# ================================
# 1. DEFINE THE PHYSICAL SYSTEM
# ================================
const argon = LennardJonesParticle(
    3.4e-10,      # m
    1.65e-21,     # J
    6.63e-26      # kg
)

# ================================
# 2. SET SIMULATION PARAMETERS
# ================================
sim = Simulation{TwoD}(
    substance = argon,
    nx = 20, ny = 20, nz = 2,
    density = 0.6,
    temp = 2.28,
    rcut = 5.0,
    dt = 0.001
)

# ================================
# 3. RUN WITH ANIMATION 
# ================================
animate_simulation!(sim, 50.0, "md_2d_argon.gif", frameskip=50, fps=20)
