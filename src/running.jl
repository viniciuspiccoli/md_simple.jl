# run the simulations

using md_simple, Plots

const argon = LennardJonesParticle(
    3.4e-10,      # m // sigma
    1.65e-21,     # J // epsilon
    6.63e-26      # kg // mass
)

sim = Simulation{TwoD}(
    substance = argon,
    nx = 20, ny = 20, nz = 2,  # nz is ignored in 2D
    density = 0.6,
    temp = 2.28,
    rcut = 5.0,
    dt = 0.001
)

results = run_md!(sim, 100.0, 100, "argon_2d.xyz")


#p1 = plot(results.time, results.temperature,
#          xlabel="Time (ps)", ylabel="Temperature (K)",
#          title="Temperature vs. Time", legend=false)

#p2 = plot(results.time, results.total_energy,
#          xlabel="Time (ps)", ylabel="Total Energy (kcal/mol)",
#          title="Total Energy vs. Time", legend=false)

#plot(p1, p2, layout=(2,1), size=(800, 600))

# Assume `trajectory`==results is the object returned by run_md!
# It should contain:
# trajectory.positions -> Vector of position matrices (Å)
# trajectory.velocities -> Vector of velocity matrices (Å/ps) -- YOU NEED TO ADD THIS
# trajectory.box -> The simulation box vector (Å)
# trajectory.total_energy -> Vector of total energy (kcal/mol)
# trajectory.temperature -> Vector of temperature (K)

# --- Calculate g(r) ---
r, g = compute_rdf("argon_2d.xyz")

using LaTeXStrings

plot(r, g, xlabel="r (Å)", ylabel="g(r)", title="Radial Distribution Function", legend=false)
savefig("gr.png")


#=
# --- Calculate Coordination & Structure ---
n_particles = size(trajectory.positions[1], 1)
volume = prod(trajectory.box[1])
ρ = n_particles / volume

n_r = running_coordination(r, g, ρ)
k, S_k = structure_factor(r, g, ρ)

plot(r, n_r, xlabel="r (Å)", ylabel="n(r)", title="Running Coordination Number", legend=false)
plot(k, S_k, xlabel="k (Å⁻¹)", ylabel="S(k)", title="Structure Factor", legend=false)


# --- Calculate Cv ---
# Make sure your trajectory object from the simulation run contains the total energy
Cv = compute_Cv(trajectory.total_energy, trajectory.temperature, n_particles)

=#
