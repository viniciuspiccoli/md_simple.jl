# Assume `trajectory` is the object returned by run_md!
# It should contain:
# trajectory.positions -> Vector of position matrices (Å)
# trajectory.velocities -> Vector of velocity matrices (Å/ps) -- YOU NEED TO ADD THIS
# trajectory.box -> The simulation box vector (Å)
# trajectory.total_energy -> Vector of total energy (kcal/mol)
# trajectory.temperature -> Vector of temperature (K)

# --- Calculate g(r) ---
r, g = compute_rdf(trajectory.positions, trajectory.box[1])
plot(r, g, xlabel="r (Å)", ylabel="g(r)", title="Radial Distribution Function", legend=false)

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
