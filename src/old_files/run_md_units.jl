
import Pkg
Pkg.activate("MD_simulation/")
using Plots, Parameters, Printf, Unitful, LaTeXStrings, Printf, LinearAlgebra, Random, StaticArrays


# Physical constants in proper units
const kB = 0.0019872041               # kcal/(mol·K)
const amu_to_kg = 1.66053906660e-27   # kg/amu
const kcal_to_J = 4184.0              # J/kcal
const angstrom_to_m = 1e-10
const ps_to_s = 1e-12
const NA = 6.02214076e23              # Avogadro's number
const J_to_kcal = 1/4184.0            # 1 J = 0.000239 kcal


########################################################## data structures ######################################

@with_kw struct AtomID
    name::String
    σ::Float64    # Å
    ε::Float64    # kcal/mol
    mass::Float64 # amu
end

@with_kw mutable struct MDinput
    box_size::Vector{Float64}  # Å
    temperature::Float64       # K
    total_time::Float64        # ps
    dt::Float64                # ps
    cutoff::Float64 = 12.0     # Å
    nsteps::Int = trunc(Int, total_time/dt)
end

@with_kw mutable struct AtomState
    id::AtomID
    r::Vector{Float64}  # Å
    v::Vector{Float64}  # Å/ps
    f::Vector{Float64}  # kcal/(mol·Å)
    f_prev::Vector{Float64} = zeros(length(r))
end

@with_kw mutable struct Simulation
    atoms::Vector{AtomState}
    input::MDinput
    current_step::Int = 0
    time::Float64 = 0.0
    trajectory::Vector{Vector{Tuple{String,Vector{Float64}}}} = []
    log_interval::Int = 100
end


##########################################################3 initial point ######################################

function setRectangular(input::MDinput, atoms::Vector{AtomID}, dims::Vector{Int}; noise=0.1)
    positions = []
    spacing = input.box_size ./ dims
    
    particle = 1
    for i in CartesianIndices(Tuple(dims))
        pos = [spacing[d] * (i[d] - 0.5) for d in 1:length(input.box_size)]
        
        # Add noise if requested
        if noise > 0
            pos .+= noise .* randn(length(pos))
        end
        
        push!(positions, (pos, atoms[mod1(particle, length(atoms))]))
        particle += 1
    end
    return positions
end

##########################################################3 minimum image ######################################
function minimum_image(rx, L)
    if rx >= L/2
        rx -= L
    elseif rx < -L/2
        rx += L
    end
    return rx
end


function apply_minimum_image(dr, box)
    for d in eachindex(dr)
        dr[d] = minimum_image(dr[d], box[d])
    end
    return dr
end



############################################################# potential energy ##############################################

function pair_lj(r, atom1::AtomID, atom2::AtomID)
  #= Lorentz-Berthelot mixing rules - it will reduce to 1 for the same atom
     even if this is terrible from the performance point of view, it will be done this way for the first version.
     the improvements of the code will be also used as a way to study more MD, since the used programs are quite well optimized, being
     complicated to understand all steps... 
  =#
  σ_mix = (atom1.σ + atom2.σ)/2
  ε_mix = sqrt(atom1.ε * atom2.ε)
  # lj calculation
  σ6 = (σ_mix)^6
  σ12 = (σ6)^2
  inv_r6 = (1/r)^6
  inv_r12 = inv_r6^2

  return 4.0 * ε_mix * (σ12 * inv_r12 - σ6 * inv_r6)
end


function utotal_lj(positions, md::MDinput)
  energy = 0.0
  n = length(positions)
  cutoff_sq = md.cutoff^2

  for i in 1:n
      (pos_i, atom_i) = positions[i]
      for j in (i+1):n
          (pos_j, atom_j) = positions[j]

          rij = pos_i .- pos_j # calculate the displacement vector
          rij = apply_minimum_image(rij, md.box_size) # apply minimum_image in all dimensions
          r2 = sum(rij.^2) # distance between two points!

          # check if r2 is smaller than the cutoff and avoid self-interaction
          if r2 < cutoff_sq && r2 > 0.0
            energy += pair_lj(r2,atom_i, atom_j)
          end

      end
  end
  return energy
end



########################################################## set velocities ######################################

function setVelocitiesNOMASS(N, T)
    vxSum = 0.0
    vySum = 0.0
    vx = zeros(N)
    vy = zeros(N) 
   
    # Random initial velocities
    for i in 1:N
        vx[i] = rand() - 0.5
        vy[i] = rand() - 0.5
        vxSum += vx[i]
        vySum += vy[i]
    end

    # Remove center-of-mass velocity
    vxcm = vxSum / N
    vycm = vySum / N
    vx .-= vxcm
    vy .-= vycm

    # Calculate current kinetic energy
    v2sum = sum(vx.^2 + vy.^2)
    current_KE_per_particle = 0.5 * v2sum / N

    # Target kinetic energy per particle (from temperature)
    kB = 0.0019872041  # kcal/(mol·K)
    target_KE_per_particle = kB * T  # 2D: KE_per_particle = kBT

    # Rescale velocities to match target temperature
    rescale = sqrt(target_KE_per_particle / current_KE_per_particle)
    vx .*= rescale
    vy .*= rescale

    return vx, vy
end


#=
function setVelocities(N, T, mass_amu)
    # Convert mass from amu to simulation units (kcal·ps²/(Å²·mol))
    amu_to_kg = 1.66053906660e-27
    J_to_kcal = 1/4184.0
    NA = 6.02214076e23
    mass = mass_amu * amu_to_kg * NA * 1e4 * J_to_kcal

    vxSum = 0.0
    vySum = 0.0
    vx = zeros(N)
    vy = zeros(N)
    
    # Random initial velocities
    for i in 1:N
        vx[i] = rand() - 0.5
        vy[i] = rand() - 0.5
        vxSum += vx[i]
        vySum += vy[i]
    end

    # Remove center-of-mass velocity
    vxcm = vxSum / N
    vycm = vySum / N
    vx .-= vxcm
    vy .-= vycm

    # Calculate current kinetic energy with mass
    v2sum = sum(vx.^2 + vy.^2)
    current_KE_per_particle = 0.5 * mass * v2sum / N

    # Target kinetic energy from temperature (2D: KE = k_B*T)
    kB = 0.0019872041  # kcal/(mol·K)
    target_KE_per_particle = kB * T

    # Rescale velocities
    rescale = sqrt(target_KE_per_particle / current_KE_per_particle)
    vx .*= rescale
    vy .*= rescale

    return vx, vy
end



# aqui precisa ajsutar para servir poara duas particulas diferetnes
function initialize_velocities!(sim::Simulation)
    T = sim.input.temperature
    N = length(sim.atoms)
    vx, vy = setVelocities(N, T, sim.atoms[1].id.mass)
    for (i,atom) in enumerate(sim.atoms)
        atom.v = [vx[i], vy[i]]
    end
end
=#





function set_velocities!(sim::Simulation)
    kb = 0.0019872041  # kcal/mol/K
    dof = 2 * length(atoms)  # 2D
    rng = MersenneTwister(42)

    velocities = [randn(rng, SVector{2, Float64}) for _ in atoms]
    for i in eachindex(atoms)
        scale = sqrt(kb * sim.input.temperature / atoms[i].id.mass)
        velocities[i] *= scale
    end

    # Remover momento linear
    total_momentum = sum([sim.atoms[i].id.mass * velocities[i] for i in eachindex(sim.atoms)])
    total_mass = sum(a.id.mass for a in atoms)
    v_cm = total_momentum / total_mass

    for i in eachindex(sim.atoms)
        sim.atoms[i].v = velocities[i] - v_cm
    end
end





########################################################## force calculation ######################################

#distance between two points
function dist(x::Vector{Real},y::Vector{Real})
    sum = 0.0
    for i in 1:length(x)
        sum = sum + (x[i] - y[i])^2
    end
    return sqrt(sum)
end



function fpair(atom1::AtomState, atom2::AtomState, dr, r)
    σ = (atom1.id.σ + atom2.id.σ)/2
    ε = sqrt(atom1.id.ε * atom2.id.ε)  # Fix typo: "aatom2" → "atom2"

    r6 = r^6
    r7 = r6 * r
    r12 = r6^2
    r13 = r12 * r

    sig6 = σ^6
    sig12 = sig6^2
    eps4 = 4 * ε

    dudr1 = -12 * (sig12 / r13)   # -12σ¹²/r¹³
    dudr2 = -6 * (sig6 / r7)      # -6σ⁶/r⁷
    dfacdr = eps4 * (dudr1 - dudr2)  # 4ε*(-12σ¹²/r¹³ + 6σ⁶/r⁷)

    # Unit vector components from corrected displacement (dr)
    drdx = dr[1] / r
    drdy = dr[2] / r

    # Force on atom1 due to atom2
    f = -dfacdr .* [drdx, drdy]
    return f
end


function compute_forces!(atoms::Vector{AtomState}, md::MDinput)

    #n = length(atoms)
    #forces = [SVector(0.0, 0.0) for _ in 1:n]
    #potential_energy = 0.0

    for atom in atoms
        atom.f .= 0.0
    end
    
    rc2 = md.cutoff^2
    n = length(atoms)

    for i in 1:n-1
        for j in (i+1):n
            # Compute displacement with minimum image convention
            r_vec = atoms[j].r  -  atoms[i].r 
            r_vec = apply_minimum_image(r_vec, md.box_size)
            r2 = dot(r_vec, r_vec)

            if r2 < rc2
                σ = (atoms[i].id.σ + atoms[j].id.σ)/2
                ε = sqrt(atoms[i].id.ε * atoms[j].id.ε) 
                
                sig2 = σ^2
                r6 = (sig2 / r2)^3
                r12 = r6^2

                f_mag = 48 * ε * (r12 - 0.5 * r6) / r2
                f_vec = f_mag * r_vec

                atoms[i].f += f_vec
                atoms[j].f -= f_vec

                # Energia corrigida com shifting
                #rc6 = (sig / md.cutoff)^6
                # rc12 = rc6^2
                #U_rc = 4 * eps * (rc12 - rc6)
                # U = 4 * eps * (r12 - r6) - U_rc
                #potential_energy += U
            end
        end
    end

end



########################################################## velocity verlet integration ######################################

function verlet_step!(sim::Simulation, md::MDinput)
    dt = sim.input.dt
    
    J_to_kcal = 1/4184.0
    amu_to_kg = 1.66053906660e-27
    NA = 6.02214076e23

    # Step 1: Update positions and apply periodic boundaries
    #previous = [atom.v for atom in sim.atoms]
    for atom in sim.atoms

        mass = atom.id.mass * amu_to_kg * NA * 1e4 * J_to_kcal

        # Update positions
        @. atom.r = atom.r + atom.v * dt + atom.f * (dt*dt / (2 * mass))

        #@. atom.r = minimum_image(atom.r, md.box_size)
        #@. atom.r = mod(atom.r, md.box_size)

        @. atom.v = atom.v +  atom.f *( dt / mass)
    end
    # Step 2: Compute new forces

    compute_forces!(sim.atoms, md)

    # Step 3: Complete velocity update

    for atom in sim.atoms
        mass = atom.id.mass  * amu_to_kg * NA * 1e4 * J_to_kcal
        atom.v .+= 0.5 .* atom.f ./ mass .* dt
    end

    sim.time += dt
    sim.current_step += 1
end




function get_positions(sim::Simulation)
    return [(atom.r, atom.id) for atom in sim.atoms]
end


function energy_calc(sim::Simulation, md::MDinput)
    J_to_kcal = 1/4184.0
    amu_to_kg = 1.66053906660e-27
    NA = 6.02214076e23

    pos = get_positions(sim)

    pe = utotal_lj(pos, md)
    
    ke = 0.0
    for atom in sim.atoms
        mass = atom.id.mass * amu_to_kg * NA * 1e4 * J_to_kcal
        ke += 0.5 * mass * sum(atom.v.^2)
    end
    
    return pe + ke, pe, ke
end



function temperature(sim::Simulation)
    kB = 0.0019872041  # kcal/(mol·K)
    amu_to_kg = 1.66053906660e-27
    J_to_kcal = 1/4184.0
    NA = 6.02214076e23

    ke = 0.0  # Initialize kinetic energy
    for atom in sim.atoms
        # Proper mass conversion (kcal·ps²/(Å²·mol))
        mass = atom.id.mass * amu_to_kg * NA * 1e4 * J_to_kcal
        ke += 0.5 * mass * sum(atom.v.^2)
    end
    
    dof = 2 * length(sim.atoms) - 2  # 2D degrees of freedom
    T = (2 * ke) / (dof * kB)
    return T
end


#=
function temperature(sim::Simulation)     # Wrap positions into the box
    sim.atoms.r .= mod.(sim.atoms.r, md.box_size)
    # Half-step velocity update
    sim.atoms.v .+= 0.5 .* sim.atoms.f ./ mass .* dt
    for atom in sim.atoms
        # Proper mass conversion (kcal·ps²/(Å²·mol))
        mass = atom.id.mass * amu_to_kg * NA * 1e4 * J_to_kcal
        ke += 0.5 * mass * sum(atom.v.^2)
    enforce_computationd
    
    dof = 2*length(sim.atoms) - 2  # 2D degrees of freedom
    (2 * ke) / (dof * kB)
end

=#

function save_trajectory!(sim::Simulation)
    frame = Tuple{String,Vector{Float64}}[]
    for atom in sim.atoms
        # Convert back to real units for storage
        pos_real = atom.r
        push!(frame, (atom.id.name, pos_real))
    end
    push!(sim.trajectory, frame)
end


function write_xyz(sim::Simulation, filename="traj.xyz")
    open(filename, "w") do io  # Changed from 'as' to 'do'
        for frame in sim.trajectory
            write(io, "$(length(frame))\n")  # Fixed missing parenthesis
            write(io, "Step $(sim.current_step), Time $(sim.time)\n")
            for (name, pos) in frame
                # Ensure 3D coordinates (add 0.0 for Z if 2D)
                x = pos[1]  
                y = pos[2] 
                write(io, "$name $x $y 0.0\n")
            end
        end
    end  # Added closing 'end' for the 'do' block
end

function run!(sim::Simulation, md::MDinput)
    # Initial velocities
    set_velocities!(sim)  

    # Initial forces
    compute_forces!(sim.atoms, md)
     
    while sim.current_step < sim.input.nsteps
        verlet_step!(sim, md)
        
        # Save trajectory and log datrun!(sim)a
        if sim.current_step % sim.log_interval == 0
            save_trajectory!(sim)
            
            T = temperature(sim)
            E_total, E_pot, E_kin = energy_calc(sim, md)
            
            @printf("Step %6d: T = %.2f K, E_total = %.3f, E_pot = %.9f, E_kin = %.3f\n",
                   sim.current_step,
                   T,  
                   E_total, E_pot, E_kin)

        end
    end
    
    write_xyz(sim)
end


## test
ar = AtomID(name="Ar", σ=3.4, ε=0.238, mass=39.948)

md = MDinput(
    box_size=[80.0, 80.0],  # 5σ box
    temperature=100.0,
    total_time=1.0,         # 10 fs
    dt=0.0001,              # 0.1 fs steps
    cutoff=8.0)

positions = setRectangular(md, [ar], [8, 8], noise=0.01)
atoms = [AtomState(id=atom, r=pos, v=zeros(2), f=zeros(2)) for (pos, atom) in positions]
sim = Simulation(atoms=atoms, input=md, log_interval=10)
run!(sim, md)