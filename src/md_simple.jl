module md_simple

    
    # packages required
    using LinearAlgebra, Random, Printf, ProgressMeter, Plots
    using StaticArrays  # performance and clarity
    
    include("properties.jl") # functions to calculate properties 
    include("xyz_parser.jl") # parser to read xyz file

    # main function to be exported
    export LennardJonesParticle, Simulation, TwoD, ThreeD, run_md!, Results, animate_simulation!


    #========================
    # 1. CORE DATA STRUCTURES
    ========================#

    # --- System Dimensionality ---
    abstract type SystemDimension end
    struct TwoD <: SystemDimension end
    struct ThreeD <: SystemDimension end

    # --- Particle Properties ---
    """
    Stores the Lennard-Jones parameters for a type of particle.
    """
    struct LennardJonesParticle
        σ::Float64       # LJ sigma (m)
        ϵ::Float64       # LJ epsilon (J)
        m::Float64       # Mass (kg)
    end




    # --- Simulation State and Parameters ---
    """
    A mutable struct to hold all information about the simulation.
    The dimension `D` is a type parameter for compile-time optimization.
    """
    mutable struct Simulation{D<:SystemDimension}
        # Particle properties
        substance::LennardJonesParticle

        # System state (dynamic)
        positions::Vector{SVector{3, Float64}}
        velocities::Vector{SVector{3, Float64}}
        forces::Vector{SVector{3, Float64}}

        # System properties (static)
        n_particles::Int
        box::SVector{3, Float64}

        # Simulation parameters
        temp_desired::Float64
        rcut::Float64
        dt::Float64

        function Simulation{D}(; 
                               substance::LennardJonesParticle,
                               nx::Int, ny::Int, nz::Int, 
                               density::Float64, 
                               temp::Float64, 
                               rcut::Float64, 
                               dt::Float64) where {D<:SystemDimension}

            dim = D()

            # Initialize positions and get system properties
            positions, box, n_particles = initialize_positions(dim, nx, ny, nz, density)

            # Initialize velocities
            velocities = Vector{SVector{3, Float64}}(undef, n_particles)
            initialize_velocities!(dim, velocities, temp, n_particles)

            # Initialize forces
            forces = zeros(SVector{3, Float64}, n_particles)
            compute_forces!(dim, forces, positions, box, rcut)

            new{D}(substance, positions, velocities, forces, n_particles, box, temp, rcut, dt)
        end
    end

    # --- Simulation Results ---
    """
    A struct to store the time-series data from the simulation.
    """
    struct Results
        time::Vector{Float64}
        temperature::Vector{Float64}
        potential_energy::Vector{Float64}
        kinetic_energy::Vector{Float64}
        total_energy::Vector{Float64}
    end


    #==============================
    # 2. INITIALIZATION FUNCTIONS
    ==============================#

    function initialize_positions(dim::ThreeD, nx, ny, nz, density)
        npart = 4 * nx * ny * nz
        vol = npart / density
        a0 = cbrt(vol / (nx * ny * nz))
        box = SVector(a0 * nx, a0 * ny, a0 * nz)

        positions = Vector{SVector{3, Float64}}(undef, npart)
        i = 0
        for iz in 1:2*nz, iy in 1:2*ny, ix in 1:2*nx
            if (ix + iy + iz) % 2 == 0
                i += 1
                x = 0.5 * a0 * (ix - 1 + 0.5 * ((iy + iz) % 2))
                y = 0.5 * a0 * (iy - 1 + 0.5 * ((ix + iz) % 2))
                z = 0.5 * a0 * (iz - 1 + 0.5 * ((ix + iy) % 2))
                positions[i] = SVector(x, y, z)
            end
        end
        return positions, box, npart
    end

    function initialize_positions(dim::TwoD, nx, ny, nz, density)
        npart = nx * ny
        area = npart / density
        a0 = sqrt(2 * area / (nx * ny * sqrt(3)))
        box = SVector(a0 * nx, a0 * ny * sqrt(3)/2, 0.0)

        positions = Vector{SVector{3, Float64}}(undef, npart)
        i = 0
        for iy in 1:ny, ix in 1:nx
            i += 1
            x = a0 * (ix - 1 + 0.5 * ((iy - 1) % 2))
            y = a0 * (iy - 1) * sqrt(3)/2
            positions[i] = SVector(x, y, 0.0)
        end
        return positions, box, npart
    end

    function initialize_velocities!(dim::D, velocities, temp, npart) where {D<:SystemDimension}
        # Maxwell-Boltzmann distribution
        sum_v = @SVector zeros(3)
        sum_v2 = 0.0
        for i in 1:npart
            v = D <: ThreeD ? randn(SVector{3, Float64}) : SVector(randn(), randn(), 0.0)
            velocities[i] = v
            sum_v += v
        end

        # Remove COM motion
        sum_v /= npart
        for i in 1:npart
            velocities[i] -= sum_v
            sum_v2 += dot(velocities[i], velocities[i])
        end

        # Scale to desired temperature
        dof = (D <: ThreeD ? 3 : 2) * npart - (D <: ThreeD ? 3 : 2)
        fs = sqrt(dof * temp / sum_v2)
        for i in 1:npart
            velocities[i] *= fs
        end
    end

    #==========================
    # 3. CORE SIMULATION LOGIC
    ==========================#

    @inline function minimum_image(dr_val, L_val)
        return dr_val - L_val * round(dr_val / L_val)
    end

    function compute_forces!(dim::D, forces, positions, box, rc) where {D<:SystemDimension}
        fill!(forces, @SVector zeros(3))
        en_pot = 0.0
        rc2 = rc^2
        npart = length(positions)

        for i in 1:npart-1
            for j in i+1:npart
                dr_vec = positions[i] - positions[j]
                dr = SVector(minimum_image(dr_vec[1], box[1]),
                             minimum_image(dr_vec[2], box[2]),
                             D <: ThreeD ? minimum_image(dr_vec[3], box[3]) : 0.0)

                r2 = dot(dr, dr)

                if r2 < rc2
                    r2i = 1.0 / r2
                    r6i = r2i^3
                    ff = 48.0 * r2i * r6i * (r6i - 0.5)

                    force_ij = ff * dr
                    forces[i] += force_ij
                    forces[j] -= force_ij

                    en_pot += 4.0 * r6i * (r6i - 1.0) # Add LJ potential tail correction for more accuracy
                end
            end
        end
        return en_pot
    end

    function integrate!(sim::Simulation{D}) where {D<:SystemDimension}
        dt2 = sim.dt / 2

        # First half-step velocity and full-step position update (mass m=1)
        for i in 1:sim.n_particles
            sim.velocities[i] += sim.forces[i] * dt2
            sim.positions[i] += sim.velocities[i] * sim.dt

            # Periodic boundary conditions
            p = sim.positions[i]
            b = sim.box
            sim.positions[i] = SVector(p[1] - b[1] * floor(p[1]/b[1] + 0.5),
                                       p[2] - b[2] * floor(p[2]/b[2] + 0.5),
                                       D <: ThreeD ? p[3] - b[3] * floor(p[3]/b[3] + 0.5) : 0.0)
        end

        # Compute new forces based on new positions
        en_pot = compute_forces!(D(), sim.forces, sim.positions, sim.box, sim.rcut)

        # Second half-step velocity update and KE calculation
        en_kin = 0.0
        for i in 1:sim.n_particles
            sim.velocities[i] += sim.forces[i] * dt2
            en_kin += 0.5 * dot(sim.velocities[i], sim.velocities[i]) # mass m=1
        end

        # Berendsen thermostat
        dof = (D <: ThreeD ? 3 : 2) * sim.n_particles - (D <: ThreeD ? 3 : 2)
        temp_inst = 2 * en_kin / dof
        lambda = sqrt(1 + (sim.dt / 0.1) * (sim.temp_desired / temp_inst - 1))

        for i in 1:sim.n_particles
            sim.velocities[i] *= lambda
        end
        en_kin *= lambda^2

        return en_pot, en_kin
    end

    #===========================#
    # 4. I/O & UNIT CONVERSIONS
    #===========================#

    # Define constants for unit conversion
    const K_B = 1.380649e-23  # Boltzmann constant, J/K
    const N_A = 6.02214076e23 # Avogadro's number, mol⁻¹

    # Conversion functions now take the substance properties as an argument
    real_time(t, substance::LennardJonesParticle) = t * sqrt(substance.m * substance.σ^2 / substance.ϵ) * 1e12 # ps
    real_temp(T, substance::LennardJonesParticle) = T * substance.ϵ / K_B # K
    real_energy(E, substance::LennardJonesParticle) = E * substance.ϵ * N_A / (4184) # kcal/mol

    function write_xyz_frame(io, sim::Simulation{D}, step::Int, time::Float64) where {D<:SystemDimension}
        σ_Å = sim.substance.σ * 1e10
        box_Å = sim.box .* σ_Å
        box_info = @sprintf("\"%.4f %.4f %.4f 90.0 90.0 90.0\"", box_Å[1], box_Å[2], box_Å[3])

        println(io, sim.n_particles)
        println(io, "Lattice=$box_info Step=$step Time=$(real_time(time, sim.substance))")

        for i in 1:sim.n_particles
            pos_Å = sim.positions[i] .* σ_Å
            @printf(io, "Ar %.4f %.4f %.4f\n", pos_Å[1], pos_Å[2], pos_Å[3])
        end
    end

    #=======================
    # 5. MAIN RUNNER FUNCTION
    =======================#

    function run_md!(sim::Simulation{D}, t_max::Float64, sample_interval::Int, traj_filename::String) where {D<:SystemDimension}
        println("Starting $(D()) MD simulation for Argon...")

        # Prepare results storage
        n_samples = floor(Int, t_max / (sim.dt * sample_interval))
        results = Results(zeros(n_samples), zeros(n_samples), zeros(n_samples), zeros(n_samples), zeros(n_samples))

        open(traj_filename, "w") do traj_file
            t = 0.0
            step = 0
            sample_count = 0

            while t < t_max
                # Run one integration step
                E_pot, E_kin = integrate!(sim)

                if step % sample_interval == 0
                    sample_count += 1

                    # 1. Write to trajectory file
                    write_xyz_frame(traj_file, sim, step, t)

                    # 2. Store thermodynamic data
                    dof = (D <: ThreeD ? 3 : 2) * sim.n_particles - (D <: ThreeD ? 3 : 2)
                    current_temp = 2 * E_kin / dof

                    results.time[sample_count] = real_time(t, sim.substance)
                    results.temperature[sample_count] = real_temp(current_temp, sim.substance)
                    results.potential_energy[sample_count] = real_energy(E_pot / sim.n_particles, sim.substance)
                    results.kinetic_energy[sample_count] = real_energy(E_kin / sim.n_particles, sim.substance)
                    results.total_energy[sample_count] = results.potential_energy[sample_count] + results.kinetic_energy[sample_count]

                    # 3. Print diagnostics to console
                    @printf("Step: %6d, Time: %7.2f ps, Temp: %6.2f K, E_tot: %8.4f kcal/mol\n",
                            step, results.time[sample_count], results.temperature[sample_count], results.total_energy[sample_count])
                end

                t += sim.dt
                step += 1
            end
        end

        println("Simulation complete. Trajectory saved to $traj_filename")
        return results
    end

    #================================
    # 6. LIVE ANIMATION FUNCTIONS
    ================================#

    """
        create_plot_layout(sim::Simulation)

    Initializes a complex plot layout for the animation.
    """
    function create_plot_layout(sim::Simulation{D}) where {D}
        # Get box size in real units (Ångströms)
        box_Å = sim.box .* sim.substance.σ .* 1e10

        # Plot 1: Simulation Box
        p1 = plot(layout=(1,1),
                  title="MD Simulation ($(D()))", xlabel="X (Å)", ylabel="Y (Å)",
                  xlims=(0, box_Å[1]), ylims=(0, box_Å[2]),
                  legend=false, framestyle=:box, aspect_ratio=:equal,
                  dpi=150)

        # Plot 2: Thermodynamics
        p2 = plot(title="System Thermodynamics",
                  xlabel="Time (ps)", ylabel="Energy (kcal/mol)",
                  legend=:bottomleft, framestyle=:box)

        # Twin axis on Plot 2 for Temperature
        p2_twin = plot!(twinx(), ylabel="Temperature (K)",
                        legend=:bottomright, framestyle=:box)

        # Combine the two plots into a single layout
        return plot(p1, p2, layout=(2,1), size=(800, 1000))
    end


    """
        update_plot!(plt, sim, step, t, history)

    Updates the plot object with the current simulation state.
    """
    function update_plot!(plt, sim::Simulation, step::Int, t::Float64, history::NamedTuple)
        # --- Update Simulation Box (Subplot 1) ---
        box_Å = sim.box .* sim.substance.σ .* 1e10

        # Using `SVector` makes this clean
        positions_Å = [p[1:2] .* sim.substance.σ * 1e10 for p in sim.positions]
        x_coords = [p[1] for p in positions_Å]
        y_coords = [p[2] for p in positions_Å]

        # Redraw particles. `clear=true` is inefficient for scatter plots.
        # It's better to just plot on top with a background.
        scatter!(plt[1], x_coords, y_coords,
                 marker=:circle, markersize=4, color=:blue,
                 xlims=(0, box_Å[1]), ylims=(0, box_Å[2]))

        # --- Update Thermodynamics (Subplot 2) ---
        plot!(plt[2], history.time, history.total_energy,
              linewidth=2, color=:dodgerblue, label="E_total")

        plot!(twinx(plt[2]), history.time, history.temperature,
              linewidth=2, color=:crimson, label="Temperature")
    end


    """
        animate_simulation!(sim, t_max, filename; frameskip=100, fps=15)

    Runs an MD simulation while generating a live animation and saving it as a GIF.
    """
    function animate_simulation!(sim::Simulation{D}, t_max::Float64, filename::String;
                                 frameskip::Int=100, fps::Int=15) where {D}

        println("Starting simulation with live animation...")

        # History is now local to this function, not global
        history = (time=Float64[], temperature=Float64[], total_energy=Float64[])

        # Setup
        steps = Int(floor(t_max / sim.dt))
        prog = Progress(steps, "Animating: ")
        plt = create_plot_layout(sim)

        # The @gif macro creates an Animation object and handles framing and saving
        anim = Plots.@gif for step in 0:steps
            # --- Run one step of MD ---
            E_pot, E_kin = integrate!(sim)

            # --- Update plot periodically ---
            if step % frameskip == 0
                # Calculate properties in real units
                dof = (D <: ThreeD ? 3 : 2) * sim.n_particles - (D <: ThreeD ? 3 : 2)
                current_temp_reduced = 2 * E_kin / dof

                T_K = real_temp(current_temp_reduced, sim.substance)
                E_kcal = real_energy((E_pot + E_kin) / sim.n_particles, sim.substance)
                t_ps = real_time(step * sim.dt, sim.substance)

                # Store history
                push!(history.time, t_ps)
                push!(history.temperature, T_K)
                push!(history.total_energy, E_kcal)

                # Update the plot with the latest data
                update_plot!(plt, sim, step, step * sim.dt, history)
            end
            next!(prog)
        end

        # The macro saves the animation once the loop is complete
        gif(anim, filename, fps=fps)
        println("Animation saved to $filename")
    end





end # end of module md_simple
