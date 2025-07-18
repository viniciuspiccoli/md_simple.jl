using Random
using LinearAlgebra

# Module to hold helper functions (equivalent to Fortran 'routines' module)
#module MD_Routines
#export read_input, read_natoms, read_positions, write_positions, write_final_positions,
#       write_statistics, randomize_velocities!, compute_engkin, pbc, check_list,
#       compute_list, compute_forces, thermostat!

# Structure to hold input parameters (similar to Fortran's subroutine arguments)
struct MDParams
    temperature::Float64
    tstep::Float64
    friction::Float64
    forcecutoff::Float64
    listcutoff::Float64
    nstep::Int
    nconfig::Int
    nstat::Int
    wrapatoms::Bool
    maxneighbours::Int
    inputfile::String
    outputfile::String
    trajfile::String
    statfile::String
    idum::Int
end

# Read input parameters from standard input
function read_input()
    # Default values
    temperature = 1.0
    tstep = 0.005
    friction = 0.0
    forcecutoff = 2.5
    listcutoff = 3.0
    nstep = 1
    nconfig = 10
    nstat = 1
    maxneighbours = 1000
    idum = 0
    wrapatoms = false
    statfile = ""
    trajfile = ""
    outputfile = ""
    inputfile = ""

    for line in eachline(stdin)
        # Remove comments and trim whitespace
        line = strip(split(line, '#')[1])
        isempty(line) && continue
        parts = split(line)
        keyword = parts[1]
        value = parts[2:end]
        
        # Parse each keyword
        if keyword == "temperature"
            temperature = parse(Float64, value[1])
        elseif keyword == "tstep"
            tstep = parse(Float64, value[1])
        elseif keyword == "friction"
            friction = parse(Float64, value[1])
        elseif keyword == "forcecutoff"
            forcecutoff = parse(Float64, value[1])
        elseif keyword == "listcutoff"
            listcutoff = parse(Float64, value[1])
        elseif keyword == "nstep"
            nstep = parse(Int, value[1])
        elseif keyword == "nconfig"
            nconfig = parse(Int, value[1])
            trajfile = value[2]
        elseif keyword == "nstat"
            nstat = parse(Int, value[1])
            statfile = value[2]
        elseif keyword == "maxneighbours"
            maxneighbours = parse(Int, value[1])
        elseif keyword == "wrapatoms"
            wrapatoms = parse(Bool, value[1])
        elseif keyword == "inputfile"
            inputfile = value[1]
        elseif keyword == "outputfile"
            outputfile = value[1]
        elseif keyword == "seed"
            idum = -parse(Int, value[1]) # Convert to negative for RNG
        else
            error("Unknown keyword: $keyword")
        end
    end

    # Validate required parameters
    isempty(inputfile) && error("Input file not specified")
    isempty(outputfile) && error("Output file not specified")
    isempty(trajfile) && error("Trajectory file not specified")
    isempty(statfile) && error("Statistics file not specified")

    MDParams(temperature, tstep, friction, forcecutoff, listcutoff, nstep, nconfig, nstat,
             wrapatoms, maxneighbours, inputfile, outputfile, trajfile, statfile, idum)
end

# Read number of atoms from XYZ file
function read_natoms(inputfile::String)
    open(inputfile) do f
        parse(Int, readline(f))
    end
end

# Read atom positions and cell from XYZ file
function read_positions(inputfile::String, natoms::Int)
    positions = zeros(3, natoms)
    cell = zeros(3)
    open(inputfile) do f
        readline(f)  # Skip natoms line
        cell .= parse.(Float64, split(readline(f)))
        for i in 1:natoms
            parts = split(readline(f))
            positions[:, i] .= parse.(Float64, parts[2:4])
        end
    end
    positions, cell
end

# Write positions to XYZ file (appends if file exists)
function write_positions(trajfile::String, natoms::Int, positions::Matrix{Float64},
                         cell::Vector{Float64}, wrapatoms::Bool)
    # Use a closure to track if first write
    first_write = Ref(true)
    mode = first_write[] ? "w" : "a"
    open(trajfile, mode) do f
        println(f, natoms)
        println(f, join(cell, " "))
        for i in 1:natoms
            pos = wrapatoms ? pbc(cell, positions[:, i]) : positions[:, i]
            println(f, "Ar ", join(pos, " "))
        end
    end
    first_write[] = false
end

# Write final positions to XYZ file
function write_final_positions(outputfile::String, natoms::Int, positions::Matrix{Float64},
                               cell::Vector{Float64}, wrapatoms::Bool)
    open(outputfile, "w") do f
        println(f, natoms)
        println(f, join(cell, " "))
        for i in 1:natoms
            pos = wrapatoms ? pbc(cell, positions[:, i]) : positions[:, i]
            println(f, "Ar ", join(pos, " "))
        end
    end
end

# Write statistics to file (appends and flushes periodically)
function write_statistics(statfile::String, istep::Int, tstep::Float64, natoms::Int,
                          engkin::Float64, engconf::Float64, engint::Float64)
    # Use closure to track last reopen time
    last_reopen = Ref(0)
    if istep - last_reopen[] > 100
        close(f) # Close and reopen to flush buffer
        open(statfile, "a") do f
            println(f, istep, " ", istep*tstep, " ", 2engkin/(3natoms), " ",
                    engconf, " ", engkin+engconf, " ", engkin+engconf+engint)
        end
        last_reopen[] = istep
    else
        open(statfile, "a") do f
            println(f, istep, " ", istep*tstep, " ", 2engkin/(3natoms), " ",
                    engconf, " ", engkin+engconf, " ", engkin+engconf+engint)
        end
    end
end

# Initialize velocities with Maxwell-Boltzmann distribution
function randomize_velocities!(velocities::Matrix{Float64}, temperature::Float64,
                               masses::Vector{Float64}, rng)
    for i in 1:size(velocities, 2)
        σ = sqrt(temperature / masses[i])
        velocities[:, i] .= σ .* randn(rng, 3)
    end
end

# Calculate kinetic energy
function compute_engkin(masses::Vector{Float64}, velocities::Matrix{Float64})
    sum(0.5 * masses[i] * sum(velocities[:, i].^2) for i in 1:size(velocities, 2))
end

# Apply periodic boundary conditions
function pbc(cell::Vector{Float64}, vin::Vector{Float64})
    vout = similar(vin)
    @. vout = vin - round(vin / cell) * cell
    vout
end

# Check if neighbor list needs updating
function check_list(positions::Matrix{Float64}, positions0::Matrix{Float64},
                    listcutoff::Float64, forcecutoff::Float64)
    δ = 0.5 * (listcutoff - forcecutoff)
    δ² = δ^2
    for i in 1:size(positions, 2)
        if sum((positions[:, i] - positions0[:, i]).^2) > δ²
            return true
        end
    end
    false
end

# Build neighbor list using Verlet algorithm
function compute_list(natoms::Int, listsize::Int, positions::Matrix{Float64},
                      cell::Vector{Float64}, listcutoff::Float64)
    listcutoff² = listcutoff^2
    list = zeros(Int, listsize)
    point = zeros(Int, natoms + 1)
    point[1] = 1
    
    for iatom in 1:natoms-1
        point[iatom+1] = point[iatom]
        for jatom in iatom+1:natoms
            Δ = pbc(cell, positions[:, iatom] - positions[:, jatom])
            if sum(Δ.^2) > listcutoff²
                continue
            end
            if point[iatom+1] > listsize
                error("Verlet list overflow")
            end
            list[point[iatom+1]] = jatom
            point[iatom+1] += 1
        end
    end
    point, list
end

# Compute forces using neighbor list
function compute_forces(natoms::Int, listsize::Int, positions::Matrix{Float64},
                        cell::Vector{Float64}, forcecutoff::Float64,
                        point::Vector{Int}, list::Vector{Int})
    forcecutoff² = forcecutoff^2
    forces = zeros(3, natoms)
    engconf = 0.0
    engcorrection = 4.0 * (1.0/forcecutoff²^6 - 1.0/forcecutoff²^3)
    
    for iatom in 1:natoms-1
        for jlist in point[iatom]:point[iatom+1]-1
            jatom = list[jlist]
            Δ = pbc(cell, positions[:, iatom] - positions[:, jatom])
            Δ² = sum(Δ.^2)
            Δ² > forcecutoff² && continue
            
            inv_Δ² = 1.0 / Δ²
            inv_Δ⁶ = inv_Δ²^3
            inv_Δ¹² = inv_Δ⁶^2
            engconf += 4.0 * (inv_Δ¹² - inv_Δ⁶) - engcorrection
            
            # Compute force
            f = 24.0 * Δ * (2.0 * inv_Δ¹² - inv_Δ⁶) * inv_Δ²
            forces[:, iatom] .+= f
            forces[:, jatom] .-= f
        end
    end
    forces, engconf
end

# Apply Langevin thermostat
function thermostat!(velocities::Matrix{Float64}, masses::Vector{Float64},
                     dt::Float64, γ::Float64, T::Float64, engint::Vector{Float64},
                     rng)
    c1 = exp(-γ * dt)
    for i in 1:size(velocities, 2)
        mass = masses[i]
        c2 = sqrt((1 - c1^2) * T / mass)
        for d in 1:3
            engint[1] += 0.5 * mass * velocities[d, i]^2
            velocities[d, i] = c1 * velocities[d, i] + c2 * randn(rng)
            engint[1] -= 0.5 * mass * velocities[d, i]^2
        end
    end
end

#end # module MD_Routines

# Main MD simulation program
function main()
    #import .MD_Routines
    
    # Read input parameters
    params = read_input()
    
    # Read system information
    natoms = read_natoms(params.inputfile)
    positions, cell = read_positions(params.inputfile, natoms)
    
    # Allocate arrays
    velocities = zeros(3, natoms)
    masses = ones(natoms)  # Masses set to 1.0 as in original code
    listsize = params.maxneighbours * natoms ÷ 2
    list = zeros(Int, listsize)
    point = zeros(Int, natoms + 1)
    positions0 = copy(positions)
    forces = zeros(3, natoms)
    engint = [0.0]  # Mutable container for energy integral
    
    # Initialize RNG
    rng = MersenneTwister(params.idum)
    
    # Randomize velocities
    randomize_velocities!(velocities, params.temperature, masses, rng)
    
    # Build initial neighbor list and compute forces
    point, list = compute_list(natoms, listsize, positions, cell, params.listcutoff)
    forces, engconf = compute_forces(natoms, listsize, positions, cell,
                                     params.forcecutoff, point, list)
    positions0 .= positions
    
    # Main MD loop
    for istep in 1:params.nstep
        # First thermostat half-step
        thermostat!(velocities, masses, 0.5params.tstep, params.friction,
                   params.temperature, engint, rng)
        
        # Velocity Verlet: update velocities (1st half) and positions
        @. velocities += forces * 0.5params.tstep / masses'
        @. positions += velocities * params.tstep
        
        # Check neighbor list and recompute if needed
        if check_list(positions, positions0, params.listcutoff, params.forcecutoff)
            point, list = compute_list(natoms, listsize, positions, cell, params.listcutoff)
            positions0 .= positions
            println("Neighbor list updated at step $istep")
        end
        
        # Compute new forces
        forces, engconf = compute_forces(natoms, listsize, positions, cell,
                                         params.forcecutoff, point, list)
        
        # Velocity Verlet: update velocities (2nd half)
        @. velocities += forces * 0.5params.tstep / masses'
        
        # Second thermostat half-step
        thermostat!(velocities, masses, 0.5params.tstep, params.friction,
                   params.temperature, engint, rng)
        
        # Write outputs
        if istep % params.nconfig == 0
            write_positions(params.trajfile, natoms, positions, cell, params.wrapatoms)
        end
        if istep % params.nstat == 0
            engkin = compute_engkin(masses, velocities)
            write_statistics(params.statfile, istep, params.tstep, natoms,
                            engkin, engconf, engint[1])
        end
    end
    
    # Write final configuration
    write_final_positions(params.outputfile, natoms, positions, cell, params.wrapatoms)
    println("Simulation completed successfully")
end

# Run the main function
main()
