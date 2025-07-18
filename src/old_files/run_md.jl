import Pkg
Pkg.activate("../")

using LinearAlgebra, Random, Printf, Plots

# ========================
# PHYSICAL SYSTEM DEFINITION
# ========================
abstract type SystemDimension end
struct TwoD <: SystemDimension end
struct ThreeD <: SystemDimension end

# User configuration - change these parameters
#const SYSTEM_TYPE = ThreeD()    # or TwoD() for 2D simulation
const SYSTEM_TYPE = TwoD()

const nx, ny, nz = 20, 20, 2     # Unit cells (nz ignored in 2D)
const density = 0.6              # Reduced density (ρ*σ^dim)
const temp = 2.28                # Reduced temperature
const temp_desired = temp
const rc = 5.0                   # Cutoff radius (σ units)
const dt = 0.001                 # Time step (reduced time units)
const tmax = 100.0               # Total simulation time
const sample_interval = 25       # Steps between samples



# ========================
# UNIT CONVERSION FACTORS (For Argon)
# ========================

#Physical constants for Argon
const σ = 3.4e-10       # LJ sigma (m)
const ϵ = 1.65e-21      # LJ epsilon (J)
const m_ar = 6.63e-26   # Mass of Argon (kg)

const σ_m = 3.4e-10        # Sigma in meters
const ϵ_J = 1.65e-21       # Epsilon in Joules
const m_kg = 6.63e-26      # Mass in kg
const N_A = 6.022e23       # Avogadro's number

# Derived conversion factors
const σ_to_Å = σ_m * 1e10                           # 3.4 Å
const ϵ_to_kcalmol = ϵ_J / (4.184 * 1e3) * N_A      # 0.238 kcal/mol
const time_to_ps = sqrt(m_kg * σ_m^2 / ϵ_J) * 1e12  # 2.15 ps
const temp_to_K = ϵ_J / (1.38e-23)                  # 120 K

# Conversion functions
real_length(x) = x * σ_to_Å          # Reduced length → Ångströms
real_energy(e) = e * ϵ_to_kcalmol    # Reduced energy → kcal/mol
real_time(t) = t * time_to_ps        # Reduced time → picoseconds
real_temp(T) = T * temp_to_K         # Reduced temp → Kelvin
real_pressure(P) = P * ϵ_J / (σ_m^3) / 1e5  # Reduced pressure → bar

# ========================
# CORE FUNCTIONS
# ========================

# Minimum image convention - works for both 2D and 3D
function minimum_image(dr, L)
    if dr < -L/2
        dr += L
    elseif dr >= L/2
        dr -= L
    end
    return dr
end

# Initialize positions (FCC lattice for 3D, triangular for 2D)
function setlat(dim::SystemDimension, nx, ny, nz, density)
    if dim isa ThreeD
        npart = 4 * nx * ny * nz  # 4 atoms per FCC unit cell
        vol = npart / density
        a0 = cbrt(vol / (nx * ny * nz))  # Lattice constant
        
        x = zeros(npart)
        y = zeros(npart)
        z = zeros(npart)
        
        i = 0
        for iz in 1:2*nz
            for iy in 1:2*ny
                for ix in 1:2*nx
                    if (ix + iy + iz) % 2 == 0
                        i += 1
                        x[i] = 0.5 * a0 * (ix - 1 + 0.5 * ((iy + iz) % 2))
                        y[i] = 0.5 * a0 * (iy - 1 + 0.5 * ((ix + iz) % 2))
                        z[i] = 0.5 * a0 * (iz - 1 + 0.5 * ((ix + iy) % 2))
                    end
                end
            end
        end
        box = (a0 * nx, a0 * ny, a0 * nz)
        return x, y, z, box, npart
    else # 2D case
        npart = nx * ny  # 1 atom per unit cell in triangular lattice
        area = npart / density
        a0 = sqrt(2area / (nx * ny * sqrt(3)))  # Lattice constant
        
        x = zeros(npart)
        y = zeros(npart)
        z = zeros(npart)
        
        i = 0
        for iy in 1:ny
            for ix in 1:nx
                i += 1
                x[i] = a0 * (ix - 1 + 0.5 * ((iy - 1) % 2))
                y[i] = a0 * (iy - 1) * sqrt(3)/2
            end
        end
        box = (a0 * nx, a0 * ny * sqrt(3)/2, 0.0)
        return x, y, z, box, npart
    end
end

# Initialize velocities - works for both 2D and 3D
function initv!(dim::SystemDimension, temp, vx, vy, vz, m)
    npart = length(vx)
    sumv = zeros(3)
    sumv2 = 0.0
    
    # Maxwell-Boltzmann distribution
    for i in 1:npart
        vx[i] = randn()
        vy[i] = randn()
        vz[i] = dim isa ThreeD ? randn() : 0.0
        sumv .+= (vx[i], vy[i], vz[i])
    end
    
    # Remove COM motion
    sumv ./= npart
    for i in 1:npart
        vx[i] -= sumv[1]
        vy[i] -= sumv[2]
        vz[i] -= sumv[3]
        sumv2 += vx[i]^2 + vy[i]^2 + (dim isa ThreeD ? vz[i]^2 : 0.0)
    end
    
    # Scale to desired temperature
    nf = (dim isa ThreeD ? 3 : 2) * npart - (dim isa ThreeD ? 3 : 2)
    fs = sqrt(temp / (sumv2/nf))
    vx .*= fs
    vy .*= fs
    if dim isa ThreeD
        vz .*= fs
    end
    
    return sumv2/nf  # Initial temperature
end


# Updated force calculation function with proper type signatures
function compute_forces!(dim::SystemDimension, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, 
                         box::Tuple{Float64,Float64,Float64}, rc::Float64,
                         fx::Vector{Float64}, fy::Vector{Float64}, fz::Vector{Float64})
    npart = length(x)
    fill!(fx, 0.0)
    fill!(fy, 0.0)
    fill!(fz, 0.0)
    en = 0.0
    rc2 = rc^2
    Lx, Ly, Lz = box
    
    for i in 1:npart-1
        for j in i+1:npart
            # Apply minimum image convention
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = dim isa ThreeD ? z[i] - z[j] : 0.0
            dx = minimum_image(dx, Lx)
            dy = minimum_image(dy, Ly)

            if dim isa ThreeD
                dz = minimum_image(dz, Lz)
            end
            
            r2 = dx^2 + dy^2 + dz^2
            
            if r2 < rc2 && r2 > 0.01  # Avoid division by zero
                r2i = 1.0 / r2
                r6i = r2i^3
                ff = 48.0 * r2i * r6i * (r6i - 0.5)
                
                fx[i] += ff * dx
                fy[i] += ff * dy
                if dim isa ThreeD
                    fz[i] += ff * dz
                    fz[j] -= ff * dz
                end
                
                fx[j] -= ff * dx
                fy[j] -= ff * dy
                
                en += 4.0 * r6i * (r6i - 1.0)
            end
        end
    end
    return en
end



function integrate!(dim::SystemDimension, x, y, z, vx, vy, vz, fx, fy, fz, dt, m, box, rc)
    npart = length(x)
    dt2 = dt/2
    Lx, Ly, Lz = box
    
    # First half-step velocity update and position update
    for i in 1:npart
        vx[i] += fx[i] * dt2 / m
        vy[i] += fy[i] * dt2 / m
        if dim isa ThreeD
            vz[i] += fz[i] * dt2 / m
        end
        
        x[i] += vx[i] * dt
        y[i] += vy[i] * dt
        if dim isa ThreeD
            z[i] += vz[i] * dt
        end
        
        # Apply periodic boundary conditions
        x[i] -= Lx * floor(x[i]/Lx + 0.5)
        y[i] -= Ly * floor(y[i]/Ly + 0.5)
        if dim isa ThreeD
            z[i] -= Lz * floor(z[i]/Lz + 0.5)
        end
    end
    
    # Compute new forces
    en = compute_forces!(dim, x, y, z, box, rc, fx, fy, fz)
    
    # Second half-step velocity update and KE calculation
    K = 0.0
    sumv = zeros(3)
    for i in 1:npart
        vx[i] += fx[i] * dt2 / m
        vy[i] += fy[i] * dt2 / m
        if dim isa ThreeD
            vz[i] += fz[i] * dt2 / m
        end
        
        K += 0.5 * m * (vx[i]^2 + vy[i]^2 + (dim isa ThreeD ? vz[i]^2 : 0.0))
        sumv .+= (vx[i], vy[i], dim isa ThreeD ? vz[i] : 0.0)
    end
    
    # Temperature and COM velocity
    nf = (dim isa ThreeD ? 3 : 2) * npart - (dim isa ThreeD ? 3 : 2)
    temp_inst = 2K / nf
    com_v = sumv ./ npart
    lambda = sqrt(1 + (dt/0.1)*(temp_desired/temp_inst - 1))
    vx .*= lambda
    vy .*= lambda
    vz .*= lambda
    K *= lambda^2
    
    return K, en, temp_inst, com_v
end


# ========================
# VISUALIZATION FUNCTIONS
# ========================


function write_xyz_frame(dim::SystemDimension, io, x, y, z, box, step, time)
    npart = length(x)
    Lx, Ly, Lz = box
    
    # Convert box dimensions to Angstroms
    Lx_Å = Lx * σ * 1e10
    Ly_Å = Ly * σ * 1e10
    Lz_Å = dim isa ThreeD ? Lz * σ * 1e10 : 0.0

    # Format box information for VMD
    if dim isa ThreeD
        # For 3D: "a b c alpha beta gamma" (cubic box)
        box_info = @sprintf("\"%.3f %.3f %.3f 90.0 90.0 90.0\"", Lx_Å, Ly_Å, Lz_Å)
    else
        # For 2D: "a b 0.0 90.0" (rectangular box in XY-plane)
        box_info = @sprintf("\"%.3f %.3f 0.0 90.0\"", Lx_Å, Ly_Å)
    end

    println(io, npart)
    println(io, "Lattice=$box_info Step=$step Time=$(real_time(time))")

    for i in 1:npart
        x_real = x[i] * σ * 1e10  # Convert to Angstroms
        y_real = y[i] * σ * 1e10
        z_real = dim isa ThreeD ? z[i] * σ * 1e10 : 0.0
        println(io, "Ar ", x_real, " ", y_real, " ", z_real)
    end
end

#=
function write_xyz_frame(dim::SystemDimension, io, x, y, z, box, step, time)
    npart = length(x)
    Lx, Ly, Lz = box
    
    # Convert box dimensions to Angstroms
    Lx_Å = Lx * σ * 1e10
    Ly_Å = Ly * σ * 1e10
    Lz_Å = dim isa ThreeD ? Lz * σ * 1e10 : 0.0

    # Format box information for VMD (origin at center)
    box_info = if dim isa ThreeD
        @sprintf("\"%.3f 0 0 0 %.3f 0 0 0 %.3f\"", Lx_Å, Ly_Å, Lz_Å)
    else
        @sprintf("\"%.3f 0 0 %.3f 0 0 0 0 1\"", Lx_Å, Ly_Å)
    end

    println(io, npart)
    println(io, "Lattice=$box_info Step=$step Time=$(real_time(time))")

    # Center coordinates in box and convert to Angstroms
    for i in 1:npart
        # Convert to Angstroms and center in box
        x_real = (x[i] * σ * 1e10 + Lx_Å/2) % Lx_Å
        y_real = (y[i] * σ   #return K, en, temp_inst, com_v
    end* 1e10 + Ly_Å/2) % Ly_Å
        z_real = dim isa ThreeD ? (z[i] * σ * 1e10 + Lz_Å/2) % Lz_Å : 0.0

        println(io, "Ar ", x_real, " ", y_real, " ", z_real)
    end
end
=#


# ========================
# MAIN SIMULATION FUNCTION
# ========================

mutable struct MDobject
    positions
    energy::Vector{Float64}
    temp::Vector{Float64}
    time::Vector{Float64}
    box
end

function run_md(dim::SystemDimension)
    # Initialize positions
    x, y, z, box, npart = setlat(dim, nx, ny, nz, density)
    
    # Initialize velocities
    vx = zeros(npart)
    vy = zeros(npart)
    vz = zeros(npart)
    init_temp = initv!(dim, temp, vx, vy, vz, 1.0)  # Using m=1 in reduced units
    
    # Initialize forces
    fx = zeros(npart)
    fy = zeros(npart)
    fz = zeros(npart)
    en = compute_forces!(dim, x, y, z, box, rc, fx, fy, fz)
    
    # Open trajectory file
    filename = dim isa ThreeD ? "argon_3d.xyz" : "argon_2d.xyz"
    traj_file = open(filename, "w")
    
    # Main MD loop
    t = 0.0
    step = 0

    traj = MDobject([], Float64[], Float64[], Float64[],Tuple{Float64,Float64,Float64}[])

    while t < tmax
        # Integration step
        K, en, temp_inst, com_v = integrate!(dim, x, y, z, vx, vy, vz, fx, fy, fz, dt, 1.0, box, rc)
        
        # Write trajectory
        if step % sample_interval == 0
            write_xyz_frame(dim, traj_file, x, y, z, box, step, t)
        end
        
        # Print diagnostics
        if step % 100 == 0
            etot = (en + K)/npart

            t_ps = real_time(t)
            E_kcal = real_energy(etot)
            if dim isa ThreeD
                # 3D: 3N-3 degrees of freedom
                T_K = real_temp(2K / (3npart - 3))
            else
                # 2D: 2N-2 degrees of freedom
                T_K = real_temp(K / (npart - 1))
            end
            @printf("Step %6d: Time %6.2f ps\n", step, t_ps)
            @printf("  T = %6.2f K  E_tot = %8.4f kcal/mol  COM_v = %6.4f Å/ps\n",
                    T_K, E_kcal, norm(com_v)*σ_to_Å/time_to_ps)

            push!(traj.positions, hcat(x.*σ_to_Å, y.*σ_to_Å))
            push!(traj.energy, E_kcal)
            push!(traj.temp, T_K)
            push!(traj.time, t_ps)
            push!(traj.box, box .* σ_to_Å)

        end
        t += dt
        step += 1
    end
    
    close(traj_file)
    println("Simulation completed. Trajectory saved to $filename")
    return traj
end

# ========================
# RUN THE SIMULATION
# ========================

# Start the simulation (change SYSTEM_TYPE at top to switch between 2D/3D)
trajectory = run_md(SYSTEM_TYPE)




# vmd commmands
#=
pbc box -on
pbc set {20.4 20.4 20.4} -all

pbc box -center origin -width 0.5  # Draw box at simulation center
pbc wrap -all                     # Ensure all atoms are in primary cell
=#
