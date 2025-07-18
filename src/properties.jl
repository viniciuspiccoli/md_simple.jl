#=
"""
    compute_rdf(positions, box; dr=0.1, rmax=nothing)

Computes the Radial Distribution Function (g(r)) from a series of particle positions.
Assumes positions and box dimensions are in Ångströms.
"""
function compute_rdf(positions::Vector, box::SVector{3, Float64}; dr=0.05)
    # System properties
    n_frames = length(positions)
    n_particles = size(positions[1], 1)
    is_3d = box[3] > 1e-6 # A small tolerance for 3D check
    
    # Set rmax to half the smallest box dimension to respect PBC
    rmax = 0.5 * minimum(is_3d ? box : box[1:2])
    
    # Histogram setup
    n_bins = Int(floor(rmax / dr))
    hist = zeros(Int, n_bins)
    
    # Density (ρ) in particles/Å³ or particles/Å²
    volume = is_3d ? (box[1] * box[2] * box[3]) : (box[1] * box[2])
    ρ = n_particles / volume

    # Loop over all frames and unique particle pairs
    for frame_positions in positions
        for i in 1:n_particles-1
            for j in i+1:n_particles
                # Minimum image convention for distance
                dr_vec = frame_positions[i,:] - frame_positions[j,:]
                dr_vec -= box .* round.(dr_vec ./ box)
                
                r = norm(dr_vec)
                
                if r < rmax
                    bin = Int(floor(r / dr)) + 1
                    if bin <= n_bins
                        hist[bin] += 1
                    end
                end
            end
        end
    end
    
    # Normalize histogram to get g(r)
    g = zeros(Float64, n_bins)
    r_vals = [(i - 0.5) * dr for i in 1:n_bins] # Bin centers
    
    # The normalization factor for the count in each bin
    # Ideal gas count = (Number of pairs) * (Prob. of being in shell) * (Number of frames)
    #                = (N*(N-1)/2) * (shell_volume / total_volume) * n_frames
    # We can rearrange this to use density ρ = N/V
    norm_const_base = (n_particles - 1) / 2 * ρ * n_frames
    
    for i in 1:n_bins
        r = r_vals[i]
        r_low, r_high = r - dr/2, r + dr/2
        
        # Exact shell volume
        shell_vol = if is_3d
            (4/3) * π * (r_high^3 - r_low^3)
        else
            π * (r_high^2 - r_low^2)
        end
        
        ideal_count = norm_const_base * shell_vol
        if ideal_count > 0
            g[i] = hist[i] / ideal_count
        end
    end
    
    return r_vals, g
end

=#
"""
    compute_rdf(filename::String; dr=0.05)

Computes the Radial Distribution Function (g(r)) from an XYZ trajectory file.
"""
function compute_rdf(filename::String; dr=0.05)
    # 1. Read trajectory data from file (no unwrapping needed for g(r))
    traj = read_xyz_trajectory(filename, unwrap=false)
    
    # 2. Extract system properties
    box = traj.box
    is_3d = box[3] > 1e-6
    rmax = 0.5 * minimum(is_3d ? box : box[1:2])
    n_bins = Int(floor(rmax / dr))
    hist = zeros(Int, n_bins)
    
    volume = is_3d ? prod(box) : box[1] * box[2]
    ρ = traj.n_particles / volume

    # 3. Loop over frames and pairs to build histogram
    for frame_positions in traj.positions
        for i in 1:traj.n_particles-1
            for j in i+1:traj.n_particles
                # Minimum image convention for distance
                dr_vec = frame_positions[i,:] - frame_positions[j,:]

                if is_3d
                    dr_vec .-= box .* round.(dr_vec ./ box) # vectorized application of MIC
                else
                    # Apply MIC only to X and Y for 2D systems
                    dr_vec[1] -= box[1] * round(dr_vec[1] / box[1])
                    dr_vec[2] -= box[2] * round(dr_vec[2] / box[2])
                    dr_vec[3] = 0.0 # Ensure z-component is explicitly zero
                end

                #dr_vec .-= box .* round.(dr_vec ./ box) # Correct MIC application
                
                r = norm(dr_vec)
                
                if r < rmax
                    bin = min(n_bins, Int(floor(r / dr)) + 1)
                    hist[bin] += 1
                end
            end
        end
    end
    
    # 4. Normalize histogram to get g(r)
    g = zeros(Float64, n_bins)
    r_vals = [(i - 0.5) * dr for i in 1:n_bins]
    norm_const_base = (traj.n_particles - 1) / 2 * ρ * traj.n_frames
    
    for i in 1:n_bins
        r_low, r_high = (i - 1) * dr, i * dr
        shell_vol = is_3d ? (4/3)π * (r_high^3 - r_low^3) : π * (r_high^2 - r_low^2)
        ideal_count = norm_const_base * shell_vol
        g[i] = ideal_count > 0 ? hist[i] / ideal_count : 0.0
    end
    
    return r_vals, g
end

export compute_rdf

# A simple trapezoidal rule for integration
function _trapz(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "Vectors must have the same length"
    integral = 0.0
    for i in 1:length(x)-1
        integral += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return integral
end

"""
    running_coordination(r, g_r, ρ)

Calculates the running coordination number n(r) = ρ ∫ g(r') dV.
"""
function running_coordination(r::Vector, g_r::Vector, ρ::Float64)
    # Integrand is g(r) * shell_volume_element
    integrand = g_r .* (4π .* r.^2) # For 3D
    n_r = similar(r)
    for i in 1:length(r)
        # Integrate up to index i
        n_r[i] = ρ * _trapz(r[1:i], integrand[1:i])
    end
    return n_r
end

export running_coordination

"""
    structure_factor(r, g_r, ρ)

Calculates the static structure factor S(k) from g(r).
"""
function structure_factor(r::Vector, g_r::Vector, ρ::Float64; k_max=20.0, dk=0.1)
    k_vals = dk:dk:k_max
    S_k = similar(k_vals)
    
    # S(k) = 1 + 4πρ ∫ [g(r) - 1] r² sin(kr)/(kr) dr
    #      = 1 + (4πρ/k) ∫ [g(r) - 1] r sin(kr) dr
    for (i, k) in enumerate(k_vals)
        if k > 0
            integrand = (g_r .- 1) .* r .* sin.(k .* r)
            integral = _trapz(r, integrand)
            S_k[i] = 1 + (4π * ρ / k) * integral
        end
    end
    
    return k_vals, S_k
end

export structure_factor

"""
    compute_vacf(velocities)

Computes the Velocity Autocorrelation Function (VACF).
Assumes `velocities` is a Vector where each element is an N_particles x 3 matrix of velocities.
"""
function compute_vacf(velocities::Vector)
    n_frames = length(velocities)
    n_particles = size(velocities[1], 1)
    vacf = zeros(n_frames)
    counts = zeros(Int, n_frames)

    for t0 in 1:n_frames
        for dt in 0:n_frames-t0
            t = t0 + dt
            # Sum dot products over all particles for this time lag dt
            vacf[dt+1] += sum(dot(velocities[t0][i,:], velocities[t][i,:]) for i in 1:n_particles)
            counts[dt+1] += n_particles
        end
    end
    
    # Normalize by the number of samples for each time lag
    vacf ./= counts
    # Normalize by the t=0 value
    vacf ./= vacf[1]

    return vacf
end

#=
"""
    compute_msd(unwrapped_positions)

Computes the Mean Squared Displacement (MSD).
Requires UNWRAPPED particle coordinates to correctly track diffusion.
"""
function compute_msd(unwrapped_positions::Vector)
    n_frames = length(unwrapped_positions)
    n_particles = size(unwrapped_positions[1], 1)
    msd = zeros(n_frames)
    counts = zeros(Int, n_frames)

    for t0 in 1:n_frames
        for dt in 0:n_frames-t0
            t = t0 + dt
            # Sum squared displacements over all particles
            msd[dt+1] += sum(norm(unwrapped_positions[t][i,:] - unwrapped_positions[t0][i,:])^2 for i in 1:n_particles)
            counts[dt+1] += n_particles
        end
    end
    
    msd ./= counts
    return msd
end
=#

"""
    compute_msd(filename::String, dt_sample::Float64)

Computes the Mean Squared Displacement (MSD) from an XYZ trajectory file.
It automatically unwraps coordinates to correctly track diffusion.
"""
function compute_msd(filename::String, dt_sample::Float64)
    # 1. Read trajectory with coordinate unwrapping ENABLED
    traj = read_xyz_trajectory(filename, unwrap=true)
    
    # 2. Calculate MSD from the unwrapped positions
    msd = zeros(traj.n_frames)
    counts = zeros(Int, traj.n_frames)

    for t0 in 1:traj.n_frames
        for dt_step in 0:traj.n_frames-t0
            t = t0 + dt_step
            # Sum squared displacements over all particles
            Δr²_sum = 0.0
            for p_idx in 1:traj.n_particles
                Δr²_sum += sum((traj.positions[t][p_idx,:] - traj.positions[t0][p_idx,:]).^2)
            end
            
            msd[dt_step+1] += Δr²_sum
            counts[dt_step+1] += traj.n_particles
        end
    end
    
    msd ./= counts
    
    # 3. Create the time axis in real units (e.g., picoseconds)
    time_axis = [dt_step * dt_sample for dt_step in 0:traj.n_frames-1]
    
    return time_axis, msd
end


export compute_msd

"""
    compute_Cv(total_energy, temperature, n_particles)

Computes the molar heat capacity (Cv) from energy fluctuations in an NVT simulation.
- total_energy: Vector of total energy per mole (e.g., in kcal/mol).
- temperature: Vector of corresponding temperatures (K).
- n_particles: Number of particles in the simulation.
"""
function compute_Cv(total_energy_kcal_mol::Vector, temperature_K::Vector, n_particles::Int)
    @assert length(total_energy_kcal_mol) == length(temperature_K) "Inputs must be same length"
    
    # Physical Constants
    R = 8.31446           # Gas constant in J/(mol·K)
    KCAL_TO_J = 4184      # Conversion factor

    # Convert energy to J/mol
    total_energy_J_mol = total_energy_kcal_mol .* KCAL_TO_J
    
    # Calculate average temperature and variance of molar energy
    T_avg = mean(temperature_K)
    var_E_molar = var(total_energy_J_mol)
    
    # Calculate molar heat capacity using the fluctuation formula
    # Cv_molar = Var(E_molar) / (R * T^2)
    # This formula gives Cv in units of R (dimensionless). To get J/(mol K), we must multiply by R.
    # The standard formula is Cv = Var(E_system) / (kB * T^2).
    # Var(E_system) = Var(E_molar / Na) = Var(E_molar) / Na^2
    # Cv = (Var(E_molar) / Na^2) / (kB * T^2) = Var(E_molar) / (Na * R * T^2)
    # This gives Cv per particle. For molar Cv, multiply by Na.
    # Cv_molar = Var(E_molar) / (R * T^2)
    
    Cv_molar = var_E_molar / (R * T_avg^2)
    
    println("Average Temperature: $(round(T_avg, digits=2)) K")
    println("Energy Variance (J/mol)²: $(round(var_E_molar, digits=2))")
    println("Molar Heat Capacity (Cv): $(round(Cv_molar, digits=2)) J/(mol·K)")
    
    return Cv_molar
end

export compute_Cv