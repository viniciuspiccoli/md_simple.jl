using Plots
using ProgressMeter

include("run_md.jl")


function animate_md(input, traj, filename; size=(800,450), fps=20)
    ENV["GKSwstype"] = "nul"
    
    # Create animation frames
    anim = @animate for i in 1:input.nsteps
        plt = plot(layout=(1,2), size=size, framestyle=:box, dpi=150)
        
        # Particle positions
        scatter!(plt[1], traj.positions[i][:,1], traj.positions[i][:,2],
                xlims=(-30, 30), ylims=(-30, 45),
                title="MD Simulation", xlabel="X (Å)", ylabel="Y (Å)",
                marker=:circle, markersize=3, color=:blue,
                legend=false, xticks=:none, yticks=:none, framestyle=:box)
        
        # Thermodynamics plot
        plot!(plt[2], traj.time[1:i], traj.energy[1:i],
              linewidth=2, color=:blue, label="Energy",
              title="System Thermodynamics",
              xlabel="Time (ps)", ylabel="Energy (kcal/mol)")
        
        # Temperature axis
        twin_ax = twinx(plt[2])
        plot!(twin_ax, traj.time[1:i], traj.temp[1:i],
              linewidth=2, color=:red, label="Temperature",
              ylabel="Temperature (K)", legend=:topright)
        
        # Annotations
        annotate!(5, 30, text(@sprintf("Step: %04d\nTime: %.2f ps", 
                 i, traj.time[i]), :left, 10, :black, :clear), subplot=1)
        
        annotate!(5, 29.80, text(@sprintf("Temp: %.2f K\nEnergy: %.2f kcal/mol", 
                 traj.temp[i], traj.energy[i]), :left, 10, :black, :clear), subplot=2)
        
        plt
    end
    
    # Save animation
    gif(anim, filename, fps=fps)
end

# ========================
# MAIN SIMULATION LOOP
# ========================

function run_md_with_animation(dim::SystemDimension)
    # Initialize system
    x, y, z, box, npart = setlat(dim, nx, ny, nz, density)
    box_Å = box .* σ_to_Å
    
    # Initialize velocities and forces
    vx, vy, vz = zeros(npart), zeros(npart), zeros(npart)
    fx, fy, fz = zeros(npart), zeros(npart), zeros(npart)
    initv!(dim, temp, vx, vy, vz, 1.0)
    en = compute_forces!(dim, x, y, z, box, rc, fx, fy, fz)
    
    # Initialize trajectory storage
    traj = (
        positions = [],
        energy = Float64[],
        temp = Float64[],
        time = Float64[]
    )
    
    steps = Int(ceil(tmax/dt))
    prog = Progress(steps, 1)
    
    # Main simulation loop
    t = 0.0
    for step in 0:steps
        # Integration step
        K, en, temp_inst, com_v = integrate!(dim, x, y, z, vx, vy, vz, fx, fy, fz, dt, 1.0, box, rc)
        
        # Store trajectory data
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
            push!(traj.positions, hcat(x.*σ_to_Å, y.*σ_to_Å))
            push!(traj.energy, E_kcal)
            push!(traj.temp, T_K)
            push!(traj.time, t_ps)
        end
        
        t = t + dt
        next!(prog)
    end
    
    # Create animation
    animate_md(
        (nsteps=length(traj.time), box=box_Å),
        traj,
        "md_animation.gif",
        fps=15
    )
end

# Run the simulation

animate_md(
      (nsteps=length(trajectory.time), box=trajectory.box),
      trajectory,
      "md_animation.gif",
      fps=15)
#run_md_with_animation(SYSTEM_TYPE)
