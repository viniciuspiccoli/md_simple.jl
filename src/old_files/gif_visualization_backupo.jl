using Plots
using ProgressMeter

include("run_md.jl")

gr()  # Set GR backend for faster rendering

#pgfplotsx()

# Initialize plot data storage
const history = (
    time = Float64[],
    temp = Float64[],
    energy = Float64[]
)

# Create initial figure layout
function create_figure(box_size)
    plt = plot(layout=(1,1), dpi=150,
               xguidefontsize=10, yguidefontsize=10,
               legendfontsize=8, titlefontsize=12)
    
    # 2D simulation view
    plot!(plt[1], xlim=(0,box_size[1]), ylim=(0,box_size[2]),
          title="MD Simulation", xlabel="X (Å)", ylabel="Y (Å)", legend=false, framestyle=:box)
    
    # Thermodynamics panel
    plot(title="System Thermodynamics", 
          xlabel="Time (ps)", ylabel="Energy (kcal/mol)",
          legend=:topleft, grid=true, framestyle=:box)
    
    # Create twin axis for temperature
    p2 = plot!(twinx(), ylabel="Temperature (K)", 
          legend=:topright, grid=true, framestyle=:box)

    hline!(p2, [0], linestyle=:dot, color=:black, label="")

    A = plot(plt,p2)
    return A
end

# Update figure for each frame
function update_plot!(plt, x, y, step, time, temp, energy, box_size)
    # Clear previous frame
    plot!(plt[1], subplot=1, clear=true)
    plot!(plt[2], subplot=2, clear=false)
    
    # Convert coordinates to Angstroms
    x_Å = x .* σ * 1e10
    y_Å = y .* σ * 1e10
    
    # Simulation view (3D scatter)
    scatter!(plt[1], x_Å, y_Å, 
            marker=:circle, markersize=2, 
            marker_z=sqrt.(x_Å.^2 + y_Å.^2),
            colorbar=false, color=:thermal)
    
    # Store thermodynamics data
    push!(history.time, real_time(time))
    push!(history.temp, real_temp(temp))
    push!(history.energy, real_energy(energy))
    
    # Plot thermodynamics
    #plot!(plt[2], history.time, history.energy, 
    #      label="Total Energy", linewidth=2, color=:blue)
    #plot!(plt[2], history.time, history.temp, 
    #      label="Temperature", linewidth=2, color=:red, yaxis=:right)

    plot!(plt[2], history.time, history.energy, 
          label=false, linewidth=2, color=:blue)
    plot!(plt[2], history.time, history.temp, 
          label=false, linewidth=2, color=:red, yaxis=:right)
 
    
    # Add current step annotation
    annotate!(plt[1], 20, 30, text("Step: $step \nTime: $(round(real_time(time), digits=2)) ps",:left, 10))
end

## MD simulation

function run_md_with_animation(dim::SystemDimension)
    # Initialize positions
    x, y, z, box, npart = setlat(dim, nx, ny, nz, density)
    box_Å = box .* σ .* 1e10

    # Create initial figure
    plt = create_figure(box_Å)
    anim = Animation()

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
    
    # Simulation parameters
    steps = Int(ceil(tmax/dt))
    prog = Progress(steps, 1)  # Progress meter

    # Main MD loop
    t = 0.0
    step = 0
    @gif for step in 0:steps
   # while t < tmax
        # Integration step
        K, en, temp_inst, com_v = integrate!(dim, x, y, z, vx, vy, vz, fx, fy, fz, dt, 1.0, box, rc)
        
        # Update plot every 10 steps
        if step % 100 == 0
            etot = (en + K)/npart
            t_ps = real_time(t)
            E_kcal = real_energy(etot)
            T_K = real_temp(2K/(2npart - 2))  # For 3D
            update_plot!(plt, x, y, step, t_ps, T_K, E_kcal, box_Å)
            frame(anim)
        end
        t += dt
        step += 1
        next!(prog)  # Update progress bar
    end every 100  # Controls GIF frame rate

    gif(anim, "md_animation.gif", fps=15)

end

run_md_with_animation(SYSTEM_TYPE)