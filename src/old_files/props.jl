# scripts to calculate cool stuff from my MD code

include("run_md.jl")


function compute_rdf(traj; dr=0.1, rmax=nothing)
	
	# Extract box dimensions from trajectory (assumes constant box)
	Lx, Ly, Lz = traj.box[1]
	is_3d = Lz > 0 # Check if 3D system
	# Set maximum distance to half the smallest box dimension
	if Lz == 0
		rmax = 0.50 * minimum([Lx,Ly])
	else Lz > 0
		rmax = 0.50 * minimum([Lx, Ly, Lz])/2
	end
	nbins = Int(ceil(rmax/dr))
	dr = rmax / nbins # Re-adjust bin width
	hist = zeros(nbins)
	N = size(traj.positions[1], 1) # Number of particles
	nframes = length(traj.positions)
	V = is_3d ? (Lx * Ly * Lz) : (Lx * Ly)
	ρ = N / V # Number density
	for frame in traj.positions
		# Extract coordinates with explicit 3D handling
		positions = is_3d ? frame : [frame[:, 1:2] zeros(size(frame, 1))]
		for i in 1:N-1
			for j in i+1:N
				# Compute distance with PBC
				dx = positions[i, 1] - positions[j, 1]
				dx -= Lx * round(dx/Lx)
				dy = positions[i, 2] - positions[j, 2]
				dy -= Ly * round(dy/Ly)
				dz = is_3d ? (positions[i, 3] - positions[j, 3] - Lz * round((positions[i, 3] - positions[j, 3])/Lz)) : 0.0
				r = sqrt(dx^2 + dy^2 + dz^2)
				if r < rmax
					bin = Int(floor(r/dr)) + 1
				if bin ≤ nbins # Prevent index overflow
					hist[bin] += 1 # Count each pair once
				end
			end
		end
	end
	
	end
	
	# Normalize to get g(r)
	r = [dr*(i-0.5) for i in 1:nbins] # Bin centers
	g = zeros(nbins)
	for i in 1:nbins
		r_val = r[i]
		#shell_vol = is_3d ? 4π*r_val^2*dr : 2π*r_val*dr
		shell_vol = is_3d ? (4/3)*π*((r_val + dr)^3 - r_val^3) : π*((r_val + dr)^2 - r_val^2)
		#ideal_count = (N * (N - 1)/2) * nframes * shell_vol * ρ * 1.e-3
		ideal_count = (N * ρ * shell_vol) * (nframes * (N - 1)/2) * 1.e-3
		g[i] = hist[i] / ideal_count
	end
	return r, g

end



# Calculate RDF from trajectory data

r, g_r = compute_rdf(trajectory, dr=0.01)


# Plot results

using Plots

plot(r, g_r,

xlabel="r (Å)", ylabel="g(r)",

title="Radial Distribution Function",

label="RDF", linewidth=2, framestyle=:box, xlims=(0,15))


savefig("gr.png")





#=

function coordination_number(r, g_r, r_max)

idx = findlast(r .< r_max)

sum(g_r[1:idx] .* 4π .* r[1:idx].^2 .* (r[2]-r[1]))

end

# First coordination shell for Argon (~3.8Å)

cn = coordination_number(r, g_r, 3.8)


function running_coordination(r, g_r)

dr = r[2] - r[1]

[sum(g_r[1:i] .* 4π .* r[1:i].^2 .* dr) for i in 1:length(r)]

end


n_r = running_coordination(r, g_r)


function structure_factor(r, g_r, ρ)

dr = r[2] - r[1]

k = 0.1:0.1:20 # Å⁻¹

S_k = similar(k)

for (i, ki) in enumerate(k)

integrand = @. (g_r - 1) * r * sin(ki*r)/(ki*r)

S_k[i] = 1 + 4π*ρ*trapz(r, integrand)

end

k, S_k

end


k, S_k = structure_factor(r, g_r, ρ)





function compute_vacf(traj)

nframes = length(traj.positions)

vacf = zeros(nframes)

velocities = [traj.positions[t+1] .- traj.positions[t] for t in 1:(nframes-1)]

@showprogress for dt in 1:nframes

for t0 in 1:(nframes - dt)

dot_sum = 0.0

for i in 1:size(velocities[t0], 1)

v0 = velocities[t0][i,:]

vt = velocities[t0+dt][i,:]

dot_sum += dot(v0, vt)

end

vacf[dt] += dot_sum/size(velocities[t0], 1)

end

end

vacf ./= (nframes - 1)

return vacf

end


function compute_msd(traj)

nframes = length(traj.positions)

msd = zeros(nframes)

box_dims = first(traj.box) # Assumes constant box size

@showprogress for dt in 1:nframes

for t0 in 1:(nframes - dt)

Δr² = 0.0

for i in 1:size(traj.positions[t0], 1)

# Use modulo for PBC wrapping

dx = mod(traj.positions[t0+dt][i,1] - traj.positions[t0][i,1], box_dims[1])

dy = mod(traj.positions[t0+dt][i,2] - traj.positions[t0][i,2], box_dims[2])

dz = mod(traj.positions[t0+dt][i,3] - traj.positions[t0][i,3], box_dims[3])

Δr² += dx^2 + dy^2 + dz^2

end

msd[dt] += Δr²/size(traj.positions[t0], 1)

end

end

msd ./= (nframes - 1)

return msd

end


compute Cv or Cp

=#
