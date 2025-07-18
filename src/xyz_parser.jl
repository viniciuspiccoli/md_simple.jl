using LinearAlgebra
using Printf

"""
    read_xyz_trajectory(filename::String; unwrap::Bool=false)

Reads an XYZ trajectory file.

Parses box dimensions from the VMD-style 'Lattice="..."' comment.
If `unwrap=true`, it calculates and returns unwrapped coordinates necessary
for diffusion calculations like MSD.

Returns a NamedTuple: `(positions, box, n_particles, n_frames)`
"""
function read_xyz_trajectory(filename::String; unwrap::Bool=false)
    positions_per_frame = Vector{Matrix{Float64}}()
    box = zeros(3)
    local n_particles

    lines = readlines(filename)
    i = 1
    frame_count = 0
    
    # For unwrapping
    previous_positions = nothing
    offsets = nothing

    while i <= length(lines)
        # --- Read Frame Header ---
        n_particles = parse(Int, lines[i])
        comment_line = lines[i+1]
        
        # On the first frame, parse box and initialize unwrapping data
        if frame_count == 0
            # Regex to find "Lattice=" followed by numbers in quotes
            m = match(r"Lattice=\"\s*([0-9\.\-eE]+)\s+([0-9\.\-eE]+)\s+([0-9\.\-eE]+)", comment_line)
            if m !== nothing
                box = [parse(Float64, s) for s in m.captures]
            else
                error("Could not parse box dimensions from the XYZ file's comment line. Expected 'Lattice=\"Lx Ly Lz ...\"'")
            end
            if unwrap
                offsets = [zeros(3) for _ in 1:n_particles]
            end
        end
        
        # --- Read Particle Coordinates ---
        frame_positions = zeros(Float64, n_particles, 3)
        for j in 1:n_particles
            parts = split(lines[i+2+j-1])
            frame_positions[j, 1] = parse(Float64, parts[2])
            frame_positions[j, 2] = parse(Float64, parts[3])
            frame_positions[j, 3] = parse(Float64, parts[4])
        end

        # --- Handle Coordinate Unwrapping for MSD ---
        if unwrap && previous_positions !== nothing
            for p_idx in 1:n_particles
                # Raw displacement from previous frame
                displacement = frame_positions[p_idx, :] - previous_positions[p_idx, :]
                
                # Check for boundary crossings in each dimension
                for dim in 1:3
                    if box[dim] > 0 # Only for periodic dimensions
                        if displacement[dim] > box[dim] / 2
                            offsets[p_idx][dim] -= 1 # Crossed negative boundary
                        elseif displacement[dim] < -box[dim] / 2
                            offsets[p_idx][dim] += 1 # Crossed positive boundary
                        end
                    end
                end
                # Apply the offset to get the unwrapped position
                frame_positions[p_idx, :] .+= offsets[p_idx] .* box
            end
        end
        
        push!(positions_per_frame, frame_positions)
        if unwrap
            previous_positions = copy(frame_positions) # Store unwrapped positions for next frame's comparison
        end
        
        # Move to the next frame
        i += n_particles + 2
        frame_count += 1
    end

    return (
        positions=positions_per_frame,
        box=box,
        n_particles=n_particles,
        n_frames=frame_count
    )
end
