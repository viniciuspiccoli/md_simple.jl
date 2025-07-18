using GLMakie  # Using GLMakie for 3D interactive plots

# FCC lattice generation function (same as before)
function setlat(nx, ny, nz, vol)
    npart = 4 * nx * ny * nz
    a0 = (vol / (nx * ny * nz))^(1/3)
    
    x0 = Vector{Float64}(undef, npart)
    y0 = similar(x0)
    z0 = similar(x0)
    
    i = 0
    xcm0 = 0.0
    ycm0 = 0.0
    zcm0 = 0.0
    
    for iz in 1:2*nz
        for iy in 1:2*ny
            for ix in 1:2*nx
                if (ix + iy + iz) % 2 == 0
                    i += 1
                    x0[i] = 0.5 * a0 * (ix - 1 + 0.5 * ((iy + iz) % 2))
                    y0[i] = 0.5 * a0 * (iy - 1 + 0.5 * ((ix + iz) % 2))
                    z0[i] = 0.5 * a0 * (iz - 1 + 0.5 * ((ix + iy) % 2))
                    
                    xcm0 += x0[i]
                    ycm0 += y0[i]
                    zcm0 += z0[i]
                end
            end
        end
    end
    
    xcm0 /= npart
    ycm0 /= npart
    zcm0 /= npart
    
    return x0, y0, z0, xcm0, ycm0, zcm0, a0
end

# Generate the lattice (2x2x2 unit cells, volume = 100)
nx, ny, nz = 2, 2, 2
vol = 100.0
x, y, z, xcm, ycm, zcm, a0 = setlat(nx, ny, nz, vol)

# Create the plot
fig = Figure(resolution = (800, 600))
ax = Axis3(fig[1, 1], title = "FCC Lattice ($nx×$ny×$nz unit cells)",
           xlabel = "X", ylabel = "Y", zlabel = "Z")

# Plot the atoms as spheres
meshscatter!(ax, x, y, z, markersize = 0.15*a0, color = :lightblue)

# (Optional) Draw unit cell edges for visualization
unit_cell_lines = [
    # Bottom face
    [0, 0, 0] => [a0, 0, 0],
    [a0, 0, 0] => [a0, a0, 0],
    [a0, a0, 0] => [0, a0, 0],
    [0, a0, 0] => [0, 0, 0],
    # Top face
    [0, 0, a0] => [a0, 0, a0],
    [a0, 0, a0] => [a0, a0, a0],
    [a0, a0, a0] => [0, a0, a0],
    [0, a0, a0] => [0, 0, a0],
    # Vertical edges
    [0, 0, 0] => [0, 0, a0],
    [a0, 0, 0] => [a0, 0, a0],
    [0, a0, 0] => [0, a0, a0],
    [a0, a0, 0] => [a0, a0, a0]
]

for (start_pt, end_pt) in unit_cell_lines
    lines!(ax, [start_pt[1], end_pt[1]], [start_pt[2], end_pt[2]], [start_pt[3], end_pt[3]],
           color = :black, linewidth = 1.5)
end

# Center the view on the lattice
limits!(ax, -0.1*a0, (2*nx)*0.5*a0 + 0.1*a0, 
        -0.1*a0, (2*ny)*0.5*a0 + 0.1*a0, 
        -0.1*a0, (2*nz)*0.5*a0 + 0.1*a0)

show(fig)
