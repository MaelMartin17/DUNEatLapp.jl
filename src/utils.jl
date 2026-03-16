"""
    fC_to_MeV(q::Vector{<:AbstractFloat}, E::Float64)

Convert charge values (in femtocoulombs, fC) to energy (in MeV) using the ArgoNEUT 2013 calibration method (JINST 8 P08005).

# Arguments
- `q`: Vector of charge values in femtocoulombs (fC).
- `E`: Electric field strength in volts per centimeter (kV/cm).

# Returns
- Vector of energy values in mega-electronvolts (MeV).

# Details
The conversion follows the recombination model:
E = (exp(q * β * W / q_e) - α) / β
where:
- β = 0.212 / 1.3849 / E (recombination parameter),
- W = 23.6e-6 MeV (work function for liquid argon),
- q_e = 1.604e-4 fC (electron charge in fC),
- α = 0.93 (recombination constant).

# References
- ArgoNEUT 2013, JINST 8 P08005: [https://doi.org/10.1088/1748-0221/8/08/P08005](https://doi.org/10.1088/1748-0221/8/08/P08005)
"""
function fC_to_MeV(q::Vector{<:AbstractFloat}, E::Float64)
    betarhoe = 0.212 / 1.3849 / E
    qe = 1.604e-4
    W = 23.6e-6
    alpha = 0.93

    return [(exp(x * betarhoe * W / qe) - alpha) / betarhoe for x in q]
end

"""
    plot_detector_planes(z_height = 338.5; alpha = 0.4)

Plot the detector planes (top, bottom, and cathode) at the specified `z_height`.

# Arguments
- `z_height`: The height of the top and bottom planes.
- `alpha`: The transparency of the planes.

# Returns
- `plt`: The plot object.
"""
function plot_detector_planes(z_height = 338.5; alpha = 0.4)
    x_range = -337.5:50:337.5
    y_range = 0:50:300
    x = [i for i in x_range, j in y_range]
    y = [j for i in x_range, j in y_range]

    top_plane = fill(z_height, length(x_range), length(y_range))
    bottom_plane = fill(-z_height, length(x_range), length(y_range))
    cathode = fill(0.0, length(x_range), length(y_range))

    plt = plot(legend=false, colorbar=false)
    surface!(plt, x, y, top_plane; color=:orange, alpha=alpha, linecolor=:orange)
    surface!(plt, x, y, bottom_plane; color=:orange, alpha=alpha, linecolor=:orange)
    surface!(plt, x, y, cathode; color=:black, alpha=alpha, linecolor=:black)

    return plt
end


"""
    fit_line_3d(x, y, z)

Fit a 3D line to the data points `(x, y, z)` using PCA.

# Arguments
- `x`, `y`, `z`: Vectors of coordinates for the points.

# Returns
- `point_on_line`: A point on the fitted line (the mean of the input points).
- `direction`: The direction vector of the fitted line.
- `rms_residual`: The root-mean-square of the perpendicular distances from the points to the line.
"""
function fit_line_3d(x, y, z)
    points = hcat(x, y, z)
    mean_point = mean(points, dims=1)
    centered_points = points .- mean_point
    cov_matrix = cov(centered_points)
    eig = eigen(cov_matrix)
    direction = eig.vectors[:, argmax(eig.values)]
    point_on_line = vec(mean(points, dims=1))

    # Compute fit quality
    residuals = [p - project_point_to_line(p, point_on_line, direction) for p in eachrow(points)]
    d_perp = [norm(r) for r in residuals]
    rms_residual = sqrt(mean(d_perp .^ 2))

    return point_on_line, direction, rms_residual
end

"""
    project_point_to_line(p, p0, v)

Project a 3D point `p` onto the line defined by `p0` and direction vector `v`.

# Arguments
- `p`: The point to project.
- `p0`: A point on the line.
- `v`: The direction vector of the line, should be unitary but it's normalized in the function in case it's not

# Returns
- The projected point on the line.
"""
function project_point_to_line(p, p0, v)
    t = dot(p - p0, v) / dot(v, v)
    return p0 + t * v
end

"""
    compute_residuals(points::AbstractMatrix, p0, v)

Compute the residuals and perpendicular distances from each point in `points` to the line defined by `p0` and direction vector `v`.

# Arguments
- `points`: An `N × 3` matrix of 3D points.
- `p0`: A point on the line.
- `v`: The direction vector of the line.

# Returns
- `residuals`: An `N × 3` matrix of residual vectors.
- `distances`: A vector of perpendicular distances from each point to the line.
"""
function compute_residuals(points::AbstractMatrix, p0, v)
    N = size(points, 1)
    residuals = zeros(N, 3)
    distances = zeros(N)

    for i in 1:N
        p = vec(points[i, :])
        p_proj = project_point_to_line(p, p0, v)
        r = p - p_proj

        residuals[i, :] .= r
        distances[i] = norm(r)
    end

    return residuals, distances
end

"""
    prepare_data_for_fit(df::AbstractDataFrame, z_value)

Prepare a DataFrame for 3D linear fitting by adjusting the `z` values and filtering rows based on a specified range.

# Arguments
- `df`: The input DataFrame containing the data.
- `z_value`: The value to add to the `z` column.

# Returns
- `filtered_df`: The filtered DataFrame where `z` values are within the range `[z_min, z_max]`.
"""
function prepare_data_for_fit(df::AbstractDataFrame, z_value)
    df.z .+= z_value
    z_max = 338.0
    z_min = 288.0
    if minimum(df.z) < -10.0
        z_min = -338.0
        z_max = -288.0
    end
    # Filter rows where z is within the specified range
    filtered_df = df[(df.z .>= z_min) .& (df.z .<= z_max), :]
    return filtered_df
end


"""
    reconstruct_low_e_hit(Q_fC)

Best for hits < 10 MeV where dx is uncertain. 
Uses the 'MIP-limit' recombination rather than local dQ/dx.
"""
function reconstruct_low_e_hit(Q_fC)
    # Constants
    W_ion = 23.6e-6  # MeV/e-
    fC_to_e = 6241.5
    
    # 1. Assume the electron is a MIP (dE/dx ~ 2.12 MeV/cm)
    # At 0.45 kV/cm, the recombination factor R is ~0.635
    # (Calculated using Birks with ICARUS params)
    #R_avg = 0.635 
    R_avg = 0.57 #From ProtoDUNE-HD + Bi-207 estimation
    
    # 2. Convert total collected charge directly to Energy
    Energy_MeV = (Q_fC * fC_to_e * W_ion) / R_avg
    
    return Energy_MeV
end

