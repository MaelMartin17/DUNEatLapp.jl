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

