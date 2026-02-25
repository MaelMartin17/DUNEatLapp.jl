"""
    check_selection_vd(track, zcorr, nHits::Int=30) -> Bool

Check if a track meets specific selection criteria based on its proximity to anodes in the drift axis.

# Arguments
- `track`: A collection of track points, where each point is expected to have at least 3 elements,
           and the third element represents the z-coordinate.
- `zcorr`: A correction value to be applied to the z-coordinate of each track point.
- `nHits::Int=30`: The minimum number of hits required for the track to be considered valid.
                   Defaults to 30.

# Returns
- `Bool`: `true` if the track meets the selection criteria, otherwise `false`.

# Details
The function checks the following conditions:
1. The track must have at least `nHits` points.
2. The track must reach near one of the anodes (z-coordinate magnitude > 330).
3. If the track reaches the positive anode side (z > 330), it must be within 10 units of the anode (338.0).
4. If the track reaches the negative anode side (z < -330), it must be within 10 units of the anode (-338.0).

# Examples
```jldoctest
track = [[0, 0, 335], [0, 0, 336], ...]  # Example track with z-coordinates
zcorr = 0.0
check_selection_vd(track, zcorr)  # Returns true if the track meets the criteria
"""

function check_selection_vd(track, zcorr, nHits::Int=30)

    N = length(track)
    N ≥ nHits || return false

    anode = 338.0

    zmin =  Inf
    zmax = -Inf
    # Find max and min in drift axis reach near one of the anodes
    @inbounds for row in track
        zc = row[3] + zcorr
        zc < zmin && (zmin = zc)
        zc > zmax && (zmax = zc)
    end
    # Must reach near one of the anodes
    max(abs(zmin), abs(zmax)) > 330 || return false

    # If it reaches +anode side, require it to be close
    if zmax > 330 && abs(zmax - anode) > 10
        return false
    end

    # If it reaches -anode side, require it to be close
    if zmin < -330 && abs(zmin + anode) > 10
        return false
    end

    return true
end

"""
    fit_line_3d_anode(track, zcorr::Float64=0.0, nFit::Int=10, zLow=288.0, zHigh=338.0) -> Union{Nothing, Tuple{NTuple{3, Float64}, Vector{Float64}, Float64}}

Fit a 3D line to a track within a specified anode window and return the mean position, direction vector, and RMS of the fit.

# Arguments
- `track`: A collection of track points, where each point is expected to have at least 3 elements (x, y, z).
- `zcorr::Float64=0.0`: A correction value to be applied to the z-coordinate of each track point.
- `nFit::Int=10`: The minimum number of hits required in the anode window for the fit to be valid.
- `zLow=288.0`: The lower bound of the anode window for the positive side.
- `zHigh=338.0`: The upper bound of the anode window for the positive side.

# Returns
- `Union{Nothing, Tuple{NTuple{3, Float64}, Vector{Float64}, Float64}}`:
  - If the fit is valid, returns a tuple containing:
    - `(μx, μy, μz)`: The mean position of the selected hits.
    - `v`: The direction vector of the fitted line (eigenvector corresponding to the largest eigenvalue of the covariance matrix).
    - `rms`: The root mean square of the residuals perpendicular to the fitted line.
  - If the fit is invalid (e.g., insufficient hits or zero residuals), returns `nothing`.

# Details
1. **Anode Window Selection**:
   - If the track's minimum z-coordinate is less than -300.0, the anode window is mirrored to the negative side.
   - Otherwise, the window is defined by `zLow` and `zHigh`.

2. **Mean Calculation**:
   - Computes the mean position `(μx, μy, μz)` of the hits within the anode window.

3. **Covariance Matrix**:
   - Computes the covariance matrix of the selected hits and extracts the direction vector `v` as the eigenvector corresponding to the largest eigenvalue.

4. **RMS Calculation**:
   - Computes the RMS of the residuals perpendicular to the fitted line.

# Examples
```jldoctest
track = [[1.0, 2.0, 300.0], [1.1, 2.1, 301.0], ...]  # Example track with x, y, z coordinates
zcorr = 0.0
result = fit_line_3d_anode(track, zcorr)
# Returns ((μx, μy, μz), v, rms) if the fit is valid, otherwise nothing
"""
function fit_line_3d_anode(track, zcorr::Float64=0.0,nFit::Int=10, zLow = 288.0, zHigh = 338.0)

    # ---- First pass: find corrected z minimum ----
    zmin =  Inf
    zmax = -Inf
    # Find max and min in drift axis reach near one of the anodes
    @inbounds for row in track
        zc = row[3] + zcorr
        zc < zmin && (zmin = zc)
        zc > zmax && (zmax = zc)
    end

    # ---- Define anode window for fit----
    if zmin < -300.0
        z_low  = -zHigh
        z_high = -zLow
    else
        z_low  = zLow
        z_high = zHigh
    end
    # ---- First pass: compute mean only on selected hits ----
    μx = 0.0; μy = 0.0; μz = 0.0
    N = 0
    #check hits inside selection
    @inbounds for row in track
        zc = row[3] + zcorr
        if z_low ≤ zc ≤ z_high
            μx += row[1]
            μy += row[2]
            μz += zc
            N += 1
        end
    end
    #Ask at least 10 hits for the fit
    N > nFit || return nothing
    

    μx /= N; μy /= N; μz /= N

    # ---- Covariance on selected hits ----
    Sxx=0.0; Sxy=0.0; Sxz=0.0
    Syy=0.0; Syz=0.0; Szz=0.0

    @inbounds for row in track
        zc = row[3] + zcorr
        if z_low ≤ zc ≤ z_high
            dx = row[1] - μx
            dy = row[2] - μy
            dz = zc - μz

            Sxx += dx*dx
            Sxy += dx*dy
            Sxz += dx*dz
            Syy += dy*dy
            Syz += dy*dz
            Szz += dz*dz
        end
    end

    cov = Symmetric([
        Sxx Sxy Sxz;
        Sxy Syy Syz;
        Sxz Syz Szz
    ])

    eigvals, eigvecs = eigen(cov)
    v = eigvecs[:, argmax(eigvals)]

    # ---- RMS only on selected hits ----
    sumsq = 0.0

    @inbounds for row in track
        zc = row[3] + zcorr
        if z_low ≤ zc ≤ z_high
            dx = row[1] - μx
            dy = row[2] - μy
            dz = zc - μz

            proj = dx*v[1] + dy*v[2] + dz*v[3]
            sumsq += dx*dx + dy*dy + dz*dz - proj^2
        end
    end
    #security check
    sumsq > 0 || return nothing

    rms = sqrt(sumsq / N)

    return (μx, μy, μz), v, rms
end

"""
    compute_true_and_residuals(track, p0, v; zcorr=0.0) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Compute the true (projected) positions and residuals of a track relative to a reference line defined by a point and a direction vector.

# Arguments
- `track`: A collection of track points, where each point is expected to have at least 3 elements (x, y, z).
- `p0`: A reference point (x₀, y₀, z₀) on the line.
- `v`: A direction vector of the line.
- `zcorr=0.0`: A correction value to be applied to the z-coordinate of each track point.

# Returns
- A tuple containing:
  - `xt`: Vector of true x-coordinates (projected onto the line).
  - `yt`: Vector of true y-coordinates (projected onto the line).
  - `zt`: Vector of true z-coordinates (projected onto the line).
  - `rx`: Vector of residuals in the x-direction.
  - `ry`: Vector of residuals in the y-direction.
  - `rz`: Vector of residuals in the z-direction.
  - `r`: Vector of total residual displacements (Euclidean norm).

# Details
For each point in `track`, the function:
1. Computes the projection of the point onto the line defined by `p0` and `v`.
2. Calculates the residual vector as the difference between the measured position and the projected position.
3. Computes the Euclidean norm of the residual vector.
# Examples
```jldoctest
track = [[1.0, 2.0, 300.0], [1.1, 2.1, 301.0], ...]  # Example track with x, y, z coordinates
p0 = [0.0, 0.0, 0.0]  # Reference point
v = [1.0, 1.0, 1.0]   # Direction vector
xt, yt, zt, rx, ry, rz, r = compute_true_and_residuals(track, p0, v)
# Returns the true positions and residuals for each point in the track
"""
function compute_true_and_residuals(track, p0, v; zcorr=0.0)

    N = length(track)
    #true values
    xt = Vector{Float64}(undef, N)
    yt = Vector{Float64}(undef, N)
    zt = Vector{Float64}(undef, N)
    #residuals values
    rx = Vector{Float64}(undef, N)
    ry = Vector{Float64}(undef, N)
    rz = Vector{Float64}(undef, N)
    #total displacement
    r  = Vector{Float64}(undef, N)

    invnorm2 = 1.0 / dot(v,v)

    @inbounds for i in 1:N
        row = track[i]

        # Measured position
        x = row[1]
        y = row[2]
        z = row[3] + zcorr

        # Relative to reference point
        dx = x - p0[1]
        dy = y - p0[2]
        dz = z - p0[3]

        # Projection parameter
        t = (dx*v[1] + dy*v[2] + dz*v[3]) * invnorm2

        # True position (projection)
        xti = p0[1] + t*v[1]
        yti = p0[2] + t*v[2]
        zti = p0[3] + t*v[3]

        xt[i] = xti
        yt[i] = yti
        zt[i] = zti

        # Residual vector
        rxi = x - xti
        ryi = y - yti
        rzi = z - zti

        rx[i] = rxi
        ry[i] = ryi
        rz[i] = rzi

        r[i]  = sqrt(rxi*rxi + ryi*ryi + rzi*rzi)
    end

    return xt, yt, zt, rx, ry, rz, r
end

"""
    compute_space_charge(filename::String) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Process track data from a Lardon output file, compute true positions and residuals for selected tracks, and return aggregated results for space charge analysis.

# Arguments
- `filename::String`: Path to the HDF5 file containing Lardon output data. The file must include `tracks3d` and `trk3d_v2` datasets.

# Returns
- A tuple containing:
  - `Xt`: Aggregated true x-coordinates (projected onto fitted lines) for all selected tracks.
  - `Yt`: Aggregated true y-coordinates (projected onto fitted lines) for all selected tracks.
  - `Zt`: Aggregated true z-coordinates (projected onto fitted lines) for all selected tracks.
  - `Xr`: Aggregated residuals in the x-direction for all selected tracks.
  - `Yr`: Aggregated residuals in the y-direction for all selected tracks.
  - `Zr`: Aggregated residuals in the z-direction for all selected tracks.
  - `R`: Aggregated total residual displacements (Euclidean norm) for all selected tracks.

# Details
1. **Data Loading**:
   - Reads `tracks3d` and `trk3d_v2` datasets from the HDF5 file.

2. **Track Selection**:
   - A track is selected if it meets all the following criteria:
     - `is_anode_crosser` is `true`.
     - `n_matched` > 1.
     - `d_match` < 2.5.
     - `t0_corr` < 9999.
     - Passes the `check_selection_vd` function with the track and its `z0_corr`.

3. **Line Fitting**:
   - For each selected track, a 3D line is fitted using `fit_line_3d_anode`.
   - Only tracks with an RMS < 0.9 are processed further.

4. **True Positions and Residuals**:
   - For each selected and well-fitted track, the true positions and residuals are computed using `compute_true_and_residuals`.
   - Results are aggregated into the output vectors.

# Notes
- Tracks that do not meet the selection criteria or have a poor fit (RMS ≥ 0.9) are skipped.
- The function assumes the HDF5 file structure matches the expected format (e.g., `tracks3d` and `trk3d_v2` datasets exist and contain the required fields).

# Examples
```jldoctest
Xt, Yt, Zt, Xr, Yr, Zr, R = compute_space_charge("path/to/lardon_output.h5")
# Returns aggregated true positions and residuals for all selected tracks in the file
"""
function compute_space_charge(filename::String)
    #read the lardon output
    dataTracks, dataHits = h5open(filename, "r") do f
        read(f["tracks3d"]), read(f["trk3d_v2"])
    end

    ntracks = length(dataTracks)
    #true and residuals
    Xt, Yt, Zt, Xr, Yr, Zr, R = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

    for itrk in 1:ntracks

        track3D = dataTracks[itrk]
        track   = dataHits[itrk]
        #selection 
        if track3D.is_anode_crosser && track3D.n_matched > 1 && track3D.d_match < 2.5 && track3D.t0_corr < 9999 &&
           check_selection_vd(track, track3D.z0_corr) 
            #linear fit if track pass selection
            res = fit_line_3d_anode(track, track3D.z0_corr )
            res === nothing && continue
            p0, v, rms = res
            #Only good fits
            if rms < 0.9

                xt, yt, zt, xr, yr, zr, r =
                    compute_true_and_residuals(track, p0, v;
                                               zcorr = track3D.z0_corr)

                append!(Xt, xt)
                append!(Yt, yt)
                append!(Zt, zt)
                append!(Xr, xr)
                append!(Yr, yr)
                append!(Zr, zr)
                append!(R,  r)
            end
        end
    end

    return Xt, Yt, Zt, Xr, Yr, Zr, R
end

"""
    process_run(dir::String) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Process all HDF5 files in a directory, compute true positions and residuals for selected tracks, and return aggregated results.

# Arguments
- `dir::String`: Path to the directory containing HDF5 files to process. Only files with the `.h5` extension are considered.

# Returns
- A tuple containing:
  - `Xt_all`: Aggregated true x-coordinates (projected onto fitted lines) for all selected tracks across all files.
  - `Yt_all`: Aggregated true y-coordinates (projected onto fitted lines) for all selected tracks across all files.
  - `Zt_all`: Aggregated true z-coordinates (projected onto fitted lines) for all selected tracks across all files.
  - `Xr_all`: Aggregated residuals in the x-direction for all selected tracks across all files.
  - `Yr_all`: Aggregated residuals in the y-direction for all selected tracks across all files.
  - `Zr_all`: Aggregated residuals in the z-direction for all selected tracks across all files.
  - `R_all`: Aggregated total residual displacements (Euclidean norm) for all selected tracks across all files.

# Details
1. **File Discovery**:
   - Lists all `.h5` files in the specified directory and sorts them alphabetically.

2. **Memory Allocation**:
   - Pre-allocates memory for the output vectors based on an estimated number of points per file (`estimate_per_file`).

3. **Processing**:
   - For each file, computes true positions and residuals using `compute_space_charge`.
   - Aggregates results into pre-allocated vectors.
   - Prints progress every 10 files.

4. **Memory Optimization**:
   - Trims unused space from the output vectors after processing all files.

# Notes
- The function assumes each file contains valid track data in the format expected by `compute_space_charge`.
- The estimated number of points per file (`estimate_per_file`) should be adjusted if the actual data size differs significantly.

# Examples
```jldoctest
Xt_all, Yt_all, Zt_all, Xr_all, Yr_all, Zr_all, R_all = process_run("/path/to/run_directory")
# Returns aggregated true positions and residuals for all selected tracks in all files
"""
function process_run(dir::String)
    
    #List all files in run directory
    files = sort(filter(f -> endswith(f, ".h5"),
                        readdir(dir, join=true)))

    nfiles = length(files)
    estimate_per_file = 35_000   # adjust if needed
    total_estimate = estimate_per_file * nfiles

    Xt_all = Vector{Float64}(undef, total_estimate)
    Yt_all = Vector{Float64}(undef, total_estimate)
    Zt_all = Vector{Float64}(undef, total_estimate)
    Xr_all = Vector{Float64}(undef, total_estimate)
    Yr_all = Vector{Float64}(undef, total_estimate)
    Zr_all = Vector{Float64}(undef, total_estimate)
    R_all  = Vector{Float64}(undef, total_estimate)

    idx = 0
    file_index = 0  # Track the file index

    for file in files
        file_index += 1  # Increment the file index
        # Print every 10 files
        if file_index % 10 == 1
            println("Processing: $(basename(file)) (file $file_index of $(length(files)))")
        end

        Xt, Yt, Zt, Xr, Yr, Zr, R = compute_space_charge(file)

        n = length(Xt)

        Xt_all[idx+1:idx+n] = Xt
        Yt_all[idx+1:idx+n] = Yt
        Zt_all[idx+1:idx+n] = Zt
        Xr_all[idx+1:idx+n] = Xr
        Yr_all[idx+1:idx+n] = Yr
        Zr_all[idx+1:idx+n] = Zr
        R_all[idx+1:idx+n]  = R

        idx += n
    end

    # Trim unused space
    resize!(Xt_all, idx)
    resize!(Yt_all, idx)
    resize!(Zt_all, idx)
    resize!(Xr_all, idx)
    resize!(Yr_all, idx)
    resize!(Zr_all, idx)
    resize!(R_all,  idx)

    return Xt_all, Yt_all, Zt_all, Xr_all, Yr_all, Zr_all, R_all
end

"""
    make_residual_heatmap(
        x::AbstractVector,
        y::AbstractVector,
        z::AbstractVector,
        dx::AbstractVector,
        dy::AbstractVector,
        dz::AbstractVector;
        axis = :y,
        y_edges,
        z_edges,
        x_edges,
        min_counts = 2,
        x_range = nothing,
        y_range = nothing
    ) -> Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Int}}

Compute binned averages of residuals (dx, dy, dz) as a function of a selected coordinate (x or y) and z, producing a residual heatmap.

# Arguments
- `x::AbstractVector`: Vector of x-coordinates.
- `y::AbstractVector`: Vector of y-coordinates.
- `z::AbstractVector`: Vector of z-coordinates.
- `dx::AbstractVector`: Vector of residuals in the x-direction.
- `dy::AbstractVector`: Vector of residuals in the y-direction.
- `dz::AbstractVector`: Vector of residuals in the z-direction.

# Keyword Arguments
- `axis = :y`: The coordinate axis to bin along. Must be either `:y` or `:x`.
- `y_edges`: Bin edges for the y-coordinate.
- `z_edges`: Bin edges for the z-coordinate.
- `x_edges`: Bin edges for the x-coordinate (used only if `axis = :x`).
- `min_counts = 2`: Minimum number of entries in a bin to compute the average. Bins with fewer entries are filled with `NaN`.
- `x_range = nothing`: Optional range `(xmin, xmax)` to filter x-coordinates. If `nothing`, no filtering is applied.
- `y_range = nothing`: Optional range `(ymin, ymax)` to filter y-coordinates. If `nothing`, no filtering is applied.

# Returns
- A tuple containing:
  - `avg_dx`: Matrix of average x-residuals for each bin.
  - `avg_dy`: Matrix of average y-residuals for each bin.
  - `avg_dz`: Matrix of average z-residuals for each bin.
  - `counts`: Matrix of counts for each bin.

# Details
1. **Binning**:
   - The function bins the data along the selected `axis` and `z`.
   - For each point, it checks if the residuals are finite and if the point falls within the specified `x_range` and `y_range` (if provided).

2. **Averaging**:
   - For each bin, the function computes the average of the residuals (`dx`, `dy`, `dz`) if the bin contains at least `min_counts` entries.
   - Bins with fewer than `min_counts` entries are filled with `NaN`.

3. **Efficiency**:
   - The function uses a single-pass algorithm to accumulate sums and counts, minimizing memory allocations.

# Examples
```jldoctest
x = rand(1000)
y = rand(1000)
z = rand(1000)
dx = randn(1000)
dy = randn(1000)
dz = randn(1000)

y_edges = range(0, 1, length=11)
z_edges = range(0, 1, length=11)

avg_dx, avg_dy, avg_dz, counts = make_residual_heatmap(x, y, z, dx, dy, dz; axis=:y, y_edges, z_edges)
# Returns binned averages of residuals and counts
"""
function make_residual_heatmap(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    dx::AbstractVector,
    dy::AbstractVector,
    dz::AbstractVector;
    axis = :y,
    y_edges,
    z_edges,
    x_edges,
    min_counts = 2,
    x_range = nothing,   # NEW
    y_range = nothing    # NEW
)

    N = length(x)

    # --- Select binning coordinate ---
    if axis == :y
        bin_edges = y_edges
        coord = y
    elseif axis == :x
        bin_edges = x_edges
        coord = x
    else
        throw(ArgumentError("axis must be :y or :x"))
    end

    nbins = length(bin_edges) - 1
    nbins_z = length(z_edges) - 1

    # --- Initialize accumulators ---
    sum_dx = zeros(nbins_z, nbins)
    sum_dy = zeros(nbins_z, nbins)
    sum_dz = zeros(nbins_z, nbins)
    counts = zeros(Int, nbins_z, nbins)

    # ==========================================================
    #                   Single-pass binning
    # ==========================================================

    @inbounds for i in 1:N

        # Skip invalid residuals
        if !isfinite(dx[i]) || !isfinite(dy[i]) || !isfinite(dz[i])
            continue
        end

        xi = x[i]
        yi = y[i]
        zi = z[i]

        # --- Apply spatial cuts (no allocations) ---
        if x_range !== nothing
            (xmin, xmax) = x_range
            if xi < xmin || xi > xmax
                continue
            end
        end

        if y_range !== nothing
            (ymin, ymax) = y_range
            if yi < ymin || yi > ymax
                continue
            end
        end

        binval = coord[i]

        ibin = searchsortedlast(bin_edges, binval)
        iz   = searchsortedlast(z_edges, zi)

        if 1 ≤ ibin ≤ nbins && 1 ≤ iz ≤ nbins_z
            sum_dx[iz, ibin] += dx[i]
            sum_dy[iz, ibin] += dy[i]
            sum_dz[iz, ibin] += dz[i]
            counts[iz, ibin] += 1
        end
    end

    # ==========================================================
    #                   Compute averages
    # ==========================================================

    avg_dx = fill(NaN, nbins_z, nbins)
    avg_dy = fill(NaN, nbins_z, nbins)
    avg_dz = fill(NaN, nbins_z, nbins)

    @inbounds for iz in 1:nbins_z, ibin in 1:nbins
        c = counts[iz, ibin]
        if c ≥ min_counts
            inv = 1.0 / c
            avg_dx[iz, ibin] = sum_dx[iz, ibin] * inv
            avg_dy[iz, ibin] = sum_dy[iz, ibin] * inv
            avg_dz[iz, ibin] = sum_dz[iz, ibin] * inv
        end
    end

    return avg_dx, avg_dy, avg_dz, counts
end
