"""
Checks that every channel in `channels` belongs to a complete pair present in `channels`.
Rejects events with "orphaned" channels from incomplete modules.
"""
function has_only_full_modules(channels, module_pairs)
    ch = Set(Int.(channels))
    
    # 1. Build a lookup mapping each channel to its paired partner
    # e.g., if (2,3) is a pair, partner_map[2] = 3 and partner_map[3] = 2
    partner_map = Dict{Int, Int}()
    for (a, b) in module_pairs
        partner_map[Int(a)] = Int(b)
        partner_map[Int(b)] = Int(a)
    end
    
    # 2. Check that for EVERY active channel, its partner is also in 'ch'
    return all(get(partner_map, c, -1) in ch for c in ch)
end


"""
    has_full_module(channels, module_pairs)

Check whether both channels of at least one module (pair) are present in `channels`.

# Arguments
- `channels`: Collection of active/detected channel IDs.
- `module_pairs`: Collection of tuples or vectors `(ch_a, ch_b)` representing paired channels for each module.

# Returns
- `true` if any pair `(a, b)` has both `a` and `b` present in `channels`, `false` otherwise.
"""
function has_full_module(channels, module_pairs)
    # Convert channel list to a Set for O(1) lookup speed
    ch_set = Set{Int}(channels)
    
    # Check if BOTH channels (a AND b) of any pair are present in the set.
    # The `any` function short-circuits as soon as a valid pair is found.
    return any(a in ch_set && b in ch_set for (a, b) in module_pairs)
end

raw"""
    find_pds_flashes_lowE(hits; dtick=1, dpeak=10.0)

Cluster low-energy PDS (Photon Detection System) hits into reconstructed optical flashes based on spatial/temporal criteria.

# Arguments
- `hits`: Vector or StructArray of hit objects.
- `dtick::Int`: Maximum allowed tick gap between the stop time of the current cluster and the start time of a candidate hit (default: `1`).
- `dpeak::Float64`: Maximum allowed difference between a hit's peak time (`max_t`) and the cluster's running mean peak time (default: `10.0`).

# Filtering Criteria
- `timestamp`: Must be between 0 and 1700 ticks.
- `glob_ch`: Global channel ID must be less than 16.
- Hit duration (`stop - start`): $\le 350$ ticks.
- Time to peak (`stop - max_t`): $\le 300$ ticks.
- `npe`: Must be $< 2200$ photoelectrons.
- Saturation: Excludes hits marked as saturated (`!h.saturates`).

# Returns
- A `Vector{PDSFlash}` containing the reconstructed clusters.
"""
function find_pds_flashes_lowE(hits; dtick::Int=1, dpeak::Float64=10.0, tmin::Real=0, tmax::Real=1700, qMax::Real = 2200, pDuration::Real=350)
    # -------------------------------------------------------------------------
    # 1. Filter Hits
    # -------------------------------------------------------------------------
    # Collect indices of hits that satisfy physical low-energy quality cuts
    valid_idx = filter(i -> begin
        h = hits[i]
        h.timestamp > tmin &&
        h.timestamp < tmax &&    
        h.glob_ch < 16 && #Only look at cathode channels
        (h.stop - h.start) <= pDuration &&
        (h.stop - h.max_t) <= 300 &&
        h.npe < qMax &&
        !h.saturates
    end, eachindex(hits))

    # Early return if no hits pass selection
    if isempty(valid_idx)
        return PDSFlash[]
    end

    # -------------------------------------------------------------------------
    # 2. Sort Hits
    # -------------------------------------------------------------------------
    # Sort remaining indices hierarchically by event ID, trigger ID, and start time
    order = sortperm(valid_idx, by = i -> (hits[i].event, hits[i].trigger, hits[i].start))
    sorted_valid_idx = valid_idx[order]

    # Output storage
    flashes = PDSFlash[]
    cluster_hits = Int[] 

    # Cluster accumulator state
    cluster_event = UInt32(0)
    cluster_trigger = UInt32(0)
    cluster_start = 0
    cluster_stop = 0

    npe_sum = 0.0
    charge_sum = 0.0
    weighted_time = 0.0
    
    # Running sum for dynamic average calculation without allocating a vector
    peak_time_sum = 0.0

    # -------------------------------------------------------------------------
    # Helper: Flush Active Cluster
    # -------------------------------------------------------------------------
    function flush_cluster()
        isempty(cluster_hits) && return

        # Retrieve original indices in `hits` corresponding to this cluster
        orig_indices = sorted_valid_idx[cluster_hits]
        ts = [hits[i].max_t for i in orig_indices]

        push!(flashes,
            PDSFlash(
                cluster_event,
                cluster_trigger,
                cluster_start,
                cluster_stop,
                minimum(ts),
                maximum(ts),
                median(ts),
                weighted_time / npe_sum, # NPE-weighted mean peak time
                npe_sum,
                charge_sum,
                length(cluster_hits),
                [hits[i].glob_ch for i in orig_indices],
                orig_indices  
            )
        )
    end

    # -------------------------------------------------------------------------
    # 3. Clustering Loop
    # -------------------------------------------------------------------------
    for i in eachindex(sorted_valid_idx)
        orig_idx = sorted_valid_idx[i]
        h = hits[orig_idx]

        # Initialize the very first cluster
        if isempty(cluster_hits)
            push!(cluster_hits, i)
            
            cluster_event = h.event
            cluster_trigger = h.trigger
            cluster_start = h.start
            cluster_stop = h.stop
            
            npe_sum = h.npe
            charge_sum = h.charge
            weighted_time = h.npe * h.max_t
            peak_time_sum = Float64(h.max_t)
            continue
        end

        # Match criteria: same readout event/trigger, continuous/overlapping time window
        same_group = (h.event == cluster_event && h.trigger == cluster_trigger)
        overlaps   = (h.start <= cluster_stop + dtick)
        
        # O(1) Dynamic check: Compare peak time to running mean using scalar sum
        current_peak_mean = peak_time_sum / length(cluster_hits)
        peaks_aligned     = abs(h.max_t - current_peak_mean) <= dpeak

        if same_group && overlaps && peaks_aligned
            # --- Expand Current Cluster ---
            push!(cluster_hits, i)
            
            cluster_start = min(cluster_start, h.start)
            cluster_stop  = max(cluster_stop, h.stop)
            npe_sum       += h.npe
            charge_sum    += h.charge
            weighted_time += h.npe * h.max_t
            peak_time_sum += h.max_t
        else
            # --- Flush Existing Cluster & Reset State ---
            flush_cluster()
            empty!(cluster_hits)
            
            push!(cluster_hits, i)
            
            cluster_event = h.event
            cluster_trigger = h.trigger
            cluster_start = h.start
            cluster_stop = h.stop
            
            npe_sum       = h.npe
            charge_sum    = h.charge
            weighted_time = h.npe * h.max_t
            peak_time_sum = Float64(h.max_t)
        end
    end

    # Flush final trailing cluster
    flush_cluster()
    return flashes
end
