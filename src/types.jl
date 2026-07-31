struct PDSFlash
    event::UInt32
    trigger::UInt32

    start::Int          # earliest pulse start
    stop::Int           # latest pulse stop

    tmin::Int           # earliest max_t
    tmax::Int           # latest max_t
    tmedian::Float64    # median peak time
    tmean::Float64      # NPE-weighted mean peak time

    npe::Float64
    charge::Float64

    nhits::Int
    channels::Vector{UInt32}
    hit_indices::Vector{Int}
end
