"""
function Single_Hits_lardon(path::String)
function to collect the data about the single hits found by Lardon.
It receives the path of the HDF5 files in your directory and it returns a DataFrame with all the single Hits and the information about
"""
function Single_Hits_lardon(path::String)
    Info = DataFrame(Id=Int32[],apa=Int32[], x=Float32[],y=Float32[],z=Float32[],timestamp=Float32[],tdc_max=Vector[],tdc_min=Vector[],charge_neg=Vector[],charge_pos=Vector[],d_bary_max=Float32[],d_track_2D=Float32[],veto=Vector[],fc_max=Float32[],tdc_start=Vector[],tdc_stop=Vector[],t_event=Float64[]);
    for (root, dirs, files) in walkdir(path)
        for i in ProgressBar(files)
            fid = h5open("$(root)/$(i)", "r")
            SHset = fid["single_hits"]
            A = read(SHset);
            Hset = fid["hits"]
            B = read(Hset);            
            Eset = fid["event"]
            C = read(Eset);   
            for i in 1:1:length(A)
                push!(Info,[A[i].ID, A[i].module,A[i].x, A[i].y, A[i].z, A[i].timestamp, A[i].tdc_max, A[i].tdc_min, A[i].charge_neg, A[i].charge_pos, A[i].d_bary_max, A[i].d_track_2D, A[i].veto,B[(A[i].hit_IDs[7])+1].fC_max,A[i].tdc_start,A[i].tdc_stop,C[(A[i].event)+1].time_s])  
            end    
        end
    end
    return Info
end
#_______________________________________________________________________________________________________________________

function Single_Hits_lardon_old(path::String)
    Info = DataFrame(Id=Int32[], x=Float32[],y=Float32[],z=Float32[],t=Float32[],tdc_max=Vector[],tdc_min=Vector[],charge_neg=Vector[],charge_pos=Vector[],d_bary_max=Float32[],d_track_2D=Float32[],veto=Vector[],fc_max=Float32[],tdc_start=Vector[],tdc_stop=Vector[]);
    for (root, dirs, files) in walkdir(path)
        for i in ProgressBar(files)
            fid = h5open("$(root)/$(i)", "r")
            SHset = fid["single_hits"]
            A = read(SHset);
            Hset = fid["hits"]
            B = read(Hset);            
            for i in 1:1:length(A)
                push!(Info,[A[i].ID,A[i].x, A[i].y, A[i].z, A[i].timestamp, A[i].tdc_max, A[i].tdc_min, A[i].charge_neg, A[i].charge_pos, A[i].d_bary_max, A[i].d_track_2D, A[i].veto,B[(A[i].ID)+1].fC_max,A[i].tdc_start,A[i].tdc_stop]) 
            end    
        end
    end
    return Info
end
