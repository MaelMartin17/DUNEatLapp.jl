"""
function Single_Hits_lardon(path::String)
function to collect the data about the single hits found by Lardon.
It receives the path of the HDF5 files in your directory and it returns a DataFrame with all the single Hits and the information about
"""
function Single_Hits_lardon(path::String)
    # Create the output Dataframe with the wanted variables
    Info = DataFrame(
    Id=Int32[],
    mod=Int32[],
    x=Float32[],
    y=Float32[],
    z=Float32[],
    nHits=Int32[],
    timestamp=Float32[],
    tdc_max=Vector[],
    tdc_min=Vector[],
    charge_neg=Vector[],
    charge_pos=Vector[],
    d_bary_max=Float32[],
    d_track_2D=Float32[],
    veto=Vector[],
    fc_max=Float32[],
    tdc_start=Vector[],
    tdc_stop=Vector[],
    t_event=Float64[],
    tns_event=Float64[]);
    
    # Run over the folder with the data files
    for (root, dirs, files) in walkdir(path)
        for i in ProgressBar(files)
            fid = h5open("$(root)/$(i)", "r")
            
            # Open the different trees in the HDF5 files
            SHset = fid["single_hits"]
            A = read(SHset);
            Hset = fid["hits"]
            B = read(Hset);            
            Eset = fid["event"]
            C = read(Eset);   
            
            for i in 1:1:length(A)

                # Add values of each event in the Dataframe
                push!(Info,[
                A[i].ID,
                A[i].module,
                A[i].x, 
                A[i].y, 
                A[i].z, 
                A[i].n_hits[3], 
                A[i].timestamp, 
                A[i].tdc_max, 
                A[i].tdc_min, 
                A[i].charge_neg, 
                A[i].charge_pos, 
                A[i].d_bary_max, 
                A[i].d_track_2D, 
                A[i].n_veto,
                B[(A[i].hit_IDs[7])+1].fC_max,
                A[i].tdc_start,
                A[i].tdc_stop,
                C[(A[i].event)+1].time_s,
                C[(A[i].event)+1].time_ns])  
            end    
        end
    end
    return Info
end
#_______________________________________________________________________________________________________________________
"""
function tracks3d_lardon(path::String)
function to collect the data about the 3d tracks found by Lardon.
It receives the path of the directory containing the HDF5 files and it returns a DataFrame with all the selected 3d tracks and the information about.
"""
function tracks3d_lardon(path::String)
    # Create the output Dataframe with the wanted variables
    Info = DataFrame(
    Id_file=Int32[],
    Id=Int32[],
    n_hits=Vector[],
    n_matched=Int32[],
    phi_end=Float32[],
    phi_ini=Float32[],
    t0_corr=Float32[],
    t_end=Int32[],
    t_ini=Int32[],
    theta_end=Float32[],
    theta_ini=Float32[],
    total_charge=Vector[],
    x_end=Float32[],
    x_ini=Float32[],
    y_end=Float32[],
    y_ini=Float32[],
    z0_corr=Float64[],
    z_end=Float32[],
    z_end_overlap=Float32[],
    z_ini=Float32[],
    z_ini_overlap=Float32[],
    exit_point=Vector[]
    )
    
    compteur_file = 0
    
    # Run over the folder with the data files
    for (root, dirs, files) in walkdir(path)
        for n in ProgressBar(files)
            fid = h5open("$(root)/$(n)", "r")
            
            compteur_file += 1
            
            # Open the tracks3d tree in the HDF5 files
            track3d_set = fid["tracks3d"]
            data_track3d = read(track3d_set)
            
            for i in 1:1:length(data_track3d)
                
                # Application of cuts to select particular tracks (cf. https://github.com/dune-lardon/cookbook/blob/main/analysis_example/ana_hit3D.py)
                 
                cut1 = data_track3d[i].t0_corr < 9999;
                cut2 = data_track3d[i].t0_corr > 0;
                cut3 = data_track3d[i].theta_ini > 100;
                cut4 = maximum(data_track3d[i].len_straight) > 20.
                cut5 = maximum(data_track3d[i].n_hits) > 10
                cut6 = data_track3d[i].n_matched > 1
                cut7 = data_track3d[i].d_match < 2.5
                cut8 = data_track3d[i].is_cathode_crosser == false
	        cut9 = data_track3d[i].is_anode_crosser == true

        	if cut1 & cut2 & cut3 & cut4 & cut5 & cut6 & cut7 & cut8 & cut9
        	    
        	    # Add values of each event in the Dataframe
                    push!(Info,[
                        compteur_file,
                        data_track3d[i].ID,
                        data_track3d[i].n_hits,
                        data_track3d[i].n_matched, 
                        data_track3d[i].phi_end, 
                        data_track3d[i].phi_ini, 
                        data_track3d[i].t0_corr, 
                        data_track3d[i].t_end, 
                        data_track3d[i].t_ini, 
                        data_track3d[i].theta_end, 
                        data_track3d[i].theta_ini, 
                        data_track3d[i].total_charge, 
                        data_track3d[i].x_end, 
                        data_track3d[i].x_ini, 
                        data_track3d[i].y_end,
                        data_track3d[i].y_ini,
                        data_track3d[i].z0_corr,
                        data_track3d[i].z_end,
                        data_track3d[i].z_end_overlap,
                        data_track3d[i].z_ini,
                        data_track3d[i].z_ini_overlap,
                        data_track3d[i].exit_point
                    ])
                end
            end
                
        end
    end
    return Info
end

#_______________________________________________________________________________________________________________________
"""
function hits_in_tracks3dNew(path::String)
function to collect the data about the hits of a 3d track.
It receives the path of the directory containing the HDF5 files and it returns a DataFrame with all the hits of the selected 3d tracks and the information about.
"""
function hits_in_tracks3dNew(path::String)
    # Create the output Dataframe with the wanted variables
    Info = DataFrame(
        Id_file = Int32[],
        Id_trk3d=Int32[], 
        x=Float32[],
        y=Float32[],
        z=Float32[],
        charge_pos=Float32[], 
        Id_hit=Int32[]
    )


    compteur_file = 0

    # Run over the folder with the data files
    for (root, dirs, files) in walkdir(path)
        for n in ProgressBar(files)
            fid = h5open("$(root)/$(n)", "r")
        
            compteur_file += 1
            
            # Open the trees in the HDF5 files
            track3d_set = fid["tracks3d"]
            data_track3d = read(track3d_set)

            trk3d_v2_set = fid["trk3d_v2"]
            data_trk3d_v2 = read(trk3d_v2_set)

            for i in 1:1:length(data_trk3d_v2)

                # Application of cuts to select particular tracks (cf. https://github.com/dune-lardon/cookbook/blob/main/analysis_example/ana_hit3D.py)                
                cut1 = data_track3d[i].t0_corr < 9999;
                cut2 = data_track3d[i].t0_corr > 0;
                cut3 = data_track3d[i].theta_ini > 100;
                cut4 = maximum(data_track3d[i].len_straight) > 20.
                cut5 = maximum(data_track3d[i].n_hits) > 10
                cut6 = data_track3d[i].n_matched > 1
                cut7 = data_track3d[i].d_match < 2.5
                cut8 = data_track3d[i].is_cathode_crosser == false
	        cut9 = data_track3d[i].is_anode_crosser == true

        	if cut1 & cut2 & cut3 & cut4 & cut5 & cut6 & cut7 & cut8 & cut9
        
                    for j in 1:1:length(data_trk3d_v2[i])

                        track3d_Id = data_track3d[i].ID
                        x = data_trk3d_v2[i][j][1]
                        y = data_trk3d_v2[i][j][2]
                        z = data_trk3d_v2[i][j][3]
                        charge_pos = data_trk3d_v2[i][j][4]
                        Hits_Id = Int(data_trk3d_v2[i][j][6]) 
                    	
                    	# Add values of each event in the Dataframe
                        push!(Info,[compteur_file, track3d_Id, x, y, z, charge_pos, Hits_Id])
                    end
         
                end
            end
        end      
    end
    return Info
end
#_______________________________________________________________________________________________________________________
