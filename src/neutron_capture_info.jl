"""
function get_capture_position(name_file::String)
function to get the capture positions of neutrons
It will take a xxx_nt_Gammas.csv to return the position of neutron captures using the first gamma of each decay
"""
function get_capture_position(name_file::String)
    data = readdlm(IOBuffer(replace(read(name_file, String), ";" => ",")), ',','\n',comments=:true)
    data = Float32.(data[:,3:6])
    data = data[ data[:,1] .==0 ,: ]
    return data[:,2:4]
end

"""
function get_rate_neutron_captures_Ar(my_file::String,name_primary::String,fidu::Real)
function to get the rate of neutrons that are captured in LAr and in a fiducial volume of LAr 
It accepts a String for my_file, a String for name_primary and a Real in centimeters for fidu. It returns two floats.
"""
function get_rate_neutron_captures_Ar(my_file::String,name_primary::String,fidu::Real=100.)
    df_neutrons = CSV.read(my_file, DataFrame,comment="#",drop=[:evt,:A,:pdg,:E,:t],header=["evt","proc","Z","A","pdg","E","x","y","z","t"])
    n_neutrons = get_n_primaries(name_primary)
    n_capture_Ar = 0
    n_capture_Ar_fidu = 0
    FD_x_size = 3100.
    FD_y_size = 755.
    FD_z_size = 700.

    for i = 1 : 1 : length(df_neutrons[!,1])
         if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18
            n_capture_Ar += 1
        end
        if df_neutrons[i,:proc] == "nCapture" && df_neutrons[i,:Z] == 18 && abs(df_neutrons[i,:x]) < (FD_x_size - fidu) && abs(df_neutrons[i,:y]) < (FD_y_size - fidu) && abs(df_neutrons[i,:z]) < (FD_z_size - fidu)
            n_capture_Ar_fidu += 1
        end
    end
    return n_capture_Ar/n_neutrons, n_capture_Ar_fidu/n_neutrons
end


"""
function get_info_neutron_captures_Ar(my_file::String,name_primary::String)
function to obtain the information of neutrons that were captured in argon
It accepts the names of the files containing the primary info and the secondaries
It return a matrix with info of the primary neutron and position of origin and capture
"""
function get_info_neutron_captures_Ar(my_file::String,name_primary::String)
    df_secondaries = CSV.read(my_file, DataFrame,comment="#",drop=[:A,:pdg,:E,:t],header=["evt","proc","Z","A","pdg","E","x","y","z","t"])
    df_primary = get_primary_vertex(name_primary)
    n_neutrons = df_primary[end,1]
    #println(n_neutrons)
    nuclei_Z = []
    nuclei_A = []
    n_evt = []
    x = []
    y = []
    z = []
    En = []
    xp = []
    yp = []
    zp = []
    for i = 1 : 1 : length(df_secondaries[!,1])
         if df_secondaries[i,:proc] == "nCapture" && df_secondaries[i,:Z] == 18
            evt = df_secondaries[i,:evt] + 1
            push!(x,df_secondaries[i,:x])
            push!(y,df_secondaries[i,:y])
            push!(z,df_secondaries[i,:z])
            push!(En,df_primary[evt,:E])
            push!(xp,df_primary[evt,:x]/10)
            push!(yp,df_primary[evt,:y]/10)
            push!(zp,df_primary[evt,:z]/10)
        end
    end
    return [x y z En xp yp zp]
end
