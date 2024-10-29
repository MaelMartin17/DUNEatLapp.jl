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
