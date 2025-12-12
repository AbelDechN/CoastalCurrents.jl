"""
This module was written for the BlueCloud2026 project, it will help you to generate surface currents
"""
module CoastalCurrents
using NCDatasets
using PhysOcean
using Dates
using Statistics
using DataStructures

# function to load all data
function loaddata(fname::AbstractString);
    ds = NCDataset(fname)
    @show fname
    time = ds["TIME"][:]
    time_qc = ds["TIME_QC"][:]

    lon = ds["LONGITUDE"][:];
    position_qc = ds["POSITION_QC"][:];
    lat = ds["LATITUDE"][:]
#    z = ds["DEPH"][:]
#    z_qc = ds["DEPH_QC"][:]

    u = ds["EWCT"][1,:]
    u_qc = ds["EWCT_QC"][1,:]

    v = ds["NSCT"][1,:]
    v_qc = ds["NSCT_QC"][1,:]

    lon = mod.(lon .+ 180,360) .- 180

    good(qc) = !ismissing(qc) && ((qc == 1) || (qc == 2))

    lon[.!good.(position_qc)] .= NaN;
    @show length(lon)
    lat[.!good.(position_qc)] .= NaN;
    @show length(lat)
#    z[.!good.(z_qc)] .= NaN;
    time[.!good.(time_qc)] .= DateTime(9999,1,1);
    u[.!good.(u_qc)] .= NaN;
    v[.!good.(v_qc)] .= NaN;

 #   lon = repeat(reshape(lon,(1,:)),size(z,1))
 #   lat = repeat(reshape(lat,(1,:)),size(z,1))
 #   time = repeat(reshape(time,(1,:)),size(z,1))
    close(ds)
    return lon[:],lat[:],time[:],u[:],v[:]
end



function loaddata(files::AbstractVector{<:AbstractString});
    data = loaddata.(files);
    # concatenate all profiles
    lon = reduce(vcat,getindex.(data,1))
    lat = reduce(vcat,getindex.(data,2))
#    z = reduce(vcat,getindex.(data,3))
    time = reduce(vcat,getindex.(data,4))
    u = reduce(vcat,getindex.(data,5))
    v = reduce(vcat,getindex.(data,6))

    return lon[:],lat[:],time[:],u[:],v[:]
end


include("hfradar.jl")
include("plotting.jl");
include("Altimetry.jl");
include("utils.jl")


end # module CoastalCurrents
