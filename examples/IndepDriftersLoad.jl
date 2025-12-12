function IndepDriftersLoad()

# Loading of the drifter's variables

londo,latdo,timedo,udo,vdo = CoastalCurrents.loaddata("/home/jovyan/tmp/BlueCloud2026/Drifter/my.cmems-du.eu/Core/INSITU_GLO_PHY_UV_DISCRETE_MY_013_044/cmems_obs-ins_glo_phy-cur_my_drifter_PT6H/history/Dr1.nc"); 

# Creation of drifter ID variable

ds_id = NCDataset("/home/jovyan/tmp/BlueCloud2026/Drifter/my.cmems-du.eu/Core/INSITU_GLO_PHY_UV_DISCRETE_MY_013_044/cmems_obs-ins_glo_phy-cur_my_drifter_PT6H/history/Dr1.nc")
drifter_id = ds_id["DRIFTER_ID"][:];
close(ds_id)

# Taking the corresponing longitude and latitude

indices_a_supprimer_long = findall(x -> !(x >= -5 && x <= 35.5), londo);
#indices_a_supprimer_long = findall(x -> !(x >= lonmin && x <= lonmax), londo);    Marche pas car taille du dataset idvalid est pas bon


deleteat!(londo, indices_a_supprimer_long);
deleteat!(latdo, indices_a_supprimer_long);
deleteat!(timedo, indices_a_supprimer_long);
deleteat!(udo, indices_a_supprimer_long);
deleteat!(vdo, indices_a_supprimer_long);
deleteat!(drifter_id, indices_a_supprimer_long);

indices_a_supprimer_lat = findall(y -> !(y >= 30 && y <= 45), latdo);
#indices_a_supprimer_lat = findall(y -> !(y >= latmin && y <= latmax), latdo);     Marche pas car taille du dataset idvalid est pas bon


deleteat!(londo, indices_a_supprimer_lat)
deleteat!(latdo, indices_a_supprimer_lat);
deleteat!(timedo, indices_a_supprimer_lat);
deleteat!(udo, indices_a_supprimer_lat);
deleteat!(vdo, indices_a_supprimer_lat);
deleteat!(drifter_id, indices_a_supprimer_lat);

# Creation of dsValid
# Repeated cell in data_prep

dsValid = NCDataset("/home/jovyan/CoastalCurrents/examples/ValidIndices.nc")

for_cv = dsValid["for_cv"][:] 

close(dsValid)

for_cvR = round.(Int, for_cv);

# Separation through indices of "for_cvR"

udov = udo[for_cvR .== 1];
vdov = vdo[for_cvR .== 1];
londov = londo[for_cvR .== 1];
latdov = latdo[for_cvR .== 1];
timedov = timedo[for_cvR .== 1];

# Retrieving NaNs

londov_ = londov[.!isnan.(udov)];
latdov_ = latdov[.!isnan.(udov)];
timedov_ = timedov[.!isnan.(udov)];
udov_ = udov[.!isnan.(udov)];
vdov_ = vdov[.!isnan.(vdov)];
    
    return(londov_,latdov_,timedov_,udov_,vdov_)
    
end