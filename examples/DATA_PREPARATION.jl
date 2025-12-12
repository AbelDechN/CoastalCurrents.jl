function DataPreparation(year)



##########################################################################################################
#                                  Altimetry
##########################################################################################################


dsa = NCDataset(altimetry_fname);
mdt = NCDatasets.loadragged(dsa["mdt"],:);
sla = NCDatasets.loadragged(dsa["slaf"],:);
lona_ = NCDatasets.loadragged(dsa["lon"],:);
lata_ = NCDatasets.loadragged(dsa["lat"],:);
timea_ = NCDatasets.loadragged(dsa["time"],:);
ida = NCDatasets.loadragged(dsa["id"],:);
close(dsa);


# Compute Absolute Dynamic topography    
    
adt = mdt .+ sla;

# To compute velocities through geostrophic currents

(lona,lata,timea,ua,va) = Altimetry.geostrophic_velocity(lona_,lata_,timea_,adt);

# Resizing into vertical vector foe easier computations

    
lona2 = reduce(vcat,lona);
lata2 = reduce(vcat,lata);
timea2 = reduce(vcat,timea);
ua2 = reduce(vcat,ua);
va2 = reduce(vcat,va);


# Main vector variable definition

n = 12  # Number of month in a year
    
lona3 = Vector{Vector{Float64}}(undef, n)
lata3 = Vector{Vector{Float64}}(undef, n)
timea3 = Vector{Vector{Float64}}(undef, n)
ua3 = Vector{Vector{Float64}}(undef, n)
va3 = Vector{Vector{Float64}}(undef, n)
ua4_ = Vector{Vector{Float64}}(undef, n)
va4_ = Vector{Vector{Float64}}(undef, n)
robsa = Vector{Vector{Float64}}(undef, n)
directionobsa = Vector{Vector{Float64}}(undef, n)
xa = Vector{Vector{Float64}}(undef, n)
ya = Vector{Vector{Float64}}(undef, n)
robsa_ = Vector{Vector{Float64}}(undef, n)
directionobsa_ = Vector{Vector{Float64}}(undef, n)
xa_ = Vector{Vector{Float64}}(undef, n)
ya_ = Vector{Vector{Float64}}(undef, n)


# Time selection
    
for i in 1:12

    if i == 2
        e = 28
    
sela = @. Dates.DateTime(year,i,1) <= timea2 < Dates.DateTime(year,i,e);

lona3[i] = lona2[sela];
lata3[i] = lata2[sela];
#timea3[1] = timea2[sela];
ua3[i] = ua2[sela];
va3[i] = va2[sela];


# Removing NaNs in velocities vectors
ua4_[i] = ua3[i][.!isnan.(ua3[i])];            
va4_[i] = va3[i][.!isnan.(va3[i])];            
            
# Robs and DirectionObs computation
robsa[i] = sqrt.(ua3[i].^2 + va3[i].^2)
directionobsa[i] = atand.(ua3[i],va3[i]);

# Grid positions
xa[i] = lona3[i];
ya[i] = lata3[i];

# Retrieving of NaNs
robsa_[i] = robsa[i][.!isnan.(robsa[i])];
ya_[i] = ya[i][.!isnan.(robsa[i])];
xa_[i] = xa[i][.!isnan.(robsa[i])];
directionobsa_[i] = directionobsa[i][.!isnan.(robsa[i])];
    
    else
        e = 30
   
# Time selection   
sela = @. Dates.DateTime(year,i,1) <= timea2 < Dates.DateTime(year,i,e);

lona3[i] = lona2[sela];
lata3[i] = lata2[sela];
#timea3[1] = timea2[sela];
ua3[i] = ua2[sela];
va3[i] = va2[sela];

            
# Removing NaNs in velocities vectors
ua4_[i] = ua3[i][.!isnan.(ua3[i])];            
va4_[i] = va3[i][.!isnan.(va3[i])];
            
# Robs and DirectionObs computation
robsa[i] = sqrt.(ua3[i].^2 + va3[i].^2)
directionobsa[i] = atand.(ua3[i],va3[i]);

# Grid position
xa[i] = lona3[i];
ya[i] = lata3[i];

# Retrieving NaNs
robsa_[i] = robsa[i][.!isnan.(robsa[i])];
ya_[i] = ya[i][.!isnan.(robsa[i])];
global xa_[i] = xa[i][.!isnan.(robsa[i])];
directionobsa_[i] = directionobsa[i][.!isnan.(robsa[i])];
     end
end
    
    

##########################################################################################################
#                                            HFRadar
##########################################################################################################


#file des fichiers .nc dans filesHFRadar

filesHFRadar = String[]; 
for (root, dirs, files2) in walkdir(hf_dir); #le directory précisé dans le common.jl dans examples
    for file in files2
        push!(filesHFRadar,joinpath(root, file)); # add path to files (names in the string)
    end
end


filesHFRadar_avg = filesHFRadar[:];



lonh,lath,timeh,rdvah,drvah = CoastalCurrents.HFRadar.loaddata(filesHFRadar_avg);



# Retrieving NaNs
maskhf = .!isnan.(rdvah)

rdvah2 = rdvah[maskhf];
drvah2 = drvah[maskhf];
timeh2 = timeh[maskhf];
lath2 = lath[maskhf];
lonh2 = lonh[maskhf];

n = 12  # Longueur du vecteur

lonhf3 = Vector{Vector{Float64}}(undef, n)
lathf3 = Vector{Vector{Float64}}(undef, n)
timehf3 = Vector{Vector{Float64}}(undef, n)
rdvahf3 = Vector{Vector{Float64}}(undef, n)
drvahf3 = Vector{Vector{Float64}}(undef, n)
robs_HF1 = Vector{Vector{Float64}}(undef, n)
directionobs_HF1 = Vector{Vector{Float64}}(undef, n)
epsilon2_HF = Vector{Vector{Int}}(undef, n)
x_HF1 = Vector{Vector{Float64}}(undef, n)
y_HF1 = Vector{Vector{Float64}}(undef, n)



for i in 1:12

# Time selection

    if i == 2
        e = 28
        
        
seldh = (timeh2 .>= Dates.DateTime(year,i,1)) .& (timeh2 .<= Dates.DateTime(year,i,e))

lonhf3[i] = lonh2[seldh];
lathf3[i] = lath2[seldh];
#timehf3[1] = timeh2[seldh];
rdvahf3[i] = rdvah2[seldh];
drvahf3[i] = drvah2[seldh];

#parametrisation 

robs_HF1[i] = (rdvahf3[i]);
#robs_HF2 = Float64.(robshc)
#robs_HF = [robs_HF1;robs_HF2]

directionobs_HF1[i] = (drvahf3[i]);
#directionobs_HF2 = atand.(uhc,vhc);
#directionobs_HF = [directionobs_HF1;directionobs_HF2]

epsilon2_HF = 1;

x_HF1[i] = (lonhf3[i]);
#x_HF2 = lonhc6;
#x_HF = [x_HF1;x_HF2]

y_HF1[i] = (lathf3[i]);
#y_HF2 = lathc6;
#y_HF = [y_HF1;y_HF2];



    else
        e = 30
    

seldh = (timeh2 .>= Dates.DateTime(year,i,1)) .& (timeh2 .<= Dates.DateTime(year,i,e))

lonhf3[i] = lonh2[seldh];
lathf3[i] = lath2[seldh];
#timehf3[1] = timeh2[seldh];
rdvahf3[i] = rdvah2[seldh];
drvahf3[i] = drvah2[seldh];

#parametrisation 

robs_HF1[i] = (rdvahf3[i]);
#robs_HF2 = Float64.(robshc)
#robs_HF = [robs_HF1;robs_HF2]

directionobs_HF1[i] = (drvahf3[i]);
#directionobs_HF2 = atand.(uhc,vhc);
#directionobs_HF = [directionobs_HF1;directionobs_HF2]

epsilon2_HF = 1;

x_HF1[i] = (lonhf3[i]);
#x_HF2 = lonhc6;
#x_HF = [x_HF1;x_HF2]

y_HF1[i] = (lathf3[i]);
#y_HF2 = lathc6;
#y_HF = [y_HF1;y_HF2];
        
    end
end

    
# U et V radar

directionobsrad_HF1 = Vector{Vector{Float64}}(undef, 12)
uRadar = Vector{Vector{Float64}}(undef, 12)
vRadar = Vector{Vector{Float64}}(undef, 12)

for i in 1:12
    

directionobsrad_HF1[i] = directionobs_HF1[i].*(pi./180);

uRadar[i] = robs_HF1[i].*cos.(directionobsrad_HF1[i]);
vRadar[i] = robs_HF1[i].*sin.(directionobsrad_HF1[i]);
    
end       

##########################################################################################################
#                                            Drifters
##########################################################################################################


londv,latdv,timedv,udv,vdv = CoastalCurrents.loaddata("/home/jovyan/tmp/BlueCloud2026/Drifter/my.cmems-du.eu/Core/INSITU_GLO_PHY_UV_DISCRETE_MY_013_044/cmems_obs-ins_glo_phy-cur_my_drifter_PT6H/history/Dr1.nc"); # ne marche pas deja ouvert

# ID for every drifters

ds_id = NCDataset("/home/jovyan/tmp/BlueCloud2026/Drifter/my.cmems-du.eu/Core/INSITU_GLO_PHY_UV_DISCRETE_MY_013_044/cmems_obs-ins_glo_phy-cur_my_drifter_PT6H/history/Dr1.nc");
drifter_id = ds_id["DRIFTER_ID"][:];
close(ds_id)


# Only from Mediterranean Sea 

indices_a_supprimer_longv = findall(x -> !(x >= -5 && x <= 35.5), londv);

deleteat!(londv, indices_a_supprimer_longv);
deleteat!(latdv, indices_a_supprimer_longv);
deleteat!(timedv, indices_a_supprimer_longv);
deleteat!(udv, indices_a_supprimer_longv);
deleteat!(vdv, indices_a_supprimer_longv);
deleteat!(drifter_id, indices_a_supprimer_longv);

indices_a_supprimer_latv = findall(y -> !(y >= 30 && y <= 45), latdv);

deleteat!(londv, indices_a_supprimer_latv)
deleteat!(latdv, indices_a_supprimer_latv);
deleteat!(timedv, indices_a_supprimer_latv);
deleteat!(udv, indices_a_supprimer_latv);
deleteat!(vdv, indices_a_supprimer_latv);
deleteat!(drifter_id, indices_a_supprimer_latv);

# Defining the subset for the cross validation  
############# test cvr ##############################################################
dsValid = NCDataset("/home/jovyan/CoastalCurrents/examples/ValidIndices.nc")

for_cv = dsValid["for_cv"][:] 

close(dsValid)

for_cvR = round.(Int, for_cv);    
##################################################################################   
    
    
udvi = udv[for_cvR .== 0];
vdvi = vdv[for_cvR .== 0];
londvi = londv[for_cvR .== 0];
latdvi = latdv[for_cvR .== 0];
timedvi = timedv[for_cvR .== 0];

#on va enlever les NaNs des valeurs de validation

londvi_ = londvi[.!isnan.(udvi)];
latdvi_ = latdvi[.!isnan.(udvi)];
timedvi_ = timedvi[.!isnan.(udvi)];
udvi_ = udvi[.!isnan.(udvi)];
vdvi_ = vdvi[.!isnan.(udvi)];



londvi_ = londv[.!isnan.(udv)];
latdvi_ = latdv[.!isnan.(udv)];
timedvi_ = timedv[.!isnan.(udv)];
udvi_ = udv[.!isnan.(udv)];
vdvi_ = vdv[.!isnan.(udv)];

n = 12  # Longueur du vecteur

#vecteurs interpolés
londvi_2 = Vector{Vector{Float64}}(undef, n)
latdvi_2 = Vector{Vector{Float64}}(undef, n)
timedvi_2 = Vector{Vector{Float64}}(undef, n)
udvi_2 = Vector{Vector{Float64}}(undef, n)
vdvi_2 = Vector{Vector{Float64}}(undef, n)

londvi_3 = Vector{Vector{Float64}}(undef, n)
latdvi_3 = Vector{Vector{Float64}}(undef, n)
udvi_3 = Vector{Vector{Float64}}(undef, n)
vdvi_3 = Vector{Vector{Float64}}(undef, n)

robsd = Vector{Vector{Float64}}(undef, n)
directionobsd = Vector{Vector{Float64}}(undef, n)
xd = Vector{Vector{Float64}}(undef, n)
yd = Vector{Vector{Float64}}(undef, n)
robsd_ = Vector{Vector{Float64}}(undef, n)
yd_ = Vector{Vector{Float64}}(undef, n)
xd_ = Vector{Vector{Float64}}(undef, n)
directionobsd_ = Vector{Vector{Float64}}(undef, n)

#vecteur validation
#londov_ = Vector{Vector{Float64}}(undef, n)
#latdov_ = Vector{Vector{Float64}}(undef, n)
#udov_ = Vector{Vector{Float64}}(undef, n)
#vdov_ = Vector{Vector{Float64}}(undef, n)
nb_elements_a_prendre = Vector{Vector{Float64}}(undef, n)
indices_a_prendre = Vector{Vector{Float64}}(undef, n)


# Per month

for i in 1:12

# Time selection

    if i == 2
        e = 28
        

seldvi = @. Dates.DateTime(year,i,1) <= timedvi_ < Dates.DateTime(year,i,e);

londvi_2[i] = londvi_[seldvi];
latdvi_2[i] = latdvi_[seldvi];
#timedvi_2[1] = timedvi_[seldvi];
udvi_2[i] = udvi_[seldvi];
vdvi_2[i] = vdvi_[seldvi];

# Inteprolation Parameters

robsd[i] = vcat(udvi_2[i],vdvi_2[i]);
robsd[i] = Float64.(nomissing(robsd[i],NaN));
directionobsd[i] = vcat(fill(90,size(udvi_2[i])), fill(0,size(vdvi_2[i])));
xd[i] = [londvi_2[i]; londvi_2[i]];
yd[i] = [latdvi_2[i]; latdvi_2[i]];


# We retrieve NaNs

robsd_[i] = robsd[i][.!isnan.(robsd[i])];
yd_[i] = yd[i][.!isnan.(robsd[i])];
xd_[i] = xd[i][.!isnan.(robsd[i])];
directionobsd_[i] = directionobsd[i][.!isnan.(robsd[i])];length(robsd_[i])
    
        
        

    else
        e = 30
   
 

# Time selection

seldvi = @. Dates.DateTime(year,i,1) <= timedvi_ < Dates.DateTime(year,i,e);

londvi_2[i] = londvi_[seldvi];
latdvi_2[i] = latdvi_[seldvi];
#timedvi_2[1] = timedvi_[seldvi];
udvi_2[i] = udvi_[seldvi];
vdvi_2[i] = vdvi_[seldvi];

# Inteprolation Parameters

robsd[i] = vcat(udvi_2[i],vdvi_2[i]);
robsd[i] = Float64.(nomissing(robsd[i],NaN));
directionobsd[i] = vcat(fill(90,size(udvi_2[i])), fill(0,size(vdvi_2[i])));
xd[i] = [londvi_2[i]; londvi_2[i]];
yd[i] = [latdvi_2[i]; latdvi_2[i]];



# We retrieve NaNs

robsd_[i] = robsd[i][.!isnan.(robsd[i])];
yd_[i] = yd[i][.!isnan.(robsd[i])];
xd_[i] = xd[i][.!isnan.(robsd[i])];
directionobsd_[i] = directionobsd[i][.!isnan.(robsd[i])];length(robsd_[i])
        

     end
end
    
    
    return(robsa_,ya_[:],xa_[:][:],directionobsa_[:],robs_HF1[:],directionobs_HF1[:],x_HF1[:],y_HF1[:],robsd_[:],yd_[:],xd_[:],directionobsd_[:],ua4_[:],va4_[:],udvi_2[:],vdvi_2[:],uRadar[:],vRadar[:],ua,va) 

    #return(robsa_,ya_[:],xa_[:][:],directionobsa_[:],robs_HF1[:],directionobs_HF1[:],x_HF1[:],y_HF1[:],robsd_[:],yd_[:],xd_[:],directionobsd_[:],ua4_[:],ua,va4_[:],va,udvi_2[:],vdvi_2[:],uRadar[:],vRadar[:])

end # End of the function














