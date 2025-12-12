using Downloads: download



# used varname
varname = "sla"

# Grid definition
lonmax = 37
lonmin = -7
latmax = 46
latmin = 30

# resolution
dlon = dlat = 0.25               # Resolution 
lonr = lonmin:dlon:lonmax        # -7 Gibraltar to 37 which is BlackSea east
latr = latmin:dlat:latmax        # 30 on Lybia to 46 onto the north of the Black Sea

product_id_altymetry = "SEALEVEL_EUR_PHY_L3_MY_008_061"        # Altimetry dataset
product_id_DRIFT_HFR = "INSITU_GLO_PHY_UV_DISCRETE_MY_014_033" # Drifters and HFRadar dataset


# Base Directory for Dataset

basedir = expanduser("~/tmp/BlueCloud2026") 


# Directories in basedir (~/tmp/BlueCloud2026)

altimetry_dir = joinpath(basedir,"Altimetry")
drifter_dir = joinpath(basedir,"Drifter")
hf_dir = joinpath(basedir,"HF") 


# File for altimetry : Variable used in "DATA_PREOARATION" to load altimetry dataset

altimetry_fname = joinpath(altimetry_dir,"all-sla.nc")
#altimetry_fname = joinpath(altimetrydir,"all-sla-subset.nc")

mkpath(basedir)        
mkpath(altimetry_dir)
mkpath(drifter_dir)
mkpath(hf_dir)



# Mask and bathymetry in basedir

bathname = joinpath(basedir,"gebco_30sec_4.nc")
bathisglobal = true

result_filename = "surface-currents.nc"

# Download of bathymetry if not present

if !isfile(bathname)
    @info "downloading $(basename(bathname))"
    download("https://dox.ulg.ac.be/index.php/s/RSwm4HPHImdZoQP/download",bathname)
end




# function to create ~/.netrc


function load_netrc(fname = expanduser("~/.netrc"))
    cred = Dict{String,Any}()
    machine = ""
    entry = Dict()

    for line in eachline(fname)
        if startswith(line,"#") || strip(line) == ""
            continue
        end

        key,value = split(line,limit=2)
        if key == "machine"
            if machine !== ""
                cred[machine] = entry
                empty!(entry)
            end
            machine = value
        elseif key in ("login","password","account","macdef")
            entry[key] = value
        end
    end

    if machine !== ""
        cred[machine] = entry
    end

    return cred
end



# Takes year and retrive DIVAnd results from dataset_DIVA
function DIVA_results(year)

filename = "dataOBSDIVA_$year.nc"

observations = NCDataset("/home/jovyan/CoastalCurrents/examples/dataset_DIVA/$filename");

# Loading DIVAnd output

# Coordinates (reversed in the dataset)
lat_obs = observations["latitude"][:];
lon_obs = observations["longitude"][:];

# Loading velocities for June
u_obs_coarse = observations["u_velocity"][:,:,:]; 
v_obs_coarse = observations["v_velocity"][:,:,:]; 

close(observations)

# Replace missing per NaN to allow computations

u_obs_coar = coalesce.(u_obs_coarse[:,:,:],NaN);
v_obs_coar = coalesce.(v_obs_coarse[:,:,:],NaN);
    
    return lon_obs,lat_obs,u_obs_coar,v_obs_coar
    
end


# Takes the speed of a sub domain for every month & year of the analysis
function speed_domain(lon_dom,lat_dom,sizelon,sizelat)

# Domain size
lonSC = sizelon
latSC = sizelat

# Define variable
S_CuTot = zeros(Float64, lonSC, latSC, 12, 9);
S_CvTot = zeros(Float64, lonSC, latSC, 12, 9);

# For every year of dataset
for year in 2013:2021 
        
    vy = year - 2012;
    #println(vy)
    lon_obs,lat_obs,u_obs_coar,v_obs_coar = DIVA_results(year)
    mask_obs,(pm_obs,pn_obs),(xi_obs,yi_obs) = DIVAnd.domain(bathname,bathisglobal,lon_obs,lat_obs) # LES DEUX INVERSES SINON CA MARCHE PAS CACA

            # For evey month
            
            for month in 1:12


# Define whoch pixel to take
xi_obs[:,1]; # Longitudes 

v = xi_obs[:,1]
indicesS_C_xi = findall(x -> x == lon_dom, v) # From longitude 9

yi_obs[1,:]; # Latitudes

w = yi_obs[1,:]
indicesS_C_yi = findall(x -> x == lat_dom, w) # From latitude 37.5

        
# Selected velocities 
S_CuTot[:,:,month,vy] = u_obs_coar[indicesS_C_xi[1]:indicesS_C_xi[1]+(lonSC-1),indicesS_C_yi[1]:indicesS_C_yi[1]+(latSC-1),month]; # Correspond à l'année définie au dessus
S_CvTot[:,:,month,vy] = v_obs_coar[indicesS_C_xi[1]:indicesS_C_xi[1]+(lonSC-1),indicesS_C_yi[1]:indicesS_C_yi[1]+(latSC-1),month]; # Correspond à l'année définie au dessus



# theta = atan(v/u)
    #theta_r_S_C = atan(mean(S_Cv21),mean(S_Cu21))
    #thetaS_C = theta_r_S_C * 180 / π

        end
    end
    
    return(S_CuTot,S_CvTot)
end



# Compute the angle defined by the maximum variation
function current_axes(u,v)

# Calcul des déviations par rapport à la moyenne
u_prim = u #.- mean(S_Cumean)  # Anomalie
v_prim = v #.- mean(S_Cvmean)  # Anomalie

# Produit élément par élément
numerator_prod = u_prim .* v_prim
numerator = 2 * mean(numerator_prod)

# Carrés élément par élément
u_primsquare = u_prim .* u_prim
v_primsquare = v_prim .* v_prim

#denominator = abs(mean(u_primsquare) - mean(v_primsquare))
denominator = mean(u_primsquare) - mean(v_primsquare)


# Angle en radians
rad_theta = 0.5 * atan(numerator / denominator)
theta = rad2deg(rad_theta)
println("Angle en degrés : ", theta)


term1 = mean(u_primsquare) + mean(v_primsquare)
term2 = (mean(u_primsquare) - mean(v_primsquare))^2
term3 = 4*(mean(u_prim .* v_prim))^2

var_maj = 0.5 * (term1 + sqrt(term2 + term3)) # valeur propre
var_min = 0.5 * (term1 - sqrt(term2+term3)) # valeur propre

#mean_u = mean(S_Cumean)
mean_u = mean(u)
a = sqrt(var_maj)

#mean_v = mean(S_Cvmean)
mean_v = mean(v)
b = sqrt(var_min)
    
    return rad_theta,var_maj,var_min,a,b
    
end



#if haskey(ENV,"CMEMS_USERNAME") && haskey(ENV,"CMEMS_PASSWORD")
#    username = ENV["CMEMS_USERNAME"]
#    password = ENV["CMEMS_PASSWORD"]
#else
#    cred = load_netrc()

    # get CMEMS credentials

    # https://help.marine.copernicus.eu/en/articles/6135460-how-to-configure-a-simple-opendap-access-directly-in-python
    # https://web.archive.org/web/20230906115443/https://help.marine.copernicus.eu/en/articles/6135460-how-to-configure-a-simple-opendap-access-directly-in-python

#    username = cred["my.cmems-du.eu"]["login"]
#    password = cred["my.cmems-du.eu"]["password"]
#end

nothing
