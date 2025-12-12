# Input : Parameters and considered year

function cosfunction(x,year)
include("common.jl")

    

    
    
#################################################################################################################################################
#                                                            Parametrization
#################################################################################################################################################    
    
    



n = 12   

lenAll = x[1] 
epsilon2A = x[2]
epsilon2D = x[3]
epsilon2H = x[4]
 
# Definition Epsilon = Error variance of the observations
    for i in 1:n
      
epsilon2_alt[i] = fill(epsilon2A, size(xa_[i]));
epsilon2_drift[i] = fill(epsilon2D, size(xd_[i]));  ### bien vers 0.005
epsilon2_HF[i] = fill(epsilon2H, size(x_HF1[i]));

# Final vector of error variance
epsilon2[i] = vcat(epsilon2_alt[i],epsilon2_drift[i],epsilon2_HF[i]); #concaténation des deux

    end


# Extra parameters 
    
len = fill(lenAll, size(h))  #longueur de correlation (valeur initiale de 50e3) valeur unique mais peut varier spatialement (scalaire ou matrice de la meme dimension que mask)
eps2_boundary_constraint = -1
eps2_div_constraint = -1
eps2_div_constraint = 1e+1 #contrainte sur la divergence (nabla u proche de 0)

    


    for i in 1:n
        
# Vector field not interpolated :
##################################################################################################
x_tot[i] = vcat(xa_[i],xd_[i],x_HF1[i]);                                                         #  
y_tot[i] = vcat(ya_[i], yd_[i],y_HF1[i]);                                                        #  
robs_tot[i] = vcat(robsa_[i],robsd_[i],robs_HF1[i]);                                             #  
directionobs_tot[i] = vcat(directionobsa_[i],directionobsd_[i],directionobs_HF1[i]);             #
##################################################################################################
    
    end
    



#################################################################################################################################################
#                                                      Definition of variables
#################################################################################################################################################   
    


    
    
# Size of the matrix for the interpolation defined from common.jl variables    
m = abs(((lonmin - lonmax)*4) - 1)
p = abs(((latmin - latmax)*4) - 1)
#m, p = 177, 65  # Dimensions de chaque matrice

    # Defining new variables

londov1_ = Vector{Vector{Float64}}(undef, n);
latdov1_ = Vector{Vector{Float64}}(undef, n);
udov1_ = Vector{Vector{Float64}}(undef, n);
vdov1_ = Vector{Vector{Float64}}(undef, n);
xo1_ = Vector{Vector{Float64}}(undef, n);
yo1_ = Vector{Vector{Float64}}(undef, n);
speedo = Vector{Vector{Float64}}(undef, n);
xval = Vector{Vector{Float64}}(undef, n);
yval = Vector{Vector{Float64}}(undef, n);
reinterpole_uri = Vector{Vector{Float64}}(undef, n);
reinterpole_vri = Vector{Vector{Float64}}(undef, n);
reinterpole_uri_ = Vector{Vector{Float64}}(undef, n);
reinterpole_vri_ = Vector{Vector{Float64}}(undef, n);
speedi_ = Vector{Vector{Float64}}(undef, n);
#itp_vri = Vector{Matrix{Float64}}(undef, n);
#itp_uri = Vector{Matrix{Float64}}(undef, n);
indicesNaNu = Vector{Vector{Bool}}(undef, n);
indicesNaNv = Vector{Vector{Bool}}(undef, n);
vecteurUisansNaN = Vector{Vector{Float64}}(undef, n);
vecteurVisansNaN = Vector{Vector{Float64}}(undef, n);
vecteurUOBSsansNaN = Vector{Vector{Float64}}(undef, n);    
vecteurVOBSsansNaN = Vector{Vector{Float64}}(undef, n);
RMSu = Vector{Float64}(undef, n);   
RMSv = Vector{Float64}(undef, n);
RMStot = Vector{Float64}(undef, n); 
sommeRMS = Vector{Float64}(undef, n); 



itp_uri = [rand(m, p) for _ in 1:n]   
itp_vri = [rand(m, p) for _ in 1:n]   

    
    
#################################################################################################################################################
#                                                      DIVAnd run : HFRAdar
################################################################################################################################################# 

    
    # Month per month : DIVAnd Interpolation
    
    for i in 1 : 12

uri[:,:,i],vri[:,:,i],ηi[:,:,i] = DIVAndrun_HFRadar(
    mask,h,(pm,pn),(xi,yi),(x_tot[i],y_tot[i]),robs_tot[i],directionobs_tot[i],len,epsilon2[i];
    eps2_boundary_constraint = eps2_boundary_constraint,
    eps2_div_constraint = eps2_div_constraint,
    # eps2_Coriolis_constraint = -1,
    # f = 0.001,
    # residual = residual,
    # g = g,
    # ratio = 100,
    # lenη = (000.0, 000.0, 24 * 60 * 60. * 10),
    # maxit = 100000,
    # tol = 1e-6,
);
        if i == 1
    println("$i st month computed")
        else 
    println("$i th month computed")
        end
    end


#################################################################################################################################################
#                                                             Observations 
################################################################################################################################################# 
    
    
    # Month per month : Observations

    for i in 1:12
     
        if i == 2
            e = 28

seldv1 = @. Dates.DateTime(year,i,1) <= timedov_ < Dates.DateTime(year,i,e);

# Coordinates of the observations
londov1_[i] = londov_[seldv1];
latdov1_[i] = latdov_[seldv1];

# Velocities of the observations                        
udov1_[i] = udov_[seldv1];
vdov1_[i] = vdov_[seldv1];

# Speed vector of the observations
speedo[i] = sqrt.(udov1_[i].^2 + vdov1_[i].^2)

# Renaming the observations coordinates
xo1_[i] = londov1_[i];
yo1_[i] = latdov1_[i];
    
            
        else
            e = 30
            
seldv1 = @. Dates.DateTime(year,i,1) <= timedov_ < Dates.DateTime(year,i,e);

            
# Coordinates of the observations
londov1_[i] = londov_[seldv1];
latdov1_[i] = latdov_[seldv1];

            
# Velocities of the observations            
udov1_[i] = udov_[seldv1];
vdov1_[i] = vdov_[seldv1];

# Speed vector of the observations
speedo[i] = sqrt.(udov1_[i].^2 + vdov1_[i].^2)

# Renaming the observations coordinates
xo1_[i] = londov1_[i];
yo1_[i] = latdov1_[i];
            
        end
    end



    
#################################################################################################################################################
#                                            Validation (comparison between observations and DIVA run
################################################################################################################################################# 
 
    
    
    


    
#                                               Interpolation of interpolated data on observations    
#################################################################################################################################################
    for i in 1:12

        
# Grid creation
itp_vri = interpolate((xi[:,1], yi[1,:]), vri[:,:,i], Gridded(Linear()));
itp_uri = interpolate((xi[:,1], yi[1,:]), uri[:,:,i], Gridded(Linear()));

# Observations coordinates
xval[i] = xo1_[i];
yval[i] = yo1_[i];



# Inteprolation of results on the observations
reinterpole_uri[i] = itp_uri.(xval[i],yval[i]);
reinterpole_vri[i] = itp_vri.(xval[i],yval[i]);
    


# Speed computation
speedi_[i] = sqrt.(reinterpole_uri[i].^2 + reinterpole_vri[i].^2)
    
  end

    
# Written output to assess that the reinterpolation is done
    
println("")    
println("reinterpolation done")
println("")     
    
println("parameters are : $x^[1] ,$x^[2] ,$x^[3] ,$x^[4]")
println("")     
    
    
#                                              Correlation between obs and interpolated data    
#################################################################################################################################################

for i in 1:12

println("$i th month of $year for the validation : ")
println("")
        
        
    if length(reinterpole_uri[i]) == 0 # Check of the availability of observations
            println("no data for this month")
            println("")
            
    else
   
# RMS computation 

    #U
    
# NaNs removal 
indicesNaNu[i] = isnan.(reinterpole_uri[i]);
vecteurUisansNaN[i] = reinterpole_uri[i][.!indicesNaNu[i]]
vecteurUOBSsansNaN[i] = udov1_[i][.!indicesNaNu[i]]

#on calcule la corrélation
println("correlation for u on the $i th month is : ",cor(vecteurUisansNaN[i],vecteurUOBSsansNaN[i]))



    #V
    
# NaNs removal 
indicesNaNv[i] = isnan.(reinterpole_vri[i]);
vecteurVisansNaN[i] = reinterpole_vri[i][.!indicesNaNv[i]]
vecteurVOBSsansNaN[i] = vdov1_[i][.!indicesNaNv[i]]

#on calcule la corrélation
println("correlation for v on the $i th month is : ",cor(vecteurVisansNaN[i],vecteurVOBSsansNaN[i]))
println("")

        
        end          
        
    end
    
    
#                                             Root Mean Squared (RMS) error computaiton    
#################################################################################################################################################
    
    for i in 1:12
        
        
println("RMS computation for the $i th month")
        
        if length(reinterpole_uri[i]) == 0 #on vérifie que le vecteur n'est pas vide
            println("no data for the $i th month")
            println("")

    
    else


# Final computation for RMS
RMSu[i] = sqrt(sum((vecteurUisansNaN[i] .- vecteurUOBSsansNaN[i]).^2)/length(vecteurUOBSsansNaN[i]));
RMSv[i] = sqrt(sum((vecteurVisansNaN[i] .- vecteurVOBSsansNaN[i]).^2)/length(vecteurVOBSsansNaN[i]));

            
println("u RMS of the $i th month is :")

println(RMSu[i])


println("v RMS of the $i th month is :")
println(RMSv[i])
println("")
        
    
println("total RMS of the $i th month is :")

            
RMStot[i] = RMSu[i] + RMSv[i];
println(RMStot[i]);            
println("")


            
            
            
            

            
        end
    end

#                                               Weighted Root Mean Squared Errpr (RMS) computation    
#################################################################################################################################################   
            
weightedRMS = (length(xo1_[1])*RMStot[1] + length(xo1_[2])*RMStot[2] + length(xo1_[3])*RMStot[3] + length(xo1_[4])*RMStot[4] + length(xo1_[5])*RMStot[5] + length(xo1_[6])*RMStot[6] + length(xo1_[7])*RMStot[7] + length(xo1_[8])*RMStot[8] + length(xo1_[9])*RMStot[9] + length(xo1_[10])*RMStot[10] + length(xo1_[11])*RMStot[11] + length(xo1_[12])*RMStot[12])/(length(xo1_[1])+length(xo1_[2])+length(xo1_[3])+length(xo1_[4])+length(xo1_[5])+length(xo1_[6])+length(xo1_[7])+length(xo1_[8])+length(xo1_[9])+length(xo1_[10])+length(xo1_[11])+length(xo1_[12])); 
    
    
    
# RMS sum on all the months computation
    
sommeRMS[1] = 0;  
    
    for i in 1:12
sommeRMS[1] += RMStot[i];  # Ajouter la valeur à la somme
       
    end
    
    
# Writing the results in a textfile
    
#f = open("fichiersortieRMS_TEST.txt", "a"); #mode écriture

    
    #println(f,"pour les paramètres $x[1:4]");
    #println(f,"La RMS par mois est de :");
    #println(f,RMStot);
            
   # println(f, "La corrélation entre le vecteur uo et ui est de :")
   # println(f, cor.(vecteurUisansNaN[3],vecteurUOBSsansNaN[3]))
   # println(f, "La corrélation entre le vecteur vo et vi est de :")
   # println(f, cor.(vecteurVisansNaN,vecteurVOBSsansNaN))
    
            
    #println(f,"La RMS additionée pour tous les mois est de  :");
    #println(f,sommeRMS[1]);
    #println(f,"");
  
    #println(f, "La weighted RMS est : ");
    #println(f, weightedRMS);
    #println(f,"");
       
#close(f)
  
    
# Output   

return (vecteurUisansNaN,
        vecteurVisansNaN,vecteurUOBSsansNaN,vecteurVOBSsansNaN,
        udov1_,vdov1_,reinterpole_uri,reinterpole_vri,xo1_,yo1_,speedo,speedi_,xval,yval,weightedRMS); # return of the function
    
    

end # End of the function









