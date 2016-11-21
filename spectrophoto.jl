#function for taking models and making into photometric points
Pkg.add("Optim")
Pkg.add("DustExtinction")
Pkg.add("Interpolations")
#Pkg.add("Gtk")
using Interpolations
using DustExtinction
using Optim
#using Gtk
include("functions.jl")
#-------------------------------START PROGRAM-----------------------------------------
#basedirectory="/home/norris/SYNTHSPEC/" #DIRECTORY OF FILES AND FOLDERS
basedirectory=input("What is the base directory for your models and data?") #DIRECTORY OF FILES AND FOLDERS
banddirectory="$basedirectory/FILTERS"  #WHERE THE TRANSMISSION CURVES ARE LOCATED
photofile="$basedirectory/photometry_azcygbandsmod.csv" #PHOTOMETRY FILE IN CSV FORMAT
modellistfile="$basedirectory/MARCS/modellist.txt"      #LIST OF MODEL FILES
modelwavelength="$basedirectory/flx_wavelengths.vac"    #WAVELENGTH FILE FOR MARCS
#read files in
photo=readdlm(photofile,',')
photosize=size(photo)
photorows=photosize[1]
photocols=photosize[2]
modellist=readdlm(modellistfile)
modeln=size(modellist,1)
modelwvvac=readdlm(modelwavelength,',')  #read in model wavelength
modelwv = modelwvvac ./ (1.0 + 2.735182E-4 + 131.4182 ./ modelwvvac.^2 + 2.76249E8./ modelwvvac.^4)
modelwvm=modelwv*1e-10
#photometry gathering
photoband,photowv,photofl,photofl_err=getphoto(photocols,photorows,photo)
photofl_err[photofl_err.==0.0]=1e-7
modelphoto=zeros(Float64,photorows)
modelphotowv=zeros(Float64,photorows)
chisquarearr=zeros(Float64,modeln)
rvs=zeros(Float64,modeln)
weights=zeros(Float64,modeln)
#get model photometric points based on transmission and bands used in obsv
for j=1:modeln
  model=modellist[j]
  println(model," ")
  modelflux="$basedirectory/MARCS/$model"
  redmodelfl_o=readdlm(modelflux,',')       #read in model flux
  chisquared=optimize(chiminphoto,[0.0,0.0])
  chisquareval=Optim.minimum(chisquared)
  chisquarearr[j]=chisquareval
  values=Optim.minimizer(chisquared)
  weights[j]=values[1]
  rvs[j]=values[2]
  println("Chi^2: ",chisquarearr[j]," Weight: ",weights[j]," RV: ",rvs[j])
end
lowchi=minimum(chisquarearr)
loweight=weights[indmin(chisquarearr)]
lowrv=rvs[indmin(chisquarearr)]
# next integrate the sprtum to get flux
bestmodel=modellist[indmin(chisquarearr)]
bolflux=boloflux(spectrum,waves)

end
