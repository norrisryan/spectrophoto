#function for taking models and making into photometric points
Pkg.add("Optim")
Pkg.add("DustExtinction")
Pkg.add("Interpolations")
using Interpolations
using DustExtinction
using Optim

#FUNCTIONS
function getphoto(cols,rows,photodata)
  band=Array{String}(rows)
  wl=Array{Float64}(rows)
  wlm=Array{Float64}(rows)
  fl=Array{Float64}(rows)
  fl_err=Array{Float64}(rows)
  for i=1:rows
    band[i]=photodata[i,1]
    wl[i]=photodata[i,2]*1e4 # in Angstroms
    wlm[i]=photodata[i,2]*1e-6 #in meters
    #fl[i]=photodata[i,3] #Keep in Jy
    #fl_err[i]=photodata[i,4] #Keep in Jy
   fl[i]=(photodata[i,3]/(wlm[i]^2))*299792458*1.0e-26 #W m-2 m
   fl_err[i]=(photodata[i,4]/(wlm[i]^2))*299792458*1.0e-26 #W m-2 m
  end
  return band,wl,fl,fl_err
end

function getbands(cols,rows,banddata)
  wl=Array{Float64}(rows)
  trans=Array{Float64}(rows)
  for i=1:rows
    wl[i]=banddata[i,1] # in AngstromsS
    trans[i]=banddata[i,2] #transmission
  end
  return wl,trans
end

function closest(arr,val)
  idx=find(arr .>= val) #indices with values greater than val
  if idx[1] .> 1 #make sure that it's not one
    above=idx[1]
    if abs(val-arr[above])<abs(val-arr[above-1])
      return above,arr[above]
    else
      return (above-1),arr[above-1]
    end
  else
    return idx[1],arr[idx[1]]
  end
end

function getline(xarr,yarr,idx)
  if idx .== 1
    delx=xarr[idx+1]-xarr[idx]
    dely=yarr[idx+1]-yarr[idx]
  else
    delx=xarr[idx]-xarr[idx-1]
    dely=yarr[idx]-yarr[idx-1]
  end
  if delx.== 0
    c1=0.0
    delx=0.00001
  else
    c1=yarr[idx]-(dely/delx)*xarr[idx]
  end
  return delx,dely,c1
end

function fill(x1,x2,y1)
  j=1
  i=1
  m=1
  x1size=size(x1,1)
  x2size=size(x2,1)
  x1start=x1[1]
  x1end=x1[x1size]
  x2start=x2[1]
  x2end=x2[x2size]
  x3=zeros(x1size+x2size)
  y3=zeros(x1size+x2size)
  xlow=zeros(x2size)
  xhigh=zeros(x2size)
  x12pos=zeros(x2size)
  for k=1:(x1size+x2size)
    if i<=x1size
      if j<=x2size
        if x1[i] < x2[j] && x1[i] != x2[j]
          x3[k]=x1[i]
          x1size=size(x1,1)
          x2size=size(x2,1)
          y3[k]=y1[i]
          i=i+1
          m=m+1
        end
        if x1[i] == x2[j]
          x3[k]=x2[j]
          y3[k]=y1[i]
          xlow[j]=k
          xhigh[j]=k+1
          x12pos[j]=k
          j=j+1
          i=i+1
          m=m+1
        end
        x1size=size(x1,1)
        x2size=size(x2,1)
        if x1[i] >x2[j] && x1[i] != x2[j]
          x3[k]=x2[j]
          xlow[j]=k
          xhigh[j]=k+1
          x12pos[j]=k
          y3[k]=(y1[i]-y1[Int(i-1)])/(x1[i]-x1[Int(i-1)])*x2[j]+(y1[Int(i-1)]-(y1[i]-y1[Int(i-1)])/(x1[i]-x1[Int(i-1)])*x1[Int(i-1)])
          j=j+1
      end
    else
      x3[k]=x1[i]
      y3[k]=y1[i]
      i=i+1
      m=m+1
    end
  end
end
x3=x3[1:m]
y3=y3[1:m]
return x3,xlow,xhigh,y3,x12pos
end


function getphotpoint(x2,y2,x3,y3,x12pos)
  x2size=size(x2,1)
  value=0
  for i=2:x2size
    xlowvalue=x12pos[Int(i-1)]
    xhighvalue=x12pos[i]
    delx2=xhighvalue-xlowvalue
    for j=1:delx2
      upper=Int(xlowvalue+j)
      lower=Int(xlowvalue+j-1)
      a=(y3[upper]-y3[lower])/(x3[upper]-x3[lower])
      b=(y2[i]-y2[Int(i-1)])/(x2[i]-x2[Int(i-1)])
      c1=y3[lower]-a*x3[lower]
      c2=y2[i]-b*x2[i]
      value=value+a*b*(x3[upper]^3-x3[lower]^3)/3.0+(a*c2+b*c1)*(x3[upper]^2-x3[lower]^2)/2.0+c1*c2*(x3[upper]-x3[lower])
    end
  end
  return value
end

#function chisquare(data,model,error)
#  chi=0.0
#  par=maximum(data)/maximum(model)
#  for i=1:length(data)
#    chi=(((data[i].-model[i]*par)./error[i]).^2.0).+chi
#  end
#  return chi
#end

function deredden(flux,wavelength,rv)
  extinct=ccm89(wave,rv)
  funred=flux*10.^(0.4*extinct)
  return extinct,funred
end

#function chisquare(alphas)
#  chi=[0.0,0.0]
#  par=maximum(photofl)/maximum(modelphoto)
#  for i=1:length(photofl)
#    chi[1,1]+=(((alphas[1]+photofl[i]-modelphoto[i]*alphas[2])/photofl_err[i])^2.0)
#  end
#  return chi[1,1]
#end

function chimin(params)
  redloc=find(modelwv .< 33333.3)
  redmodelwv=modelwv[redloc]
  redmodelfl=redmodelfl_o.*9.999999994059551e-14.*params[1]
  redmodelflloc=redmodelfl[redloc]
  extinct=ccm89(redmodelwv,params[2])
  modelfl=redmodelflloc.*10.^(0.4.*extinct)
  for i=1:photorows
    bandname=photoband[i]
    currentband=readdlm("$banddirectory/$bandname.csv",',')  #open transmission file
    sband=size(currentband)
    bandrows=sband[1]
    bandcols=sband[2]
    bandwv,bandtrans=getbands(bandcols,bandrows,currentband)
    bandbegin=bandwv[1]
    bandend=bandwv[bandrows]
    value=0.0
    #interpolate and multiply together and integrate
    x1=redmodelwv
    x2=bandwv
    y1=modelfl
    y2=bandtrans
    x1size=size(x1,1)
    x2size=size(x2,1)
    x1end=x1[x1size]
    x2end=x2[x2size]
    if x1[1]>x2[1]
      #println("Lowest wavelength of filter of out Bounds! Exiting.")
      break
    end
    if x1end<x2end
      #println("Highest wavelength of filter of out Bounds! Exiting.")
      break
    end
    x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
    photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
    modelphoto[i]=photopoint
    modelphotowv[i]=photowv[i]
  end
  chi=0.0
  for i=1:length(photofl)
    chi+=(((photofl[i].-modelphoto[i]).^2.0)./photofl_err[i])
  end
  return chi
  #lowest=minimum(chisquarearr)
  #chi_min=chisquarearr[indmin(chisquarearr)]
  #modelmin=modellist[indmin(chisquarearr)]
end

function findalpha(alpha)
  redmodelfl=redmodelfl_o*9.999999994059551e-14*alpha
  #extinct=ccm89(modelwv,params[2])
  #modelfl=redmodelfl*10.^(0.4*extinct)
  modelfl=redmodelfl
  for i=1:photorows
    bandname=photoband[i]
    currentband=readdlm("$banddirectory/$bandname.csv",',')  #open transmission file
    sband=size(currentband)
    bandrows=sband[1]
    bandcols=sband[2]
    bandwv,bandtrans=getbands(bandcols,bandrows,currentband)
    bandbegin=bandwv[1]
    bandend=bandwv[bandrows]
    value=0.0
    #interpolate and multiply together and integrate
    x1=modelwv
    x2=bandwv
    y1=modelfl
    y2=bandtrans
    x1size=size(x1,1)
    x2size=size(x2,1)
    x1end=x1[x1size]
    x2end=x2[x2size]
    if x1[1]>x2[1]
      #println("Lowest wavelength of filter of out Bounds! Exiting.")
      break
    end
    if x1end<x2end
      #println("Highest wavelength of filter of out Bounds! Exiting.")
      break
    end
    x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
    photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
    modelphoto[i]=photopoint
    modelphotowv[i]=photowv[i]
  end
  chi=0.0
  for i=1:length(photofl)
    chi+=(((photofl[i].-modelphoto[i]).^2.0)/photofl_err[i])
  end
  return chi
  #lowest=minimum(chisquarearr)
  #chi_min=chisquarearr[indmin(chisquarearr)]
  #modelmin=modellist[indmin(chisquarearr)]
end

#-------------------------------START PROGRAM-----------------------------------------
basedirectory="/home/norris/SYNTHSPEC/"
banddirectory="$basedirectory/FILTERS"
photofile="$basedirectory/photometry_azcygbandsmod.csv"
modellistfile="$basedirectory/MARCS/modellist.txt"
modelwavelength="$basedirectory/flx_wavelengths.vac"
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
  #alphas=zeros(Float64,100)
  #alphaarr=zeros(Float64,100)
  #betas=zeros(Float64,100)
  #numberparm=length(alphas)
#  for i=1:length(alphas)
#      alpha=alphas[i]+i*0.05
#      alphaarr[i]=alpha
#      chisquarearr[i]=findalpha(alpha)
#  end
  chisquared=optimize(chimin,[0.0,0.0])
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
end































  #interpolate
  #for i=2:(size(allwv,1))
  #  value=0.0
  #  idxbegin1,beginx1=closest(modelwv,allwv[(i-1.0),1])
  #  idxend1,endx1=closest(modelwv,allwv[i,1])
  #  idxbegin2,beginx2=closest(bandwv,allwv[(i-1.0),1])
  #  idxend2,endx2=closest(bandwv,allwv[i,1])
  #  dely1=y1[idxend1]-y1[idxbegin1]
  #  dely2=y2[idxend2]-y2[idxbegin2]
  #  delw1=endx1-beginx1
  #  delw2=endx2-beginx2
  #  alpha=dely1/delw1
  #  beta=dely2/delw2
  #  c1=y1[idxbegin1]-alpha*delw1
  #  c2=y2[idxbegin2]-beta*delw2
  #  value=value+(1.0/3.0)*(endx1^3.0-beginx1^3.0)*alpha*beta+(1.0/2.0)*(endx1^2.0-beginx1^2.0)*(alpha*c2+beta*c1)+delw1*c1*c2
  #end
#if bandrows .< bandflxsize #if more points in the flux file
#  if bandrows .< bandflxsize #if more points in the flux file
#    for i=1:(bandrows-1.0)
      #if i .== bandrows
      #  rangeidxbegin,rangevalbegin=closest(modelwv,bandwv[i-1.0,1])
      #  rangeidxend,rangevalend=closest(modelwv,bandwv[i,1])
      #else
#        rangeidxbegin,rangevalbegin=closest(modelwv,bandwv[i,1])
#        rangeidxend,rangevalend=closest(modelwv,bandwv[(i+1.0),1])
#        rangeflxsize=rangeidxend-rangeidxbegin                     #number of model points in this range
#      for j=1:(rangeflxsize-1)
#        delmodelfl=modelfl[rangeidxbegin+j+1]-modelfl[rangeidxbegin+(j)]
#        delbandtrans=bandtrans[i+1]-bandtrans[i]
#        delmodelwv=modelwv[rangeidxbegin+j+1]-modelwv[rangeidxbegin+(j)]
#        delbandwv=bandwv[i+1]-bandwv[i]
#        c1=modelfl[idxbegin+j+1]-(delmodelfl/delmodelwv)*modelwv[idxbegin+j+1]
#        c2=bandtrans[i]-(delbandtrans/delbandwv)*bandwv[i+1]
#        interpmodelflx[m]=(delmodfl/delmodelwv)*allwv[idxbegin+i+j]+c1
        #possible problem site
#        interpmodelwv[m]=allwv[idxbegin+i+j]
#        interpbandflx[m]=(delbandtrans/delbandwv)*allwv[idxbegin+i+j]+c2
#        interpbandwv[m]=allwv[idxbegin+i+j]
#        m=m+1.0
#      end
#    end
#  end
#  if bandrows .>= bandflxsize #if more points in the band file
#    for i=1:bandflxsize
#      if i .== bandflxsize
#        rangeidxbegin,rangevalbegin=closest(bandwv,modelwv[idxbegin+(i-1),1])
#        rangeidxend,rangevalend=closest(bandwv,modelwv[(idxbegin+i),1])
#      else
#        rangeidxbegin,rangevalbegin=closest(modelwv,bandwv[i,1])
#        rangeidxend,rangevalend=closest(modelwv,bandwv[(i+1),1])
#      end
#      rangeflxsize=rangeidxend-rangeidxbegin                     #number of model points in this range
#      for j=1:rangeflxsize
#        delmodfl=modelfl[idxbegin+i]-modelfl[idxbegin+(i-1)]
#        deljbandtrans=bandtrans[rangeidxbegin+j]-bandtrans[rangeidxbegin+j-1]
#        delmodelwv=modelwv[idxbegin+i]-modelwv[idxbegin+(i-1)]
#        delbandwv=bandwv[rangeidxbegin+j]-bandwv[rangeidxbegin+j-1]
#        c1=modelfl[idxbegin+i]-(delmodelfl/delmodelwv)*modelwv[idxbegin+i]
#        c2=bandtrans[rangeidxbegin+i]-(delbandtrans/delbandwv)*bandwv[rangeidxbegin+i]
#        interpmodelflx[i]=(delmodfl/delmodelwv)*allwv[idxbegin+i+j]+c1
        #possible problem site
#        interpmodelwv[i]=allwv[idxbegin+i+j]
#        interpbandflx[i]=(delbandtrans[i]/delbandwv[i])*allwv[idxbegin+i+j]+c2
#        interpbandwv[i]=allwv[idxbegin+i+j]
#      end
#    end
#  end
#end
#value=0.0
#y1=interpmodelflx
#y2=interpbandflx
#x1=interpbandwv
#for i=2:size(y2,1)
#           dely1=y1[i]-y1[i-1]
#           dely2=y2[i]-y2[i-1]
#           delwv=x1[i]-x1[i-1]
#           alpha=dely1/delwv
#           beta=dely2/delwv
#           c1=y1[i-1]-alpha*x1[i-1]
#           c2=y2[i-1]-beta*x2[i-1]
#           value=value+(1.0/3.0)*(x1[i]^3.0-x1[i-1]^3.0)*alpha*beta+(1.0/2.0)*(x1[i]^2.0-x1[i-1]^2.0)*(alpha*c2+beta*c1)+delwv*c1*c2
#end



#end
