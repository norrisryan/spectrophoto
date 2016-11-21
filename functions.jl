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

function getspectro(cols,rows,spectrodata)
  wl=Array{Float64}(rows)
  wlm=Array{Float64}(rows)
  fl=Array{Float64}(rows)
  fl_err=Array{Float64}(rows)
  for i=1:rows
    wl[i]=photodata[i,1]*1e4 # in Angstroms
    wlm[i]=photodata[i,1]*1e-6 #in meters
   fl[i]=(spectrodata[i,3]/(wlm[i]^2))*299792458*1.0e-26 #W m-2 m
   fl_err[i]=(spectrodata[i,4]/(wlm[i]^2))*299792458*1.0e-26 #W m-2 m
  end
  return wl,fl,fl_err
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

function chiminspectro(params)
  redloc=find(modelwv .< 33333.3)
  redmodelwv=modelwv[redloc]
  redmodelfl=redmodelfl_o.*9.999999994059551e-14.*params[1]
  redmodelflloc=redmodelfl[redloc]
  extinct=ccm89(redmodelwv,params[2])
  modelfl=redmodelflloc.*10.^(0.4.*extinct)
  for i=1:spectrorows
    value=0.0
    #interpolate and multiply together and integrate
    x1=redmodelwv
    x2=spectrowv
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


function chiminphoto(params)
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

function input(prompt::AbstractString="")
  print(prompt)
  return chomp(readline())
end

function boloflux(spectrum,waves)
    n=size((model),1)
    bflux=0.0
    range=find(waves .< 33333.3)
    wvl=waves[range]
    bolofluxmodel=spectrum[range]
    #probably need to put onto same scale
    for i=1:n
      bflux=bolofluxmodel[i]*(wvl[i+1]-wvl[i])+blfux
    return bflux
end

function parameters(flux,angdiameter,temp)
  #F=4*pi*sigma*T^4/R^2
  #angdiameter=202605*2*R/D
  #what exactly do I need to find here?

end
