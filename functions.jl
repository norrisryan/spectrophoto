

function fill(x1,x2,y1)
  j=1
  i=1
  m=1
  tol=1e-20
  x1size=size(x1,1)
  x2size=size(x2,1)
  x1start=x1[1]
  x1end=x1[x1size]
  x2start=x2[1]
  x2end=x2[x2size]
  if x1[1]>x2[1]
    error("Lowest wavelength of filter of out Bounds! Exiting.")
  end
  if x1end<x2end
    error("Highest wavelength of filter of out Bounds! Exiting.")
  end
  x3=zeros(x1size+x2size)
  y3=zeros(x1size+x2size)
  xlow=zeros(x2size)
  xhigh=zeros(x2size)
  x12pos=zeros(x2size)
  for k=1:(x1size+x2size)
    if i<=x1size
      if j<=x2size
        if x1[i] < x2[j] && abs((x1[i]-x2[j])) > tol
          x3[k]=x1[i]
          y3[k]=y1[i]
          i=i+1
          m=m+1
        end
        if abs((x1[i]-x2[j])) < tol
          x3[k]=x2[j]
          y3[k]=y1[i]
          xlow[j]=k
          xhigh[j]=k+1
          x12pos[j]=k
          j=j+1
          i=i+1
          m=m+1
        end
        if x1[i] >x2[j] && abs((x1[i]-x2[j])) > tol
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

x1=linspace(-pi,pi,1000)
y1=sin(x1)
x2=linspace(0,3/4*pi,10)
y2=cos(x2)
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
println("Equal Spacing: Integral of Cos(x)*Sin(x) from 0 to 3/4*pi: ",photopoint)

x1=[-pi,-1.5,0,0.2,0.3,0.45,0.47,0.48,.5,.75,1.0,2.0,pi]
y1=sin(x1)
x2=[0,0.1,0.5,1.0,1.5,3/4*pi]
y2=cos(x2)
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
println("Unequal Spacing: Integral of Cos(x)*Sin(x) from 0 to 3/4*pi: ",photopoint)

x1=linspace(0,100,400)
y1=x1*0+3
x2=linspace(20,40,10)
y2=x2*0+2
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
println("Equal Spacing: Integral of 3*2 from 20 to 40: ",photopoint)

x1=[0,1,2,3,4,5,10,20,25,35,37,38,39,41,42,50,60]
y1=x1*0+3
x2=[20,25,27,28,33]
y2=x2*0+2
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
println("Unequal Spacing: Integral of 3*2 from 20 to 40: ",photopoint)

println("Testing Error: Filter bound too low")
x1=[0,1,2,3,4,5,10,20,25,35,37,38,39,41,42,50,60]
y1=x1*0+3
x2=[-10,20,25,27,28,33]
y2=x2*0+2
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)

println("Testing Error: Filter bound too high")
x1=[0,1,2,3,4,5,10,20,25,35,37,38,39,41,42,50,60]
y1=x1*0+3
x2=[10,20,25,27,28,33,75]
y2=x2*0+2
x3,xlow,xhigh,y3,x12pos=fill(x1,x2,y1)
photopoint=getphotpoint(x2,y2,x3,y3,x12pos)
end
