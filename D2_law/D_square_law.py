from numpy import *
from vtk import *
from vtk.util import numpy_support as VN
import time as myTime
from matplotlib.pyplot import *
from scipy.optimize import curve_fit


fig=figure(1)
ax = fig.add_subplot(111)
ax.set_xlabel("Dimensionless time t*")
ax.set_ylabel("(D/D0)²")

        
#VTK files

D=[]
t=[]

def func(x,a,b):      #fiiting curve for kinematic phase
    return a*x+b

startTime = myTime.time()
limit=0.018
nstart=0
niter=3000
nfreq=100
nu=0.1/3
drho=0.431-0.00149  #parameters for Tr=0.55
gamma=8.871e-3
rhol=0.431 
do=60
u=(gamma*20/(drho*do))**0.5
Re=do*u/nu
reader = vtkStructuredPointsReader()
ycenter=[]
xcenter=[]
time=[]
#diameter=[] 
dsquare=[]
for i in range(nstart,niter+1,nfreq):    #loop to read all the VTK files
    reader.SetFileName('static0.6/VTKResults/t'+str(i)+".vtk") #Need to perform a simulation and to call it "static0.6", and use the VTK files. See input.command
    print('static0.6/VTKResults/t'+str(i)+".vtk")
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    n = data.GetNumberOfPoints()    #number of points
    dim = data.GetDimensions()      #dimensions of the domain
    nx=dim[0]
#Find the center position of the droplet
    
    yres=xres=c=0
    for p in range(n):
        pos=data.GetPoint(p)    #to get the position of the point
        x=pos[0]
        y=pos[1]
        d=data.GetPointData()
        array=d.GetArray('rho')
        density=array.GetValue(p)   #to get the value of the point density
        if density>=0.2:      #to find the center of the droplet
            c+=1
            yres+=y
            xres+=x
    yc=yres/c
    xc=xres/c
    xcenter.append(xc)
    ycenter.append(yc)
    
#Find the diameter of the droplet
    dmax=dmin=0
    k=int(yc*nx)
    z=k+nx
    while dmax==0:
        posk=data.GetPoint(k)
        dk=data.GetPointData()
        arrayk=dk.GetArray('rho')
        rhok=arrayk.GetValue(k)   
        if rhok>=0.2:      
            dmax=posk[0]
        k+=1
    while dmin==0:
        posz=data.GetPoint(z)
        dz=data.GetPointData()
        arrayz=dz.GetArray('rho')
        rhoz=arrayz.GetValue(z)   
        if rhoz>=0.2:     
            dmin=posz[0]
        z-=1
    #diameter.append((dmin-dmax))
    """
    if i*nu/do**2>limit:
        dsquare.append(((dmin-dmax)/do)**2)
        time.append(i*nu/do**2-limit)
    """
    dsquare.append(((dmin-dmax)/do)**2)
    time.append(i)
ax.plot(time,dsquare,'ro',label='lambda=2/3 ; cv=5')  
  
res,error=curve_fit(func,time,dsquare)
print('a=',res[0],'; b=',res[1])
t=linspace(min(time),max(time),1000)
D=[]
for i in t:
    D.append(res[0]*i+res[1])
    
ax.plot(t,D,'-r',label="fit curve")
legend()
savefig("D_square_law_0.6.eps")
savefig("D_square_law_0.6.png")
show()
print('Total execution time = {:5.2f}'.format((myTime.time()-startTime)/60.),"min")
