from numpy import *
from vtk import *
from vtk.util import numpy_support as VN
import time as myTime
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

testWe=[10,30,50,'water'] #Different We numbers

kinematic=True  #True if we want to test kinematic phase, False if we want all the phases.

### For kinematic phase study, see the input.command file for an example ###

if kinematic==True:
    testD=[40,50,60]    #Different diameters, then different Re numbers
else:
    testD=[40]

if kinematic==True:
    niter=3000
    nfreq=10
    nstart=0
else:
    niter=30000
    nfreq=100
    nstart=0
        
fig=figure(1)
ax = fig.add_subplot(111)
ax.set_xlabel("Dimensionless time t*")
ax.set_ylabel("D/D0")
        
#VTK files
nu=0.1/3
drho=0.431-0.00149  #parameters for Tr=0.55
gamma=8.871e-3
rhol=0.431 
D=[]
t=[]

def func(x,a):      #fiiting curve for kinematic phase
    return a*x**0.5

startTime = myTime.time()

for we in testWe:
    for do in testD:
        if we=='water' and do==40:
            print(we)
            nstart=9000
            niter=11000
            do=52
            u=(gamma*20/(drho*do))**0.5
        elif we!='water':
            print('we'+str(we)+'d'+str(do))
            u=(gamma*we/(drho*do))**0.5
        Re=do*u/nu
        reader = vtkStructuredPointsReader()
        ycenter=[]
        xcenter=[]
        time=[]
        #diameter=[] 
        dratio=[]
        if (we=='water' and do==52) or we!='water':
            for i in range(nstart,niter+1,nfreq):    #loop to read all the VTK files
                if kinematic==True:
                    if we=='water':
                        reader.SetFileName('kinematic/water/VTKResults/t'+str(i)+".vtk")    #Need to create a simulation file with the corresponding VTK files
                    else:
                        reader.SetFileName('kinematic/we'+str(we)+'d'+str(do)+"/VTKResults/t"+str(i)+".vtk")
                else:
                    reader.SetFileName('time evolution spreading/we'+str(we)+'d'+str(do)+"/VTKResults/t"+str(i)+".vtk")
                #print('we'+str(we)+'d'+str(do)+"/VTKResults/t"+str(i)+".vtk")
                reader.ReadAllScalarsOn()
                reader.Update()
                data = reader.GetOutput()
                n = data.GetNumberOfPoints()    #number of points
                dim = data.GetDimensions()      #dimensions of the domain
                nx=dim[0]
    #Find the center position of the droplet
                """
                yres=xres=c=0
                for p in range(n):
                    pos=data.GetPoint(p)    #to get the position of the point
                    x=pos[0]
                    y=pos[1]
                    d=data.GetPointData()
                    array=d.GetArray('rho')
                    density=array.GetValue(p)   #to get the value of the point density
                    if density>=rhol*0.8:      #to find the center of the droplet
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
                    if rhok>=rhol*0.9:      
                        dmax=posk[0]
                    k+=1
                while dmin==0:
                    posz=data.GetPoint(z)
                    dz=data.GetPointData()
                    arrayz=dz.GetArray('rho')
                    rhoz=arrayz.GetValue(z)   
                    if rhoz>=rhol*0.9:     
                        dmin=posz[0]
                    z-=1
                if dmin-dmax>0:
                    diameter.append((dmin-dmax))
                    dratio.append((dmin-dmax)/do)
                    time.append(i*u/do)
                """
                dmax=dmin=0
                k=z=int(5*nx+nx/2)
                pos=data.GetPoint(k)
                d=data.GetPointData()
                array=d.GetArray('rho')
                rho=array.GetValue(k)
                if rho>0.9*rhol:
                    while dmax==0:
                        posk=data.GetPoint(k)
                        dk=data.GetPointData()
                        arrayk=dk.GetArray('rho')
                        rhok=arrayk.GetValue(k)   
                        if rhok<=rhol*0.9:      
                            dmax=posk[0]
                        k-=1
                    while dmin==0:
                        posz=data.GetPoint(z)
                        dz=data.GetPointData()
                        arrayz=dz.GetArray('rho')
                        rhoz=arrayz.GetValue(z)   
                        if rhoz<=rhol*0.9:     
                            dmin=posz[0]
                        z+=1
                
                    #if dmin-dmax>0:
                    #diameter.append((dmin-dmax))
                    dratio.append((dmin-dmax)/do)
                    time.append(i*u/do)
    
            if kinematic==True:
                k=0
                t0=0
                dlim=[]
                tlim=[]
                while dratio[k+1]<=dratio[k]:
                    k+=1
                dlim.append(dratio[k])
                t0=time[k]
                tlim.append(time[k]-t0)
                k+=1
                while time[k]<t0+0.1:
                    dlim.append(dratio[k])
                    tlim.append(time[k]-t0)
                    k+=1
                D.extend(dlim) #D and t for the fitting curve D=a*t**0.5
                t.extend(tlim)   
                del dlim[0]
                del tlim[0]
                            
    
            #tight_layout()
            if kinematic==True:
                if we==10:
                    ax.plot(tlim,dlim,'^',label='We='+str(we)+' ; Re='+str(round(Re,1)))
                if we==30:
                    ax.plot(tlim,dlim,'o',label='We='+str(we)+' ; Re='+str(round(Re,1)))
                if we==50:
                    ax.plot(tlim,dlim,'s',label='We='+str(we)+' ; Re='+str(round(Re,1)))
                if we=='water' and do==52:
                    ax.plot(tlim,dlim,'v',label='We= 20 ; Re='+str(round(Re,1)))
            elif kinematic==False:
                if we==10:
                    ax.plot(time,dratio,'^',label='We='+str(we)+' ; Re='+str(round(Re,1)))
                if we==30:
                    ax.plot(time,dratio,'o',label='We='+str(we)+' ; Re='+str(round(Re,1)))
                if we==50:
                    ax.plot(time,dratio,'s',label='We='+str(we)+' ; Re='+str(round(Re,1)))

if kinematic==True:
    res,error=curve_fit(func,t,D)
    t=linspace(0,max(t),1000)
    D=[]
    for i in t:
        D.append(res[0]*i**0.5)
    ax.plot(t,D,'-k',label="fit curve: a="+str(round(res[0],2)))
    lgd=ax.legend(loc='upper left', bbox_to_anchor=(0, -0.14,1,0), ncol=2, mode="expand")
#    savefig("Spread_factor_kinematic.eps", bbox_extra_artists=(lgd,), bbox_inches='tight')
#    savefig("Spread_factor_kinematic.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
#    savefig("Spread_factor_kinematic.eps")
#    savefig("Spread_factor_kinematic.png")
if kinematic==False:
    legend()
    semilogx()
    xlim(1,20)
    savefig("Spread_factor.eps")
    savefig("Spread_factor.png")  
show()
print('Total execution time = {:5.2f}'.format((myTime.time()-startTime)/60.),"min")
