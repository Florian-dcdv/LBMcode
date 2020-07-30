import numpy
from vtk import *
from vtk.util import numpy_support as VN
import time as myTime
from matplotlib.pyplot import *
import csv

test='water'

#Bertola graph, height of the center of the droplet in respect of time.

wat=csv.reader(open('Bertola_water.csv','r'),delimiter=',')
peo=csv.reader(open('Bertola_PEO.csv','r'),delimiter=',')
xg=csv.reader(open('Bertola_XG.csv','r'),delimiter=',')

if test=='water':
    t_w=[]
    h_w=[]
    k=0
    for i in wat:
        k+=1
        if k>6:
            t_w.append(float(i[0]))
            h_w.append(float(i[1]))
    D=48            #Droplet diameter in lattice units
    Dp=3.09e-3
    g=1.2661e-5     #gravity in lattice units, found with 'Dimensionless numbers.py'
    gp=9.81
elif test=='peo':        
    t_peo=[]
    h_peo=[]
    k=0
    for i in peo:
        k+=1
        if k>6:
            t_peo.append(float(i[0]))
            h_peo.append(float(i[1]))
    D=44
    Dp=2.93e-3
    g=1.3935e-5
    gp=9.81
elif test=='xg':
    t_xg=[]
    h_xg=[]
    k=0
    for i in xg:
        k+=1
        if k>6:
            t_xg.append(float(i[0]))
            h_xg.append(float(i[1]))
    D=1
    Dp=3.12e-3
    g=1
    gp=9.81

#VTK files
        
startTime = myTime.time()

nfreq=500
niter=25000

reader = vtkStructuredPointsReader()
xcenter=[]
ycenter=[]
time=[]
height=[]

dx=Dp/D                 #physical lattice constant in m
dt=(g*dx/gp)**0.5*1000  #physical time step in ms

rhol=0.455
xmin=10000

for i in range(0,niter+1,nfreq):    #loop to read all the VTK files
    reader.SetFileName(test+"/VTKResults/t"+str(i)+".vtk")  #Need to perform a simulation of a bouncing droplet with the parameters calculated from 'Dimensionless numbers.py'.
    print(test+"/VTKResults/t"+str(i)+".vtk")
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    n = data.GetNumberOfPoints()    #number of points
    dim = data.GetDimensions()      #dimensions of the domain
    nx=dim[0]
    ny=dim[1]
    yres=xres=c=0
    for p in range(n):
        pos=data.GetPoint(p)    #to get the position of the point
        x=pos[0]
        y=pos[1]
        d=data.GetPointData()
        array=d.GetArray('rho')
        density=array.GetValue(p)   #to get the value of the point density
        if density>=rhol*0.95:      #to find the center of the droplet
            c+=1
            yres+=y
            xres+=x
    ycenter.append(yres/c)
    xcenter.append(xres/c)
    time.append((i-xmin)*dt)
    height.append(yres/c/D)

print('Total execution time = {:5.2f}'.format((myTime.time()-startTime)/60.),"min")
print()
print("dx=",dx,"m")
print("dt=",dt,"ms")

figure(1)
xlabel("Time (ms)")
ylabel("H/D")
tight_layout()
if test=='water':
    plot(time,height,'g^',markerfacecolor='none',label="LBM Water")
    plot(t_w,h_w,'k-',label="Bertola Water")
elif test=='peo':
    plot(time,height,'g^',markerfacecolor='none',label="LBM PEO")
    plot(t_peo,h_peo,'ro',markerfacecolor='none',label="Bertola PEO")
elif test=='xg':
    plot(time,height,'g^',markerfacecolor='none',label="LBM XG")
    plot(t_xg,h_xg,'bd',markerfacecolor='none',label="Bertola XG")
xlim(0,140)
ylim(0,2.5)
legend()
savefig("Droplet_Height_"+test+".eps")
savefig("Droplet_Height_"+test+".png")
show()
