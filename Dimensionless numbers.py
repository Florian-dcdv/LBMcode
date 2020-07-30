from numpy import *

#non dimensional study of the Black & Bertola (2013) case
test='water'
a=0.50
Tc=a/(2.65018*4)
print("a=",a,"| Tc=",Tc)
#Physical parameters at 20°C

if test=='water':
    gamma_p=0.072
    D_p=3.09e-3
    D=50
elif test=='peo':
    gamma_p=0.070
    D_p=2.93e-3
    D=44
elif test=='xg':
    gamma_p=0.071
    D_p=3.12e-3

rhol_p=1000
rhog_p=1.2
mug_p=1.85e-5
nug_p=mug_p/rhog_p
g_p=9.81
drho_p=rhol_p-rhog_p

print("Dimensionless numbers:")
We=20
print("We =",We)
u_p=((We*gamma_p)/(drho_p*D_p))**0.5
Re=(D_p*u_p)/nug_p
print("Re =",Re)
Oh=We**0.5/Re
print("Oh =",Oh)
H_D=(We**3*gamma_p**3)/(Re**4*2*g_p*drho_p**3*nug_p**4)+1
print("H/D =",round(H_D,1))

print()
#Lattice parameters

print("Lattice parameters:")
Tr=0.55              #reduced temperature
rhol=0.431
rhog=0.00149
print("Tr =",Tr)
print("rhol=",rhol," |rhog=",rhog)
print("D =",D)
print()

print("Results:")

H=H_D*D     #height of the center of the droplet
#H=250
print("H=",round(H,0))

#drho=0.455-6.698e-4     #density difference
drho=rhol-rhog
#gamma = 1.265e-2*a + 3.892e-3        #surface tension
gamma=8.871e-3
print("gamma=",gamma)

u=((We*gamma)/(drho*D))**0.5    #velocity
print("u=",u)
g=u**2/(2*(H-D))    #gravity
print("g=",g)
Fb=g*drho           #body force
print("Fb=",Fb)
nu=D*u/Re               #kinematic viscosity
print("viscosity=",nu)
tauv=3*nu+0.5           #relaxation time
print("tauv=",tauv)
print()
print("Conversion coefficient:")
dx=D_p/D                 #physical lattice constant in m
dt=(g*dx/g_p)**0.5       #physical time step in ms
dm=drho_p/drho*dx**3     #physical mass in kg
print("dx=",dx,"m")
print("dt=",dt,"s")
print("dm=",dm,"kg")
#print("gamma_p=",gamma*dm/dt**2)
print()
print("Heated wall:")
cv=5
print("specific heat=",cv)

cond=2/3
print("conductivity=",cond)
T_top=Tr*Tc
T_bottom=Tc
print("Temperature: Top=",T_top," Bottom=",T_bottom)
