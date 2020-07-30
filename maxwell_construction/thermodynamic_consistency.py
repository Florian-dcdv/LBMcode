from numpy import *
from matplotlib.pyplot import *
from scipy import integrate

test=False
#Initialisation
EOS='CS'
ratioT=0.8

if (EOS=='CS') or (EOS=='NDCS'):       # Input parameters for Carnahan-Starling equation of state (Li 2013)
    a=0.5                            #NDCS: Non Dimensionalised CS
    b=4
    #a=0.4963*R**2*Tc**2/pc
    #b=0.18727*R*Tc/pc
    R=1 
    Tc=a/(2.65018*R*b)
    T=ratioT*Tc
    rhoc=0.5218/b
    pc=0.070663*a/b**2
elif EOS=='PR':     #Input parameters for Peng-Robinson equation of state (Li 2015)
    a=3/49
    b=2/21
    #a=0.45724*R**2*Tc**2/pc
    #b=0.0778*R*Tc/pc
    R=1    
    Tc=a/(5.87712*R*b)
    T=ratioT*Tc
    omega=0.344
    phi=(1+(0.37464+1.54226*omega-0.26992*omega**2)*(1-(ratioT)**0.5))**2
    pc=0.0778*R*Tc/b
elif EOS=='VdW':        #Input parameters for Van-der-Waals equation of state (Li 2015)
    a=9/49
    b=2/21
    #a=(27*R**2*Tc**2)/(64*pc)
    #b=(R*Tc)/(8*pc)
    R=1
    rhoc=1/(3*b)
    Tc=(8*a)/(27*b*R)
    T=ratioT*Tc
    pc=a/(27*b**2)

#Equation of State (PR, CS)
def peos(rho,a,b,T,EOS,R):
    if EOS=='CS':
        peos=rho*R*T*(1+b*rho/4+(b*rho/4)**2-(b*rho/4)**3)/(1-b*rho/4)**3-a*rho**2
    if EOS=='PR':
        peos=(rho*R*T)/(1-b*rho)-(a*phi*rho**2)/(1+2*b*rho-b**2*rho**2)
    if EOS=='NDCS':
        peos=2.786*pc*(rho/rhoc)*((ratioT*(1+0.13045*(rho/rhoc)+(0.13045*(rho/rhoc))**2 \
                                 +(0.13045*(rho/rhoc))**3)/(1-0.13045*(rho/rhoc))**3)-1.3829*(rho/rhoc))
    if EOS=='VdW':
        peos=(rho/rhoc)*R*ratioT/(1-b*(rho/rhoc))-a*(rho/rhoc)**2
    return peos

#Expression of the function for the Maxwell construction
def maxwell(rho,p0):
    return (p0-peos(rho,a,b,T,EOS,R))/rho**2

#Integral of the above function for the Equal Area Rule

def integral(rhog,rhol,p0):
    (res,err)=integrate.quad(maxwell,rhog,rhol,args=(p0,))
    return res,err

#To find the max and min of the possible rhog and rhol
def extrema(acc):    
    rho=2*acc
    imax=imin=0
    pmin=pmax=0
    if EOS=='VdW':
        condition=rhoc/b
    elif EOS=='CS':
        condition=4/b
    elif EOS=='NDCS':
        condition=rhoc/0.13045
    elif EOS=='PR':
        condition=1/b
    while ((imax==0) or (imin==0)) and (rho<condition):
        if sign(peos(rho+acc,a,b,T,EOS,R)-peos(rho,a,b,T,EOS,R))!=sign(peos(rho,a,b,T,EOS,R)-peos(rho-acc,a,b,T,EOS,R)):
            if sign(peos(rho+acc,a,b,T,EOS,R)-peos(rho,a,b,T,EOS,R))<0:
                pmax=peos(rho,a,b,T,EOS,R)
                imax=rho
            else:
                pmin=peos(rho,a,b,T,EOS,R)
                imin=rho
        rho+=acc
    return pmin,pmax,imin,imax

#(pmin,pmax,imin,imax)=extrema(0.001)
#print("pmin=",pmin,"| pmax=",pmax)
#print("imin=",imin,"| imax=",imax)

#To define the limits of the domain where the maxwell construction will be applied
def limit(acc):
    (pmin,pmax,imin,imax)=extrema(acc)
    if (imin==0) or (imax==0):
        raise SystemExit("It does not respect the Maxwell's rule")
    rhomin=rhomax=-1
    rho=acc
    while ((rhomin==-1) or (rhomax==-1)) and (peos(rho,a,b,T,EOS,R)<pmax+acc):
        if (peos(rho,a,b,T,EOS,R)>pmax-acc) and (peos(rho,a,b,T,EOS,R)<pmax+acc) and (rho>imin):
            rhomax=rho
        elif (peos(rho,a,b,T,EOS,R)>pmin-acc) and (peos(rho,a,b,T,EOS,R)<pmin+acc) and (rho<imax):
            rhomin=rho
        rho+=acc
    if rhomin==-1:
        rhomin=acc
    return imin,imax,pmin,pmax,rhomin,rhomax

#(imin,imax,rhomin,rhomax)=limit(0.0001)

#give a list of gas and liquid densities (rhog,rhol) to have p0=peos(rhog)=peos(rhol)
def pgl(acc):
    pgl=[]
    (imin,imax,pmin,pmax,rhomin,rhomax)=limit(acc)
    i=rhomin
    j=rhomax
    while i<imax:
        while j>imin:
            if (abs(peos(i,a,b,T,EOS,R)-peos(j,a,b,T,EOS,R))<(pmax-pmin)*acc):
                pgl.append((round(i,6),round(j,6)))
            j-=acc
        j=rhomax
        i+=acc
    return rhomin,rhomax,pgl

#calculate the densities to respect the Maxwell rule
def densities(acc):
    (rhomin,rhomax,listp0)=pgl(acc)
    #print(len(listp0))
    zero=err=1000
    rhog=rhol=pzero=p0=-1
    for k in listp0:
        p0=peos(k[0],a,b,T,EOS,R)
        (res,error)=integral(k[0],k[1],p0)
        #res2=integral2(k[0],k[1],int(1/acc),maxwell)
        #print(p0,k[0],k[1],res)
        if abs(res)<abs(zero):
            zero=res
            err=error
            (rhog,rhol)=k
            pzero=p0
    return pzero,rhomin,rhomax,rhog,rhol,zero,err

(p0,rhomin,rhomax,rhog,rhol,zero,err)=densities(0.0001)

print()
print("p0=",round(p0,10))
print("rhog=",rhog,"| rhol=",rhol,"| Density ratio rhol/rhog:",round(rhol/rhog,2))
print("EoS:",EOS,"a=",a,"b=",b,"R=",R)
print("Temperature ratio:",ratioT,"| Tc=",round(Tc,4),"| pc=",round(pc,4))
print("Integral from the Maxwell rule:",zero,"| error:",err)
print("min=",rhomin,"| max=",rhomax)

#Graph

rho=linspace(rhomin,rhomax,10000)
listp0=[]
t=linspace(rhog,rhol,1000)
listfill=[]
for i in rho:
    listp0.append(p0)
for j in t:
    listfill.append(p0)
    
figure(1)
#title("p_EoS")
xlabel("Density")
ylabel("Pressure")
plot(rho,peos(rho,a,b,T,EOS,R),'k',label='peos')
plot(rho,listp0,'r',label='p0')
plot(rhog,peos(rhog,a,b,T,EOS,R),'ro')
plot(rhol,peos(rhol,a,b,T,EOS,R),'ro')
#fill_between(t,listfill,peos(t,a,b,T,EOS,R))
legend()
tight_layout()
savefig('pressure_graph.eps')
savefig('pressure_graph.png')
figure(2)
#title("(p0-peos)/rho²")
xlabel("Density")
ylabel("Maxwell function")
plot(t,maxwell(t,p0),'b')
fill_between(t,maxwell(t,p0),0)
grid(True)
tight_layout()
savefig('maxwell_rule.png')
savefig('maxwell_rule.eps')
show()
