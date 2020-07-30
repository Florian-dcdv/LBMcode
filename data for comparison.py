from numpy import *
from matplotlib.pyplot import *
from scipy import integrate
    
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
    return imin,imax,rhomin,rhomax

#give a list of gas and liquid densities (rhog,rhol) to have p0=peos(rhog)=peos(rhol)
def pgl(acc):
    pgl=[]
    (imin,imax,rhomin,rhomax)=limit(acc)
    i=rhomin
    j=rhomax
    pmax=peos(imax,a,b,T,EOS,R)
    pmin=peos(imin,a,b,T,EOS,R)
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
    #print(listp0)
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



EOS='PR'
LBM='_LBM'
#LBM=''
Li=False
MC=False
spec=''
save=True

tempRatio=[]
pressure=[]
gasDensity=[]
liquidDensity=[]
acc=0.0001
 
if MC==True:
    t1=0.5
    t2=0.95
    deltaT=0.05   
    ratioT=t2
    while ratioT>=t1:
        if (EOS=='CS') or (EOS=='NDCS'):       # Input parameters for Carnahan-Starling equation of state (Li 2013)
            a=0.5                             #NDCS: Non Dimensionalised CS
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
        
        if ratioT<=0.5:
            (p0,rhomin,rhomax,rhog,rhol,zero,err)=densities(acc/100)
        else:
            (p0,rhomin,rhomax,rhog,rhol,zero,err)=densities(acc)
        pressure.append(p0)
        gasDensity.append(rhog)
        liquidDensity.append(rhol)
        tempRatio.append(ratioT)
        print(ratioT)
        ratioT=round(ratioT-deltaT,3)
    
    print("T/Tc=",tempRatio)
    print("rhog=",gasDensity)
    print("rhol=",liquidDensity)

#Results from pyLBM
if LBM=='_LBM':
    if EOS=='CS':
        tempLBM=[0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5]
        gasLBM=[0.066,0.0453,0.0316,0.0217,0.0144,0.00914,0.00543,0.00298,0.00149,0.000668]
        liquidLBM=[0.214,0.251,0.282,0.309,0.335,0.360,0.385,0.409,0.431,0.455]
    elif EOS=='PR':
        tempLBM=[0.95,0.9,0.85,0.8,0.75]
        gasLBM=[1.041,0.6078,0.3679,0.2231,0.135]
        liquidLBM=[4.9439,5.9217,6.6408,7.2136,7.6935]
        
#Reads Q.Li values from the fig.1 in (Li et al.,2013)    
if Li==True:
    T=[]
    R=[]
    file=open("QiResultsRhof.dat","r")
    line=file.readline()
    while line:
        rho=line.split()
        T.append(float(rho[1]))
        R.append(float(rho[0]))
        line=file.readline()
    file.close()
    T.reverse()
    R.reverse()
    T2=[]
    R2=[]
    file=open("QiResultsRhog.dat","r")
    line=file.readline()
    while line:
        rho=line.split()
        T2.append(float(rho[1]))
        R2.append(float(rho[0]))
        line=file.readline()
    file.close()
    T2.reverse()
    R2.reverse()

#Reads previous values from Maxwell construction
if MC==False:
    file=open(EOS+"_"+str(acc)+".txt","r")
    line=file.readline()
    line=file.readline()
    for i in range(10):
        rho=line.split()
        tempRatio.append(float(rho[0]))
        gasDensity.append(float(rho[1]))
        liquidDensity.append(float(rho[2]))
        line=file.readline()
    file.close()
    
#Graph and saves
figure(1)
#title("Thermodynamic consistency")
xlabel("Density")
ylabel("Reduced temperature")
plot(gasDensity,tempRatio,'ro',label='Maxwell construction',markerfacecolor="none")
plot(liquidDensity,tempRatio,'ro',markerfacecolor="none")
if LBM=='_LBM':
    plot(gasLBM,tempLBM,'bx',label='IMRT')
    plot(liquidLBM,tempLBM,'bx')
if Li==True:
    plot(R,T,'g+',label='Q.Li results')
    plot(R2,T2,'g+')
semilogx()
legend()
if save==True:
    savefig(EOS+"_"+str(acc)+LBM+spec+".png")
    savefig(EOS+"_"+str(acc)+LBM+spec+".eps")
    results=open(EOS+"_"+str(acc)+LBM+spec+".txt",'w')
    results.write("#Tr \t rho_g \t \t rho_l \n")
    for i in range(0,size(tempRatio)):
        if len(str(gasDensity[i]))>5:
            results.write('{} \t {} \t {} \n'.format(tempRatio[i],gasDensity[i],liquidDensity[i]))
        else:
            results.write('{} \t {} \t\t {} \n'.format(tempRatio[i],gasDensity[i],liquidDensity[i]))
    if LBM=='_LBM':
        results.write('\n')
        results.write("#Tr \t rho_g_LBM \t rho_l_LBM \n")
        for i in range(0,size(tempLBM)):
            if len(str(gasLBM[i]))>5:
                results.write('{} \t {} \t {} \n'.format(tempLBM[i],gasLBM[i],liquidLBM[i]))
            else:
                results.write('{} \t {} \t \t {} \n'.format(tempLBM[i],gasLBM[i],liquidLBM[i]))
    results.close()
show()