from numpy import *; from numpy.linalg import *
import numpy
from vtk.util import numpy_support
from globalVars import GlobalVars
from interface import Interface
from localVars import LocalVars
from lattice import *
from bcs import BCs
from heatEquation import HeatEquationSolver
from ioClass import ioLBM
import time as myTime

# Global variable class
globalVars = GlobalVars()
globalVars.init()

# Create interface class and read input file
case = Interface()
case.readInputVars() 

# Init variable class
localVars = LocalVars()
localVars.init(globalVars,case)

# Define lattice class
lat = Lattice()
lat.init(globalVars,case)

# Define bcs
bcs = BCs()
bcs.define(globalVars,case,lat)

# Init body force    
localVars.bodyForce(globalVars,case)

# Initial density field
localVars.Density(globalVars,case,lat,localVars.feq,localVars.fin)

# Heat equation solver class
heatEq = HeatEquationSolver()
heatEq.init(globalVars,case,localVars)

# Create IO class
ioCase = ioLBM() 

# CPU-start time
startTime = myTime.time()
#rho=ones((case.nx,case.ny))*case.rhow
#notobs=invert(bcs.obstacleb)
# Main Time Loop:
for time in range(case.restartIter+1,case.maxIter+1):
    if(time%case.nFreq==0): print("Timestep = ", time)
    
    # Apply wall bc to rho
    #rho[bcs.obstacleb] = case.rhow  
    #print(rho)      
    # Density
    #rho[notobs]= lat.sumpop(localVars.fin[:,notobs])   
    #print(rho)    
    
     
    # Density
    rho = lat.sumpop(localVars.fin)   
    #print(rho)      
    # Apply wall bc to rho
    rho[bcs.obstacleb] = case.rhow  
    #print("rho=",rho[87,131])      

    if case.FlowConfiguration=='MultiPhase':       
        ## Equation of State:
        if case.EoS=="SC":
            localVars.phi = case.rho0*(1-exp(-rho/case.rho0))
            localVars.Peos = rho/3. + case.cZero*case.g*localVars.phi**2./2. # For visualisation purposes
        elif case.EoS=="PR":    
            ksi = (1+(0.37464+1.54226*case.constAsc-0.26992*case.constAsc**2.)*(1-(localVars.Tim/case.Tc)**0.5))**2 
            localVars.Peos = (rho*globalVars.R*localVars.Tim)/(1.-case.constB*rho) - (case.constA*ksi*rho**2)  \
                            /(1+2*case.constB*rho-(case.constB**2)*(rho**2))
            localVars.phi = (2.*(localVars.Peos[:,:]-(rho[:,:]/3.))/(case.cZero*case.g))**0.5
        elif case.EoS=="CS":
            localVars.Peos = rho*globalVars.R*localVars.Tim*(1.+case.constB*rho/4.+(case.constB*rho/4.)**2.-(case.constB*rho/4.)**3.)\
                            /(1.-case.constB*rho/4.)**3.-case.constA*rho**2.
            localVars.phi = (2.*(localVars.Peos[:,:]-(rho[:,:]/3.))/(case.cZero*case.g))**0.5  
        #print("Peos=",localVars.Peos[87,131])
        #print("phi=",localVars.phi[87,131])
        
        #########  pseudopotential for neighbour lattices   #####
        for i in range(globalVars.q):       # Pseudopotential
            localVars.rh[i,:,:] = roll(roll(localVars.phi[:,:],lat.c[i,0],axis=0),lat.c[i,1],axis=1)
    
        ###### Shan-Chen force ########                   
        for i in range(globalVars.q):
            if case.latticeWeighting==True:
                localVars.fc[i,0,:,:]=lat.t[i]*localVars.rh[i,:,:]*lat.c[i,0]
                localVars.fc[i,1,:,:]=lat.t[i]*localVars.rh[i,:,:]*lat.c[i,1]
            else:
                localVars.fc[i,0,:,:] = lat.tPP[i]*localVars.rh[i,:,:]*lat.c[i,0]
                localVars.fc[i,1,:,:] = lat.tPP[i]*localVars.rh[i,:,:]*lat.c[i,1]     
        localVars.fsc[0,:,:] = localVars.phi[:,:]*case.g*lat.sumpop(localVars.fc[:,0,:,:])
        localVars.fsc[1,:,:] = localVars.phi[:,:]*case.g*lat.sumpop(localVars.fc[:,1,:,:])
        
        # Interaction force + Gravity force
        if case.rhol != case.rhog: # To run multiphase code for single phase
            ftot = localVars.fsc + localVars.fb*(rho-case.rhog)/(case.rhol-case.rhog)
        else:
            ftot = localVars.fsc + localVars.fb
    else:
        # Body force for single phase flow
        ftot = localVars.fb    
        
    # Velocity 
    u = dot(lat.c.transpose(), localVars.fin.transpose((1,0,2)))/rho + ftot/(2.*rho)
    #print("ftot=",ftot[:,87,131])
    #print("u=",u[:,87,131])
    #print("ux=",u[0,87,131])
    #print("uy=",u[1,87,131])
    
    # Initial velocity for heavier fluid       
    if (time==0 and case.boolInitVel==True): localVars.Velocity(case,rho,u)

    # Compute localVars.feq
    localVars.feq = lat.equilibrium(globalVars,case,lat,rho,u)
    # Updating relaxation Parameter for non-Newtonian fluid
    if (time>=0) and (case.newtonN!=1):
        fneq = localVars.fin - localVars.feq   # Non-equilibrium part of the distribution function
        for i in range (globalVars.q): e1 = fneq[i,:,:]*lat.c[i,0]*lat.c[i,1]   
        exy=(3*case.omega/2)*lat.sumpop(e1)
        uLB2=lat.nulb*(2*exy*exy)**((case.newtonN-1)/2)
        case.omega=1/(3*uLB2+0.5)
        print("There is a problem with non-Newtonian fluid law coding in line-65:pyLBM.py",case.omega)
    
    ## Solve heat equation:
    if case.solveTemp==True: 
        if(time%case.nFreq==0): print("\t Solving heat equation...")
        if case.HeatEqSolverType=='Explicit':
            heatEq.RK4(case,localVars,rho,u)
        elif case.HeatEqSolverType=='Implicit':
            heatEq.Implicit(case,localVars,rho,u)    
        #print("Tim=",localVars.Tim[87,131])                       
    #### Source Term Si #####    
    if case.relaxScheme=="SRT": 
        f2 = lat.SRT(globalVars,case,lat,ftot,u)                   
        # Collision step
        fout = localVars.fin - case.omega * (localVars.fin - localVars.feq) + f2
    #MRT Force Scheme 
    elif case.relaxScheme=="MRT" or case.relaxScheme=="IMRT":
        # Equation (9) and (25) in Reference [1].
        fout = lat.MRT(globalVars,case,localVars,ftot,u)          

    # Apply boundary conditions
    bcs.apply(globalVars,case,fout,localVars.fin)    
    
    # Streaming step.
    for i in range(globalVars.q):    
        localVars.fin[i,:,:] = roll(roll(fout[i,:,:],lat.c[i,0],axis=0),lat.c[i,1],axis=1)
    #print("fin=",localVars.fin[:,87,131]) 
                          
    # Output
    if ( (time%case.nFreq==0) or (time==case.maxIter) ):
        ioCase.output(case,localVars,time,rho,u)
        
    # Restart output
    if (time>0) and (time%case.RestartNFreq)==0:
        ioCase.restart(globalVars,case,localVars,time)    

# Print out total CPU time
print('Total execution time = {:5.2f}'.format((myTime.time()-startTime)/60.),"min")
