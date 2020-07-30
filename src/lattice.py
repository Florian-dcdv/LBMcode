from numpy import *; from numpy.linalg import *
import numpy

class Lattice:
    
    def init(self,globalVars,case): 
        # Defining lattice constants and nulb
        self.nulb = ((1./case.omega)-0.5)/3.
        #self.c = array([(x,y) for x in [0,-1,1] for y in [0,-1,1]])                 # Lattice velocities.
        #self.t = 1./36. * ones(globalVars.q)                                        # Lattice weights.
        #self.t[asarray([norm(ci)<1.1 for ci in self.c])] = 1./9.; self.t[0] = 4./9.
        self.c = array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]) # Lattice velocities.
        self.t = array([4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.])                      # Lattice weights.
        self.tPP = array([0, 1./3.,1./3.,1./3.,1./3.,1./12.,1./12.,1./12.,1./12.])  # Lattice weights for pseudopotential 
        
    # Density computation
    def sumpop(self,fin): 
        tot = sum(fin,axis=0) 
        return tot
    
    # Equilibrium distribution function.     
    def equilibrium(self,globalVars,case,lat,rho,u):             
        cu   = 3.0 * dot(lat.c,u.transpose(1,0,2))
        usqr = 3./2.*(u[0]**2+u[1]**2)
        feq = zeros((globalVars.q,case.nx,case.ny))
        for i in range(globalVars.q):
            feq[i,:,:] = rho*lat.t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
        return feq
        
    # SRT 
    def SRT(self,globalVars,case,lat,ftot,u):
        e = zeros((globalVars.q,2,case.nx,case.ny))
        p = zeros((globalVars.q,case.nx,case.ny))
        p2 = zeros((globalVars.q,2,case.nx,case.ny))
        f2 = zeros((globalVars.q,case.nx,case.ny))
        for i in range(globalVars.q): 
            e[i, 0, :, :] = subtract(lat.c[i, 0], u[0, :, :])
            e[i, 1, :, :] = subtract(lat.c[i, 1], u[1, :, :])
            p[i, :, :] = lat.c[i, 0] * u[0, :, :] + lat.c[i, 1] * u[1, :, :]
            p2[i, 0, :, :] = p[i, :, :] * lat.c[i, 0]
            p2[i, 1, :, :] = p[i, :, :] * lat.c[i, 1]
            z2 = e + 3.*p2
            for i in range(globalVars.q): 
                f2[i, :, :] = (1-0.5*case.omega)*lat.t[i]*3.*(z2[i, 0, :, :] * ftot[0, :, :] + z2[i, 1, :, :] * ftot[1, :, :])
        return f2
    
    # MRT & IMRT
    def MRT(self,globalVars,case,localVars,ftot,u):
        # MRT parameters
        #tp=1.;te=1.1;tt=1.1;tj=1.;tq=1.1;tv=case.omega   # relaxation time
        # Relaxation time s
        tv=1./case.omega 
        tp=1.;tj=1.     
        te=1./1.1;tt=1./1.1;tq=1./1.1
        #rel=array([tp,te,tt,tj,tq,tj,tq,tv,tv]) # relaxation time , eq.(2) in reference [1].
        rel=array([1./tp,1./te,1./tt,1./tj,1./tq,1./tj,1./tq,1./tv,1./tv]) # relaxation time , eq.(2) in reference [1].
        A = diag(rel)                 #diagonal matrix, eq.(2) in reference [1].

        ## Transformation matrix eq (10.30) in Reference [3] ,Page 420.
        M = array([[1,1,1,1,1,1,1,1,1],
                   [-4,-1,-1,-1,-1,2,2,2,2],
                   [4,-2,-2,-2,-2,1,1,1,1],
                   [0,1,0,-1,0,1,-1,-1,1],
                   [0,-2,0,2,0,1,-1,-1,1],
                   [0,0,1,0,-1,1,1,-1,-1],
                   [0,0,-2,0,2,1,1,-1,-1],
                   [0,1,-1,1,-1,0,0,0,0],
                   [0,0,0,0,0,1,-1,1,-1]])
     
        I = identity(9) #Unit Tensor    
        
        # Compute forces
        f3 = zeros((globalVars.q,case.nx,case.ny))   
        f3[1, :, :] = 6.*(u[0,:,:]*ftot[0,:,:]+u[1,:,:]*ftot[1,:,:])
        f3[2, :, :] = -6.*(u[0,:,:]*ftot[0,:,:]+u[1,:,:]*ftot[1,:,:])
        f3[3, :, :] = ftot[0,:,:]
        f3[4, :, :] = -ftot[0,:,:]
        f3[5, :, :] = ftot[1,:,:]
        f3[6, :, :] = -ftot[1,:,:]
        f3[7, :, :] = 2.*(u[0,:,:]*ftot[0,:,:]-u[1,:,:]*ftot[1,:,:])
        #f3[8, :, :] = (u[0,:,:]*ftot[0,:,:]+u[1,:,:]*ftot[1,:,:])
        f3[8, :, :] = (u[0,:,:]*ftot[1,:,:]+u[1,:,:]*ftot[0,:,:])
        
        if case.relaxScheme=="IMRT":
            f3[1, :, :] += (12.*case.sigma*(ftot[0,:,:]**2+ftot[1,:,:]**2))/(localVars.phi[:,:]**2.*(te-0.5))
            f3[2, :, :] -= (12.*case.sigma*(ftot[0,:,:]**2+ftot[1,:,:]**2))/(localVars.phi[:,:]**2.*(tt-0.5))

        ### Equation (4) in Refrence [1]:
        mi = dot(M,localVars.fin.transpose(1,0,2))
        meq = dot(M,localVars.feq.transpose(1,0,2))
        mto = mi-dot(A,(mi-meq).transpose(1,0,2))+dot((I-A/2),f3.transpose(1,0,2))
            
        # Collision step
        fout = zeros((globalVars.q,case.nx,case.ny))
        fout = dot(inv(M),mto.transpose(1,0,2))
        return fout              
                
                    
