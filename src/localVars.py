from numpy import *; from numpy.linalg import *
import numpy
import os 
import sys
import pickle

class LocalVars:

    def init(self,globalVars,case):
        # fs
        self.fin = zeros((globalVars.q,case.nx,case.ny))
        self.feq = zeros((globalVars.q,case.nx,case.ny))
        # Temperature
        if case.EoS=="SC":
            self.Tim = ones((case.nx,case.ny))
        else:    
            self.Tim = case.ratioT*case.Tc*ones((case.nx,case.ny))
        if case.restart==True:
            if not os.path.exists("./Restarts/"+case.restartFilename):
                print("File:./Restarts/"+case.restartFilename+" does not exist")
                sys.exit()
            restartFile = open("./Restarts/"+case.restartFilename,"rb")
            print("Restart file:./Restarts/"+case.restartFilename+" is opening...")
            if case.solveTemp==True:
                [case.restartIter, self.fin, self.Tim] = pickle.load(restartFile)
            else:
                [case.restartIter, self.fin] = pickle.load(restartFile)
        
        #print(case.restartIter)        
    
        #self.fout = zeros((globalVars.q,case.nx,case.ny))
        if case.FlowConfiguration=='MultiPhase':
            #Shan-Chen force
            self.fsc = zeros((2,case.nx,case.ny))   
            self.fc = zeros((globalVars.q,2,case.nx,case.ny))  # sum part of Shan-Chen force
            #rho for neighbours lattice.    
            self.rh = zeros((globalVars.q,case.nx,case.ny))
            self.phi = zeros((case.nx,case.ny))
            self.Peos = zeros((case.nx,case.ny))
            #
                
    def Density(self,globalVars,case,lat,feq,fin):
        if case.boolInitDensity==True:                    
            if case.shape=='sphere':
                # 
                # Li et al. 2013, Eq. 31
                #for i in range (globalVars.q):
                #    for j in range(0,case.ny):
                #        for k in range(0,case.nx):
                #            r = sqrt((k-case.x0)**2.+(j-case.y0)**2.) 
                #            feq[i,k,j] = ((case.rhol+case.rhog)/2.-(case.rhol-case.rhog)/2.*tanh(2.*(r-case.r0)/5.))*lat.t[i]
                #            fin[i,k,j] = feq[i,k,j]
                for i in range (globalVars.q):
                    for j in range(0,case.ny):
                        for k in range(0,case.nx):
                             feq[i,k,j] = ((case.rhol+case.rhog)/2.-((case.rhol-case.rhog)/2.) \
                                          *tanh((2.*(((k-case.x0)**2.+(j-case.y0)**2.)**0.5-case.r0)/5.)))*lat.t[i]
                             fin[i,k,j] = feq[i,k,j] 
            elif case.shape=='rectangular':
                for i in range (globalVars.q):
                    for j in range(0,case.ny):
                        for k in range(0,case.nx):
                             feq[i,k,j] = case.rhog*lat.t[i]
                             if j>(case.y0-case.yLen/2-1) and j<(case.y0+case.yLen/2):
                                 if k>(case.x0-case.xLen/2-1) and k<(case.x0+case.xLen/2):
                                     feq[i,k,j] = case.rhol*lat.t[i]
                             fin[i,k,j] = feq[i,k,j]                                                              
        else:
            if case.restart==True:
                pass # Already read from restart file; lines 20-23
            elif case.restart==False:
                for i in range (globalVars.q):
                    feq[i,:,:] = case.rhol*lat.t[i]
                    fin[i,:,:] = feq[i,:,:]                        

        #print(fin)        
                
    def Velocity(self,case,rho,u):
        if case.boolInitVel==True:
            for j in range(case.ny):
                for i in range(case.nx):    
                    if (rho[i,j] > case.rhog):  
                        u[0,i,j] = case.initVx
                        u[1,i,j] = case.initVy
                        
    # Body force
    def bodyForce(self,globalVars,case):
        # Initial values
        self.fb = zeros((2,case.nx,case.ny))
        self.fb[0,:,:] = case.fbx*ones((case.nx,case.ny))
        self.fb[1,:,:] = case.fby*ones((case.nx,case.ny))
    
                        
