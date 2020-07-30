from numpy import *; from numpy.linalg import *
import matplotlib.pyplot as plt; from matplotlib import cm
import numpy
import sys
import scipy.sparse
import scipy.sparse.linalg

class HeatEquationSolver:
    
    def init(self,globalVars,case,localVars):
        # Gas constant
        self.Rgas = globalVars.R
        # Initial temperature field 
        #self.Tim = case.ratioT*case.Tc*ones((case.nx,case.ny))
        #self.Tim = case.ratioT*case.Tc*ones((case.nx,case.ny))
        #if case.FlowConfiguration=='MultiPhase':
        if case.boolInitTemperature==True:
            # If two-phase sphere
            if case.shapeTemp=='sphere':
                #obstacle = fromfunction(lambda x,y: (x-case.x0Temp)**2.+(y-case.y0Temp)**2.<case.r0Temp**2., (case.nx,case.ny))
                #self.Tim[obstacle] = case.initT0
                #localVars.Tim[obstacle] = case.initT0 
                for j in range(0,case.ny):
                        for k in range(0,case.nx):
                             localVars.Tim[k,j] = ((case.InitT0+case.Tc)/2.-((case.InitT0-case.Tc)/2.) \
                                          *tanh((2.*(((k-case.x0Temp)**2.+(j-case.y0Temp)**2.)**0.5-case.r0Temp)/5.)))
                
        if case.solveTemp==False:
            pass
        elif case.solveTemp==True:        
                #    
                # Apply Temp bcs (periodic is default: do nothing)
                # Left
                if case.leftBcTemp=='Adiabatic':
                    #obstacleb = fromfunction(lambda x,y: x==0, (case.nx,case.ny))        
                    #self.Tim[obstacleb] = case.leftBcTempValue
                    #self.Tim[:,0] = case.leftBcTempValue
                    localVars.Tim[0,:] = case.leftBcTempValue
                elif case.leftBcTemp=='Periodic':
                    pass 
                else:
                    print("Not implemented")
                    sys.exit()
                # Right                                
                if case.rightBcTemp=='Adiabatic':        
                    #obstacleb = fromfunction(lambda x,y: x==case.nx-1, (case.nx,case.ny))
                    #self.Tim[obstacleb] = case.rightBcTempValue
                    #self.Tim[:,case.ny-1] = case.rightBcTempValue
                    localVars.Tim[case.nx-1,:] = case.rightBcTempValue
                elif case.rightBcTemp=='Periodic':
                    pass 
                else:
                    print("Not implemented")
                    sys.exit()
                # Bottom
                if case.bottomBcTemp=='Adiabatic':        
                    #obstacleb = fromfunction(lambda x,y: y==case.ny-1, (case.nx,case.ny))
                    #self.Tim[obstacleb] = case.bottomBcTempValue
                    #self.Tim[case.nx-1,:] = case.bottomBcTempValue
                    localVars.Tim[:,0] = case.bottomBcTempValue
                elif case.bottomBcTemp=='Periodic':
                    pass
                else:
                    print("Not implemented")
                    sys.exit()
                # Top                               
                if case.topBcTemp=='Adiabatic':        
                    #obstacleb = fromfunction(lambda x,y: y==0, (case.nx,case.ny))
                    #self.Tim[obstacleb] = case.topBcTempValue
                    #self.Tim[0,:] = case.topBcTempValue
                    localVars.Tim[:,case.ny-1] = case.topBcTempValue
                elif case.topBcTemp=='Periodic':
                    pass    
                else:
                    print("Not implemented")
                    sys.exit()
        else:
            print("Wrong solveTemp key:",case.solveTemp)
            sys.exit()         
                                                
        
    # Heat equation function.
    def rks(self,case,rho,u,Tk):                  
   
        dTx1=(roll(Tk,-1,1) - roll(Tk,1,1))/2.
        dTy1=(roll(Tk,1,0) - roll(Tk,-1,0))/2.
        dTx2=(roll(Tk,-1,1)-2.*Tk +roll(Tk,1,1))
        dTy2=roll(Tk,1,0)-2.*Tk +roll(Tk,-1,0)
        dux1=(roll(u[1,:,:],-1,1) - roll(u[1,:,:],1,1))/2.
        duy1=(roll(u[0,:,:],1,0) - roll(u[0,:,:],-1,0))/2.
        dPeos=rho*self.Rgas/(1.-case.constB*rho) 
        #dPeos = rho*self.Rgas*(1.+case.constB*rho/4.+(case.constB*rho/4.)**2.-(case.constB*rho/4.)**3.)/(1.-case.constB*rho/4.)**3.-case.constA*rho**2.
        
        DT=-1.*(u[0,:,:]*dTx1+u[1,:,:]*dTy1)+(case.lam/(rho*case.cv))*(dTx2 + dTy2)-(Tk/(rho*case.cv))*(dPeos)*(dux1-duy1)
        
        # Apply Temp bc for DT (recomputing bc in thre previous avoided for simpler coding)
        # Left
        if case.leftBcTemp=='Adiabatic':        
            DT[0,:] = 0.
        # Right                                
        if case.rightBcTemp=='Adiabatic':        
            DT[case.nx-1,:] = 0.
        # Bottom
        if case.bottomBcTemp=='Adiabatic':        
            DT[:,0] = 0.
        # Top                                  
        if case.topBcTemp=='Adiabatic':        
            DT[:,case.ny-1] = 0.

        return DT
        
    # 4th order Runge-Kutta scheme
    def RK4(self,case,localVars,rho,u):
        Tk = localVars.Tim
        h1 = self.rks(case,rho,u,Tk)

        Tk = localVars.Tim + h1/2.
        h2 = self.rks(case,rho,u,Tk)

        Tk = localVars.Tim + h2/2.
        h3 = self.rks(case,rho,u,Tk)

        Tk = localVars.Tim + h3/2.
        h4 = self.rks(case,rho,u,Tk)

        localVars.Tim += (h1 + 2.*h2 + 2.*h3+ h4)/6.
        
    # Implicit solver
    ''''
    def Implicit(self,case,localVars,rho,u):
        #Tn = localVars.Tim.reshape(case.nx*case.ny,1)
        Tn = zeros(case.nx*case.ny)
        for j in range(0,case.ny):
            for i in range(0,case.nx):
                Tn[case.nx*j+i] = localVars.Tim[i,j]
        X = case.lam/(rho*case.cv)
        deltaVel = zeros((case.nx,case.ny))
        for j in range(1,case.ny-1):
            for i in range(1,case.nx-1):
                drhodx = rho[i,j+1] - rho[i,j]
                drhody = rho[i+1,j] - rho[i,j]        
                duxdx = (u[0,i,j+1]-u[0,i,j-1])/2.
                if drhodx==0: duxdx = u[0,i,j]-u[0,i,j-1]
                duydy = -(u[1,i+1,j]-u[1,i-1,j])/2. 
                if drhody==0: duydy = -u[1,i,j]-u[1,i-1,j]    
                deltaVel[i,j] = duxdx + duydy
        # Based on EOS
        #deltaVel = 0.01*ones((case.nx,case.ny))
        #if case.EoS=="SC":
        #    Y = 1./(rho*case.cv)*(rho*self.Rgas)*deltaVel
        #elif case.EoS=="PR":    
        #    Y = 1./(rho*case.cv)*(rho*self.Rgas)/(1.-case.constB*rho)*deltaVel
        #elif case.EoS=="CS":
        #    Y = 1./(rho*case.cv)*rho*self.Rgas*(1.+case.constB*rho/4.+(case.constB*rho/4.)**2.-(case.constB*rho/4.)**3.)/((1.-case.constB*rho/4.)**3.)*deltaVel
        #print(amax(Y),amin(Y))
        Y = 1./(rho*case.cv)*(rho*self.Rgas)*deltaVel
        #Y = zeros((case.nx,case.ny))
        D = 1.+Y+4.*X
        E = u[0,:,:]/2.-X
        F = -(u[0,:,:]/2.+X)
        G = u[1,:,:]/2.-X
        H = -(u[1,:,:]/2.+X)
        LHS = identity(case.nx*case.ny)         
        for j in range(1,case.ny-1):
            for i in range(1,case.nx-1):
                LHS[case.nx*j+i,case.nx*j+i] = D[i,j]
                LHS[case.nx*j+i,case.nx*j+i+1] = E[i,j]
                LHS[case.nx*j+i,case.nx*j+i-1] = F[i,j]
                LHS[case.nx*j+i,case.nx*(j-1)+i] = H[i,j]
                LHS[case.nx*j+i,case.nx*(j+1)+i] = G[i,j]                
        mm = scipy.sparse.csr_matrix(LHS)
        Tn1 = scipy.sparse.linalg.spsolve(mm,Tn)
        #localVars.Tim = Tn1.reshape(case.nx,case.ny)         
        for j in range(0,case.ny):
            for i in range(0,case.nx):
                localVars.Tim[i,j] = Tn1[case.nx*j+i] 

    '''
    # Implicit solver
    def Implicit(self,case,localVars,rho,u):                
        Tn = zeros(case.nx*case.ny)
        Xn = zeros(case.nx*case.ny)
        Yn = zeros(case.nx*case.ny)
        uxn = zeros(case.nx*case.ny)
        uyn = zeros(case.nx*case.ny)
        for j in range(0,case.ny):
            for i in range(0,case.nx):
                Tn[case.nx*j+i] = localVars.Tim[i,j]
                Xn[case.nx*j+i] = case.lam/(rho[i,j]*case.cv) 
                deltaVel = 0
                rhoI = rho[i,j]
                if (i>0) and (j>0) and (i<case.nx-1) and (j<case.ny-1): 
                    drhodx = rho[i,j+1] - rho[i,j]
                    drhody = rho[i+1,j] - rho[i,j]        
                    duxdx = (u[0,i,j+1]-u[0,i,j-1])/2.
                    if drhodx==0: duxdx = u[0,i,j]-u[0,i,j-1]
                    duydy = -(u[1,i+1,j]-u[1,i-1,j])/2. 
                    if drhody==0: duydy = -u[1,i,j]-u[1,i-1,j]    
                    deltaVel = duxdx + duydy
                    #deltaVel = 0.5*(u[0,i,j+1]-u[0,i,j-1]-u[1,i+1,j]-u[1,i-1,j])
                #print(deltaVel)
                deltaVel = 0              
                if case.EoS=="SC":
                    Yn[case.nx*j+i] = 1./(rhoI*case.cv)*(rhoI*self.Rgas)*deltaVel
                elif case.EoS=="PR":    
                    Yn[case.nx*j+i] = 1./(rhoI*case.cv)*(rhoI*self.Rgas)/(1.-case.constB*rhoI)*deltaVel
                elif case.EoS=="CS":
                    Yn[case.nx*j+i] = 1./(rhoI*case.cv)*rhoI*self.Rgas* \
                     (1.+case.constB*rhoI/4.+(case.constB*rhoI/4.)**2.-(case.constB*rhoI/4.)**3.) \
                    /((1.-case.constB*rhoI/4.)**3.)*deltaVel
                #print(amax(Y),amin(Y))
                #Yn[case.nx*j+i] = self.Rgas/case.cv/(1.-case.constB*rho[i,j])*deltaVel
                uxn[case.nx*j+i] = u[0,i,j]
                uyn[case.nx*j+i] = u[1,i,j]
    
        D = 1.+Yn+4.*Xn       
        E = uxn/2.-Xn
        F = -(uxn/2.+Xn)
        G = uyn/2.-Xn
        H = -(uyn/2.+Xn)
        
        # LHS
        LHS = scipy.sparse.diags([D, E, F[1:],G,H[case.nx:]], [0, 1, -1, case.nx, -case.nx], format='csc')
    
        # Update boundaries
        i=j=0   
        LHS[case.nx*j+i,case.nx*j+i] = 1
        LHS[case.nx*j+i,case.nx*j+i+1] = 0
        #LHS[nx*j+i,nx*j+i-1] = 0
        #LHS[nx*j+i,nx*(j-1)+i] = 0
        LHS[case.nx*j+i,case.nx*(j+1)+i] = 0
    
        for i in range(1,case.nx):
            j = 0 
            LHS[case.nx*j+i,case.nx*j+i] = 1
            LHS[case.nx*j+i,case.nx*j+i+1] = 0
            LHS[case.nx*j+i,case.nx*j+i-1] = 0
            #LHS[nx*j+i,nx*(j-1)+i] = 0
            LHS[case.nx*j+i,case.nx*(j+1)+i] = 0
    
        i=case.nx-1
        j=case.ny-1
        LHS[case.nx*j+i,case.nx*j+i] = 1
        #LHS[nx*j+i,nx*j+i+1] = 0
        LHS[case.nx*j+i,case.nx*j+i-1] = 0
        LHS[case.nx*j+i,case.nx*(j-1)+i] = 0
        #LHS[nx*j+i,nx*(j+1)+i] = 0
    
        for i in range(0,case.nx-1):
            j = case.ny-1 
            LHS[case.nx*j+i,case.nx*j+i] = 1
            LHS[case.nx*j+i,case.nx*j+i+1] = 0
            LHS[case.nx*j+i,case.nx*j+i-1] = 0
            LHS[case.nx*j+i,case.nx*(j-1)+i] = 0
            #LHS[nx*j+i,nx*(j+1)+i] = 0
    
        for j in range(1,case.ny-1):
            i = 0
            LHS[case.nx*j+i,case.nx*j+i] = 1
            LHS[case.nx*j+i,case.nx*j+i+1] = 0
            LHS[case.nx*j+i,case.nx*j+i-1] = 0
            LHS[case.nx*j+i,case.nx*(j-1)+i] = 0
            LHS[case.nx*j+i,case.nx*(j+1)+i] = 0

        for j in range(1,case.ny-1):
            i = case.nx-1
            LHS[case.nx*j+i,case.nx*j+i] = 1
            LHS[case.nx*j+i,case.nx*j+i+1] = 0
            LHS[case.nx*j+i,case.nx*j+i-1] = 0
            LHS[case.nx*j+i,case.nx*(j-1)+i] = 0
            LHS[case.nx*j+i,case.nx*(j+1)+i] = 0
                            
        # Solve system   
        Tn1 = scipy.sparse.linalg.spsolve(LHS,Tn)
        for j in range(0,case.ny):
            for i in range(0,case.nx):
                localVars.Tim[i,j] = Tn1[case.nx*j+i]                
       

    
  
