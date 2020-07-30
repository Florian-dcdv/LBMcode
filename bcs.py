from numpy import *; from numpy.linalg import *
import numpy
import sys

class BCs:
    
    def define(self,globalVars,case,lat): 
        if case.flowBcs['Format']=='Standard':
            self.noslip = [lat.c.tolist().index((-lat.c[i]).tolist()) for i in range(globalVars.q)]
            self.outlet = arange(globalVars.q)[asarray([ci[0]<0  for ci in lat.c])]
            self.obstacleb = False
            for key, bc in case.flowBcs.items():
                if key!='Format' and key!='WallDensity':
                    if bc=='Periodic':
                        pass
                    elif bc=='NoSlip':
                        if key=='Left':
                            self.obstacleb = fromfunction(lambda x,y: x==0, (case.nx,case.ny))
                        elif key=='Right':
                            self.obstacleb += fromfunction(lambda x,y: x==case.nx-1, (case.nx,case.ny))
                        elif key=='Bottom':
                            self.obstacleb += fromfunction(lambda x,y: y==case.ny-1, (case.nx,case.ny))
                        elif key=='Top':
                            self.obstacleb += fromfunction(lambda x,y: y==0, (case.nx,case.ny))
                    else:
                        print("Bc:",bc,"Not implemented...(bcs.py)") 
                                
        elif case.flowBcs['Format']=='Image':
            print("Not implemented...(bcs.py)") 
            sys.exit()
            
        #print("obstacleb",self.obstacleb)    
        
    def apply(self,globalVars,case,fout,fin):                               
        for i in range(globalVars.q): 
            fout[i,self.obstacleb] = fin[self.noslip[i],self.obstacleb]
