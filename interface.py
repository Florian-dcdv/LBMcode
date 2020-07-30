from numpy import *; from numpy.linalg import *
import numpy
import os
import vtk
from vtk.util import numpy_support
import ast
import sys
class Interface:
    def readInputVars(self): 
    
        # Dictionaries & default variables
        self.maxIter=self.nFreq=1
        self.RestartNFreq=0
        self.restart=False;self.restartFilename='t0.dat';self.restartIter=-1
        solverOptions = {}

        self.nx=self.ny=self.nz=1
        domain = {}

        self.fbx=self.fby=0
        bodyForces = {}

        self.initVx=self.initVy=0
        boolInitVel=False
        initialVelocities = {}

        self.FlowConfiguration='SinglePhase'
        self.rho0=self.rhol=self.rhog=self.newtonN=1
        self.g=-1
        fluidProperties = {}

        self.EoS='SC'
        self.rho0=self.Tc=self.ratioT=1
        self.constA=self.constB=1
        self.sigma=0.11
        self.cZero=1./3.
        latticeWeighting=True
        EoSParams={}
            
        self.solveTemp=False
        self.HeatEqSolverType='Explicit'
        self.cv=self.lam=1
        self.leftBcTemp=self.rightBcTemp=self.bottomBcTemp=self.topBcTemp='Periodic'
        self.leftBcTempValue=self.rightBcTempValue=self.bottomBcTempValue=self.topBcTempValue=self.Tc
        heatEquationParams = {}
    
        self.relaxScheme='SRT'
        self.omega=1
        relaxationScheme = {}
     
        self.outputF='screen'
        self.saveSnapshot=False
        self.lineXmin=self.lineYmin=self.lineXmax=self.lineYmax=1
        outputFormat={}

        self.boolInitDensity=False
        self.shape='sphere'
        self.x0=self.y0=self.r0=1
        self.xLen=self.yLen=1
        initialDensities={} 

        self.boolInitTemperature=False
        self.shapeTemp='sphere'
        self.x0Temp=self.y0Temp=self.r0Temp=1
        self.xLenTemp=self.yLenTemp=1
        initialTemperatures={} 
               
        self.formatBc='Standard'
        self.leftBc=self.rightBc=self.bottomBc=self.topBc='Periodic'
        self.flowBcs={}
        self.rhow=0.
        
        # Open input file
        finput = open('input.command','r')
        while True:
            line = finput.readline()
            if not line.startswith('#') and line.strip():
                if line.startswith('solverOptions',0,13):
                    solverOptions = ast.literal_eval(line[14:])
                    # Max. iterations
                    self.maxIter = solverOptions['iteration']
                    # Output freq.
                    self.nFreq = solverOptions['nFreq']
                    self.RestartNFreq = solverOptions['restartOutputFreq']
                    self.restart = solverOptions['Restart']
                    if self.restart==True:
                        self.restartFilename=solverOptions['RestartFilename']
                elif line.startswith('domain',0,6): 
                    domain = ast.literal_eval(line[7:])
                    # nx,ny,nz
                    self.nx = domain['nx']
                    self.ny = domain['ny']
                    self.nz = domain['nz']
                elif line.startswith('bodyForces',0,10): 
                    bodyForces = ast.literal_eval(line[11:])
                    # fbx,fby
                    self.fbx = bodyForces['bodyForcex']
                    self.fby = bodyForces['bodyForcey']
                elif line.startswith('initialVelocities',0,17): 
                    initialVelocities = ast.literal_eval(line[18:])
                    # initVx,initVy
                    self.boolInitVel = initialVelocities['initVelocity']
                    if self.boolInitVel==True:
                        self.initVx = initialVelocities['initVx']
                        self.initVy = initialVelocities['initVy']   
                elif line.startswith('fluidProperties',0,15): 
                    fluidProperties = ast.literal_eval(line[16:])
                    # rhol, rhog, n, g
                    self.newtonN = fluidProperties['newtonN'] 
                    self.rhol = fluidProperties['liquidDensity']
                    self.FlowConfiguration = fluidProperties['FlowConfiguration'] 
                    if self.FlowConfiguration=='MultiPhase':
                        self.rhog = fluidProperties['liquidGas'] 
                        self.g = fluidProperties['cohesionStrength']  
                elif line.startswith('heatEquationParams',0,18): 
                    heatEquationParams = ast.literal_eval(line[19:])
                    # solveTemp, Tc, cv, lam
                    self.solveTemp = heatEquationParams['solveTemperature']
                    if self.solveTemp==True:
                        self.HeatEqSolverType=heatEquationParams['SolverType']
                        self.cv = heatEquationParams['specificHeat']  
                        self.lam = heatEquationParams['thermalConductivity'] 
                        self.leftBcTemp = heatEquationParams['Left']
                        self.rightBcTemp = heatEquationParams['Right']
                        self.bottomBcTemp = heatEquationParams['Bottom']
                        self.topBcTemp = heatEquationParams['Top']
                        if self.leftBcTemp=='Adiabatic':
                            self.leftBcTempValue = heatEquationParams['LeftValue']
                        if self.rightBcTemp=='Adiabatic':
                            self.rightBcTempValue = heatEquationParams['RightValue']
                        if self.bottomBcTemp=='Adiabatic':
                            self.bottomBcTempValue = heatEquationParams['BottomValue'] 
                        if self.topBcTemp=='Adiabatic':
                            self.topBcTempValue = heatEquationParams['TopValue']
                        initialTemperatures = heatEquationParams['initialTemperatures']
                        #print(initialTemperatures)     
                        self.boolInitTemperature = initialTemperatures['initTemperature']
                        if self.boolInitTemperature==True:
                            self.shapeTemp = initialTemperatures['shape']
                            if self.shapeTemp=='sphere':
                                self.x0Temp = initialTemperatures['x0']
                                self.y0Temp = initialTemperatures['y0']
                                self.r0Temp = initialTemperatures['r0'] 
                                self.InitT0 = initialTemperatures['Value'] 
                            elif self.shapeTemp=='rectangular':
                                self.x0Temp = initialTemperatures['x0']
                                self.y0Temp = initialTemperatures['y0']
                                self.xLenTemp = initialTemperatures['xLen']
                                self.yLenTemp = initialTemperatures['yLen'] 
                                self.InitT0 = initialTemperatures['Value']                             
                            else:
                                sys.exit()                                                                                                               
                elif line.startswith('relaxationScheme',0,16): 
                    relaxationScheme = ast.literal_eval(line[17:])
                    self.relaxScheme = relaxationScheme['scheme']
                    # Omega
                    self.omega=1./relaxationScheme['tauv']
                    if self.relaxScheme=="IMRT":
                        self.sigma = relaxationScheme['sigma']
                elif line.startswith('EoSParams',0,9): 
                    EoSParams = ast.literal_eval(line[10:])
                    self.EoS = EoSParams['EoS']
                    self.latticeWeighting = EoSParams['latticeWeighting']
                    if self.EoS=='SC':
                        self.rho0 = EoSParams['refDensity']                     
                    elif self.EoS=='PR':
                        self.Tc = EoSParams['criticalTemperature']
                        self.ratioT = EoSParams['temperatureRatio']
                        self.constAsc = EoSParams['constAsc']
                        self.constA = EoSParams['constA']
                        self.constB = EoSParams['constB']
                        self.cZero = EoSParams['cZero']
                    elif self.EoS=='CS':
                        self.Tc = EoSParams['criticalTemperature']
                        self.ratioT = EoSParams['temperatureRatio']
                        self.constA = EoSParams['constA']
                        self.constB = EoSParams['constB'] 
                        self.cZero = EoSParams['cZero'] 
                elif line.startswith('outputFormat',0,12):
                    outputFormat = ast.literal_eval(line[13:])
                    self.outputF = outputFormat['format']
                    if  self.outputF=='screen':
                        self.saveSnapshot = outputFormat['saveSnapshot']
                    elif self.outputF=='line':
                        self.lineXmin = outputFormat['xMin']
                        self.lineYmin = outputFormat['yMin']                                                               
                        self.lineXmax = outputFormat['xMax']
                        self.lineYmax = outputFormat['yMax']
                elif line.startswith('initialDensities',0,16): 
                    initialDensities = ast.literal_eval(line[17:])
                    self.boolInitDensity = initialDensities['initDensity']
                    if self.boolInitDensity==True:
                        self.shape = initialDensities['shape']
                        if self.shape=='sphere':
                            self.x0 = initialDensities['x0']
                            self.y0 = initialDensities['y0']
                            self.r0 = initialDensities['r0'] 
                        elif self.shape=='rectangular':
                            self.x0 = initialDensities['x0']
                            self.y0 = initialDensities['y0']
                            self.xLen = initialDensities['xLen']
                            self.yLen = initialDensities['yLen']                             
                        else:
                            sys.exit()                             
                elif line.startswith('flowBoundaryConditions',0,22): 
                    self.flowBcs = ast.literal_eval(line[23:])
                    self.formatBc = self.flowBcs['Format']
                    if self.formatBc=='Standard': 
                        self.leftBc = self.flowBcs['Left']
                        self.rightBc = self.flowBcs['Right']
                        self.bottomBc = self.flowBcs['Bottom']
                        self.topBc = self.flowBcs['Top']
                    if self.leftBc == 'NoSlip' or self.rightBc == 'NoSlip' or \
                    self.bottomBc == 'NoSlip' or self.topBc == 'NoSlip':
                        self.rhow = self.flowBcs['WallDensity']                      
                    elif self.formatBc=='Image':
                        print("Not implemented...")
            if not line:
                break
        finput.close()
        
        # Update some variables
        if self.restart==True:
            self.boolInitDensity = False
            self.boolInitVel = False
            self.boolInitTemperature = False
            
        print("\n")
        print("---------------------------------------------------")
        print("Domain:",self.nx,"x",self.ny,"x",self.nz)
        print("Relaxation scheme =",self.relaxScheme)
        print("\t Omega =",self.omega)
        if self.relaxScheme=="IMRT":
            print("\t Sigma =",self.sigma)
        print("Max. iteration =",self.maxIter)
        print("\t Freq. output =",self.nFreq)
        print("\t Restart freq. =",self.RestartNFreq)
        if self.restart==True:
            print("\t Restarting from file:",self.restartFilename)    
        print("Body forces = [",self.fbx,",",self.fby,"]")
        if self.FlowConfiguration=='SinglePhase':
            print("Liquid density =",self.rhol)
            print("\t Shear law coeff. =",self.newtonN)            
        elif self.FlowConfiguration=='MultiPhase':
            print("Liquid density =",self.rhol,"Gas density=",self.rhog)
            print("\t Shear law coeff. =",self.newtonN)
            print("\t Cohesive strength =",self.g)
            print("\t Lattice weightings for PP =",self.latticeWeighting)
            print("Equation of state =",self.EoS)
            print("\t Critical Temp =",self.Tc,"T/Tc =",self.ratioT)
            if self.EoS=='SC':
                print("\t Ref. density =",self.rho0)
            elif self.EoS=='PR':
                print("Constants asc. =",self.constAsc,"a =",self.constA,"b =",self.constB, "cZero =",self.cZero)
            elif self.EoS=='CS':
                print("Constants: a =",self.constA,"b =",self.constB,"sigma =",self.sigma, "cZero =",self.cZero)
        print("Flow boundary conditions:",self.formatBc)
        if self.formatBc=='Standard': 
            print("\t Left:",self.leftBc,"Right:",self.rightBc)
            print("\t Bottom:",self.bottomBc,"Top:",self.topBc)
            if self.leftBc == 'NoSlip' or self.rightBc == 'NoSlip' or \
                self.bottomBc == 'NoSlip' or self.topBc == 'NoSlip':
                print("\t Wall density:",self.rhow)     
        elif self.formatBc=='Image':
            print("Not implemented...")
        print("Solving heat equation =",self.solveTemp)
        if self.solveTemp==True: 
            print("\t Solver type =",self.HeatEqSolverType)
            print("\t Specific heat cap. =",self.cv,"Thermal cond. =",self.lam)
            if self.boolInitTemperature==True:
                print("Initial temperature shape =",self.shape)
                if self.shape=='sphere': print("\t x0 =",self.x0,"y0 =",self.y0,"r0 =",self.r0, "Value =",self.InitT0)
                if self.shape=='rectangular': print("\t x0 =",self.x0,"y0 =",self.y0,"xLen =",self.xLenTemp,"yLen =",self.yLenTemp,"Value =",self.initT0) 
            print("\t Bcs for temperature:")
            print("\t \t Left surface:",self.leftBcTemp)
            if self.leftBcTemp=='Adiabatic': print("\t \t Value =",self.leftBcTempValue)
            print("\t \t Right surface:",self.rightBcTemp)
            if self.rightBcTemp=='Adiabatic': print("\t \t Value =",self.rightBcTempValue)
            print("\t \t Bottom surface:",self.bottomBcTemp)
            if self.bottomBcTemp=='Adiabatic': print("\t \t Value =",self.bottomBcTempValue)
            print("\t \t Top surface:",self.topBcTemp)
            if self.topBcTemp=='Adiabatic': print("\t \t Value =",self.topBcTempValue)                      
        print("Initial vel. = [",self.initVx,",",self.initVy,"]")
        if self.boolInitDensity==True:
            print("Initial heavier fluid shape =",self.shape)
            if self.shape=='sphere': print("\t x0 =",self.x0,"y0 =",self.y0,"r0 =",self.r0)
            if self.shape=='rectangular': print("\t x0 =",self.x0,"y0 =",self.y0,"xLen =",self.xLen,"yLen =",self.yLen)   
        print("Output format:",self.outputF)
        if self.outputF=="line":
            print("\t Line = [",self.lineXmin,":",self.lineXmax,",",self.lineYmin,":",self.lineYmax,"]")          
        print("---------------------------------------------------")
        print("\n")

        # Variables for VTK outputs (lattices are StructuredPoints)
        if self.outputF=='VTK':
            self.domain = vtk.vtkStructuredPoints()
            self.domain.SetDimensions(self.nx, self.ny, self.nz)
            self.domain.SetOrigin(0,0,0)
            self.domain.SetSpacing(1.0, 1.0, 1.0)
            self.sg = vtk.vtkStructuredPointsWriter()
            # Output folder
            if not os.path.exists("./VTKResults"):
                os.makedirs("./VTKResults")
        elif self.outputF=='screen':
            if self.saveSnapshot==True:
                if not os.path.exists("./Snapshots"):
                    os.makedirs("./Snapshots")                
        elif self.outputF=='line':
            if not os.path.exists("./Lines"):
                os.makedirs("./Lines")
                
        # Restart folder
        if not os.path.exists("./Restarts"):
            os.makedirs("./Restarts") 
                                                           
#case = Interface()
#case.readInputVars()         
