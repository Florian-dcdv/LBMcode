# Flow Domain
domain={'nx': 200,'ny': 200, 'nz': 1}  

# Relaxation schemes: SRT, MRT, IMRT
#relaxationScheme={'scheme': 'MRT', 'tauv': 0.6}   
relaxationScheme={'scheme': 'IMRT', 'tauv': 0.6, 'sigma': 0.11}

# Solution options
#solverOptions={'iteration': 300000, 'nFreq': 1000, 'restartOutputFreq': 1000, 'Restart': False}     
solverOptions={'iteration': 1000, 'nFreq': 100, 'restartOutputFreq': 500, 'Restart': False, 'RestartFilename': 't500.dat'}                   

# Fluid properties
fluidProperties={'FlowConfiguration': 'MultiPhase', 'liquidDensity': 0.431, 'liquidGas': 0.00149, 'newtonN': 1, 'cohesionStrength': -1}
#fluidProperties={'FlowConfiguration': 'SinglePhase', 'liquidDensity': 1, 'newtonN': 1}

# Body forces
bodyForces={'bodyForcex': 0, 'bodyForcey': 0}

# Parameters For Heat Equation
#heatEquationParams={'solveTemperature': False} 
#heatEquationParams={'solveTemperature': True, 'SolverType': 'Implicit', 'specificHeat': 6, 'thermalConductivity': 0.6666, 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'Adiabatic', 'Top': 'Periodic', 'LeftValue': 0.047, 'RightValue': 0.047, 'BottomValue': 0.047, 'TopValue': 0.047, 'initialTemperatures': {'initTemperature': False, 'shape': 'sphere', 'x0': 100, 'y0': 100, 'r0': 30, 'Value': 0.04042}}
heatEquationParams={'solveTemperature': True, 'SolverType': 'Implicit', 'specificHeat': 5, 'thermalConductivity': 0.666667, 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'Adiabatic', 'Top': 'Adiabatic', 'BottomValue': 0.047, 'TopValue': 0.02585, 'initialTemperatures': {'initTemperature': False}}
                                   
# Equation of state parameters
#EoSParams={'EoS': 'SC', 'latticeWeighting': False, 'refDensity': 1}
#EoSParams={'EoS': 'PR', 'criticalTemperature': 0.1094, 'temperatureRatio': 0.85, 'latticeWeighting': False, 'constAsc': 0.344, 'constA': 0.06122, 'constB': 0.09524, 'cZero': 1.}
EoSParams={'EoS': 'CS', 'criticalTemperature': 0.047, 'temperatureRatio': 0.55, 'latticeWeighting': False, 'constA': 0.5, 'constB': 4, 'cZero': 1}

# Boundary conditions
flowBoundaryConditions={'Format': 'Standard', 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'NoSlip', 'Top': 'NoSlip', 'WallDensity': 1.95e-2}

# Initial velocities
initialVelocities={'initVelocity': False}
#initialVelocities={'initVelocity': True, 'initVx': 0, 'initVy': 0}

# Initial densities
initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 100, 'y0': 80, 'r0':30}
#initialDensities={'initDensity': False}

# Output format: VTK or screen
#outputFormat={'format': 'line', 'xMin': 0, 'yMin': 32, 'xMax': 63, 'yMax': 32}
#outputFormat={'format': 'screen', 'saveSnapshot': False}
outputFormat={'format': 'VTK'}
