# Flow Domain
domain={'nx': 100,'ny': 100, 'nz': 1}  

# Relaxation schemes: SRT, MRT, IMRT
relaxationScheme={'scheme': 'SRT', 'tauv': 0.6}   
#relaxationScheme={'scheme': 'IMRT', 'tauv': 0.6, 'sigma': 0.11}

# Solution options
solverOptions={'iteration': 130000, 'nFreq': 1000, 'restartOutputFreq': 130000, 'Restart': False}                   
#solverOptions={'iteration': 10, 'nFreq': 1,'restartOutputFreq': 10, 'Restart': False }                   

# Fluid properties
#fluidProperties={'FlowConfiguration': 'MultiPhase', 'liquidDensity': 0.431, 'liquidGas': 0.00149, 'newtonN': 1, 'cohesionStrength': -1}
fluidProperties={'FlowConfiguration': 'SinglePhase', 'liquidDensity': 1, 'newtonN': 1}

# Body forces
bodyForces={'bodyForcex': 1e-7, 'bodyForcey': 0}

# Parameters For Heat Equation
heatEquationParams={'solveTemperature': False} 
#heatEquationParams={'solveTemperature': True, 'SolverType': 'Implicit', 'specificHeat': 5, 'thermalConductivity': 0.666667, 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'Adiabatic', 'Top': 'Periodic', 'BottomValue': 0.047, 'initialTemperatures': {'initTemperature': False}}
                                   
# Equation of state parameters
EoSParams={'EoS': 'SC', 'latticeWeighting': False, 'refDensity': 1}
#EoSParams={'EoS': 'PR', 'criticalTemperature': 0.0729, 'temperatureRatio': 0.86, 'latticeWeighting': False, 'constAsc': 0.344, 'constA': 0.0408, 'constB': 0.09523, 'cZero': 1.}
#EoSParams={'EoS': 'CS', 'criticalTemperature': 0.047, 'temperatureRatio': 0.55, 'latticeWeighting': False, 'constA': 0.5, 'constB': 4, 'cZero': 1}

# Boundary conditions
flowBoundaryConditions={'Format': 'Standard', 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'NoSlip', 'Top': 'NoSlip', 'WallDensity': 1}

# Initial velocities
initialVelocities={'initVelocity': False}
#initialVelocities={'initVelocity': True, 'initVx': 0, 'initVy': 0}

# Initial densities
#initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 50, 'y0': 50, 'r0':20}
initialDensities={'initDensity': False}

# Output format: VTK or screen
#outputFormat={'format': 'line', 'xMin': 0, 'yMin': 32, 'xMax': 63, 'yMax': 32}
#outputFormat={'format': 'screen', 'saveSnapshot': False}
outputFormat={'format': 'VTK'}
