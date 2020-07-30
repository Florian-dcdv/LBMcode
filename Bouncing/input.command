# Flow Domain
domain={'nx': 400,'ny': 500, 'nz': 1}  

# Relaxation schemes: SRT, MRT, IMRT
#relaxationScheme={'scheme': 'SRT', 'tauv': 0.8}   
relaxationScheme={'scheme': 'IMRT', 'tauv': 0.6, 'sigma': 0.11}

# Solution options
solverOptions={'iteration': 20000, 'nFreq': 500, 'restartOutputFreq': 20000, 'Restart': False}                   
#solverOptions={'iteration': 1, 'nFreq': 1}                   

# Fluid properties
fluidProperties={'FlowConfiguration': 'MultiPhase', 'liquidDensity': 0.455, 'liquidGas': 6.698e-4, 'newtonN': 1, 'cohesionStrength': -1}

# Body forces
bodyForces={'bodyForcex': 0, 'bodyForcey': -5.7524e-6}

# Parameters For Heat Equation
heatEquationParams={'solveTemperature': False} 
                                   
# Equation of state parameters
#EoSParams={'EoS': 'SC', 'latticeWeighting': True, 'refDensity': 1}
#EoSParams={'EoS': 'PR', 'criticalTemperature': 0.0729, 'temperatureRatio': 0.86, 'latticeWeighting': False, 'constAsc': 0.344, 'constA': 0.0408, 'constB': 0.09523, 'cZero': 1.}
EoSParams={'EoS': 'CS', 'criticalTemperature': 0.047, 'temperatureRatio': 0.50, 'latticeWeighting': False, 'constA': 0.5, 'constB': 4, 'cZero': 1}

# Boundary conditions
flowBoundaryConditions={'Format': 'Standard', 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'NoSlip', 'Top': 'NoSlip'}

# Initial velocities
initialVelocities={'initVelocity': False}
#initialVelocities={'initVelocity': True, 'initVx': 0, 'initVy': 0}

# Initial densities
initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 200, 'y0': 417, 'r0':24}

# Output format: VTK or screen
#outputFormat={'format': 'line', 'xMin': 0, 'yMin': 32, 'xMax': 63, 'yMax': 32}
#outputFormat={'format': 'screen', 'saveSnapshot': False}
outputFormat={'format': 'VTK'}
