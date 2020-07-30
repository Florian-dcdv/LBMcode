# Flow Domain
domain={'nx': 200,'ny': 200, 'nz': 1}  

# Relaxation schemes: SRT, MRT, IMRT
#relaxationScheme={'scheme': 'SRT', 'tauv': 0.8}   
relaxationScheme={'scheme': 'IMRT', 'tauv': 0.6, 'sigma': 0.11}

# Solution options
solverOptions={'iteration': 35000, 'nFreq': 1000, 'restartOutputFreq': 35000, 'Restart': False}                   
#solverOptions={'iteration': 1, 'nFreq': 1}                   

# Fluid properties
fluidProperties={'FlowConfiguration': 'MultiPhase', 'liquidDensity': 0.455, 'liquidGas': 6.698e-4, 'newtonN': 1, 'cohesionStrength': -1}

# Body forces
bodyForces={'bodyForcex': 0, 'bodyForcey': 0}

# Parameters For Heat Equation
heatEquationParams={'solveTemperature': False} 
                                   
# Equation of state parameters
#EoSParams={'EoS': 'SC', 'latticeWeighting': True, 'refDensity': 1}
#EoSParams={'EoS': 'PR', 'criticalTemperature': 0.0729, 'temperatureRatio': 0.86, 'latticeWeighting': False, 'constAsc': 0.344, 'constA': 0.0408, 'constB': 0.09523, 'cZero': 1.}
EoSParams={'EoS': 'CS', 'criticalTemperature': 0.0377, 'temperatureRatio': 0.50, 'latticeWeighting': False, 'constA': 0.5, 'constB': 4, 'cZero': 1}

# Boundary conditions
flowBoundaryConditions={'Format': 'Standard', 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'Periodic', 'Top': 'Periodic'}

# Initial velocities
initialVelocities={'initVelocity': False}
#initialVelocities={'initVelocity': True, 'initVx': 0, 'initVy': 0}

# Initial densities
#initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 100, 'y0': 100, 'r0':50}
initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 100, 'y0': 100, 'r0':CHARRAD}

# Output format: VTK or screen
outputFormat={'format': 'line', 'xMin': 0, 'yMin': 100, 'xMax': 199, 'yMax': 100}
#outputFormat={'format': 'screen', 'saveSnapshot': False}
#outputFormat={'format': 'VTK'}
