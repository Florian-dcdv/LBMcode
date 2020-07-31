# Flow Domain
domain={'nx': 300,'ny': 200, 'nz': 1}  

# Relaxation schemes: SRT, MRT, IMRT
#relaxationScheme={'scheme': 'SRT', 'tauv': 0.8}   
relaxationScheme={'scheme': 'IMRT', 'tauv': 0.6, 'sigma': 0.11}

# Solution options
solverOptions={'iteration': 3000, 'nFreq': 10, 'restartOutputFreq': 30000, 'Restart': False}                   
#solverOptions={'iteration': 1, 'nFreq': 1}                   

# Fluid properties
fluidProperties={'FlowConfiguration': 'MultiPhase', 'liquidDensity': 0.431, 'liquidGas': 0.00149, 'newtonN': 1, 'cohesionStrength': -1}

# Body forces
bodyForces={'bodyForcex': 0, 'bodyForcey': -2.7722e-5}

# Parameters For Heat Equation
heatEquationParams={'solveTemperature': False} 
                                   
# Equation of state parameters
#EoSParams={'EoS': 'SC', 'latticeWeighting': True, 'refDensity': 1}
#EoSParams={'EoS': 'PR', 'criticalTemperature': 0.0729, 'temperatureRatio': 0.86, 'latticeWeighting': False, 'constAsc': 0.344, 'constA': 0.0408, 'constB': 0.09523, 'cZero': 1.}
EoSParams={'EoS': 'CS', 'criticalTemperature': 0.047, 'temperatureRatio': 0.55, 'latticeWeighting': False, 'constA': 0.5, 'constB': 4, 'cZero': 1}

# Boundary conditions
flowBoundaryConditions={'Format': 'Standard', 'Left': 'Periodic', 'Right': 'Periodic', 'Bottom': 'NoSlip', 'Top': 'NoSlip', 'WallDensity': 1.95e-2}

# Initial velocities
initialVelocities={'initVelocity': False}
#initialVelocities={'initVelocity': True, 'initVx': 0, 'initVy': 0}

# Initial densities
initialDensities={'initDensity': True, 'shape': 'sphere', 'x0': 150, 'y0': 80, 'r0':20}

# Output format: VTK or screen
#outputFormat={'format': 'line', 'xMin': 0, 'yMin': 32, 'xMax': 63, 'yMax': 32}
#outputFormat={'format': 'screen', 'saveSnapshot': False}
outputFormat={'format': 'VTK'}
