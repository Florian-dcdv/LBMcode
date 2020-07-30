# LBMcode
These python codes can be used to simulate multiphase flows with Lattice Boltzmann Method (LBM). 
Each file contain different codes and input file for different cases:

- src: LBM code (pyLBM.py) with its classes.

- maxwell_construction: 
  - coexistence_curve.py: code to plot the coexistence curve
  - thermodynamic_consistency.py: code to plot graphs for a better understanding of the maxwell construction
  
- Poiseuille: input file to test the Poiseuille flow and a python code to compare simualtion results with analytical ones.

- Laplace: code to perform and to plot Laplace test and to test surface tension adjustment.

- Bouncing: 
  - graph_height_time.py: code to compare LBM results with Bertola results for Leidenfrost effect.
  - Dimensionless numbers.py: code to find lattice units to fit with real case with physical units.

- Spreading: code to study the spreading factor of a droplet during the kinematic phase.

- D2_law: code to test the D^2 law for heat transfer.
