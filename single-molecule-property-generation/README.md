Contains

Files:

*`manipulateparameters.py`:script for testing allowable accurate jumps in parameter space for reweighting. Carries out reweighting with MBAR from a 
                           prescribed state and compares reweighted observables (like average bond length or its variance) to values from full MD simulations. 
                           Results from this are less critical now that we are considering surrogate modeling as the major efficiency booster in the posterior 
                           sampling process. Still useful for getting accurate local approximations of observables in order to construct better surrogates for 
                           MCMC sampling. Still a work in progress, but we have some rough idea of how far in parameter space you can jump and still have 
                           accurate observable estimates (at least for bonds and angles). Can utilize multiple trajectories in order to reweight to unsampled 
                           states using the Python implementation of MBAR and calculate various single molecule observables including:   

  * Bond equilibrium lengths
  * Bond equilibrium standard deviations
  * Angle equilibrium angles
  * Angle equilibrium standard deviations
  * Torsion Fourier representations (0, 1, 2, 3, 6 + phase angle) (determined by least square fit) (still a WIP)

*`run_molecule.py`:script used to run simulations whose trajectories are utilized in the reweighting script

*`least_squares_fit_example.py`:Script being used to develop fits for torsion probability histograms. The fitting form is currently a fourier series 
                                approximation with sin^2 terms. Parameters for the fit are fourier coefficients A_i, periodicities P_i and phase angles phi_i. 
                                Still a work in progress, but fourier series fits seem like a promising avenue for describing torsion distributions. 
                                (still a WIP) 
 
*`construct_posterior.py`:Script for developing the posterior sampling process for bonds, angles and torsions parameterization. Still in early infancy and 
                          working on fitting data to surrogate models in the meantime. This mainly represents my efforts in writing an MCMC algorithm and 
                          convincing myself of the process outline. Integration with MBAR and MCMC with a surrogate model to come soon and will appear on 
                          my fork. (still a WIP)

Directories:

*`Mol2_files`:file of .mol2 files representing most of AlkEthOH. Used in `manipulateparameter.py` and `run_molecule.py`

*`mbar_analyses`:output reweighting statistics from `manipulateparameters.py` in the form of .pkl and .csv  

*`torsion_histograms`:output images of torsion probability distributions from `least_squares_fit_example.py`

*`traj`:long (50 ns) trajectory files being used in `least_squares_fit_example.py`

*`traj4ns`:4 ns trajectory files being used in `manipulateparameters.py`

*`traj20ns`:20 ns trajectory files being used in `manipulateparameters.py`


Author: Bryce Manubay
Affiliation: CU Boulder 
