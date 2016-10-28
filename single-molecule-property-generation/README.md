Contains

* manipulateparameters.py: script for testing allowable accurate jumps in parameter space for reweighting. Can utilize multiple trajectories in order to reweight to unsampled states using the Python implementation of MBAR and calculate various single molecule observables including:   

  * Bond equilibrium lengths
  * Bond equilibrium standard deviations
  * Angle equilibrium angles
  * Angle equilibrium standard deviations
  * Torsion Fourier representations (0, 1, 2, 3, 6 + phase angle) (determined by least square fit) (still a WIP)

Simulation trajectories of the "unsampled" states are also used as input in order to assesss the accuracy of the reweighting. Current criteria of an accurate reweighted estimate are:

  * 1) Estimated observable is within 2 standard deviations away of the true sampled value
  * 2) Error in estimated uncertainties is less than 20%

* run_molecule.py: script used to run simulations whose trajectories are utilized in the reweighting script

* least_squares_fit_example.py: script being used to develop eventual fitting process for torsion representations 
 
* oe_license.txt: the license for using the OpenEye tools necessary for the reweighting and simulation processes

