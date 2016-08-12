In this directory, I'm doing some exploration jumping off of SampleParameters_Energy in the directory above this one in the tree.

Manifest:
* tools.py: Utility functions used by other modules here.
* generate_reference_data.py: Generate gas phase reference data for a set of molecules -- compute energies after minimization with parm@frosst and store for later use.
* get_uncertainties.py: Use conformer expansion for a few molecules to get an idea of how much the energy varies across conformers so we can use this as the uncertainty in the observed energies.
* sample_parameters.py: MC sample parameters
