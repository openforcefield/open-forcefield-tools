#!/bin/env python

from tools import *
import pickle

# Load pre-generated data if already generated, otherwise re-generate
if os.path.isfile('ref_energies.pickle') and os.path.isfile('topologies.pickle') and os.path.isfile('oemols_molnames.pickle'):
    file = open('ref_energies.pickle', 'r')
    ref_energies = pickle.load(file)
    file.close()
    file = open('topologies.pickle', 'r')
    topologies = pickle.load(file)
    file = open('oemols_molnames.pickle', 'r')
    oemols, mol_names = pickle.load(file)
    file.close()
    print("Loaded reference data.")
else:
    print("Re-generating reference data.")
    commands.getoutput('python generate_reference_data.py')
    file = open('ref_energies.pickle', 'r')
    ref_energies = pickle.load(file)
    file.close()
    file = open('topologies.pickle', 'r')
    topologies = pickle.load(file)
    file = open('oemols_molnames.pickle', 'r')
    oemols, mol_names = pickle.load(file)
    file.close()
    print("Loaded reference data.")

if os.path.isfile('s_E_mean.pickle'):
    file = open('s_E_mean.pickle', 'r')
    s_E, mean_unc = pickle.load(file)
    file.close()
else:
    print("Re-generating uncertainty data.")
    commands.getoutput('python get_uncertainties.py')
    file = open('s_E_mean.pickle', 'r')
    s_E, mean_unc = pickle.load(file)
    file.close()

# Utility function to determine which parameters are in which molecules
def param_in_mols(labels, smirks):
    """Return list of True/False values as to whether specified SMIRKS pattern is in the molecules
    for which labels are provided. Labels should be as output by ForceField.labelMolecules"""
    smirks_in_mol = []
    for mol_entry in range(len(labels)):
        found = False
        for force in labels[mol_entry].keys():
            for (atom_indices, pid, s) in labels[mol_entry][force]:
                if s==smirks and not found:
                    smirks_in_mol.append(True)
                    found = True
        if not found:
            smirks_in_mol.append(False)

    return smirks_in_mol

def params_in_mols( labels, smirkslist ):
    """Return list of True/False values as to whether any of the specified SMIRKS patterns are in each of the specified molecules. Arguments are labels, as output by ForceField.labelMolecules, and smirkslist, list of SMIRKS of interest."""
    params_is_in_mols = []
    for smirks in smirkslist:
        # Handle case where this is the first smirks
        if len(params_is_in_mols)==0:
            param_is_in_mols == param_in_mols(labels, smirks)
        # If this is a subsequent smirks, just update where needed
        else:
            tmp = param_in_mols(labels, smirks)
            # Make updates
            for idx in range(len(tmp)):
                if tmp[idx]:
                    params_is_in_mols[idx] = True
    return params_is_in_mols

# Define a uniform prior
def uniform_prior( theta, prange, parameter_key):
    """Defines a uniform prior

    Parameters
    ----------
    theta : float
        Parameter value to consider
    prange : dict
        Dictionary, keyed by parameter type, defining the parameter range to consider.
    parameter_key : str
        key from the dictionary which applies to this parameter
    """
    bounds = prange[parameter_key]
    if theta < bounds[0]:
        return 0
    elif theta > bounds[1]:
        return 0
    else: return 1

# Initialization function
def initialize( Nmols, smirkslist, modify_key_by_smirks, ff, prange, oemols, measE, mean_unc):
    """Initialize for sampling.

    Parameters
    ----------
    Nmols : int
        Number of molecules to use in sampling
    smirkslist : list of strs
        SMIRKS strings of parameters we will sample
    modify_key_by_smirks : dict
        Dictionary, keyed by SMIRKS, of parameters to modify
    ff : ForceField
        ForceField object to start with
    prange : dict
        Dictionary, keyed by parameter name, of what parameter values to consider
    oemols : list of OEMols
        List of molecules we are studying/using or reference
    measE : list of floats
        List of measured (reference) energies for molecules considered
    mean_unc : float
        Uncertainty in reference energies
    startparms : str
        Choice of starting parameter values; use 'ff' to start with existing forcefield parameters; use 'uniform' to pick random starting point within range.

    Returns
    -------
    molidx_list : list of ints
        List of indexes of molecules we will use, of length Nmols
    moveff : ForceField
        ForceField object we will be using for sampling
    measE : list of floats
        Measured energies for reference molecules
    meas_s_E : list of floats
        Uncertainties in measured energies for reference molecules
    start_param_vals : dict
        Dictionary of starting parameter values, keyed by SMIRKS and parameter key
    """

    # Find molecules for which this parameter changes uncertainty - list of True false falues
    labels = ff.labelMolecules(oemols)
    params_is_in_mols = params_in_mols( labels, smirkslist )

    # Retrieve Nmols molecules containing at least one of those params
    molidx_list = []
    ct=0
    while len(molidx_list)< Nmols and ct < len(molidx_list):
        if params_is_in_mols[ct]:
            molidx_list.append(ct)
        ct+=1
    if ct==len(molidx_list):
        raise Exception("Error: Didn't find %s molecules containing that parameter." % Nmols)

    # Determine starting parameter values
    start_param_vals = {}
    for smirks in smirkslist:
        start_param_vals[smirks] = {}
        for param_key in modify_key_by_smirks[smirks]:
            thisrange = prange[param_key]
            if startparams=='uniform':
                start_param_vals[smirks][param_key] = random.uniform( thisrange[0], thisrange[1])
            elif startparams=='ff':
                start_param_vals[smirks][param_key] = ff.getParameter(smirks=smirks)[param_key]
            else:
                raise Exception("Error: Please choose starting parameter values.")

    # Retrieve measured and uncertainties
    measE = [ ref_energies[idx] for idx in molidx_list ]
    meas_s_E = np.array([mean_unc]*Nmols)

    # Copy forcefield
    moveff = copy.deepcopy(ff)

    # Initialize forcefield
    for smirks in smirkslist:
        params = moveff.getParameter(smirks=smirks)
    for param_key in modify_key_by_smirks[smirks]:
        params[param_key]=str(start_param_vals[smirks][param_key])
    moveff.setParameter(params, smirks=smirks)


    # Return stuff
    return molidx_list, moveff, measE, meas_s_E, start_param_vals


# THEN WRITE A FUNCTION FOR PARAMETER SAMPLING

def

# THEN TAKE USER INPUTS AND APPLY FUNCTIONS

