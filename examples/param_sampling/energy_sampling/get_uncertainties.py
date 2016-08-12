#!/bin/env python

"""
Determine uncertainties via conformer gneneration/minimization of multiple conformers.
"""

from tools import *
import pickle

# Load oemols
file = open('oemols_molnames.pickle', 'r')
oemols, mol_names = pickle.load(file)
file.close()
file = open('topologies.pickle', 'r')
topologies = pickle.load(file)
file.close()

# Input info
nmols = 25 #Number of multi-conformer molecules to examine
verbose = True

# Load forcefield file
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')
ff = ForceField(ffxml)

# Storage
s_E = np.zeros((nmols), np.float32)

# Configure Omega
maxConfs = 25
omega = oeomega.OEOmega()
omega.SetMaxConfs(maxConfs)
omega.SetStrictStereo(False) #Go ahead and generate even if stereocenter is unspecified
omega.SetSampleHydrogens(True)
omega.SetEnumNitrogen(oeomega.OENitrogenEnumeration_All)

# Pick molecules at random til we've processed at least nmols multi-conformer molecules
nprocessed = 0
used_idxs=[]
while nprocessed < nmols:
    found = False

    # Pick random molecule
    molidx = random.randint(0, len(oemols)-1)
    if molidx in used_idxs:
        continue

    if verbose: print("Attempting to process %s as molecule number %s..." % (molidx, nprocessed+1))

    # Generate conformers
    mol = OEMol(oemols[molidx])
    if not omega(mol):
        print("Conformer generation failed for molecule %s..." % molidx)
        continue

    if not mol.NumConfs() > 1:
        continue

    # Reformat coordinates of minima into format OpenMM will like
    cpositions=reformat_oemol_coordinates(mol)
    topology = topologies[molidx]
    # Create initial system
    system = ff.createSystem(topology, [mol], verbose=False)
    nconfs = mol.NumConfs()
    minimized_energies = np.zeros((nconfs), np.float32)
    for idx in range(nconfs):
        # To do: this constructs a simulation internally which is inefficient, should just change conformations
        en = get_minimum_energy(system, topology, cpositions[:,:,idx])
        # Store minimized energy
        minimized_energies[idx]=en
    if verbose:
        print("After minimization, conformer energies are (kcal/mol):", minimized_energies)
        print("After minimization, energy (of first conformer) is %.4g kcal/mol" % minimized_energies[0])

    # Compute uncertainty in average
    unc = minimized_energies.std()/np.sqrt(float(nconfs))

    # Store
    s_E[nprocessed] = unc
    if verbose: print("Uncertainty %.4g kcal/mol" % unc)

    nprocessed += 1

# Average
print(s_E)
mean_unc = s_E.mean()
print(mean_unc)

# Store
file = open('s_E_mean.pickle', 'w')
pickle.dump( (s_E, mean_unc), file)
file.close()

