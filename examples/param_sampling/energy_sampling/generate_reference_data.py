#!/bin/env python

"""
Generate gas phase reference data for molecules from parm@frosst set.
This differes from SampleParameters_Energy.ipynb in that it uses minimization rather than just omega conformers (to ensure we get variation in bond lengths, angles, etc. that omega might not have).
"""

from tools import *
import pickle


# Specify where to pull data from
datapath = '/Users/dmobley/GoogleDrive/research/github/open-forcefield-data/Model-Systems/AlkEthOH_distrib'
datafiles = ['AlkEthOH_chain_filt1.oeb', 'AlkEthOH_rings_filt1.oeb']
maxmolnr = 1500 # Limit number of molecules for testing purposes; set to 1500 to process all

# Create storage for energies, other things. Also save actual oemols for reuse later
ref_energies = np.zeros((maxmolnr), float)
mols_found = 0
mol_names = []
oemols = []
topologies = [] #Save so we don't have to re-generate again and again later

# Load forcefield file
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')
ff = ForceField(ffxml)


# Infrastructure for topology generation
from smarty.forcefield import generateTopologyFromOEMol

# Process molecules, evaluate energies and store
mol = oechem.OEMol()
for ifile in datafiles:
    istream = oemolistream( os.path.join(datapath, ifile))
    print("Processing %s..." % ifile)
    while OEReadMolecule(istream, mol):
        if mols_found%10==0:
            print("   Molecule %d processed..." % mols_found)
        if mols_found >= maxmolnr:
            break
        title = mol.GetTitle()
        mol_names.append(title)


        # These OEB files don't bring in names, and we require atom names, so assign them.
        oechem.OETriposAtomNames(mol)
        oemols.append(OEMol(mol))

        # Pull coordinates
        cpositions=reformat_oemol_coordinates(mol)

        # Generate topology
        topology = generateTopologyFromOEMol(mol)
        topologies.append(topology)

        # Create initial system
        system = ff.createSystem(topology, [mol], verbose=False)

        # Get energy and store reference data
        ref_energies[mols_found] = get_minimum_energy(system, topology, cpositions[:,:,0])
        mols_found+=1

    istream.close()

    if mols_found >= maxmolnr:
        break


# Reduce storage to the number of molecules actually processed
ref_energies = ref_energies[0:mols_found]

print("Minimum and maximum energies found (kcal/mol) %.4g and %.4g" % (ref_energies.min(), ref_energies.max()))



# Store stuff to pickle files
# reference energies
file = open('ref_energies.pickle', 'w')
pickle.dump(ref_energies, file)
file.close()
# topologies
file = open('topologies.pickle', 'w')
pickle.dump(topologies,file)
file.close()

# oemols, mol names
file = open('oemols_molnames.pickle', 'w')
pickle.dump((oemols, mol_names), file)
file.close()

