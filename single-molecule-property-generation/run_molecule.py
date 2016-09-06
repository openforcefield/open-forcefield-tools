#!/bin/env python

import time
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *
import numpy as np
from mdtraj.reporters import NetCDFReporter
from smarty import *

#Define what molecule to work on, and a few simulation parameters
molname = 'AlkEthOH_r51'
mol_filename =  molname+'.mol2'
time_step = 2 #Femtoseconds
temperature = 300 #kelvin
friction = 1 # per picosecond
num_steps = 100000 
trj_freq = 1000 #steps
data_freq = 1000 #steps

# Load OEMol
mol = oechem.OEGraphMol()
ifs = oechem.oemolistream(mol_filename)
flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol )
oechem.OETriposAtomNames(mol)

# Get positions
coordinates = mol.GetCoords()
natoms = len(coordinates)
positions = np.zeros([natoms,3], np.float32)
for index in range(natoms):
    (x,y,z) = coordinates[index]
    positions[index,0] = x
    positions[index,1] = y
    positions[index,2] = z
positions = Quantity(positions, unit.angstroms)

# Load forcefield
forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

# Define system
topology = generateTopologyFromOEMol(mol)
system = forcefield.createSystem(topology, [mol])

# Choose SMIRKS string to apply k parameter change to
s = '[a,A:1]-[#6X4:2]-[a,A:3]'

# Alter forcefield and re-run simulation 
param = forcefield.getParameter(smirks=s)
param['k'] = str(10) # choosing to simulate adjacent value of parameter in manipulateparameters.py script
forcefield.setParameter(param, smirks = s)
system = forcefield.createSystem(topology, [mol])


#Do simulation
integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)
netcdf_reporter = NetCDFReporter('AlkEthOH_r51_k10.nc', trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(app.StateDataReporter('data_k10.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

print("Starting simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")
