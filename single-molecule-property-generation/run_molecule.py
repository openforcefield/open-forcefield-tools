#!/bin/env python

import time
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *
import numpy as np
from mdtraj.reporters import NetCDFReporter
from smarty import *
import sys
import numpy as np

#Define what molecule to work on, and a few simulation parameters
molname = 'AlkEthOH_r48'
mol_filename = 'Mol2_files/'+molname+'.mol2'
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
positions = np.zeros([natoms,3], np.float64)
for index in range(natoms):
    (x,y,z) = coordinates[index]
    positions[index,0] = x
    positions[index,1] = y
    positions[index,2] = z
positions = Quantity(positions, unit.angstroms)

# Load forcefield
forcefield = ForceField(get_data_filename('forcefield/smirff99Frosst.ffxml'))

# Define system
topology = generateTopologyFromOEMol(mol)
params = forcefield.getParameter(smirks='[#1:1]-[#8]')
params['rmin_half']='0.01'
params['epsilon']='0.01'
forcefield.setParameter(params, smirks='[#1:1]-[#8]')
system = forcefield.createSystem(topology, [mol])

paramlist = np.arange(1107,1140,1)

param = forcefield.getParameter(smirks='[#8:1]-[#1:2]')
for i in paramlist:
    param['k'] = str(i)
    forcefield.setParameter(param, smirks='[#8:1]-[#1:2]')
    system = forcefield.createSystem(topology, [mol])


    #Do simulation
    integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
    platform = mm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature*kelvin)
    netcdf_reporter = NetCDFReporter('traj/'+molname+'_[#8:1]-[#1:2]_kbond'+str(i)+'.nc', trj_freq)
    simulation.reporters.append(netcdf_reporter)
    simulation.reporters.append(app.StateDataReporter('StateData/data_'+molname+'_[#8:1]-[#1:2]_kbond'+str(i)+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

    print("Starting simulation")
    start = time.clock()
    simulation.step(num_steps)
    end = time.clock()

    print("Elapsed time %.2f seconds" % (end-start))
    netcdf_reporter.close()
    print("Done!")

#Do simulation
#integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
#platform = mm.Platform.getPlatformByName('Reference')
#simulation = app.Simulation(topology, system, integrator)
#simulation.context.setPositions(positions)
#simulation.context.setVelocitiesToTemperature(temperature*kelvin)
#netcdf_reporter = NetCDFReporter('traj/AlkEthOH_r51.nc', trj_freq)
#simulation.reporters.append(netcdf_reporter)
#simulation.reporters.append(app.StateDataReporter('StateData/data_r51.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

#print("Starting simulation")
#start = time.clock()
#simulation.step(num_steps)
#end = time.clock()

#print("Elapsed time %.2f seconds" % (end-start))
#netcdf_reporter.close()
#print("Done!")
