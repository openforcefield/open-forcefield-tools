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
molname = ['AlkEthOH_r0','AlkEthOH_r48','AlkEthOH_r51','AlkEthOH_c581','AlkEthOH_c100','AlkEthOH_c1161','AlkEthOH_c1266','AlkEthOH_c38','AlkEthOH_r118','AlkEthOH_r12']
mol_filename = ['Mol2_files/'+m+'.mol2' for m in molname]
time_step = 2 #Femtoseconds
temperature = 300 #kelvin
friction = 1 # per picosecond
num_steps = 100000 
trj_freq = 1000 #steps
data_freq = 1000 #steps

# Load OEMol
for ind,j in enumerate(mol_filename):
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(j)
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

    paramlist = np.arange(100,200,20)
    smirkseries = '[#6X4:1]-[#1:2]'
    paramtype = 'k'


    param = forcefield.getParameter(smirks=smirkseries)
    for i in paramlist:
        param[paramtype] = str(i)
        forcefield.setParameter(param, smirks=smirkseries)
        system = forcefield.createSystem(topology, [mol])


        #Do simulation
        integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(temperature*kelvin)
        netcdf_reporter = NetCDFReporter('traj/'+molname[ind]+'_'+smirkseries+'_'+paramtype+str(i)+'.nc', trj_freq)
        simulation.reporters.append(netcdf_reporter)
        simulation.reporters.append(app.StateDataReporter('StateData/data_'+molname[ind]+'_'+smirkseries+'_'+paramtype+str(i)+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

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
#netcdf_reporter = NetCDFReporter('traj/'+molname+'.nc', trj_freq)
#simulation.reporters.append(netcdf_reporter)
#simulation.reporters.append(app.StateDataReporter('StateData/data_'+molname+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

#print("Starting simulation")
#start = time.clock()
#simulation.step(num_steps)
#end = time.clock()

#print("Elapsed time %.2f seconds" % (end-start))
#netcdf_reporter.close()
#print("Done!")
