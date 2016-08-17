# imports needed
from smarty import ForceField
import openeye
from openeye import oechem
import smarty
from smarty.utils import get_data_filename
from simtk import openmm
from simtk import unit
import numpy as np
import netCDF4 as netcdf
import collections as cl
import pandas as pd
import pymbar
from pymbar import timeseries

#######################################################################
# Constants
########################################################################

kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)

#######################################################################
# Utility Functions
#######################################################################

def read_col(filename,colname,frames):
    """Reads in columns from .csv outputs of OpenMM StateDataReporter 
    ARGUMENTS
	filename (string) - the path to the folder of the csv
	colname (string) - the column you wish to extract from the csv
	frames (integer) - the number of frames you wish to extract		
    """

    print "--Reading %s from %s/..." % (colname,filename)

    # Read in file output as pandas df
    df = pd.read_csv(filename, sep= ',')
	
    # Read values direct from column into numpy array
    dat = df.as_matrix(columns = colname)
    dat = dat[-frames:]


    return dat

def readtraj(ncfiles):

    """
    Take multiple .nc files and read in coordinates in order to re-valuate energies based on parameter changes

    ARGUMENTS
    ncfiles - a list of trajectories in netcdf format
    """

    for fname in ncfiles:
        data = netcdf.Dataset(fname)
        xyz = data.variables['coordinates']

    return data, xyz 


def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    return energy

#######################################################################

# Load simple OEMol
verbose = False
# Load one of the provided files
ifs = oechem.oemolistream(get_data_filename('molecules/AlkEthOH_r51.mol2'))
mol = oechem.OEMol()
# This uses parm@frosst atom types, so make sure to use the forcefield-flavor reader
flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol )
# Perceive tripos types
oechem.OETriposAtomNames(mol)

#Get positions for use below
data, xyz = readtraj(['traj/trajectory.nc'])
xyzn = unit.Quantity(xyz[:], unit.angstroms)

# Load forcefield file
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')
ff = ForceField(ffxml)

# Generate a topology
from smarty.forcefield import generateTopologyFromOEMol
topology = generateTopologyFromOEMol(mol)

# Create initial system
system = ff.createSystem(topology, [mol], verbose=verbose)

# Get params for an angle
params = ff.getParameter(smirks='[a,A:1]-[#6X4:2]-[a,A:3]')

# Make parameter moves getting new energy of system for each move
frames = 100 
k_vals = np.arange(0,600,10)
energies = np.zeros([len(k_vals),frames],np.float64)
for ind,val in enumerate(k_vals):
    temp = np.zeros(frames,np.float64)
    params['k'] = str(val)
    ff.setParameter(params, smirks = '[a,A:1]-[#6X4:2]-[a,A:3]')
    system = ff.createSystem(topology, [mol], verbose=verbose)
    for i,a in enumerate(xyzn):  
        e = np.float(get_energy(system, a))
        temp[i] = e
    energies[ind] = temp

# Post process energy distributions to find expectation values, analytical uncertainties and bootstrapped uncertainties
T_from_file = read_col('StateData/data.csv',["Temperature (K)"],100)
K = len(energies)
N_k = np.zeros(K, np.int32)
Temp_k = T_from_file
g = np.zeros(K, np.float64)
nBoots = 100

for k in range(K):
    g[k] = timeseries.statisticalInefficiency(energies[k])
    indices = np.array(timeseries.subsampleCorrelatedData(energies[k],g=g[k])) # indices of uncorrelated samples 
    N_k[k] = len(indices) # number of uncorrelated samples
    energies[k,0:N_k[k]] = energies[k,indices]

beta_k = 1 / (kB*Temp_k)


E_kn = np.zeros([K,frames],np.float64)

for k in range(K):
    E_kn[k,0:N_k[k]] = energies[k,0:N_k[k]]

#################################################################
# Compute reduced potentials
#################################################################

print "--Computing reduced potentials..."

u_kln = np.zeros([K,K,frames], np.float64)
E_kn_samp = np.zeros([K,frames], np.float64)

nBoots_work = nBoots + 1

allE_expect = np.zeros([K,nBoots_work], np.float64)
allE2_expect = np.zeros([K,nBoots_work], np.float64)
dE_expect = np.zeros([K], np.float64)

E_kn = energies

for n in range(nBoots_work):
    if (n > 0):
        print "Bootstrap: %d/%d" % (n,nBoots)
    for k in range(K):
        if N_k[k] > 0:
	    if (n == 0):
    		booti = np.array(range(N_k[k]))
	    else:
		booti = np.random.randint(N_k[k], size = N_k[k])
 	    E_kn_samp[k,0:N_k[k]] = E_kn[k,booti]            

    for k in range(K):
        for l in range(K):
            u_kln[k,l,0:N_k[k]] = beta_k[l] * E_kn_samp[k,0:N_k[k]] 

#########################################################################
# Initialize MBAR
############################################################################

# Initialize MBAR with Newton-Raphson
    if (n==0):  # only print this information the first time		
	print ""
	print "Initializing MBAR:"
	print "--K = number of parameter values with data = %d" % (K)
	print "--L = number of total parameter values tested = %d" % (K) 
	print "--N = number of Energies per parameter value = %d" % (np.max(N_k))

        # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
    if (n==0):
	initial_f_k = None # start from zero 
    else:
	initial_f_k = mbar.f_k # start from the previous final free energies to speed convergence
		
    mbar = pymbar.MBAR(u_kln, N_k, verbose=False, relative_tolerance=1e-12, initial_f_k=initial_f_k)

    #-----------------------------------------------------------------------
    # Compute Expectations for E_kt and E2_kt as E_expect and E2_expect
    #------------------------------------------------------------------------


    print ""
    print "Computing Expectations for E..."
    E_kln = u_kln  # not a copy, we are going to write over it, but we don't need it any more.
    for k in range(K):
        E_kln[:,k,:]*=beta_k[k]**(-1)  # get the 'unreduced' potential -- we can't take differences of reduced potentials because the beta is different.
    (E_expect, dE_expect) = mbar.computeExpectations(E_kln, state_dependent = True)
    allE_expect[:,n] = E_expect[:]
    # expectations for the differences, which we need for numerical derivatives  
    (DeltaE_expect, dDeltaE_expect) = mbar.computeExpectations(E_kln,output='differences', state_dependent = True)
    
    print "Computing Expectations for E^2..."
    (E2_expect, dE2_expect) = mbar.computeExpectations(E_kln**2, state_dependent = True)
    allE2_expect[:,n] = E2_expect[:]

if nBoots > 0:
    dE_boot = np.zeros([K])

    for k in range(K):
	dE_boot[k] = np.std(allE_expect[k,1:nBoots_work])
    print "E_expect: %s  dE_Expect: %s  dE_boot: %s" % (E_expect,dE_expect,dE_boot)
   
########################################################################

# Load forcefield file
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')
ff = ForceField(ffxml)

# Get a parameter by parameter id
param = ff.getParameter(paramID='b0001')
print(param)

# Get a parameter with a search restricted to a particular section, by smirks
param = ff.getParameter(smirks='[$([#1]-C):1]', force_type='NonbondedForce')
print(param)
