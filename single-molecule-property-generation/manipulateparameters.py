# imports needed
from smarty import ForceField
import openeye
from openeye import oechem
import smarty
from smarty.forcefield_labeler import *
from smarty.utils import get_data_filename
from simtk import openmm
from simtk import unit
import numpy as np
import netCDF4 as netcdf
import collections as cl
import pandas as pd
import pymbar
from pymbar import timeseries
import glob
import sys

np.set_printoptions(threshold=np.inf)

#######################################################################
# Constants
########################################################################

kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)

#######################################################################
# Utility Functions
#######################################################################
def constructDataFrame(mol_files):
    """ 
    Construct a pandas dataframe to be populated with computed single molecule properties. Each unique bond, angle and torsion has it's own column for a value
    and uncertainty.
    inputs: a list of mol2 files from which we determine connectivity using OpenEye Tools and construct the dataframe using Pandas.
    """    
    
    molnames = []
    for i in mol_files:
        molname = i.replace(' ', '')[:-5]
        molname = molname.replace(' ' ,'')[13:]
        molnames.append(molname)

    OEMols=[]
    for i in mol_files:
        mol = oechem.OEGraphMol()
        ifs = oechem.oemolistream(i)
        flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
        ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
        oechem.OEReadMolecule(ifs, mol)
        oechem.OETriposAtomNames(mol)
        OEMols.append(mol)

    labeler = ForceField_labeler( get_data_filename('/data/forcefield/Frosst_AlkEtOH.ffxml') )

    labels = []
    lst0 = []
    lst1 = []
    lst2 = []
    lst00 = [[] for i in molnames]
    lst11 = [[] for i in molnames]
    lst22 = [[] for i in molnames] 
    
    for ind, val in enumerate(OEMols):
        label = labeler.labelMolecules([val], verbose = False) 
        for entry in range(len(label)):
            for bond in label[entry]['HarmonicBondForce']:
                lst0.extend([str(bond[0])])
	        lst00[ind].extend([str(bond[0])])
	    for angle in label[entry]['HarmonicAngleForce']:
	        lst1.extend([str(angle[0])])
	        lst11[ind].extend([str(angle[0])])
	    for torsion in label[entry]['PeriodicTorsionForce']:  
                lst2.extend([str(torsion[0])])
	        lst22[ind].extend([str(torsion[0])])

    # Return unique strings from lst0
    cols0 = set()
    for x in lst0:
	cols0.add(x)
    cols0 = list(cols0)


    # Generate data lists to populate dataframe
    data0 = [[] for i in range(len(lst00))]
    for val in cols0:
	for ind,item in enumerate(lst00):
	    if val in item:
		data0[ind].append(1)
	    else: 
		data0[ind].append(0)

    # Return unique strings from lst1
    cols1 = set()
    for x in lst1:
	cols1.add(x)
    cols1 = list(cols1)

    # Generate data lists to populate frame (1 means val in lst11 was in cols1, 0 means it wasn't)
    data1 = [[] for i in range(len(lst11))]
    for val in cols1:
	for ind,item in enumerate(lst11):
	    if val in item:
		data1[ind].append(1)
	    else: 
	        data1[ind].append(0)

    # Return unique strings from lst2
    cols2 = set()
    for x in lst2:
	cols2.add(x)
    cols2 = list(cols2)   
    
    # Generate data lists to populate frame (1 means val in lst22 was in cols2, 0 means it wasn't)
    data2 = [[] for i in range(len(lst22))]
    for val in cols2:
	for ind,item in enumerate(lst22):
	    if val in item:
		data2[ind].append(1)
	    else: 
		data2[ind].append(0)

    # Clean up clarity of column headers and molecule names
    cols0t = ["BondEquilibriumLength " + i for i in cols0]
    cols0temp = ["BondEquilibriumLength_std " + i for i in cols0]
    cols0 = cols0t + cols0temp

    cols1t = ["AngleEquilibriumAngle " + i for i in cols1]
    cols1temp = ["AngleEquilibriumAngle_std " + i for i in cols1]
    cols1 = cols1t + cols1temp

    cols2t = ["TorsionFourier1 " + i for i in cols2]
    cols2temp = ["TorsionFourier1_std " + i for i in cols2]
    cols2 = cols2t + cols2temp

    data0 = [i+i for i in data0]
    data1 = [i+i for i in data1]
    data2 = [i+i for i in data2]

    # Construct dataframes
    df0 = pd.DataFrame(data = data0, index = molnames, columns = cols0)
    df0['molecule'] = df0.index
    df1 = pd.DataFrame(data = data1, index = molnames, columns = cols1)
    df1['molecule'] = df1.index
    df2 = pd.DataFrame(data = data2, index = molnames, columns = cols2)
    df2['molecule'] = df2.index

    dftemp = pd.merge(df0, df1, how = 'outer', on = 'molecule')
    df = pd.merge(dftemp, df2, how = 'outer', on = 'molecule')

    return df


def ComputeBondsAnglesTorsions(xyz, bonds, angles, torsions):
    """ 
    compute a 3 2D arrays of bond lengths for each frame: bond lengths in rows, angle lengths in columns
    inputs: the xyz files, an array of length-2 arrays.
    we calculate all three together since the torsions and angles
    require the bond vectors to be calculated anyway.
    """

    niterations = xyz.shape[0] # no. of frames
    natoms = xyz.shape[1]

    nbonds = np.shape(bonds)[0]
    nangles = np.shape(angles)[0]
    ntorsions = np.shape(torsions)[0] 
    bond_dist = np.zeros([niterations,nbonds])
    angle_dist = np.zeros([niterations,nangles])
    torsion_dist = np.zeros([niterations,ntorsions])

    for n in range(niterations):
        xyzn = xyz[n] # coordinates this iteration
        bond_vectors = np.zeros([nbonds,3])
	for i, bond in enumerate(bonds):
	    bond_vectors[i,:] = xyzn[bond[0]-1] - xyzn[bond[1]-1]  # calculate the length of the vector
            bond_dist[n,i] = np.linalg.norm(bond_vectors[i]) # calculate the bond distance

        # we COULD reuse the bond vectors and avoid subtractions, but would involve a lot of bookkeeping
        # for now, just recalculate

        bond_vector1 = np.zeros(3)
        bond_vector2 = np.zeros(3)
        bond_vector3 = np.zeros(3)

        for i, angle in enumerate(angles):
            bond_vector1 = xyzn[angle[0]-1] - xyzn[angle[1]-1]  # calculate the length of the vector
            bond_vector2 = xyzn[angle[1]-1] - xyzn[angle[2]-1]  # calculate the length of the vector
            dot = np.dot(bond_vector1,bond_vector2)
            len1 = np.linalg.norm(bond_vector1)
            len2 = np.linalg.norm(bond_vector2)
            angle_dist[n,i] = np.arccos(dot/(len1*len2))  # angle in radians

        for i, torsion in enumerate(torsions):
            # algebra from http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates, Daniel's answer
            bond_vector1 = xyzn[torsion[0]-1] - xyzn[torsion[1]-1]  # calculate the length of the vector
            bond_vector2 = xyzn[torsion[1]-1] - xyzn[torsion[2]-1]  # calculate the length of the vector
            bond_vector3 = xyzn[torsion[2]-1] - xyzn[torsion[3]-1]  # calculate the length of the vector
            bond_vector1 /= np.linalg.norm(bond_vector1)
            bond_vector2 /= np.linalg.norm(bond_vector2)
            bond_vector3 /= np.linalg.norm(bond_vector3)
            n1 = np.cross(bond_vector1,bond_vector2)
            n2 = np.cross(bond_vector2,bond_vector3)
            m = np.cross(n1,bond_vector2)
            x = np.dot(n1,n2)
            y = np.dot(m,n2)
            torsion_dist[n,i] = np.arctan2(y,x)  # angle in radians

    return bond_dist, angle_dist, torsion_dist


def calculateBondsAnglesTorsionsStatistics(properties, bond_dist, angle_dist, torsion_dist, bonds, angles, torsions, torsionbool):

    """Inputs:
    properties: A list of property strings we want value for
    bond_dist: a Niterations x nbonds list of bond lengths
    angle_dist: a Niterations x nbonds list of angle angles (in radians)
    torsion_dist: a Niterations x nbonds list of dihedral angles (in radians)
    bonds: a list of bonds (ntorsions x 2)
    angles: a list of angles (ntorsions x 3)
    torsions: a list of torsion atoms (ntorsions x 4)

    # we assume the bond_dist / bonds , angle_dist / angles, torsion_dist / torsion were constucted in the same order.
    """
    PropertyDict = dict()
    nbonds = np.shape(bonds)[0]
    nangles = np.shape(angles)[0]
    ntorsions = np.shape(torsions)[0]
    
    nsamp = np.shape(bond_dist)[0]-1 #WARNING: assumes data points uncorrelated!
    for p in properties:        
        AtomList = p.split(' ', 1)[1:]  # figure out which bond this is: 
	AtomList = [i.lstrip('[').rstrip(']') for i in AtomList]  # we assume bond_dist /bond is in the same order.
	for i in AtomList:
            AtomList = i.strip().split(',')
        AtomList = map(int, AtomList) 

        if 'BondEquilibriumLength' in p:
            for i in range(nbonds):
                if np.array_equal(AtomList, bonds[i]): 
                    value = np.mean(bond_dist[:,i])
                    uncertainty = np.std(bond_dist[:,i])/np.sqrt(nsamp)
                    PropertyDict[p] = [value,uncertainty]

        if 'BondEquilibriumLength_std' in p:
            for i in range(nbonds):
        	if np.array_equal(AtomList, bonds[i]): 
                    value = np.std(bond_dist[:,i])
                    uncertainty = np.std(bond_dist[:,i])**2/np.sqrt(nsamp/2)
                    PropertyDict[p] = [value,uncertainty]

	if 'AngleEquilibriumAngle' in p:
       	    for i in range(nangles):
                if np.array_equal(AtomList, angles[i]): 
                    value = np.mean(angle_dist[:,i])
                    uncertainty = np.std(angle_dist[:,i])/np.sqrt(nsamp)
                    PropertyDict[p] = [value,uncertainty]

        if torsionbool==True:
	    if 'TorsionFourier1' in p:
                for i in range(ntorsions):
                    if np.array_equal(AtomList, torsions[i]): 
                    	value = np.mean(torsion_dist[:,i])
                    	uncertainty = np.std(torsion_dist[:,i])/np.sqrt(nsamp)
                    	PropertyDict[p] = [value,uncertainty]

	    if 'TorsionFourier1_std' in p:
	    	    for i in range(ntorsions):
	                if np.array_equal(AtomList, torsions[i]):
	            	    value = np.std(torsion_dist[:,i])
		    	    uncertainty = np.std(torsion_dist[:,i])**2/np.sqrt(nsamp/2)
		    	    PropertyDict[p] = [value,uncertainty]

	# Circular distribution alternate for torsion calculation
        
	    if 'TorsionFourier1' in p:
		for i in range(ntorsions):
		    if np.array_equal(AtomList, torsions[i]):
		        value = np.array([])
			for j in range(nsamp):
			    val = np.real((np.exp(cmath.sqrt(-1)*torsion_dist[:,i]))**j)
			    value = np.append(value, val)
			    value = (1/nsamp)*np.sum(value)
			    uncertainty = np.std(torsion_dist[:,i])/np.sqrt(nsamp)
			    PropertyDict[p] = [value, uncertainty]

	    if 'TorsionFourier1_std' in p:
		for i in range(ntorsions):
                    if np.array_equal(AtomList, torsions[i]):
                        value = np.std(torsion_dist[:,i])
                        uncertainty = np.std(torsion_dist[:,i])**2/np.sqrt(nsamp/2)
                        PropertyDict[p] = [value,uncertainty]
	else:
	    pass
                 
    return PropertyDict


def get_properties_from_trajectory(ncfiles, torsionbool=True):

    """take multiple .nc files with identifier names and a pandas dataframe with property 
    names for single atom bonded properties (including the atom numbers) and populate 
    those property pandas dataframe.
    ARGUMENTS dataframe (pandas object) - name of the pandas object
       that contains the properties we want to extract.  ncfile
       (netcdf file) - a list of trajectories in netcdf format.  Names
       should correspond to the identifiers in the pandas dataframe.
    """

    PropertiesPerMolecule = dict()

    # here's code that generate list of properties to calculate for each molecule and 
    # populate PropertiesPerMolecule
     
    mol_files = glob.glob('./Mol2_files/AlkEthOH_*.mol2')
 
    df = constructDataFrame(mol_files)
    MoleculeNames = df.molecule.tolist()
    properties = df.columns.values.tolist()
 
    for ind, val in enumerate(MoleculeNames):
        defined_properties  = list()
        for p in properties:
            if (p is not 'molecule') and ('_std' not in p):
                if df.iloc[ind][p] != 0:
		    defined_properties.append(p)
                PropertiesPerMolecule[val] = defined_properties

   
    AtomDict = dict()
    AtomDict['MolName'] = list()
    for fname in ncfiles:
        MoleculeName = fname.split('.')[0]
        AtomDict['MolName'].append(MoleculeName)
         	
        # extract the xyz coordinate for each frame
     
	data = netcdf.Dataset(fname)
        xyz = data.variables['coordinates']
	

        # what is the property list for this molecule
        PropertyNames = PropertiesPerMolecule[MoleculeName]

	# extract the bond/angle/torsion lists
        AtomDict['Bond'] = list()
        AtomDict['Angle'] = list()
        AtomDict['Torsion'] = list()

        # which properties will we use to construct the bond list
        ReferenceProperties = ['BondEquilibriumLength','AngleEquilibriumAngle','TorsionFourier1']
	for p in PropertyNames:
            PropertyName = p.split(' ', 1)[0]
            AtomList = p.split(' ', 1)[1:]
	    AtomList = [i.lstrip('[').rstrip(']') for i in AtomList]
	    for i in AtomList:
                AtomList = i.strip().split(',')
            AtomList = map(int, AtomList) 
            if any(rp in p for rp in ReferenceProperties):
                if 'Bond' in p:
                    AtomDict['Bond'].append(AtomList)
                if 'Angle' in p:
                    AtomDict['Angle'].append(AtomList)
                if 'Torsion' in p:
                    AtomDict['Torsion'].append(AtomList)
         

        bond_dist, angle_dist, torsion_dist = computeBondsAnglesTorsions(xyz,
                                                                         AtomDict['Bond'],
                                                                         AtomDict['Angle'],
                                                                         AtomDict['Torsion'])
		

        Properties = calculateBondsAnglesTorsionsStatistics(PropertyNames,
                                                            bond_dist, angle_dist, torsion_dist,
                                                            AtomDict['Bond'], AtomDict['Angle'], AtomDict['Torsion'], torsionbool)

        #Put properties back in dataframe and return

    return [bond_dist, angle_dist, torsion_dist, Properties]


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
# Parameters
#######################################################################

K_k = np.array([600, 500, 450, 100, 10]) # angular force constants for sampled states
N_k = np.array([0, 0, 0, 100, 0]) # no. of samples from each state (it's okay if some are zero)
K_extra = np.array([550, 300, 250, 60, 0]) # unsampled force constants
observables = ['angle', 'potential energy']
ncfiles = ['AlkEthOH_r51.nc']

seed = None

seed = 0
np.random.seed(seed)

######################################################################
# System setup
######################################################################
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
data, xyz = readtraj(['AlkEthOH_r51.nc'])
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

#######################################################################
# Main
#######################################################################

# Determine number of simulations
K = np.size(N_k)
if np.shape(K_k) != np.shape(N_k): raise "K_k and N_k must have same dimensions"


# Determine max number of samples to be drawn from any state
N_max = np.max(N_k)

# Calculate energies and the sampled angle distribution
energies = np.zeros([len(K_k),N_max],np.float64)
for ind,val in enumerate(K_k):
    temp = np.zeros(N_max,np.float64)
    params['k'] = str(val)
    ff.setParameter(params, smirks = '[a,A:1]-[#6X4:2]-[a,A:3]')
    system = ff.createSystem(topology, [mol], verbose=verbose)
    for i,a in enumerate(xyzn):  
        e = np.float(get_energy(system, a))
        temp[i] = e
	#b = get_properties_from_trajectory(
    energies[ind] = temp

energiesnew = np.zeros([len(K_extra),N_max],np.float64)
for ind,val in enumerate(K_extra):
    temp = np.zeros(N_max,np.float64)
    params['k'] = str(val)
    ff.setParameter(params, smirks = '[a,A:1]-[#6X4:2]-[a,A:3]')
    system = ff.createSystem(topology, [mol], verbose=verbose)
    for i,a in enumerate(xyzn):  
        e = np.float(get_energy(system, a))
        temp[i] = e
	#b = get_properties_from_trajectory(
    energiesnew[ind] = temp



# here's code that generate list of properties to calculate for each molecule and 
# populate PropertiesPerMolecule
PropertiesPerMolecule = dict()    
mol_files = glob.glob('./Mol2_files/AlkEthOH_*.mol2')
 
df = constructDataFrame(mol_files)
MoleculeNames = df.molecule.tolist()
properties = df.columns.values.tolist()
 
for ind, val in enumerate(MoleculeNames):
    defined_properties  = list()
    for p in properties:
        if (p is not 'molecule') and ('_std' not in p):
            if df.iloc[ind][p] != 0:
		defined_properties.append(p)
            PropertiesPerMolecule[val] = defined_properties

   
AtomDict = dict()
AtomDict['MolName'] = list()
for fname in ncfiles:
    MoleculeName = fname.split('.')[0]
    AtomDict['MolName'].append(MoleculeName)
         	
    # extract the xyz coordinate for each frame
    data = netcdf.Dataset(fname)
    xyz = data.variables['coordinates']
	

    # what is the property list for this molecule
    PropertyNames = PropertiesPerMolecule[MoleculeName]

    # extract the bond/angle/torsion lists
    AtomDict['Bond'] = list()
    AtomDict['Angle'] = list()
    AtomDict['Torsion'] = list()

    # which properties will we use to construct the bond list
    ReferenceProperties = ['BondEquilibriumLength','AngleEquilibriumAngle','TorsionFourier1']
    for p in PropertyNames:
        PropertyName = p.split(' ', 1)[0]
        AtomList = p.split(' ', 1)[1:]
        AtomList = [i.lstrip('[').rstrip(']') for i in AtomList]
	for i in AtomList:
            AtomList = i.strip().split(',')
        AtomList = map(int, AtomList) 
        if any(rp in p for rp in ReferenceProperties):
            if 'Bond' in p:
                AtomDict['Bond'].append(AtomList)
            if 'Angle' in p:
                AtomDict['Angle'].append(AtomList)
            if 'Torsion' in p:
                AtomDict['Torsion'].append(AtomList)

# Compute bond lengths and angles and return array of angles
a = ComputeBondsAnglesTorsions(xyzn,AtomDict['Bond'],AtomDict['Angle'],AtomDict['Torsion'])

angles = a[1] # for K_k = 100, angles in radians
numatom = len(a[1][1])

angtimeser = [angles[:,ind] for ind in range(numatom)]

# Post process energy distributions to find expectation values, analytical uncertainties and bootstrapped uncertainties
T_from_file = read_col('StateData/data.csv',["Temperature (K)"],100)
#K = len(energies)
Temp_k = T_from_file
nBoots = 100


g = np.array([timeseries.statisticalInefficiency(en) for en in energies])
N_k = np.array([len(timeseries.subsampleCorrelatedData(en,g=b)) for en, b in zip(energies,g)])
N = N_k.sum()
E_in = np.array([en[timeseries.subsampleCorrelatedData(en,g=b)] for en, b in zip(energies,g)])
E_n = np.zeros([N],np.float64)
for k in range(0,K): E_n[N_k[0:k].sum():N_k[0:k+1].sum()] = E_in[k]

gnew = np.array([timeseries.statisticalInefficiency(en) for en in energiesnew])
N_knew = np.array([len(timeseries.subsampleCorrelatedData(en,g=b)) for en, b in zip(energiesnew,gnew)])
Nnew = N_knew.sum()
E_innew = np.array([en[timeseries.subsampleCorrelatedData(en,g=b)] for en, b in zip(energiesnew,gnew)])
E_nnew = np.zeros([Nnew],np.float64)
for k in range(0,K): E_nnew[N_knew[0:k].sum():N_knew[0:k+1].sum()] = E_innew[k]

gang = np.array([timeseries.statisticalInefficiency(ang) for ang in angtimeser])
N_kang = np.array([len(timeseries.subsampleCorrelatedData(ang,g=b)) for ang, b in zip(angtimeser,gang)])
Nang = N_kang.sum()
Ang_in = np.array([ang[timeseries.subsampleCorrelatedData(ang,g=b)] for ang, b in zip(angtimeser,gang)])
Ang_n = np.zeros([Nang],np.float64)
for k in range(0,numatom): Ang_n[N_kang[0:k].sum():N_kang[0:k+1].sum()] = Ang_in[k]

beta_k = 1 / (kB*Temp_k)

E_kn = np.zeros([K,N_max],np.float64)
E_knnew = np.zeros([K,N_max],np.float64)
Ang_kn = np.zeros([numatom,N_max],np.float64)

for k,j in zip(range(K),range(numatom)):
    E_kn[k,0:N_k[k]] = E_in[k][0:N_k[k]]
    E_knnew[k,0:N_knew[k]] = E_innew[k][0:N_knew[k]]
    Ang_kn[j,0:N_knew[j]] = Ang_kn[j][0:N_knew[j]]

#################################################################
# Compute reduced potentials
#################################################################

print "--Computing reduced potentials..."

L = np.size(K_k)
Lnew = np.size(K_extra)

u_kln = np.zeros([K,L,N_max], np.float64)
E_kn_samp = np.zeros([K,N_max], np.float64)
u_klnnew = np.zeros([K,Lnew,N_max], np.float64)
E_kn_sampnew = np.zeros([K,N_max], np.float64)
Ang_kn_samp = np.zeros([numatom,N_max], np.float64)

nBoots_work = nBoots + 1

allE_expect = np.zeros([K,nBoots_work], np.float64)
allAng_expect = np.zeros([numatom,nBoots_work],np.float64)
allE2_expect = np.zeros([K,nBoots_work], np.float64)
dE_expect = np.zeros([K], np.float64)
allE_expectnew = np.zeros([K,nBoots_work], np.float64)
allE2_expectnew = np.zeros([K,nBoots_work], np.float64)
dE_expectnew = np.zeros([K], np.float64)
dAng_expect = np.zeros([numatom],np.float64)

E_kn = E_in
E_knnew = E_innew

for n in range(nBoots_work):
    if (n > 0):
        print "Bootstrap: %d/%d" % (n,nBoots)
    for k in range(K):
        if N_k[k] > 0:
	    if (n == 0):
    		booti = np.array(range(N_k[k]))
	    else:
		booti = np.random.randint(N_k[k], size = N_k[k])
 	    E_kn_samp[k,0:N_k[k]] = E_kn[k][booti]    

        if N_knew[k] > 0:
	    if (n ==0):
		bootnewi = np.array(range(N_knew[k]))
	    else:
		bootnewi = np.random.randint(N_knew[k], size = N_knew[k])        
            E_kn_sampnew[k,0:N_knew[k]] = E_knnew[k][bootnewi]
    for k in range(numatom):
        if N_kang[k] > 0:
            if (n == 0):
		bootangi = np.array(range(N_kang[k]))
	    else:
		bootangi = np.random.randint(N_kang[k],size=N_kang[k])
	    Ang_kn_samp[k,0:N_kang[k]] = Ang_kn[k][bootangi]
    for k in range(K):
        for l in range(L):
            u_kln[k,l,0:N_k[k]] = beta_k[l] * E_kn_samp[k,0:N_k[k]] 
	    u_klnnew[k,l,0:N_knew[k]] = beta_k[l] * E_kn_sampnew[k,0:N_knew[k]]      
   #for k in range(K):
       # for l in range(numatom):
           # ang_k[k,l,0:N_kang[k]] = Ang_kn_samp[k  
#########################################################################
# Initialize MBAR
############################################################################

# Initialize MBAR with Newton-Raphson
    if (n==0):  # only print this information the first time		
	print ""
	print "Initializing MBAR:"
	print "--K = number of parameter values with data = %d" % (K)
	print "--L = number of total parameter values tested = %d" % (L) 
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
    E_klnnew = u_klnnew
    for k in range(K):
        E_kln[:,k,:]*=beta_k[k]**(-1)  # get the 'unreduced' potential -- we can't take differences of reduced potentials because the beta is different.
	E_klnnew[:,k,:]*=beta_k[k]**(-1)
    (E_expect, dE_expect) = mbar.computeExpectations(E_kln, state_dependent = True)
    (E_expectnew, dE_expectnew) = mbar.computeExpectations(E_klnnew, state_dependent = True)
    (Ang_expect, dAng_expect) = mbar.computeExpectations(Ang_kn_samp,state_dependent=True)
   #    print Angles_expectnew
    allE_expect[:,n] = E_expect[:]
    allE_expectnew[:,n] = E_expectnew[:]
    allAng_expect[:,n] = Ang_expect[:]
    # expectations for the differences, which we need for numerical derivatives  
    (DeltaE_expect, dDeltaE_expect) = mbar.computeExpectations(E_kln,output='differences', state_dependent = True)
    (DeltaE_expectnew, dDeltaE_expectnew) = mbar.computeExpectations(E_klnnew,output='differences', state_dependent = True)

    #print "Computing Expectations for E^2..."
    #(E2_expect, dE2_expect) = mbar.computeExpectations(E_kln**2, state_dependent = True)
    #allE2_expect[:,n] = E2_expect[:]



if nBoots > 0:
    dE_boot = np.zeros([K])
    dE_bootnew = np.zeros([K])
    dAng_boot = np.zeros([numatom])
    for k,j in zip(range(K),range(numatom)):
	dE_boot[k] = np.std(allE_expect[k,1:nBoots_work])
        dE_bootnew[k] = np.std(allE_expectnew[k,1:nBoots_work])
        dAng_boot[j] = np.std(allAng_expect[j,1:nBoots_work])
    print "E_expect: %s  dE_Expect: %s  dE_boot: %s" % (E_expect,dE_expect,dE_boot)
    print "E_expectnew: %s  dE_Expectnew: %s  dE_bootnew: %s" % (E_expectnew,dE_expectnew,dE_bootnew)
    print "Ang_expect: %s  dAng_Expect: %s  dAng_boot: %s" % (Ang_expect,dAng_expect,dAng_boot)
sys.exit()

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
