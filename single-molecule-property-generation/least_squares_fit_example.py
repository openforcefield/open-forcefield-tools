# Ben Coscia
# Homework 10
# Part II Question 3
import matplotlib as mpl

mpl.use('Agg')

import math
import scipy.optimize as sci
import matplotlib.pyplot as plt
import numpy as np
import openeye
from openeye import oechem
import glob 
from smarty import *
from smarty.utils import *
from smarty.forcefield import *
from simtk import unit
import pandas as pd
import netCDF4 as netcdf
import sys


# Given Data

I = [50, 80, 130, 200, 250, 350, 450, 550, 700]
P = [99, 177, 202, 248, 229, 219, 173, 142, 72]

# Fit the given data to the mathematical model:
# P = Pm*(I/Isat)*exp(-(I/Isat) + 1)   where
#   P == Photosynthesis Rate (mg / m^3 d)
#   Pm == Maximum Photosynthesis Rate (mg / m^3 d)
#   I == Solar Radiation (mu E /m^4 s)
#   Isat == Optimal Solar Radiation (mu E /m^4 s)

# The simplest way to do it is to minimize the sum of the residuals so that is what I am going to do
# The residual at each point is calculated as (y - f(a0, a1, x))**2
# for this specific case it will be (P - Pm*(I/Isat)*exp(-(I/Isat) + 1))**2

PmIsat = [100, 10]  # initial guesses at Pm and Isat. Format: [Pm, Isat]

def chi(PmIsat):
    sum_chi = 0
    for i in range(0, len(I)):
        sum_chi += (P[i] - PmIsat[0]*(I[i]/PmIsat[1])*math.exp(-(I[i]/PmIsat[1]) + 1))**2
    return sum_chi

A = sci.minimize(chi, PmIsat)


def goodness_of_fit(P, I, A):
    St = 0
    Sr = 0
    mean = np.mean(P)
    for i in range(0, len(P)):
        Sr += (P[i] - A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1))**2
        St += (P[i] - mean)**2
    s = math.sqrt(Sr/(len(P) - len(A.x)))
    R_squared = 1 - (Sr/St)
    return float(R_squared), s

def plot_fit(A, I):
    y = np.zeros((len(I)))
    for i in range(0, len(I)):
        y[i] = A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1)
    return y
#----------------------------------------------------------------------
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

    ff = ForceField(get_data_filename('/data/forcefield/Frosst_AlkEtOH.ffxml'))

    labels = []
    lst0 = []
    lst1 = []
    lst2 = []
    lst00 = [[] for i in molnames]
    lst11 = [[] for i in molnames]
    lst22 = [[] for i in molnames] 
    
    for ind, val in enumerate(OEMols):
        label = ff.labelMolecules([val], verbose = False) 
        for entry in range(len(label)):
            for bond in label[entry]['HarmonicBondGenerator']:
                lst0.extend([str(bond[0])])
	        lst00[ind].extend([str(bond[0])])
	    for angle in label[entry]['HarmonicAngleGenerator']:
	        lst1.extend([str(angle[0])])
	        lst11[ind].extend([str(angle[0])])
	    for torsion in label[entry]['PeriodicTorsionGenerator']:  
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

#------------------------------------------------------------------

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

#------------------------------------------------------------------
def get_small_mol_dict(mol2, traj):
    """
    Return dictionary specifying the bond, angle and torsion indices to feed to ComputeBondsAnglesTorsions()

    Parameters
    ----------
    mol2: mol2 file associated with molecule of interest used to determine atom labels
    traj: trajectory from the simulation ran on the given molecule
     
    Returns
    -------
    AtomDict: a dictionary of the bond, angle and torsion indices for the given molecule

    """
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
    for fname in traj:
        MoleculeName = fname.split('.')[0][:13]
        AtomDict['MolName'].append(MoleculeName)
         	
        
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

    return AtomDict
#------------------------------------------------------------------

def readtraj(ncfiles):

    """
    Take multiple .nc files and read in coordinates in order to re-valuate energies based on parameter changes

    ARGUMENTS
    ncfiles - a list of trajectories in netcdf format
    """

    for fname in ncfiles:
        data = netcdf.Dataset(fname)
        xyz = data.variables['coordinates']
        time = data.variables['time']

    return data, xyz, time 

#------------------------------------------------------------------
# Fit the given data to a 5 term fourier series:
# S_6(x) = (a0/2) + a1*np.sin(2*math.pi*x/P + phi1) + a2*np.sin(4*math.pi*x/P + phi2) + a3*np.sin(6*math.pi*x/P + phi3) + a6*np.sin(12*math.pi*x/P + phi6)
#   S_6(x)  == Torsion angle in radians
#   Pm == Maximum Photosynthesis Rate (mg / m^3 d)
#   I == Solar Radiation (mu E /m^4 s)
#   Isat == Optimal Solar Radiation (mu E /m^4 s)

# The simplest way to do it is to minimize the sum of the residuals so that is what I am going to do
# The residual at each point is calculated as (y - f(a0, a1, x))**2
# for this specific case it will be (P - Pm*(I/Isat)*exp(-(I/Isat) + 1))**2

PmIsat = [100, 10]  # initial guesses at Pm and Isat. Format: [Pm, Isat]

def chi(PmIsat):
    sum_chi = 0
    for i in range(0, len(I)):
        sum_chi += (P[i] - PmIsat[0]*(I[i]/PmIsat[1])*math.exp(-(I[i]/PmIsat[1]) + 1))**2
    return sum_chi

A = sci.minimize(chi, PmIsat)


def goodness_of_fit(P, I, A):
    St = 0
    Sr = 0
    mean = np.mean(P)
    for i in range(0, len(P)):
        Sr += (P[i] - A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1))**2
        St += (P[i] - mean)**2
    s = math.sqrt(Sr/(len(P) - len(A.x)))
    R_squared = 1 - (Sr/St)
    return float(R_squared), s

def plot_fit(A, I):
    y = np.zeros((len(I)))
    for i in range(0, len(I)):
        y[i] = A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1)
    return y

mol2= 'molecules/AlkEthOH_c581.mol2'
traj = ['AlkEthOH_c581_50ns.nc']

AtomDict = get_small_mol_dict(mol2, traj)

data, xyz, time = readtraj(traj)
xyzn = unit.Quantity(xyz[:], unit.angstroms)
time = unit.Quantity(time[:], unit.picoseconds)

# Compute bond lengths and angles and return array of angles
a = ComputeBondsAnglesTorsions(xyzn,AtomDict['Bond'],AtomDict['Angle'],AtomDict['Torsion'])

# Pull out torsion data
torsions = a[2]

# Get number of angles in molecule
numatom = len(torsions[0])

# Re-organize data into timeseries
torstimeser = [torsions[:,ind] for ind in range(numatom)]

# Using the angle at index 0 for test case
#torsion = np.zeros([1,100],np.float64)

torsion = torstimeser[0]

step = 0.05
bins = np.arange(np.around(np.min(torsion),decimals=1),np.around(np.max(torsion),decimals=1)+step,step)

binplace = np.digitize(torsion, bins)


likelihood = (np.bincount(binplace))/100.



R_squared, s = goodness_of_fit(P, I, A)
P_fit = plot_fit(A, I)
plt.plot(I, P)
plt.plot(I, P_fit)
plt.title('P versus I, R-squared = %s, s = %s' %(R_squared, s))
plt.xlabel('Solar Radiation (mu E /m^4 s)')
plt.ylabel('Photosynthesis Rate (mg / m^3 d)')
plt.show()

plt.savefig("some figure.png")

num_bins = 100

plt.figure()
plt.hist(torsion, num_bins)
plt.ylabel('Likelihood that configuration is sampled')
plt.xlabel('Torsion angle (radians)')
plt.savefig('Torsion_likelihood_c581.png')
