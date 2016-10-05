# imports needed
from smarty.forcefield import *
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
import glob
import sys
from smarty.forcefield import generateTopologyFromOEMol
import pdb


#np.set_printoptions(threshold=np.inf)

#----------------------------------------------------------------------
# CONSTANTS
#----------------------------------------------------------------------

kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)

#----------------------------------------------------------------------
# UTILITY FUNCTIONS
#----------------------------------------------------------------------
def constructDataFrame(mol_files):
    """ 
    Construct a pandas dataframe to be populated with computed single molecule properties. Each unique bond, angle and torsion has it's own column for a value
    and uncertainty.
    
    Parameters
    -----------
    mol_files -  a list of mol2 files from which we determine connectivity using OpenEye Tools and construct the dataframe using Pandas.
    
    Returns
    -----------
    df - data frame in form molecules x property id that indicates if a specific property exists for a molecule (1 in cell if yes, 0 if no)
    """    
    
    molnames = []
    for i in mol_files:
        molname = i.replace(' ', '')[:-5]
        molname = molname.rsplit('/' ,1)[1]
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

    ff = ForceField(get_data_filename('/data/forcefield/smirff99Frosst.ffxml'))

    labels = []
    lst0 = []
    lst1 = []
    lst2 = []
    lst00 = [[] for i in molnames]
    lst11 = [[] for i in molnames]
    lst22 = [[] for i in molnames] 
    lst_0 = [[] for i in molnames]
    lst_1 = [[] for i in molnames]
    lst_2 = [[] for i in molnames] 
  

 
    for ind, val in enumerate(OEMols):
        label = ff.labelMolecules([val], verbose = False) 
        for entry in range(len(label)):
            for bond in label[entry]['HarmonicBondGenerator']:
                lst0.extend([str(bond[0])])
	        lst00[ind].extend([str(bond[0])])
                lst_0[ind].append([str(bond[0]),str(bond[2])])
	    for angle in label[entry]['HarmonicAngleGenerator']:
	        lst1.extend([str(angle[0])])
	        lst11[ind].extend([str(angle[0])])
                lst_1[ind].append((str(angle[0]),str(angle[2])))
	    for torsion in label[entry]['PeriodicTorsionGenerator']:  
                lst2.extend([str(torsion[0])])
	        lst22[ind].extend([str(torsion[0])])
                lst_2[ind].append([str(torsion[0]),str(torsion[2])])

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

    return df, lst_0, lst_1, lst_2

#------------------------------------------------------------------

def ComputeBondsAnglesTorsions(xyz, bonds, angles, torsions):
    """ 
    compute a 3 2D arrays of bond lengths for each frame: bond lengths in rows, angle lengths in columns
    
    Parameters 
    -----------
    xyz - xyz files, an array of length-2 arrays
    bonds, angles, torsions - numbered atom indices tuples associated with all unqiue bonds, angles and torsions
 
    Returns
    ----------
    bond_dist, angle_dist, torsion_dist - computed bonds, angles and torsions across the provided time series
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

def calculateBondsAnglesTorsionsStatistics(properties, bond_dist, angle_dist, torsion_dist, bonds, angles, torsions, torsionbool):

    """
    Parameters
    -----------
    properties: A list of property strings we want value for
    bond_dist: a Niterations x nbonds list of bond lengths
    angle_dist: a Niterations x nbonds list of angle angles (in radians)
    torsion_dist: a Niterations x nbonds list of dihedral angles (in radians)
    bonds: a list of bonds (ntorsions x 2)
    angles: a list of angles (ntorsions x 3)
    torsions: a list of torsion atoms (ntorsions x 4)
    torsionbool: boolean which suppresses torsion statistical analysis if False

    # we assume the bond_dist / bonds , angle_dist / angles, torsion_dist / torsion were constucted in the same order.

    PropertyDict - dictionary of average value of bond, angle or torsion across time series with associated uncertainty in mean and uncertainty in uncertainty
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

#------------------------------------------------------------------

def get_properties_from_trajectory(mol2, ncfiles, torsionbool=True):

    """
    take multiple .nc files with identifier names and a pandas dataframe with property 
    names for single atom bonded properties (including the atom numbers) and populate 
    those property pandas dataframe.
    
    Parameters
    -----------
    mol2 - mol2 files used to identify and index molecules  
    ncfiles -  a list of trajectories in netcdf format. Names should correspond to the identifiers in the pandas dataframe.
    torsionbool - boolean value passed to computeBondsAnglesTorsionsStatistics() to supress torsion statistics analysis. Default set to True (torsion calculatio                  n not supressed). 
    
    Returns
    ----------
    bond_dist - calculated bond distribution across trajectory
    angle_dist - calculated angle distribution across trajectory
    torsion_dist - calculated torsion distribution across trajectory
    Properties - dictionary of an average value of bond, angle or torsion across time series with associated uncertainty in mean and uncertainty in uncertainty
    """

    PropertiesPerMolecule = dict()

    # here's code that generate list of properties to calculate for each molecule and 
    # populate PropertiesPerMolecule
     
    mol_files = get_data_filename(mol2)
 
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

#------------------------------------------------------------------

def read_col(filename,colname,frames):
    """
    Reads in columns from .csv outputs of OpenMM StateDataReporter 
    
    Parameters
    -----------
    filename (string) - the path to the folder of the csv
    colname (string) - the column you wish to extract from the csv
    frames (integer) - the number of frames you wish to extract		
    
    Returns
    ----------
    dat - the pandas column series written as a matrix
    """

    print "--Reading %s from %s/..." % (colname,filename)

    # Read in file output as pandas df
    df = pd.read_csv(filename, sep= ',')
	
    # Read values direct from column into numpy array
    dat = df.as_matrix(columns = colname)
    dat = dat[-frames:]


    return dat

#------------------------------------------------------------------

def readtraj(ncfiles):
    """
    Take multiple .nc files and read in coordinates in order to re-valuate energies based on parameter changes

    Parameters
    ----------- 
    ncfiles - a list of trajectories in netcdf format

    Returns
    ----------
    data - all of the data contained in the netcdf file
    xyzn - the coordinates from the netcdf in angstroms
    """
    
    data = netcdf.Dataset(ncfiles)
    xyz = data.variables['coordinates']
    xyzn = unit.Quantity(xyz[:], unit.angstroms)
   
    return data, xyzn 

#------------------------------------------------------------------

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

#------------------------------------------------------------------

def new_param_energy(mol2, traj, smirkss, N_k, params, paramtype, samps, *coords):
    """
    Return potential energies associated with specified parameter perturbations.

    Parameters
    ----------
    mol2: mol2 file associated with molecule of interest used to construct OEMol object
    traj: trajectory from the simulation ran on the given molecule
    smirkss: list of smirks strings we wish to apply parameter changes to (Only changing 1 type of string at a time now. All bonds, all angles or all torsions)
    N_k: numpy array of number of samples per state
    params: a numpy array of the parameter values we wish to test
    paramtype: the type of ff param being edited
        BONDS - k (bond force constant), length (equilibrium bond length) 
        ANGLES - k (angle force constant), angle (equilibrium bond angle)
        TORSIONS - k{i} (torsion force constant), idivf{i} (torsional barrier multiplier), periodicity{i} (periodicity of the torsional barrier), phase{i} 
                   (phase offset of the torsion)
        NONBONDED - epsilon and rmin_half (where epsilon is the LJ parameter epsilon and rmin_half is half of the LJ parameter rmin)
    samps: samples per energy calculation

    Returns
    -------
    energies: a list of the energies associated with the forcfield parameters used as input

    """
    #-------------------
    # PARAMETERS
    #-------------------
    params = params
    N_k = N_k
    ncfiles = traj
    

    # Determine number of simulations
    K = np.size(N_k)
    #if np.shape(params) != np.shape(N_k): raise "K_k and N_k must have same dimensions"


    # Determine max number of samples to be drawn from any state

    #-------------
    # SYSTEM SETUP
    #-------------
    verbose = False # suppress echos from OEtoolkit functions
    ifs = oechem.oemolistream(get_data_filename(mol2))
    mol = oechem.OEMol()
    # This uses parm@frosst atom types, so make sure to use the forcefield-flavor reader
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    # Perceive tripos types
    oechem.OETriposAtomNames(mol)

    # Get positions for use below
    if not coords:
        data, xyz = readtraj(traj)
        xyzn = unit.Quantity(xyz[:], unit.angstroms)
    else:
        xyzn = coords
        
    # Load forcefield file
    ffxml = get_data_filename('forcefield/smirff99Frosst.ffxml')
    ff = ForceField(ffxml)

    # Generate a topology
    from smarty.forcefield import generateTopologyFromOEMol
    topology = generateTopologyFromOEMol(mol)

    #-----------------
    # MAIN
    #-----------------

    # Calculate energies 
    
    energies = np.zeros([len(smirkss),len(params),samps],np.float64)
    for inds,s in enumerate(smirkss):
        temp0 = np.zeros([len(params),samps],np.float64)
        param = ff.getParameter(smirks=s)
        for ind,val in enumerate(params):
            for p in paramtype:
                temp1 = np.zeros(samps,np.float64)
                for a,b in zip(val,p):
                    param[b] = str(a)       
                ff.setParameter(param, smirks = s)
                system = ff.createSystem(topology, [mol], verbose=verbose)
                for i,a in enumerate(xyzn):
                    e = np.float(get_energy(system, a)) * 4.184 #(kcal to kJ)
                    energies[inds,ind,i] = e
   
    return energies, xyzn, system

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
    mol_files = [] 
    for i in mol2:
        temp = get_data_filename(i)
        mol_files.append(temp)
    df,lst_0,lst_1,lst_2 = constructDataFrame(mol_files)
    MoleculeNames = df.molecule.tolist()
    properties = df.columns.values.tolist()
    #print MoleculeNames 
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
        MoleculeName = fname.split('.')[0][5:]
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

    return AtomDict,lst_0,lst_1,lst_2

#------------------------------------------------------------------

def subsampletimeseries(timeser,xyzn,N_k):
    """
    Return a subsampled timeseries based on statistical inefficiency calculations.

    Parameters
    ----------
    timeser: the timeseries to be subsampled
    xyzn: the coordinates associated with each frame of the timeseries to be subsampled
    N_k: original # of samples in each timeseries
    
    Returns
    ---------
    N_k_sub: new number of samples per timeseries
    ts_sub: the subsampled timeseries
    xyz_sub: the subsampled configuration series
    """
    # Make a copy of the timeseries and make sure is numpy array of floats
    ts = timeser
    xyz = xyzn

    # initialize array of statistical inefficiencies
    g = np.zeros(len(ts),np.float64)    


    for i,t in enumerate(ts):
        if np.count_nonzero(t)==0:
	    g[i] = np.float(1.)
            print "WARNING FLAG"
        else:
            g[i] = timeseries.statisticalInefficiency(t)
     
    N_k_sub = np.array([len(timeseries.subsampleCorrelatedData(t,g=b)) for t, b in zip(ts,g)])
    ind = [timeseries.subsampleCorrelatedData(t,g=b) for t,b in zip(ts,g)]
    
    #xyz_sub = np.array([unit.Quantity(c[i], unit.angstroms) for c,i in zip(xyz,ind)])
    if (N_k_sub == N_k).all():
        ts_sub = ts
        xyz_sub = xyz
        print "No sub-sampling occurred"
    else:
        print "Sub-sampling..." 
        ts_sub = np.array([t[timeseries.subsampleCorrelatedData(t,g=b)] for t,b in zip(ts,g)])
        for c in xyz:
            xyz_sub = [c[timeseries.subsampleCorrelatedData(t,g=b)] for t,b in zip(ts,g)]
    return ts_sub, N_k_sub, xyz_sub, ind

#------------------------------------------------------------------

# MAIN

#-----------------------------------------------------------------
# PARAMETERS
#-----------------------------------------------------------------
mol2 = ['molecules/AlkEthOH_r48.mol2']
mol2en = 'molecules/AlkEthOH_r48.mol2'
traj = 'traj/AlkEthOH_r48.nc'
smirkss = ['[#8:1]-[#1:2]']

#N_k = np.array([100, 100, 100, 100, 100])
#N_k = np.array([100,100])
N_k= np.array([100])

K = np.size(N_k)
print K

N_max = np.max(N_k)

#K_k = np.array([[106], [104], [102], [100], [98]])
#K_k = np.array([[104.],[100.]])
K_k = np.array([[1106.]])

#K_extra = np.array([[96], [99], [103], [105], [108]]) # unsampled force constants
#K_extra = np.array([[110.],[98.]])
K_extra = np.array([[1050.]])

paramtype = [['k']]

# Calculate energies at various parameters of interest
energies, xyzn, system = new_param_energy(mol2en,traj, smirkss, N_k, K_k, paramtype, N_max)
energiesnew, xyznnew, systemnew = new_param_energy(mol2en, traj, smirkss, N_k, K_extra, paramtype, N_max)

# Return AtomDict needed to feed to ComputeBondsAnglesTorsions()
AtomDict,lst_0,lst_1,lst_2 = get_small_mol_dict(mol2, [traj])
#mylist = [i[1] for i in lst_0[0]]
#myset = set(mylist)
#for i in myset:
#    print "%s occurs %s times" %(i, mylist.count(i))
#print lst_0
#print AtomDict['Bond']
#sys.exit()

# Read in coordinate data 
# Working on functionalizing this whole process of organizing the single molecule property data
#trajs = ['traj/AlkEthOH_r51_k106.nc','traj/AlkEthOH_r51_k104.nc','traj/AlkEthOH_r51_k102.nc','traj/AlkEthOH_r51.nc','traj/AlkEthOH_r51_k98.nc']
#trajs = ['traj/AlkEthOH_r51_k104.nc','traj/AlkEthOH_r51.nc']
trajs = ['traj/AlkEthOH_r48.nc']#,'traj/AlkEthOH_r48_[#8:1]-[#1:2]_kbond1120.nc']
#trajstest = ['traj/AlkEthOH_r51_k110.nc','traj/AlkEthOH_r51_k98.nc']
trajstest = ['traj/AlkEthOH_r48_[#8:1]-[#1:2]_kbond1050.nc']#,'traj/AlkEthOH_r48_[#8:1]-[#1:2]_kbond1150.nc']

xyznsampled = [[] for i in trajs] 
A = np.zeros([K,N_max],np.float64)
for i,x in enumerate(trajs):
    coord = readtraj(x)[1]
    #coord = unit.Quantity(coord[:], unit.angstroms)
    xyznsampled[i] = coord   
    obs = ComputeBondsAnglesTorsions(coord,AtomDict['Bond'],AtomDict['Angle'],AtomDict['Torsion'])[0]# Compute angles and return array of angles
    numatom = len(obs[0]) # get number of unique angles in molecule
    timeser = [obs[:,ind] for ind in range(numatom)] # re-organize data into timeseries
    A[i] = timeser[11] # pull out single angle in molecule for test case



xyznnewtest = [[] for i in trajstest] 
Anewtest = np.zeros([K,N_max],np.float64)
for i,x in enumerate(trajstest):
    coordtest = readtraj(x)[1]
    #coord = unit.Quantity(coord[:], unit.angstroms)
    xyznnewtest[i] = coordtest   
    obstest = ComputeBondsAnglesTorsions(coordtest,AtomDict['Bond'],AtomDict['Angle'],AtomDict['Torsion'])[0]# Compute angles and return array of angles
    numatomtest = len(obstest[0]) # get number of unique angles in molecule
    timesertest = [obstest[:,ind] for ind in range(numatomtest)] # re-organize data into timeseries
    Anewtest[i] = timesertest[11] # pull out single angle in molecule for test case

#print anglesnewtest 
# Subsample timeseries and return new number of samples per state
A_sub, N_kA, xyzn_A_sub, indA  = subsampletimeseries(A, xyznsampled, N_k)
En_sub, N_kEn, xyzn_En_sub, indEn = subsampletimeseries(energies[0], xyznsampled, N_k) 
Ennew_sub, N_kEnnew, xyzn_Ennew_sub, indEnnew = subsampletimeseries(energiesnew[0], xyznsampled, N_k)
A_sub_test,N_kA_test,xyzn_A_test,indAtest = subsampletimeseries(Anewtest,xyznnewtest,N_k)


A_kn = np.zeros([sum(N_kA)],np.float64)
count = 0
for x in A_sub:
    for y in x:
        A_kn[count] = y
        count += 1
#--------------------------------------------------------------
# Re-evaluate potenitals at all subsampled coord and parameters
#--------------------------------------------------------------
verbose = False # suppress echos from OEtoolkit functions
ifs = oechem.oemolistream(get_data_filename(mol2en))
mol = oechem.OEMol()
# This uses parm@frosst atom types, so make sure to use the forcefield-flavor reader
flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol )
# Perceive tripos types
oechem.OETriposAtomNames(mol)

# Load forcefield file
ffxml = get_data_filename('forcefield/smirff99Frosst.ffxml')
ff = ForceField(ffxml)

# Generate a topology
from smarty.forcefield import generateTopologyFromOEMol
topology = generateTopologyFromOEMol(mol)


# Re-calculate energies     
E_kn = np.zeros([len(K_k),sum(N_kA)],np.float64)
for inds,s in enumerate(smirkss):
    param = ff.getParameter(smirks=s)
    for ind,val in enumerate(K_k):
        count = 0
        for p in paramtype:
            for a,b in zip(val,p):
                param[b] = str(a)
            ff.setParameter(param, smirks = s)
            system = ff.createSystem(topology, [mol], verbose=verbose)  
            while count < sum(N_kA):
                for k_ind, pos in enumerate(xyzn_A_sub):
                    for i,a in enumerate(pos):
                        e = np.float(get_energy(system, a)) * 4.184 #(kcal to kJ)
                        E_kn[ind,count] = e
                        count += 1 
                       
E_knnew = np.zeros([len(K_extra),sum(N_kA)],np.float64)
for inds,s in enumerate(smirkss):
    param = ff.getParameter(smirks=s)
    for ind,val in enumerate(K_extra):
        count = 0
        for p in paramtype:
            for a,b in zip(val,p):    
                param[b] = str(a)
            ff.setParameter(param, smirks = s)
            system = ff.createSystem(topology, [mol], verbose=verbose)  
            while count < sum(N_kA):
                for k_ind, pos in enumerate(xyzn_A_sub):
                    for i,a in enumerate(pos):
                        e = np.float(get_energy(system, a)) * 4.184 #(kcal to kJ)
                        E_knnew[ind,count] = e
                        count += 1

# Post process energy distributions to find expectation values, analytical uncertainties and bootstrapped uncertainties
T_from_file = read_col('StateData/data.csv',["Temperature (K)"],100)
Temp_k = T_from_file
T_av = np.average(Temp_k)

nBoots = 100

beta_k = 1 / (kB*T_av)
bbeta_k = 1 / (kB*Temp_k)

#################################################################
# Compute reduced potentials
#################################################################

print "--Computing reduced potentials..."

# Initialize matrices for u_kn/observables matrices and expected value/uncertainty matrices
u_kn = np.zeros([K, sum(N_kA)], dtype=np.float64)
E_kn_samp = np.zeros([K,sum(N_kA)],np.float64)
u_knnew = np.zeros([K,sum(N_kA)], np.float64)
E_knnew_samp = np.zeros([K,sum(N_kA)], np.float64)
A_kn_samp = np.zeros([sum(N_kA)],np.float64)
A2_kn = np.zeros([sum(N_kA)],np.float64)


nBoots_work = nBoots + 1

allE_expect = np.zeros([K,nBoots_work], np.float64)
allA_expect = np.zeros([K,nBoots_work],np.float64)
allE2_expect = np.zeros([K,nBoots_work], np.float64)
dE_expect = np.zeros([K], np.float64)
allE_expectnew = np.zeros([K,nBoots_work], np.float64)
allE2_expectnew = np.zeros([K,nBoots_work], np.float64)
dE_expectnew = np.zeros([K], np.float64)
dA_expect = np.zeros([K],np.float64)
dA_expect_unsamp = np.zeros([K],np.float64)
allA_expect_unsamp = np.zeros([K,nBoots_work],np.float64)


# Begin bootstrapping loop
for n in range(nBoots_work):
    if (n > 0):
        print "Bootstrap: %d/%d" % (n,nBoots)
    for k in range(K):        
        if N_kA[k] > 0:
	    if (n == 0):
                booti = np.array(range(N_kA[k]))
	    else:
	        booti = np.random.randint(N_kA[k], size = N_kA[k])
           
            E_kn_samp[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])] = E_kn[:,booti]
            E_knnew_samp[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])] = E_knnew[:,booti]
            A_kn_samp[sum(N_kA[0:k]):sum(N_kA[0:k+1])] = A_kn[booti] 
            
#        if N_kang[k] > 0:
#	    if (n ==0):
#		bootnewi = np.array(range(N_kang[k]))
#	    else:
#		bootnewi = np.random.randint(N_kang[k], size = N_kang[k])        
#            if (k==0):
#                E_knnew_samp[k,0:N_kang[k]] = E_knnew[k,bootnewi]
#	    else:
#                E_knnew_samp[k,N_kang[k]:N_kang[k+1]] = E_knnew[k,bootnewi]
#    for k in range(K):
#        if N_kang[k] > 0:
#            if (n == 0):
#		bootangi = np.array(range(sum(N_kang)))
#	    else:
#		bootangi = np.random.randint(sum(N_kang), size=sum(N_kang))
#	    Ang_kn_samp[0:sum(N_kang)] = Ang_kn[bootangi]
        
    for k in range(K): 
        u_kn[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])] = beta_k * E_kn_samp[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])]     
        u_knnew[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])] = beta_k * E_knnew_samp[:,sum(N_kA[0:k]):sum(N_kA[0:k+1])]
	A2_kn[sum(N_kA[0:k]):sum(N_kA[0:k+1])] = A_kn_samp[sum(N_kA[0:k]):sum(N_kA[0:k+1])]
    print repr(u_kn)
############################################################################
# Initialize MBAR
############################################################################

# Initialize MBAR with Newton-Raphson
    if (n==0):  # only print this information the first time		
	print ""
	print "Initializing MBAR:"
	print "--K = number of parameter values with data = %d" % (K)
	print "--L = number of unsampled parameter values tested = %d" % (len(K_extra)) 
	print "--N = number of Energies per parameter value = %d" % (np.max(N_k))

        # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
#    if (n==0):
#	initial_f_k = None # start from zero 
#    else:
#	initial_f_k = mbar.f_k # start from the previous final free energies to speed convergence
		
    mbar = pymbar.MBAR(u_kn, N_kA, verbose=False, relative_tolerance=1e-12 )

    #------------------------------------------------------------------------
    # Compute Expectations for energy and angle distributions
    #------------------------------------------------------------------------
    if (n==40):
        pdb.set_trace()
    
    print ""
    print "Computing Expectations for E..."
    E_kn2 = u_kn  # not a copy, we are going to write over it, but we don't need it any more.
    E_knnew2 = u_knnew
    for k in range(K):
        E_kn2[k,:]*=beta_k**(-1)  # get the 'unreduced' potential -- we can't take differences of reduced potentials because the beta is different.
	E_knnew2[k,:]*=beta_k**(-1)
    (E_expect, dE_expect) = mbar.computeExpectations(E_kn2,state_dependent = True)
    (E_expectnew, dE_expectnew) = mbar.computeExpectations(E_knnew2,state_dependent = True)
    (A_expect, dA_expect) = mbar.computeExpectations(A2_kn,state_dependent = False) 


    allE_expect[:,n] = E_expect[:]
    allE_expectnew[:,n] = E_expectnew[:]
    allA_expect[:,n] = A_expect[:]
    
    # expectations for the differences, which we need for numerical derivatives  
    # To be used once the energy expectations are fixed
    (DeltaE_expect, dDeltaE_expect) = mbar.computeExpectations(E_kn2,output='differences', state_dependent = False)
    (DeltaE_expectnew, dDeltaE_expectnew) = mbar.computeExpectations(E_knnew2,output='differences', state_dependent = False)

    print "Computing Expectations for E^2..."
    (E2_expect, dE2_expect) = mbar.computeExpectations(E_kn2**2, state_dependent = True)
    allE2_expect[:,n] = E2_expect[:]

    (A_expect_unsamp, dA_expect_unsamp) = mbar.computeExpectations(A2_kn,u_knnew,state_dependent=False)
    allA_expect_unsamp[:,n] = A_expect_unsamp[:]
   
if nBoots > 0:
    dE_boot = np.zeros([K])
    dE_bootnew = np.zeros([K])
    dA_boot = np.zeros([K])
    dA_boot_unsamp = np.zeros([K])
    for k in range(K):
	dE_boot[k] = np.std(allE_expect[k,1:nBoots_work])
        dE_bootnew[k] = np.std(allE_expectnew[k,1:nBoots_work])
        dA_boot[k] = np.std(allA_expect[k,1:nBoots_work])
        dA_boot_unsamp[k] = np.std(allA_expect_unsamp[k,1:nBoots_work])
    print "E_expect: %s  dE_expect: %s  dE_boot: %s" % (E_expect,dE_expect,dE_boot)
    print "E_expectnew: %s  dE_expectnew: %s  dE_bootnew: %s" % (E_expectnew,dE_expectnew,dE_bootnew)
    print "delta_E_expect: %s  percent_delta_E_expect: %s" % (E_expectnew-E_expect, 100.*(E_expectnew-E_expect)/E_expect)
    print "A_expect: %s  dA_expect: %s  dA_boot: %s" % (A_expect,dA_expect,dA_boot)
    print "A_expect_unsamp: %s  dA_expect_unsamp: %s  dA_boot_unsamp: %s" % (A_expect_unsamp,dA_expect_unsamp,dA_boot_unsamp)
    print "The mean of the sampled series = %s" % ([np.average(A) for A in A_sub])
    print "The true sampled mean of the observable we're reweighting to = %s" % ([np.average(A) for A in A_sub_test])
    print "The mean of the energies corresponding to the sampled series = %s" % ([np.average(E) for E in E_kn]) 
    print "The mean of the energies corresponding to the unsampled series = %s" % ([np.average(E) for E in E_knnew])

    # calculate standard deviations away that estimate is from sampled value
    E_mean_samp = np.array([np.average(E) for E in E_kn])
    E_mean_unsamp = np.array([np.average(E) for E in E_knnew]) 
    A_mean_samp = np.array([np.average(A) for A in A_sub])
    A_mean_test = np.array([np.average(A) for A in A_sub_test])
    
    E_samp_stddevaway = np.zeros([K],np.float64)
    E_unsamp_stddevaway = np.zeros([K],np.float64)
    A_samp_stddevaway = np.zeros([K],np.float64)
    A_test_stddevaway = np.zeros([K],np.float64) 
    for k in range(K):  
        for i in allE_expect[k,1:]:
            E_samp_stddevaway[k] = np.array([np.average(np.abs(i-j)/u) for j,u in zip(E_mean_samp,dE_boot)])
        for i in allE_expectnew[k,1:]:
            E_unsamp_stddevaway[k] = np.array([np.average(np.abs(i-j)/u) for j,u in zip(E_mean_unsamp,dE_boot)])
        for i in allA_expect[k,1:]:
            A_samp_stddevaway[k] = np.array([np.average(np.abs(i-j)/u) for j,u in zip(A_mean_samp,dA_boot)])
        for i in allA_expect_unsamp[k,1:]:
            A_test_stddevaway[k] = np.array([np.average(np.abs(i-j)/u) for j,u in zip(A_mean_test,dA_boot_unsamp)])
    
    print "Standard deviations away from true sampled observables for E_expect: %s  E_expectnew: %s  A_expect: %s A_expect_unsamp: %s" % (E_samp_stddevaway,E_unsamp_stddevaway,A_samp_stddevaway,A_test_stddevaway) 
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
