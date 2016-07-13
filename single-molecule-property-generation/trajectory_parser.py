import pdb 
import os, sys, getopt
import numpy as np
import math
import netCDF4 as netcdf
import pandas as pd

def computeBondsAnglesTorsions(xyz, bonds, angles, torsions):
    """ 
    compute a 3 2D arrays of bond lengths for each frame: bond lengths in rows, angle lengths in columns

    inputs: the xyz files, an array of length-2 arrays.

    we calculate all three together since the torsions and angles
    require the bond vectors to be calculated anyway.

    """

    niterations = xyz.shape[0]
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
    
def calculateBondsAnglesTorsionsStatistics(properties, bond_dist, angle_dist, torsion_dist, bonds, angles, torsions):

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
        AtomList = np.array(p.split('Atom')[1:], dtype=int) # figure out which bond this is: 
                                                            # we assume bond_dist /bond is in the same order.
        if 'BondEquilibriumLength' in p:
            for i in range(nbonds):
                if np.array_equal(AtomList, bonds[i]): 
                    value = np.mean(bond_dist[:,i])
                    uncertainty = np.std(bond_dist[:,i])/np.sqrt(nsamp)
                    PropertyDict[p] = [value,uncertainty]

        if 'BondEquilibriumStd' in p:
            for i in range(nbonds):
                if np.array_equal(AtomList, bonds[i]): 
                    value = np.std(bond_dist[:,i])
                    uncertainty = np.std(bond_dist[:,i])**2/np.sqrt(nsamp/2)
                    PropertyDict[p] = [value,uncertainty]

    return PropertyDict


def get_properties_from_trajectory(dataframe, ncfiles):

    """take multiple .nc files with identifier names and a pandas dataframe with property 
    names for single a tom bonded properties (including the atom numbers) and populate 
    those property pandas dataframe.

    ARGUMENTS dataframe (pandas object) - name of the pandas object
       that contains the properties we want to extract.  ncfile
       (netcdf file) - a list of trajectories in netcdf format.  Names
       should correspond to the identifiers in the pandas dataframe.
    """

    PropertiesPerMolecule = dict()

    # here's code that generate list of properties to calculate for each molecule and 
    # populate PropertiesPerMolecule

    if 0:
        df = dataframe.set_index('molecule')
        MoleculeNames = df.molecule.tolist()
        properties = df.columns.values.tolist()

        for m in MoleculeNames:
            defined_properties  = list()
            for p in properties:
                if (p is not 'molecule') and ('_std' not in p):
                    if df.iloc[m][p] != np.nan:
                        defined_properties.append(p)
                    PropertiesPerMolecule[m] = defined_properties
    else:

        # hard coded properties
        PropertyNames = list()

        # bond properties
        PropertyNames.append('BondEquilibriumLengthAtom1Atom6')
        PropertyNames.append('BondEquilibriumLengthAtom2Atom6')
        PropertyNames.append('BondEquilibriumLengthAtom2Atom13')
        PropertyNames.append('BondEquilibriumLengthAtom2Atom14')

        PropertyNames.append('BondEquilibriumStdAtom1Atom6')
        PropertyNames.append('BondEquilibriumStdAtom2Atom6')
        PropertyNames.append('BondEquilibriumStdAtom2Atom13')
        PropertyNames.append('BondEquilibriumStdAtom2Atom14')

        #angle properties
        PropertyNames.append('AngleEquilibriumAngleAtom1Atom6Atom2')
        PropertyNames.append('AngleEquilibriumAngleAtom6Atom2Atom13')
        PropertyNames.append('AngleEquilibriumAngleAtom6Atom2Atom14')

        PropertyNames.append('AngleEquilibriumStdAtom1Atom6Atom2')
        PropertyNames.append('AngleEquilibriumStdAtom6Atom2Atom13')
        PropertyNames.append('AngleEquilibriumStdAtom6Atom2Atom14')
        
        # torsion properties
        PropertyNames.append('TorsionFourier1Atom1Atom6Atom2Atom13')
        PropertyNames.append('TorsionFourier1Atom1Atom6Atom2Atom14')

        PropertyNames.append('TorsionFourier2Atom1Atom6Atom2Atom13')
        PropertyNames.append('TorsionFourier2Atom1Atom6Atom2Atom14')

        PropertyNames.append('TorsionFourier3Atom1Atom6Atom2Atom13')
        PropertyNames.append('TorsionFourier3Atom1Atom6Atom2Atom14')

        PropertyNames.append('TorsionFourier6Atom1Atom6Atom2Atom13')
        PropertyNames.append('TorsionFourier6Atom1Atom6Atom2Atom14')

        PropertyNames.append('TorsionFourierPhaseAtom1Atom6Atom2Atom13')
        PropertyNames.append('TorsionFourierPhaseAtom1Atom6Atom2Atom14')

        PropertiesPerMolecule['AlkEthOH_c581'] = PropertyNames

    for fname in ncfiles:
        MoleculeName = fname.split('.')[0]

        # extract the xyz coordinate for each frame
        data = netcdf.Dataset(fname)
        xyz = data.variables['coordinates']

        # what is the property list for this molecule
        PropertyNames = PropertiesPerMolecule[MoleculeName]
        # extract the bond/angle/torsion lists
        AtomDict = dict() 
        AtomDict['Bond'] = list()
        AtomDict['Angle'] = list()
        AtomDict['Torsion'] = list()

        # which properties will we use to construct the bond list
        ReferenceProperties = ['BondEquilibriumLength','AngleEquilibriumAngle','TorsionFourier1']
        for p in PropertyNames:
            PropertyName = p.split('Atom')[0]
            AtomList = np.array(p.split('Atom')[1:],dtype=int)
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
                                                            AtomDict['Bond'], AtomDict['Angle'], AtomDict['Torsion'])

