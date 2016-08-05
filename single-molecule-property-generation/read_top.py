from parmed.amber import *
import numpy as np
import glob
import pandas as pd
from smarty import *
from smarty.forcefield_labeler import *

mol_files = glob.glob('./Mol2_files/AlkEthOH_*.mol2')
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

tops = []
for i in OEMols:
    top = generateTopologyFromOEMol(i)
    tops.append(top)

labeler = ForceField_labeler( get_data_filename('/data/forcefield/Frosst_AlkEtOH.ffxml') )


labels = []
lstt0 = []
lstt1 = []
lstt2 = []
lst00 = []
lst11 = []
lst22 = [] 
for ind, val in enumerate(OEMols):
    print('Molecule %s') % molnames[ind]
    label = labeler.labelMolecules([val], verbose = False) 
    for entry in range(len(label)):
        for bond in label[entry]['HarmonicBondForce']:
            lstt0.extend(bond[2])
	    lst00.append(bond[2])
	for angle in label[entry]['HarmonicAngleForce']:
	    lstt1.extend(angle[2])
	    lst11.append(angle[2])
	for torsion in label[entry]['PeriodicTorsionForce']:
	    lstt2.extend(torsion[2])
	    lst22.append(torsion[2])
	  # for (atom_indices, pid, smirks) in label[entry][force]:
		#print [atom_indices, pid, smirks]
	#	atomstr = ''
	#	for idx in atom_indices:
	#	    atomstr += '%6s' % idx 
		#rint("%s : %s \t smirks %s" % (atomstr, pid, smirks))
print lst00, lst11, lst22
files = glob.glob('./amb_tops/AlkEthOH_*.top')

def drop(mylist, m, n):
	mylist = list(mylist)
	del mylist[m::n]	
	return mylist

# Reading in and cleaning up atoms involved in bonds
lst0name = []
lstt0 = []
lst00 = []
print("PRINTING BOND PAIRS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	a1 = fin.parm_data['BONDS_INC_HYDROGEN']
	a2 = fin.parm_data['BONDS_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the bond length
	a1 = drop(a1,2,3)
	a2 = drop(a2,2,3)
	
	# Don't need to distinguish between bonded to H or not
	a1.extend(a2)
	
	# Return true atom numbers based on AMBER documentation
	a1 = np.array(a1)/3 + 1

	# Subdivide array into those of length 2 to make assigning column titles easier later
	#a2 = np.array_split(a1, len(a1)/2)

	# Need to create multiple lists for this to work
	# lst0name and lst0 allow me to keep the bond pairs indexed with the molecule
	# lst00 will allow me to create the column names after finding the unique pairs 
	lst0name.append(FileName)
	lstt0.append(a1)
	lst00.extend(a1)

# Convert lst00 into list of strings
lstt0 = [map(str,i) for i in lstt0]
lst00 = map(str, lst00)

print lstt0

# Join every two entries into space delimited string
lst0 = []
for sublst in lstt0:
	temp  = [i+' '+j for i,j in zip(sublst[::2], sublst[1::2])]
	lst0.append(temp)
lst00 = [i+' '+j for i,j in zip(lst00[::2], lst00[1::2])]

print lst00

# Return unique strings from lst00
cols0 = set()
for x in lst00:
	cols0.add(x)
cols0 = list(cols0)

# Generate data lists to populate dataframe
data0 = [[] for i in range(len(lst0))]
for val in cols0:
	for ind,item in enumerate(lst0):
		if val in item:
			data0[ind].append(1)
		else: 
			data0[ind].append(0)

# Reading in and cleaning up atoms involved in angles
lst1name = []
lstt1 = []
lst11 = []
print("PRINTING ANGLE TRIPLETS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	b1 = fin.parm_data['ANGLES_INC_HYDROGEN']
	b2 = fin.parm_data['ANGLES_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the angles
	b1 = drop(b1,3,4)
	b2 = drop(b2,3,4)
	
	# Don't need to distinguish between angles including H or not
	b1.extend(b2)
	
	# Return true atom numbers based on AMBER documentation
	b1 = np.array(b1)/3 + 1
	
	# Need to create multiple lists for this to work
	# lst1name and lst1 allow me to keep the angle trios indexed with the molecule
	# lst11 will allow me to create the column names after finding the unique trios 
	lst1name.append(FileName)
	lstt1.append(b1)
	lst11.extend(b1)

# Convert lstt1 and lst11 into list of strings
lstt1 = [map(str, i) for i in lstt1]
lst11 = map(str, lst11)

# Join every three entries into space delimited string
lst1 = []
for sublst in lstt1:
	temp  = [i+' '+j+' '+k for i,j,k in zip(sublst[::3], sublst[1::3], sublst[2::3])]
	lst1.append(temp)
lst11 = [i+' '+j+' '+k for i,j,k in zip(lst11[::3], lst11[1::3], lst11[2::3])]

# Return unique strings from lst11
cols1 = set()
for x in lst11:
	cols1.add(x)
cols1 = list(cols1)

# Generate data lists to populate frame (1 means val in lst1 was in cols1, 0 means it wasn't)
data1 = [[] for i in range(len(lst1))]
for val in cols1:
	for ind,item in enumerate(lst1):
		if val in item:
			data1[ind].append(1)
		else: 
			data1[ind].append(0)

# Reading in and cleaning up atoms involved in dihedrals
lstt2 = [] 
lst2name = []
lst22 = []
print("PRINTING DIHEDRAL QUARTETS...")
for FileName in files:
	# read in AMBER prmtop
	fin = AmberFormat(FileName)
 
	# pull out specified parm data
	c1 = fin.parm_data['DIHEDRALS_INC_HYDROGEN']
	c2 = fin.parm_data['DIHEDRALS_WITHOUT_HYDROGEN']
	
	# Get rid of the index identifier for the value of the torsions
	c1 = drop(c1,4,5)
	c2 = drop(c2,4,5)
	
	# Don't need to distinguish between torsions including H or not
	c1.extend(c2)
	
	# Return true atom numbers based on AMBER documentation
	for i in range(len(c1)):
		if c1[i] >= 0:
			c1[i] = np.array(c1[i])/3 + 1
		else:
			c1[i] = -(abs(np.array(c1[i]))/3 + 1)
				
	# Need to create multiple lists for this to work
	# lst2name and lst2 allow me to keep the torsion quartets indexed with the molecule
	# lst22 will allow me to create the column names after finding the unique quartets 
	lst2name.append(FileName)
	lstt2.append(c1)
	lst22.extend(c1)

# Convert lstt2 and lst22 into list of strings
lstt2 = [map(str,i) for i in lstt2]
lst22 = map(str, lst22)

# Join every four entries into space delimited string
lst2 = []
for sublst in lstt2:
	temp  = [i+' '+j+' '+k+' '+l for i,j,k,l in zip(sublst[::4], sublst[1::4], sublst[2::4], sublst[3::4])]
	lst2.append(temp)
lst22 = [i+' '+j+' '+k+' '+l for i,j,k,l in zip(lst22[::4], lst22[1::4], lst22[2::4], lst22[3::4])]

# Return unique strings from lst11
cols2 = set()
for x in lst22:
	cols2.add(x)
cols2 = list(cols2)

# Generate data lists to populate frame (1 means val in lst2 was in cols2, 0 means it wasn't)
data2 = [[] for i in range(len(lst2))]
for val in cols2:
	for ind,item in enumerate(lst2):
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
df0 = pd.DataFrame(data = data0, index = lst0name, columns = cols0)
df0['molecule'] = df0.index
df1 = pd.DataFrame(data = data1, index = lst1name, columns = cols1)
df1['molecule'] = df1.index
df2 = pd.DataFrame(data = data2, index = lst2name, columns = cols2)
df2['molecule'] = df2.index

dftemp = pd.merge(df0, df1, how = 'outer', on = 'molecule')
df = pd.merge(dftemp, df2, how = 'outer', on = 'molecule')

df['molecule'] = df.molecule.map(lambda x: x.replace(' ', '')[11:])
df['molecule'] = df.molecule.map(lambda x: x.replace(' ', '')[:-4])

print(df.molecule)

df.to_csv("check.csv")
