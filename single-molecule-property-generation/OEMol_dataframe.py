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

print molnames

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
lst0 = []
lst1 = []
lst2 = []
lst00 = [[] for i in molnames]
lst11 = [[] for i in molnames]
lst22 = [[] for i in molnames] 

for ind, val in enumerate(OEMols):
    print('Molecule %s') % molnames[ind]
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

#f['molecule'] = df.molecule.map(lambda x: x.replace(' ', '')[11:])
#f['molecule'] = df.molecule.map(lambda x: x.replace(' ', '')[:-4])

print df

df.to_csv("check.csv")
