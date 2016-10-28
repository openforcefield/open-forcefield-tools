import pdb
import trajectory_parser
import numpy as np
import pandas as pd


np.set_printoptions(threshold=np.inf)

# don't actually pass in a pandas dataset yet: isn't accessed yet anyway

a = trajectory_parser.get_properties_from_trajectory(['AlkEthOH_r51.nc'],torsionbool=False)  


bonds = a[0]
angles = a[1]


for ind,val in enumerate(angles):
    print "Printing Angles of frame %s : %s " % (ind,val)	
