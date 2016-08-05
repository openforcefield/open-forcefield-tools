import pdb
import trajectory_parser
import numpy as np
import pandas as pd


# don't actually pass in a pandas dataset yet: isn't accessed yet anyway

d = trajectory_parser.get_properties_from_trajectory(['AlkEthOH_c581.nc', 'AlkEthOH_r13.nc', 'AlkEthOH_r48.nc'])  

print(d)
	
