import pdb
import trajectory_parser
import numpy as np
import pandas as pd


# don't actually pass in a pandas dataset yet: isn't accessed yet anyway
df = pd.read_csv('check.csv', sep = ',')

trajectory_parser.get_properties_from_trajectory(df,['AlkEthOH_c581.nc'])  



