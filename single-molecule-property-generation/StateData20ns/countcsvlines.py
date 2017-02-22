import csv
import glob

fnames = glob.glob('data_AlkEthOH*')

for filename in fnames:
    with open(filename,'r') as f:
        reader = csv.reader(f,delimiter=';')
        data = list(reader)
        row_count = len(data)
        if row_count-1 < 10000:
            print row_count,filename 
     
