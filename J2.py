#print("Hello")
import numpy as np
import os
import pandas as pd
import time
start_time = time.time()

bed_file="JT_MmSc_IgGi_0120.bed"

 
#a = np.loadtxt(bed_file, dtype=np.str, skiprows=1, usecols=(0,1,2))
b=pd.read_table(bed_file, delim_whitespace=True, skiprows=0, usecols=(0,1,2),engine='c',nrows=10)
b.columns = ['Chr', 'Start', 'End']


print(b.shape)
print("--- %s seconds ---" % (time.time() - start_time))

convert_dict = {'Chr': str, 
                'Start': int,
                'End' : int
               } 
b.sort_values(by=['Start'])
b = b.astype(convert_dict) 
print(b.dtypes)
print(b.Start)
 








