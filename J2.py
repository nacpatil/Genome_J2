from dict import chrdict
import numpy as np
import os
import pandas as pd
import time
import math
from sys import getsizeof
import matplotlib.pyplot as plt
import os
import psutil


start_time = time.time()

rows = 6.5e6
skip = np.arange(rows)
skip = np.delete(skip, np.arange(0, rows, 1))


bed_file="JT_MmSc_IgGi_0120.bed"
window=40

b=pd.read_table(bed_file, delim_whitespace=True, skiprows=skip, usecols=(0,1,2),engine='c')#   ,nrows=10)
b.columns = ['Chr', 'Start', 'End']
b.Chr=b['Chr'].map(chrdict)       # note: if the dictionary does not exhaustively map all

print("--- %s seconds ---" % (time.time() - start_time))

convert_dict = {'Chr': int, 
                'Start': int,
                'End' : int
               } 
b = b.astype(convert_dict) 


b.sort_values(by=['Chr','Start'])
n_rows=len(b.index)
max_end=b['End'].max()

num_bins=int(np.ceil(max_end/window))

bin_enrich=np.zeros([num_bins,],dtype=int)
enrich=np.zeros([max_end*22,1],dtype=int)


a=b.Start.to_numpy()
a1=b.End.to_numpy()
chr_i=b.Chr.to_numpy()

print("--- %s seconds ---" % (time.time() - start_time))




print(b.Chr.unique())
for i in range(0,n_rows):#n_rows loop
	pad=(chr_i[i]-1)*max_end
	st=a[i]+pad
	en=a1[i]+pad
	enrich[st:en,0]+=1
	

print("--- %s seconds ---" % (time.time() - start_time))

#plt.plot(enrich)
#plt.show()


