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
#=====INPUTS===========================================================================================

bed_file="JT_MmSc_IgGi_0120.bed"  # input bed file
window=40 
read_n_th_line=1 # reading all data or n th line data, 1 means will read all lines
plot_hist=0 # plots histogram of input enrichment site widths
#======READING=============================================================================================
start_time = time.time()
rows = 6.5e6
skip = np.arange(rows)
skip = np.delete(skip, np.arange(0, rows, read_n_th_line))

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
num_bins=int(np.ceil(max_end/window))+1

bin_enrich=np.zeros([22,num_bins],dtype=int)
chr_i=b.Chr.to_numpy()

if plot_hist:
	b['diff']=b.End-b.Start
	new = b['diff'].copy()
	ax = new.hist(bins=window, alpha=0.8)
	print(new.min())
	plt.show()

#======com_data_for_binning==============================================================
#enrich=np.zeros([max_end,],dtype=int)

com_data=np.zeros([n_rows,5],dtype=int)
a=b.Start.to_numpy()
com_data[:,0]=b.Start.to_numpy()
com_data[:,1]=b.End.to_numpy()
del b

a_temp=np.floor(com_data[:,0]/window)
a_temp.astype(np.int)
com_data[:,2]=a_temp
a_temp=np.floor(com_data[:,1]/window)
a_temp.astype(np.int)
com_data[:,3]=a_temp
com_data[:,4]=chr_i
del a_temp
del chr_i

print("--- %s seconds ---" % (time.time() - start_time))
print("size_of_genome",max_end, num_bins)

#====BINNING_ENRICHMENTS=================================================================
for i in range(0,n_rows):#n_rows loop
        inds=com_data[i,:]
	st=inds[0]
        en=inds[1]
        bst=inds[2]
        ben=inds[3]
	chr_num=inds[4]-1

        if (ben-bst)>1:
                bin_enrich[chr_num,bst]+=(bst+1)*window-st+1
                bin_enrich[chr_num,ben]+=en-(ben)*window
                bin_enrich[chr_num,bst+1:ben]+=window
        elif (ben-bst)==1:
                mid=(bst+1)*window
		bin_enrich[chr_num,bst]+= mid-st+1
                bin_enrich[chr_num,ben]+=en-mid
        else:
		bin_enrich[chr_num,bst]+=1+en-st

print("--- %s seconds ---" % (time.time() - start_time))
#======END_Computation===============================================

plt.plot(bin_enrich[0,:])
plt.show()

