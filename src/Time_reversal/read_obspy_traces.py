# from obspy import read
import numpy as np
from math import sin, cos
from tqdm import tqdm
import h5py
import os as os
import obspy

f= open("earthquake.txt", "r")
earthquake= f.read()
f.close()

# a= os.listdir("/mnt/datawave3_new/eq_data/"+earthquake+"/")

# trace_idx= -1
# for i in range(len(a)):
#     fname= a[i]
#     if(fname[-2:]== ".a"):
#         trace_idx= i

# st= []
# if(trace_idx!=-1):
#     st= obspy.read("/mnt/datawave3_new/eq_data/"+earthquake+"/good_traces/*")
#     nr= len(st)

st= obspy.read("/mnt/datawave3_new/"+earthquake+"/*.a/good_traces/*")
nr= len(st)
    
for tr in st:
    tr.trim(tr.stats.starttime+1000+tr.stats.sac.a-200, tr.stats.starttime+1000+tr.stats.sac.a+200)
    tr.filter("highpass", freq=0.02, zerophase=True)
    d=tr.data
    d=d-np.mean(d)
    d=d/np.std(d)
    tr.data=d
    # tr.id=eq_name+' '+tr.id

data= np.empty((nr, len(st[1].data)))
for i in range(nr):
    data[i,:]= st[i].data

with h5py.File('/home/abhinav/workspace/Bp/processed_traces/'+ earthquake+'.h5', 'w') as f:
    f.create_dataset('traces', data= data);
    
#     print("Processed traces for "+ earthquake+ " loaded and stored for future use.")

# else:
#     print(".a file for this earthquake does not exist.")