# from obspy import read
import numpy as np
from math import sin, cos
from tqdm import tqdm
import h5py
import obspy

from obspy.taup import TauPyModel
model = TauPyModel(model="prem")

f= open("earthquake.txt", "r")
earthquake= f.read()
f.close()



st= obspy.read("/mnt/data2/isha/unique_pwin33_traces/"+earthquake+"*");
# st= obspy.read("/mnt/datawave3_new/"+earthquake+"/*.a/good_traces/*")
               # abhinav/workspace/Bp/traces/"+earthquake+"/*"
nr= len(st)

Re= 6371; # radius of earth in km (as used in obspy)

recx= np.zeros(len(st))
recy= np.zeros_like(recx)
recz= np.zeros_like(recx)

ids= np.empty((nr), dtype= object)

theta= np.empty((nr))
phi= np.empty((nr))

evlat, evlon, evdep= [st[0].stats.sac.evla, st[0].stats.sac.evlo, st[0].stats.sac.evdp]
travel_times= np.empty((nr))
network= np.empty((nr), dtype= object)
for i in range(0, nr):
    st_lat, st_lon, st_dep= st[i].stats.sac.stla, st[i].stats.sac.stlo, st[i].stats.sac.stdp


    epd_deg=obspy.geodetics.locations2degrees(evlat, evlon,
                                              st[i].stats.sac.stla, st[i].stats.sac.stlo)

    tt= model.get_travel_times(source_depth_in_km= evdep, distance_in_degree=epd_deg, 
                                phase_list= ["ttp"], receiver_depth_in_km=st_dep)

    theta[i]= obspy.geodetics.gps2dist_azimuth(evlat, evlon, st_lat, st_lon, a= Re*1000, f= 0)[1]
    if theta[i]>180:
        theta[i]= theta[i]- 360
    phi[i]= tt[0].takeoff_angle
    travel_times[i]= tt[0].time
    network[i]= st[i].stats.network;
    ids[i]= st[i].id;

# "~/workspace/Bp/utilities/"
with h5py.File("/Bp_eq_data_processed/"+ earthquake+".h5", 'w') as f:
    f.create_dataset('center', data= [st[0].stats.sac.evla, st[0].stats.sac.evlo, st[0].stats.sac.evdp])
    f.create_dataset('theta', data= theta);
    f.create_dataset('phi', data= phi);
    f.create_dataset("travel_times", data= travel_times)
    f.create_dataset("network", data= network)
    f.create_dataset("ids", data= ids)

# print("Process for "+ earthquake+ " completed")