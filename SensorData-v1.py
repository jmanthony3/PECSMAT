'''
v 1 JMC 7/23/2024
This iteration of the code is just reviewing how it works, and making the code more readable. 
It appears the code will work on specific data from Sam Scott's 304 data, but the Copper data has different 
data which is confusing this code here
'''


#Import necessary libraries
from turtle import back
import numpy as np
import scipy.signal as sps
import scipy.integrate as spi
import pandas as pd
import matplotlib.pyplot as plt
import time
#import json
## Code here is used to import from CSV files that were saved for MATLAB purposes

#Declare vars
#time step of the sensors, (should also be called baud)
ts = np.float64(1/44100) 
#Base name of the file set we are looking at
import_filename = "04-02-25_14_30_Fat-Hex1mm1bar_" 
#Save every how many impacts, this shows up in the Config .csv
save_every = 10
#Total number of impacts, this shows up in the Config .csv
num_loc = 2
#If multiple impacts per location
num_per_loc = 1
# ms #I dunno
cut_time = 2 
#how many samples backwards from the impact should we look into?
back_sample = 5
#The force tolerance which is required for the data to be from an actual impact
tolerance = 5 # in lbf? Looks like it...
# Filter parameters
Freq = 44100
C_freq = 1500
freq_step = 100
start_freq = 2500
num_freq = 1
order = 1
#plotrange = [0,38000]

#Declare functions
def fa_convert(f,a):
    ## Converts data from DAQ into actual units
    # FA_CONVERT converts force and acceleration inputs from the USB DAQ to
    # proper values. USB DAQ is PCB piezo 485B39, sn 1558. Force transducer is PCB piezo
    # 208C04, sn LW57182, accelerometer is PCB piezo 355B02, sn 51922 
    # Outputs are in lbf and in/s^2, and acceleration is multiplied by -1 for convention

    # DAQ Calibration Parameters
    FreqAdjust = 1.00329051
    A_sense = 806990 * FreqAdjust
    B_sense = 808064 * FreqAdjust
    ForceSens = 0.005153 # 208C04 transducer sensetivity, V/bar
    AccelSens = 0.01021  # 355B02 accelerometer sensetivity, V/g   
    #these should return as bar and g for the force and acceleration respectively
    return(f[:]*(np.float64(2**23))/(A_sense*ForceSens), -1*a[:]*(np.float64(2**23))*386.08858267717/(B_sense*AccelSens))

def fa_allign(f,a,ts,tol = 25, ct = 10, back_samples = 5):
    ## Alligns impact data to beginning of impact as decided by force tol
    # ts is timestep, tol is impact force START in lbf, and ct is cutoff time in ms
    # Back samples is samples BEHIND tolerance to consider

    for cn in range(0,f.shape[0]):
        #Find index of first force entry to go above the tolerance value "tol"
        if (f[cn] >= tol):
            ID = cn-back_samples
            print("Tolerance found")
            break

        #Find index of first force entry to go above the tolerance value "tol"?
        #But are these different tol values?
    for cn in range(0,a.shape[0]):
        if (a[cn] >= tol):
            AID = cn-back_samples
            break

    cut_samples = np.int32(ct*44100/1000)
    CID = np.int16(ID+back_samples+cut_samples)
    return(ID,CID,AID)

def fa_integrate(f,a,ID,CID,ts):
    ## This function integrates to obtain velocity, FV, energy, position

    va = spi.cumtrapz(a,initial=0)*ts
    sa = spi.cumtrapz(va,initial=0)*ts
    fv = f*va
    en = spi.cumtrapz(fv[ID:(CID+1)],initial=0)*ts
    return(va,fv,sa,en)

k = 0
for cn in range(num_loc):
    f_index = np.int32(cn*save_every)

    for cn2 in range(num_per_loc):
        fafile = import_filename + str(f_index) + "_" + str(cn2) + ".csv"
        fa = np.loadtxt(fafile, delimiter=",", dtype="str")[1:3,:].astype(dtype=np.float64)
        # print(fa)

        a = fa[1,:] #Read in the acceleration data
        f = fa[0,:] #Read in the force data
        # print("Force data")
        # print(f)
        # print("Acceleration data")
        # print(a)
        fp, ap = fa_convert(f,a)
        #Convert to SI if desired
        # fp = fp* 4.44822 #Convert to Newtons
        # ap = ap* 0.0254 #Convert to m/s^2
        
        ID, CID, AID = fa_allign(fp,ap,ts,tolerance,cut_time,back_sample)
        # plt.figure()
        # plt.plot(ap[ID:CID])

        # plt.figure()
        # plt.plot(fp[ID:CID])
        
        print(ID,CID)
        VL, FV, PS, EN = fa_integrate(f,a,ID,CID,ts)
        
# Plotting 
        for cn2 in range(num_freq):
            C_freq = start_freq + cn2*freq_step
            fb, fa = sps.butter(order,C_freq,'lowpass',analog=False,fs=Freq)
            fig = plt.figure()
            gs = fig.add_gridspec(2,hspace=0.55)
            axs = gs.subplots()
            fig.suptitle("Pre and Post filter, freq = "+str(C_freq))
            axs[0].plot(ap[ID:CID])
            a = sps.lfilter(fb, fa, ap)
            axs[0].plot(a[ID:CID])
            axs[0].set(title="Acceleration",ylabel="in/s^2")
            axs[1].plot(fp[ID:CID])
            f = sps.lfilter(fb, fa, fp)
            axs[1].plot(f[ID:CID])
            axs[1].set(title="Force",ylabel="lbf")                                                
                  
        #plt.show()

        fname = "FF" +str(k) + ".csv"
        aname = "AA" +str(k) + ".csv"
        iname = "II" +str(k) + ".csv"
        vname = "VV" +str(k) + ".csv"
        ename = "EE" +str(k) + ".csv"
        k = k+1
        np.savetxt(fname,fp,delimiter=",")
        np.savetxt(aname,ap,delimiter=",")
        np.savetxt(iname,[ID,CID,AID],delimiter=",")
        np.savetxt(vname,VL,delimiter=",")
        np.savetxt(ename,EN,delimiter=",")

plt.show()
# plt.plot(f[ID:CID])
# plt.show()
# plt.plot(PS[0:CID])
# print(max(PS[0:CID]))

# fig = plt.figure()
# gs = fig.add_gridspec(4,hspace=0.55)
# axs = gs.subplots()
# fig.suptitle("Initial Run 1 Data")
# axs[0].plot(a)
# axs[0].set(title="Acceleration",ylabel="in/s^2")
# axs[1].plot(VL)
# axs[1].set(title="Velocity",ylabel="in/s")
# axs[2].plot(f)
# axs[2].set(title="Force",ylabel="lbf")
# axs[3].plot(EN)
# axs[3].set(title="Energy",ylabel="lbf-in")

# plt.show()