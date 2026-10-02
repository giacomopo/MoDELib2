# /opt/local/bin/python3.14 plotLoops.py
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib
sys.path.append("../../python/")
from modlibUtils import *

# 0. Collect preliminary information
simulationDir=os.path.abspath(".")
materialFile=getStringInFile('inputFiles/polycrystal.txt','materialFile')
F,Flabels=readFfile(simulationDir+'/F')
ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
b_SI=ddBase.poly.b_SI
cs_SI=ddBase.poly.cs_SI
configIO=pyMoDELib.DDconfigIO(simulationDir+'/evl')
defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
dislocationNetwork=defectiveCrystal.dislocationNetwork()
runIDs=getFarray(F,Flabels,'runID')
times=getFarray(F,Flabels,'time [b/cs]')*b_SI/cs_SI/3600 # time in hours

# 1. Collect data from simulations
loopRaii=[]
loopTimes=[]
loopPlanes=[] # basal or prismatic
loopTypes=[] # vacancy or interstitial
dislocationDensity=[]
for k in range(0,len(runIDs)):
    runID=int(runIDs[k])
    time=times[k]
    configIO.read(runID)
    defectiveCrystal.initializeConfiguration(configIO)
    networkLengthTuple=dislocationNetwork.networkLength()
    dislocationDensity.append((networkLengthTuple[0]+networkLengthTuple[1]+networkLengthTuple[2]+networkLengthTuple[3])/ddBase.mesh.volume()/b_SI/b_SI)
    for loopID in dislocationNetwork.loops():
        loop=dislocationNetwork.loops().getRef(loopID)
        planeID=loop.glidePlane.planeBaseID() # 0=basal, 1,2,3=prismatic
        loopTypes.append(np.trace(loop.averagePlasticDistortion())<0.0)
        r=np.sqrt(loop.slippedArea()/np.pi)*b_SI*1.0e9 # average radius of current loop
        loopRaii.append(r)
        loopTimes.append(time)
        loopPlanes.append(planeID)
#        loopTypes.append(1)

# Convert to numPy
loopTimes = np.array(loopTimes)
loopRaii = np.array(loopRaii)
loopPlanes = np.array(loopPlanes)
loopTypes = np.array(loopTypes)

# 2. Create mappping for distinguishing loop types
mapping ={
    (0, 1): "basal vacancy",
    (0, 0): "basal interstitial",
    (1, 1): "prismatic vacancy", # (<a>1)
    (1, 0): "prismatic interstial", # (<a>1)
    (2, 1): "prismatic vacancy", # (<a>2)
    (2, 0): "prismatic interstial", # (<a>2)
    (3, 1): "prismatic vacancy", # (<a>3)
    (3, 0): "prismatic interstial", # (<a>3)
}

# 3. Get UNIQUE strings and assign a unique color to each string
unique_strings = sorted(list(set(mapping.values())))
cmap = plt.colormaps['plasma']
string_colors = {string: cmap(i / max(1, len(unique_strings) - 1))
                 for i, string in enumerate(unique_strings)}

# 4. Plot by unique mapping values
fig, axs = plt.subplots(1,1)
# 4. Group data by the final String using an array mask
for string_s in unique_strings:
    # Build a mask that finds all (p, t) pairs mapping to this string S
    string_mask = np.zeros(len(loopRaii), dtype=bool)
    
    for (pk, tk), mapped_string in mapping.items():
        if mapped_string == string_s:
            # Accumulate all points matching this specific (p, t) pair
            string_mask |= (loopPlanes == pk) & (loopTypes == tk)
            
    # Plot all data points for this string category in one clean batch
    axs.scatter(loopTimes[string_mask], loopRaii[string_mask],
                color=string_colors[string_s],
                label=string_s,
                s=5, alpha=0.6)

plt.legend(loc="upper right", markerscale=2)
fig.savefig("loopRadii.pdf", bbox_inches='tight', format='pdf')

#fig, axs = plt.subplots(1,1)
#axs.plot(times,rho,0.1,color='k')
#axs.set_xlabel('time [hours]')
#axs.set_ylabel('dislocation density [1/m^2]')
#fig.savefig("density.pdf", bbox_inches='tight', format='pdf')

fig, axs = plt.subplots(1,1)
axs.plot(times,dislocationDensity,0.1,color='k')
axs.set_xlabel('time [hours]')
axs.set_ylabel('dislocation density [1/m^2]')
fig.savefig("density.pdf", bbox_inches='tight', format='pdf')
