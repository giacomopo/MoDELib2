# /opt/local/bin/python3.14 plotLoops.py
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib
sys.path.append("../../python/")
from modlibUtils import *

simulationDir=os.path.abspath(".")
materialFile=getStringInFile('inputFiles/polycrystal.txt','materialFile')
#b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')
#cs_SI=getValueInFile('inputFiles/'+materialFile,'cs_SI')

F,Flabels=readFfile(simulationDir+'/F')
ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
b_SI=ddBase.poly.b_SI
cs_SI=ddBase.poly.cs_SI
configIO=pyMoDELib.DDconfigIO(simulationDir+'/evl')
defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
dislocationNetwork=defectiveCrystal.dislocationNetwork()
runIDs=getFarray(F,Flabels,'runID')
times=getFarray(F,Flabels,'time [b/cs]')*b_SI/cs_SI/3600 # time in hours
#rhoG=getFarray(F,Flabels,'glissile density [m^-2]')
#rhoS=getFarray(F,Flabels,'sessile density [m^-2]')
#rho=rhoG+rhoS


loopRaii=[]
loopTimes=[]
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
        loopRaii.append(np.sqrt(loop.slippedArea()/np.pi)*b_SI*1.0e9)
        loopTimes.append(time)
#        loopRunIDs.append(runID)

fig, axs = plt.subplots(1,1)
axs.scatter(loopTimes,loopRaii,0.1,color='k')
axs.set_xlabel('time [hours]')
axs.set_ylabel('average loop radii [nm]')
#fig.show()
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
