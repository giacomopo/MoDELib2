# sudo /usr/local/bin/python3 -m pip install PyQt5
# sudo /usr/local/bin/python3 -m pip install matplotlib
# sudo /usr/local/bin/python3 -m pip install numpy
import sys, string, os
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
sys.path.insert(0, '../../../../python')
from modlibUtils import *

font = {'family' : 'normal',
        'size'   : 16}
plt.rc('font', **font)


def getLoopRadius(folderName):
    F=np.loadtxt(folderName+'/F/F_0.txt');
    print(folderName+'/F/F_0.txt has size ' + str(np.shape(F)));
    R=np.empty(np.size(F[:,0]))
    n=0;
    for runID in F[:,0]:
        #print(int(runID))
        evlFile=folderName+'/evl/evl_'+str(int(runID))
        print('reading '+ evlFile);
        evl=readEVLtxt(evlFile)
        c=np.mean(evl.nodes, axis=0)
        C=np.tile(c,[evl.nodes.shape[0],1])
        nodesC=evl.nodes-C;
        nodesR=np.sqrt(np.sum(np.square(nodesC),axis=1))
        R[n]=np.mean(nodesR);
        n=n+1;
    return F[:,1], R;

# main code
materialFile='inputFiles/'+getStringInFile('inputFiles/polycrystal.txt','materialFile')
print('materialFile='+materialFile)
mu_SI=getValueInFile(materialFile,'mu0_SI');     #[Pa]
rho_SI=getValueInFile(materialFile,'rho_SI');     #[kg/m^3]
b_SI=getValueInFile(materialFile,'b_SI');     #[m]
v_dd2SI=np.sqrt(mu_SI/rho_SI);
t_dd2SI=b_SI/v_dd2SI;


# simulation data
t,R=getLoopRadius('.')
F,Flabels=readFfile('./F')
trBp=getFarray(F,Flabels,'tr(betaP)')


fig1 = plt.figure()
ax1=plt.subplot(1,1,1)
ax1.plot(t*t_dd2SI, R*b_SI*1e10,color='b')
ax1.grid()
#ax1.legend()
ax1.set_xlabel('time [sec]')
ax1.set_ylabel('loop radius [$\AA$]',color='b')

ax2=ax1.twinx()
ax2.plot(t*t_dd2SI, trBp,color='r')
ax2.set_ylabel('trBp [-]',color='r')

fig1.autofmt_xdate()
#plt.show()
fig1.savefig("fig1.pdf", bbox_inches='tight')




