# sudo /usr/local/bin/python3 -m pip install PyQt5
# sudo /usr/local/bin/python3 -m pip install matplotlib
# sudo /usr/local/bin/python3 -m pip install numpy
import sys, string, os
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, '../../../../python')
from modlibUtils import *
plt.rcParams['text.usetex'] = True

# main code
materialFile='inputFiles/'+getStringInFile('inputFiles/polycrystal.txt','materialFile')
print('materialFile='+materialFile)
mu_SI=getValueInFile(materialFile,'mu0_SI')/1e6
print('mu='+ str(mu_SI) + ' MPa')
nu=getValueInFile(materialFile,'nu')
print('nu='+ str(nu))
E_SI=2*mu_SI*(1.0+nu)
print('E='+ str(E_SI) + ' MPa')

rho_SI=getValueInFile(materialFile,'rho_SI');    #[kg/m^3]
b_SI=getValueInFile(materialFile,'b_SI');  #[m]
v_dd2SI=np.sqrt(mu_SI/rho_SI);
t_dd2SI=b_SI/v_dd2SI;
T=getValueInFile('inputFiles/polycrystal.txt','absoluteTemperature')

F,Flabels=readFfile('./F')

runDI=getFarray(F,Flabels,'runID')
trBp=getFarray(F,Flabels,'tr(betaP)')
fig1 = plt.figure()
ax11=plt.subplot(1,1,1)
ax11.plot(runDI,trBp)
ax11.grid()
plt.show()
#fig1.savefig("fig1.pdf", bbox_inches='tight')




