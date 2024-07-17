# /opt/local/bin/python3.12 modelibPy11.py
import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("../build/tools/pyMoDELib")
import pyMoDELib
#sys.path.append("../build/")
#import libMoDELib

simulationDir="../tutorials/DislocationDynamics/periodicDomains/uniformLoadController/"
ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)

points=np.array([[0.0,0.0,0.0],[1.0,0.0,0.0],[2.0,0.0,0.0]])
print(points)
disp=defectiveCrystal.displacement(points)

DN=defectiveCrystal.dislocationNetwork()
print(len(DN.loops()))
print(len(DN.loopNodes()))

for loopID in DN.loops():
    loop=DN.loops().getRef(loopID)
    #print(type(loop))
    #point=np.array([0.0,0.0,0.0])
    #loop.meshed(100)
    meshedLoopVector=loop.meshed(100)
#    print(type(meshedLoopVector[0]))
#    for meshedLoop in meshedLoopVector:
#        a=1
#    loop.solidAngle(point)
#    print(loop.solidAngle(point))


#for node in DN.loopNodes():
#    print(node)

rp=ddBase.poly.grain(1)
#print(rp)

#mesh=pyMoDELib.SimplicialMesh()
#print(mesh.xMin())

#
## create the DefectsFieldsExtractor object
#dfe=DefectsFields.DefectsFieldsExtractor(simulationDir)
#
## Read a pre-existing configuration file (e.g. readConfiguration(X) reads file simulationDir/evl/evl_X.txt)
#dfe.readConfiguration(0)
#
## Alternatively, generate a new configuration using simulationDir/inputFiles/initialMicrostructure.txt
##dfe.readMicrostructure()
##dfe.writeConfiguration(0) # Optional. Write the generated configuration to file (writeConfiguration(X) writes file simulationDir/evl/evl_X.txt)
#
## grab the domain corners
#ldc=dfe.lowerDomainCorner()
#udc=dfe.upperDomainCorner()
#
## Extracting grids of values on a a y-z plane: e.g. solid angle and sigma_11
#n=200
#x=np.linspace(4*ldc[0], 4*udc[0], num=n) # grid x-range
#z=np.linspace(4*ldc[2], 4*udc[2], num=n) # grid z-range
#y=0.5*(ldc[1]+udc[1]) # grid position in y
#solidAngle=np.empty([n, n]) # grid of solid angle values
#s13=np.empty([n, n]) # grid of sigma_11 values
#for i in range(0,z.size):
#    for j in range(0,x.size):
#        solidAngle[j,i]=dfe.solidAngle(x[j],y,z[i])
#        stress=dfe.dislocationStress(x[j],y,z[i])
#        s13[i,j]=stress[0,2]
#
#fig=plt.figure()
#plt.imshow(solidAngle, origin='lower',cmap='jet')
#plt.colorbar()
#plt.show()
#
#fig=plt.figure()
#plt.imshow(s13, origin='lower',cmap='jet',vmin = -0.01,vmax = 0.01)
#plt.colorbar()
#plt.show()

