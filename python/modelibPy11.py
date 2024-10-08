# /opt/local/bin/python3.12 modelibPy11.py
import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("../build/tools/pyMoDELib")
import pyMoDELib

simulationDir="../tutorials/DislocationDynamics/periodicDomains/climb/"
ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)

# Microstructure Generation
microstructureGenerator=pyMoDELib.MicrostructureGenerator(ddBase)

spec0=pyMoDELib.ShearLoopDensitySpecification()
spec0.targetDensity=5.0e12
spec0.radiusDistributionMean=1.0e-07
spec0.radiusDistributionStd=0.0e-8
spec0.numberOfSides=20
#microstructureGenerator.addShearLoopDensity(spec0)

spec1=pyMoDELib.ShearLoopIndividualSpecification()
spec1.slipSystemIDs=[0,-1]
spec1.loopRadii=[27.0e-8,27.0e-8]
spec1.loopCenters=np.array([[200.0,0.0,0.0],[0.0,0.0,0.0]])
spec1.loopSides=[10,10]
#microstructureGenerator.addShearLoopIndividual(spec1)

spec2=pyMoDELib.PeriodicDipoleDensitySpecification()
spec2.targetDensity=5.0e13
#microstructureGenerator.addPeriodicDipoleDensity(spec2)

spec3=pyMoDELib.PeriodicDipoleIndividualSpecification()
spec3.slipSystemIDs=[0,-1]
spec3.exitFaceIDs=[1,0]
spec3.dipoleCenters=np.array([[200.0,0.0,0.0],[0.0,0.0,0.0]])
spec3.dipoleHeights=[100.,100.]
spec3.nodesPerLine=[4,4]
spec3.glideSteps=[10.,10.]
#microstructureGenerator.addPeriodicDipoleIndividual(spec3)

spec4=pyMoDELib.PrismaticLoopDensitySpecification()
spec4.targetDensity=5.0e13
spec4.radiusDistributionMean=3e-08
spec4.radiusDistributionStd=0e-08
#microstructureGenerator.addPrismaticLoopDensity(spec4)

spec5=pyMoDELib.PrismaticLoopIndividualSpecification()
spec5.slipSystemIDs=[0,7,13]
spec5.loopRadii=[5e-8,2e-7,2e-8]
spec5.loopCenters=np.array([[0,0,0],[500,600,500],[500,500,500]])
spec5.glideSteps=[10.,10.,300.]
#microstructureGenerator.addPrismaticLoopIndividual(spec5)

spec6=pyMoDELib.FrankLoopsDensitySpecification()
spec6.targetDensity=5.0e12
spec6.radiusDistributionMean=1.0e-07
spec6.radiusDistributionStd=0.0e-8
spec6.numberOfSides=20
spec6.areVacancyLoops=1
#microstructureGenerator.addFrankLoopsDensity(spec6)

spec7=pyMoDELib.FrankLoopsIndividualSpecification()
spec7.planeIDs=[0,-1]
spec7.loopRadii=[27.0e-8,27.0e-8]
spec7.loopCenters=np.array([[200.0,0.0,0.0],[0.0,0.0,0.0]])
spec7.loopSides=[10,10]
spec7.isVacancyLoop=[1,0]
#microstructureGenerator.addFrankLoopsIndividual(spec7)

spec8=pyMoDELib.StackingFaultTetrahedraDensitySpecification()
spec8.targetDensity=1.0e21
spec8.sizeDistributionMean=2.5e-08
spec8.sizeDistributionStd=1.0e-9
#microstructureGenerator.addStackingFaultTetrahedraDensity(spec8)

spec9=pyMoDELib.StackingFaultTetrahedraIndividualSpecification()
spec9.planeIDs=[0,1]
spec9.areInverted=[0,0]
spec9.sizes=[2.5e-08,2.5e-08]
spec9.basePoints=np.array([[0.0,0.0,0.0],[100.0,0.0,0.0]])
microstructureGenerator.addStackingFaultTetrahedraIndividual(spec9)


microstructureGenerator.writeConfigFiles(0) # write evel_0.txt (optional)

# DefectiveCrystal
defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)

points=np.array([[0.0,0.0,0.0],[1.0,0.0,0.0],[2.0,0.0,0.0]])
#
disp=defectiveCrystal.displacement(points)
print(disp)

# DislocationNetwork
DN=defectiveCrystal.dislocationNetwork()
print(len(DN.loops()))
print(len(DN.loopNodes()))

meshSize=100
for loopID in DN.loops():
    loop=DN.loops().getRef(loopID)
    meshedLoopVector=loop.meshed(meshSize)
    for meshedLoop in meshedLoopVector:
        disp=meshedLoop.plasticDisplacement(points)
        print(disp)
#
#rp=ddBase.poly.grain(1)





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

