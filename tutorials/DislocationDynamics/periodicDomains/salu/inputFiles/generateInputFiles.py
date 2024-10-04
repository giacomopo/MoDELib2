import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

# Make a local copy of simulation parameters file and modify that copy if necessary
DDfile='DD.txt'
shutil.copy2('../../'+DDfile, '.')
setInputVariable(DDfile,'useFEM','0')
setInputVariable(DDfile,'useDislocations','1')
setInputVariable(DDfile,'useInclusions','1')
setInputVariable(DDfile,'useElasticDeformation','1')
setInputVariable(DDfile,'useClusterDynamics','0')
setInputVariable(DDfile,'Nsteps','100000')
setInputVariable(DDfile,'timeSteppingMethod','adaptive')
setInputVariable(DDfile,'glideSolverType','Galerkin')
setInputVariable(DDfile,'climSolverType','none')
setInputVariable(DDfile,'outputFrequency','10')

# Generate polycrystal file
pf=PolyCrystalFile('../../../MaterialsLibrary/W.txt');
pf.absoluteTemperature=1300;
pf.dislocationMobilityType='default'
#pf.meshFile='../../../MeshLibrary/n10-id1.msh'
#pf.meshFile='../../../MeshLibrary/bicrystal_8188.msh'
pf.meshFile='../../../MeshLibrary/unitCube24.msh'
#pf.grain1globalX1=np.array([1,2,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([1,1,-3])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.alignToSlipSystem0=1
#pf.boxEdges=np.array([[1,2,1],[2,1,1],[1,1,-3]]) # i-throw is the direction of i-th box edge
#pf.boxScaling=np.array([1500,1500,1500]) # must be a vector of integers
pf.boxScaling=np.array([400,400,400]) # must be a vector of integers
pf.X0=np.array([0.0,0.0,0.0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])

pf.write()

# make a local copy of first microstructure file, and modify that copy if necessary
microstructureTemplate='periodicDipoleIndividual.txt';
shutil.copy2('../../../MicrostructureLibrary/'+microstructureTemplate, '.') # target filename is /dst/dir/file.ext
setInputVector(microstructureTemplate,'slipSystemIDs',np.array([0]),'slip system IDs for each dipole')
setInputVector(microstructureTemplate,'exitFaceIDs',np.array([2]),'4 is for edge, 2 for screw')
setInputVector(microstructureTemplate,'dipoleCenters',np.array([0.0,0.0,0.0]),'center of the dipole dipole')
setInputVector(microstructureTemplate,'nodesPerLine',np.array([10]),'number of extra nodes on each dipole')
setInputVector(microstructureTemplate,'dipoleHeights',np.array([50]),'height of each dipole, in number of planes')
setInputVector(microstructureTemplate,'glideSteps',np.array([10.0]),'step of each dipole in the glide plane')

# make a local copy of second microstructure file, and modify that copy if necessary
microstructureTemplate='sphericalInclusionsIndividual.txt';
shutil.copy2('../../../MicrostructureLibrary/'+microstructureTemplate, '.') # target filename is /dst/dir/file.ext
setInputVector(microstructureTemplate,'radii_SI',np.array([4.3e-9]),'inclusion radius')
setInputVector(microstructureTemplate,'eigenDistortions',np.array([0.03,0,0,0,0.03,0,0,0,0.03]),'1 is for edge, 0 for screw')
setInputVector(microstructureTemplate,'centers',np.array([0,0,0]),'center of the inclusion')
setInputVector(microstructureTemplate,'phaseIDs',np.array([1]),'number of extra nodes on each dipole')

