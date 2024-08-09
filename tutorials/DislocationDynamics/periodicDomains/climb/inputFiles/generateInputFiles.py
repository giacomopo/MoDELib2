import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

# Make a local copy of simulation parameters file and modify that copy if necessary
DDfile='DD.txt'
shutil.copy2('../../'+DDfile, '.') 

pf=PolyCrystalFile('../../../MaterialsLibrary/Al.txt');
pf.absoluteTemperature=470;
pf.dislocationMobilityType='default'
#pf.meshFile='../../../MeshLibrary/n10-id1.msh'
#pf.meshFile='../../../MeshLibrary/bicrystal_8188.msh'
pf.meshFile='../../../MeshLibrary/unitCube24.msh'
#pf.grain1globalX1=np.array([1,2,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([1,1,-3])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.alignToSlipSystem0=0
#pf.boxEdges=np.array([[1,2,1],[2,1,1],[1,1,-3]]) # i-throw is the direction of i-th box edge
#pf.boxScaling=np.array([1500,1500,1500]) # must be a vector of integers
pf.boxScaling=np.array([1061,1061,1061]) # must be a vector of integers
pf.X0=np.array([0.0,0.0,0.0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])

pf.write()

#print(pf.A)
