/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleGenerator_H_
#define model_PeriodicDipoleGenerator_H_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
//#include <MicrostructureGeneratorBase.h>
#include <PeriodicDipoleDensitySpecification.h>

namespace model
{

class PeriodicDipoleGenerator //: public MicrostructureGeneratorBase
{
    
    typedef Eigen::Matrix<double,3,1> VectorDimD;

    
    static void generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& dipolePoint,const int& exitFaceID,const int& dipoleHeight,const int& dipoleNodes,double glideStep);
    
    
//    static void insertJunctionLoop(MicrostructureGenerator& mg,
//                            std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>>& uniqueNetworkNodeMap,
//                            const std::vector<VectorDimD>& loopNodePos,
//                            const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
//                            const VectorDimD& b,
//                            const VectorDimD& unitNormal,
//                            const VectorDimD& P0,
//                            const size_t& grainID,
//                            const DislocationLoopIO<dim>::DislocationLoopType& loopType);
    
public:
    
//    PeriodicDipoleGenerator(const std::string& fileName);
    PeriodicDipoleGenerator(const PeriodicDipoleDensitySpecification& spec,MicrostructureGenerator& mg);
    PeriodicDipoleGenerator(const PeriodicDipoleIndividualSpecification& spec,MicrostructureGenerator& mg);

    
    
//    void generate(MicrostructureGenerator& mg) override;
//    void generateIndividual(MicrostructureGenerator& mg) override;
//    void generateDensity(MicrostructureGenerator& mg) override;
    
    
};

}
#endif
