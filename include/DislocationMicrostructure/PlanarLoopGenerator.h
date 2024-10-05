/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarLoopGenerator_H_
#define model_PlanarLoopGenerator_H_


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
#include <PlanarLoopIndividualSpecification.h>

namespace model
{

    class PlanarLoopGenerator 
//: public MicrostructureGeneratorBase
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;

        
        static void generateSingle(MicrostructureGenerator& mg,const VectorDimD& b,const VectorDimD& n,const std::vector<VectorDimD>& loopPoints);

//        static void insertJunctionLoop(MicrostructureGenerator& mg,
//                                std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>>& uniqueNetworkNodeMap,
//                                const std::vector<VectorDimD>& loopNodePos,
//                                const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
//                                const VectorDimD& b,
//                                const VectorDimD& unitNormal,
//                                const VectorDimD& P0,
//                                const size_t& grainID,
//                                const DislocationLoopIO<dim>::DislocationLoopType& loopType);
        
    public:
        
//        PlanarLoopGenerator(const PlanarLoopIndividualSpecification& spec,const std::string& fileName);
        PlanarLoopGenerator(const PlanarLoopIndividualSpecification& spec,MicrostructureGenerator& mg);

//        void generate(MicrostructureGenerator& mg) override;
//        void generateIndividual(MicrostructureGenerator& mg) override;
//        void generateDensity(MicrostructureGenerator& mg) override;

    };

}
#endif
