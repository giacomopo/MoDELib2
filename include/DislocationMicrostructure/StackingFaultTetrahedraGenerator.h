/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StackingFaultTetrahedraGenerator_H_
#define model_StackingFaultTetrahedraGenerator_H_


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
#include <StackingFaultTetrahedraDensitySpecification.h>
#include <StackingFaultTetrahedraIndividualSpecification.h>

namespace model
{

    class StackingFaultTetrahedraGenerator //: public MicrostructureGeneratorBase
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        
        static bool generateSingleSFT(MicrostructureGenerator& mg,const int& planeID,const VectorDimD& basePoint,const bool& inverted,const int& sftSize);
        void generateIndividualFCC(const StackingFaultTetrahedraIndividualSpecification&,MicrostructureGenerator& mg);
        void generateDensityFCC(const StackingFaultTetrahedraDensitySpecification&,MicrostructureGenerator& mg);

    public:
        
//        StackingFaultTetrahedraGenerator(const std::string& fileName);
        StackingFaultTetrahedraGenerator(const StackingFaultTetrahedraDensitySpecification&,MicrostructureGenerator& mg);
        StackingFaultTetrahedraGenerator(const StackingFaultTetrahedraIndividualSpecification&,MicrostructureGenerator& mg);

        
        
//        void generateIndividual(MicrostructureGenerator& mg) override;
//        void generateDensity(MicrostructureGenerator& mg) override;
    };

}
#endif
