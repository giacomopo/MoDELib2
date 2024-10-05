/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionsGenerator_H_
#define model_SphericalInclusionsGenerator_H_


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
#include <SphericalInclusionDensitySpecification.h>
#include <SphericalInclusionIndividualSpecification.h>


namespace model
{

class SphericalInclusionsGenerator 
//: public MicrostructureGeneratorBase
{
    
    static constexpr int dim=3;
    typedef Eigen::Matrix<double,dim,1> VectorDimD;

    
    bool generateSingle(MicrostructureGenerator& mg,const VectorDimD& C,const double& R, const Eigen::Matrix<double,1,dim*dim>& eT, const double& vrc,const int&type, const bool& allowOutside,const bool& allowOverlap);
    

    
public:
    

//    const bool allowOverlap;
//    const bool allowOutside;

  
    SphericalInclusionsGenerator(const SphericalInclusionDensitySpecification& spec,MicrostructureGenerator& mg);
    SphericalInclusionsGenerator(const SphericalInclusionIndividualSpecification& spec,MicrostructureGenerator& mg);

    
//    SphericalInclusionsGenerator(const std::string& fileName);
    
//    void generate(MicrostructureGenerator& mg) override;
//    void generateIndividual(MicrostructureGenerator& mg) override;
//    void generateDensity(MicrostructureGenerator& mg) override;
    
    
};

}
#endif
