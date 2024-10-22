/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_H_
#define model_BCClattice_H_

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <SlipSystem.h>
#include <PolycrystallineMaterialBase.h>
#include <DislocationMobilityBCC.h>
#include <RationalLatticeDirection.h>
#include <SingleCrystalBase.h>
#include <DislocationMobilitySelector.h>
//#include <DislocationMobilityBCC.h>
//#include <SecondPhase.h>

namespace model
{
    
    template<int dim>
    struct BCClattice
    {
        
    };
    
    template<>
    struct BCClattice<3> : public SingleCrystalBase<3>
    /*                  */,private SingleCrystalBase<3>::PlaneNormalContainerType
    /*                  */,private SingleCrystalBase<3>::SlipSystemContainerType
    /*                  */,private SingleCrystalBase<3>::SecondPhaseContainerType
    {
//        static constexpr auto name="BCC";
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename SingleCrystalBase<dim>::MatrixDim MatrixDim;
        typedef typename SingleCrystalBase<dim>::PlaneNormalContainerType PlaneNormalContainerType;
        typedef typename SingleCrystalBase<dim>::SlipSystemContainerType SlipSystemContainerType;
        typedef typename SingleCrystalBase<dim>::SecondPhaseContainerType SecondPhaseContainerType;

                
        BCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile);
        static Eigen::Matrix<double,dim,dim> getLatticeBasis();
        std::vector<std::shared_ptr<GlidePlaneBase>> getPlaneNormals(const PolycrystallineMaterialBase& material,const std::string& polyFile) const;
        std::vector<std::shared_ptr<SlipSystem>> getSlipSystems(const PolycrystallineMaterialBase& material,const PlaneNormalContainerType& plN) const;
        SecondPhaseContainerType getSecondPhases(const PolycrystallineMaterialBase& material,const PlaneNormalContainerType& plN);
        
        const PlaneNormalContainerType& planeNormals() const override;
        const SlipSystemContainerType& slipSystems() const override;
        const SecondPhaseContainerType& secondPhases() const override;
        
    };
    
} // namespace model
#endif

