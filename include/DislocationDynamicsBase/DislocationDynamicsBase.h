/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * Additional contributors:
 *   Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationDynamicsBase_H_
#define model_DislocationDynamicsBase_H_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <FiniteElement.h>
#include <GlidePlaneFactory.h>
#include <PeriodicGlidePlaneFactory.h>
#include <VoigtTraits.h>

namespace model
{
    template <int _dim>
    struct DislocationDynamicsBase
    {
        // desired members:
        //  those common to MicrostructureGenerator and DefectiveCrystal
        //  so that their members can be references to those in this class
        static constexpr int dim=_dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef DislocationDynamicsBase<dim> DislocationDynamicsBaseType;
        typedef LagrangeElement<3,2> ElementType;
        typedef typename ElementType::BaryType BaryType;

        DislocationDynamicsBase( const std::string& folderName);
        DefectiveCrystalParameters simulationParameters;
        const SymmetricVoigtTraits<dim> voigtTraits;
        const SimplicialMesh<dim> mesh;

        const bool isPeriodicDomain;
        const std::vector<int> periodicImageSize;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeBasis;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeReciprocalBasis;
        const std::vector<VectorDim> periodicShifts;
        const Polycrystal<dim> poly;
        std::shared_ptr<FiniteElement<ElementType>> fe;
        GlidePlaneFactory<3> glidePlaneFactory;
        PeriodicGlidePlaneFactory<3> periodicGlidePlaneFactory;
        const double EwaldLength;

        Eigen::VectorXd periodicCoordinates(const VectorDim& x) const;

    };
} // namespace model
#endif
