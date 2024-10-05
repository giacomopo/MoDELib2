/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationClimbSolverFactory_h_
#define model_DislocationClimbSolverFactory_h_

#include <DislocationVelocitySolverBase.h>
#include <memory>
#include <Eigen/Dense>
#include <ClusterDynamics.h>
#include <TypeTraits.h>
#include <ClusterDynamicsParameters.h>


namespace model
{
    
    template <typename DislocationNetworkType>
    struct DislocationClimbSolverBase : public DislocationVelocitySolverBase<DislocationNetworkType>
    /*                               */,private std::vector<Eigen::Array<double,1,ClusterDynamicsParameters<TypeTraits<DislocationNetworkType>::dim>::mSize>>

    {
//        const DislocationNetworkType& DN;
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef std::vector<Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize>> ScalarVelocitiesContainerType;
        const double glideEquilibriumRate;
        const ClusterDynamics<dim>* const CD;
        const double vClimbRef;
        
        DislocationClimbSolverBase(const DislocationNetworkType&,const ClusterDynamics<dim>* const);
//        const ClusterDynamics<dim>*  getCD() const;
        double getVclimbRef() const;
        virtual void computeClimbScalarVelocities() =0;

        const ScalarVelocitiesContainerType& scalarVelocities() const;
        ScalarVelocitiesContainerType& scalarVelocities();

        
  //      virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };

    template <typename DislocationNetworkType>
    struct DislocationClimbSolverFactory
    {
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;

//        static const ClusterDynamics<dim>*  getCD(const DislocationNetworkType& DN);

        static std::shared_ptr<DislocationClimbSolverBase<DislocationNetworkType>> getClimbSolver(const DislocationNetworkType& DN, const std::string& solverType);
    };
    
}
#endif
