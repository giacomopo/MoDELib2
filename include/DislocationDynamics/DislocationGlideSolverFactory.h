/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationGlideSolverFactory_h_
#define model_DislocationGlideSolverFactory_h_

#include <DislocationVelocitySolverBase.h>
#include <memory>
#include <Eigen/Dense>

namespace model
{
    
    template <typename DislocationNetworkType>
    struct DislocationGlideSolverBase : public DislocationVelocitySolverBase<DislocationNetworkType>
    {
//        const DislocationNetworkType& DN;
        DislocationGlideSolverBase(const DislocationNetworkType& );
//        virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };

    template <typename DislocationNetworkType>
    struct DislocationGlideSolverFactory
    {
        static std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>> getGlideSolver(const DislocationNetworkType& DN, const std::string& solverType);
    };
    
}
#endif
