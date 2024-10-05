/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PyGlideSolver_h_
#define model_PyGlideSolver_h_

//#include <PlanarDislocationNode.h>


#ifndef NDEBUG
#define VerbosePyGlideSolver(N,x) if(verbosePyGlideSolver>=N){std::cout<<x;}
#else
#define VerbosePyGlideSolver(N,x)
#endif

#include <iostream>

#include <DislocationDynamicsModule.h>
#include <DislocationGlideSolverFactory.h>

#ifdef _MODEL_PYBIND11_
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#endif

namespace model
{
    
    template <typename DislocationNetworkType>
    class PyGlideSolver : public DislocationGlideSolverBase<DislocationNetworkType>
    {
        typedef typename DislocationNetworkType::VectorDim   VectorDim;
        typedef typename DislocationNetworkType::MatrixDim   MatrixDim;
        typedef typename DislocationNetworkType::NetworkNodeType   NetworkNodeType;
        typedef typename DislocationNetworkType::LoopNodeType   LoopNodeType;
        typedef typename DislocationNetworkType::LoopType   LoopType;
        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;

        public:
        
        const std::string pyModuleName;
        
#ifdef _MODEL_PYBIND11_
        pybind11::module pyModule;
#endif

        Eigen::MatrixXd updateUnconnectedNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& removedNodeIDmap) const;
        Eigen::MatrixXd updateBoundaryNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& removedNodeIDmap) const;

        PyGlideSolver(const DislocationNetworkType&);
        Eigen::VectorXd getNodeVelocities() const override;
        
    };
    
}
#endif
