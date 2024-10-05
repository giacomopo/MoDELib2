/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinClimbSolver_h_
#define model_GalerkinClimbSolver_h_

//#include <PlanarDislocationNode.h>


#ifndef NDEBUG
#define VerboseGalerkinClimbSolver(N,x) if(verboseGalerkinClimbSolver>=N){std::cout<<x;}
#else
#define VerboseGalerkinClimbSolver(N,x)
#endif

#include <vector>

#include <DislocationDynamicsModule.h>
#include <DislocationClimbSolverFactory.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    class GalerkinClimbSolver : public DislocationClimbSolverBase<DislocationNetworkType>
    {
        typedef typename DislocationNetworkType::VectorDim   VectorDim;
        typedef typename DislocationNetworkType::MatrixDim   MatrixDim;
        typedef typename DislocationNetworkType::NetworkNodeType   NetworkNodeType;
        typedef typename DislocationNetworkType::NetworkLinkType NetworkLinkType;
        typedef typename DislocationNetworkType::LoopNodeType   LoopNodeType;
        typedef typename DislocationNetworkType::LoopType   LoopType;
        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;     // Cv, Ci, C2i, C3i
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;  // Nc, Na1, Na2, Na3, rc, ra1, ra2, ra3
        typedef Eigen::Matrix<double,2,mSize> ForceVectorMatrixType;
        typedef Eigen::Matrix<double,2*mSize,2> StiffnessMatrixType;
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        typedef std::deque<Eigen::Triplet<double> > TripletContainerType;

        ForceVectorMatrixType clusterForceVector(const NetworkLinkType& networkLink) const;
        ForceVectorMatrixType clusterForceKernel(const int& k,const NetworkLinkType& networkLink) const;
        StiffnessMatrixType clusterStiffnessKernel(const int& k,const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const;
        StiffnessMatrixType clusterStiffnessMatrix(const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const;
        
        
        Eigen::VectorXd getNodeVelocitiesBulk() const;
        Eigen::VectorXd getNodeVelocitiesPipe() const;

        void computeClimbScalarVelocitiesBulk();

        
        public:
        
        GalerkinClimbSolver(const DislocationNetworkType&,const ClusterDynamics<dim>* const );
        Eigen::VectorXd getNodeVelocities() const override;
        void computeClimbScalarVelocities() override;
        
    };
    
}
#endif
