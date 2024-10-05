/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_H_
#define model_ClusterDynamicsFEM_H_

#ifdef MODELIB_CHOLMOD // SuiteSparse Cholmod module
#include <Eigen/CholmodSupport>
#endif

#ifdef MODELIB_UMFPACK // SuiteSparse UMFPACK module
#include <Eigen/UmfPackSupport>
#endif

#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <EvalFunction.h>
#include <FixedDirichletSolver.h>
#include <DDconfigIO.h>
#include <MicrostructureBase.h>
#include <SecondOrderReaction.h>
#include <MicrostructureContainer.h>

namespace model
{

    template <int dim>
    struct FluxMatrix : public EvalFunction<FluxMatrix<dim>>
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        constexpr static int mSize=ClusterDynamicsParameters<dim>::mSize;
        constexpr static int rows=dim*mSize;
        constexpr static int cols=rows;
        typedef Eigen::Matrix<double,rows,cols> MatrixType;

        const ClusterDynamicsParameters<dim>& cdp;
        
        FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in);
        const MatrixType operator() (const ElementType& elem, const BaryType& bary) const;
    };

    template <int dim>
    struct InvDscaling : public EvalFunction<InvDscaling<dim>>
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        constexpr static int mSize=ClusterDynamicsParameters<dim>::mSize;
        constexpr static int rows=mSize;
        constexpr static int cols=rows;
        typedef Eigen::Matrix<double,rows,cols> MatrixType;

        const ClusterDynamicsParameters<dim>& cdp;
        
        InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in);
        const MatrixType operator() (const ElementType& elem, const BaryType& bary) const;
    };

    template<int dim>
    struct ClusterDynamicsFEM
    {
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
        typedef FiniteElement<ElementType> FiniteElementType;
        static constexpr int dVorder=4;
        typedef IntegrationDomain<FiniteElementType,0,dVorder,GaussLegendre> VolumeIntegrationDomainType;
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;
        typedef TrialFunction<'m',mSize,FiniteElementType> MobileTrialType;
        typedef TrialFunction<'i',iSize,FiniteElementType> ImmobileTrialType;
//        typedef TrialFunction<'z',dim,FiniteElementType> DiffusiveTrialType;
        typedef TrialGrad<MobileTrialType> MobileGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileGradType> MobileFluxType;
        
        typedef TrialProd<InvDscaling<dim>,MobileTrialType> MobileTestType;
        typedef TrialGrad<MobileTestType> MobileTestGradType;
        typedef BilinearForm<MobileTestGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;

//        typedef BilinearForm<MobileGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;
        typedef BilinearWeakForm<MobileBilinearFormType,VolumeIntegrationDomainType> MobileBilinearWeakFormType;

        typedef TrialFunction<'d',mSize,FiniteElementType> MobileIncrementTrialType;
        typedef TrialProd<InvDscaling<dim>,MobileIncrementTrialType> MobileIncrementTestType;
        typedef TrialGrad<MobileIncrementTestType> MobileIncrementTestGradType;

        typedef TrialGrad<MobileIncrementTrialType> MobileIncrementGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileIncrementGradType> MobileIncrementFluxType;
        typedef BilinearForm<MobileIncrementTestGradType,TrialProd<Constant<double,1,1>,MobileIncrementFluxType>> MobileIncrementBilinearFormType;
        typedef BilinearWeakForm<MobileIncrementBilinearFormType,VolumeIntegrationDomainType> MobileIncrementBilinearWeakFormType;
        typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMatrixType;
#ifdef CHOLMOD_H // SuiteSparse Cholmod (LLT) module
    typedef Eigen::CholmodSupernodalLLT<SparseMatrixType> LltSolverType;
#else
    typedef Eigen::SimplicialLLT<SparseMatrixType> LltSolverType;
#endif
        
#ifdef UMFPACK_H // SuiteSparse Umfpack (LU) module
    typedef Eigen::UmfPackLU<SparseMatrixType> LuSolverType;
#else
    typedef Eigen::SparseLU<SparseMatrixType> LuSolverType;
#endif
                
        typedef FixedDirichletSolver<LltSolverType,Eigen::ConjugateGradient<SparseMatrixType>> MobileSolverType;
        typedef FixedDirichletSolver<LuSolverType,Eigen::BiCGSTAB<SparseMatrixType>> MobileReactionSolverType;

        const DislocationDynamicsBase<dim>& ddBase;
        const ClusterDynamicsParameters<dim>& cdp;
        const InvDscaling<dim> iDs;
        
        const Eigen::Matrix<double,mSize,mSize> invTrD;
        
        MobileTrialType mobileClusters;
        MobileGradType mobileGrad;
        MobileFluxType mobileFlux;
        ImmobileTrialType immobileClusters;

        const int nodeListInternalExternal;
        MobileIncrementTrialType mobileClustersIncrement;
        VolumeIntegrationDomainType dV;
        MobileBilinearWeakFormType mBWF;
        MobileIncrementBilinearWeakFormType dmBWF;

        MobileSolverType mSolver;
        bool solverInitialized;

        const Eigen::VectorXd cascadeGlobalProduction;

        ClusterDynamicsFEM(const DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
        void solveMobileClusters();
        void solveImmobileClusters();
        void solve();
//        void applyBoundaryConditions();
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
        void initializeSolver();
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

    };
    
}
#endif

