/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamics_H_
#define model_ClusterDynamics_H_

#include <MicrostructureContainer.h>
#include <ClusterDynamicsParameters.h>
#include <ClusterDynamicsFEM.h>
//#include <SimpleNullSpaceSolver.h>
#include <UniformController.h>
//#include <SecondOrderReaction.h>

namespace model
{

    template<int dim>
    struct ClusterDynamics : public MicrostructureBase<dim>
//    /*                   */, public ClusterDynamicsFEM<dim>
    {
        
//        typedef Eigen::SparseMatrix<double> SparseMatrixType;
//    #ifdef _MODEL_PARDISO_SOLVER_
//        typedef Eigen::PardisoLLT<SparseMatrixType> DirectSPDSolverType;
//        typedef Eigen::PardisoLU<SparseMatrixType> DirectSquareSolverType;
//    #else
//        typedef Eigen::SimplicialLLT<SparseMatrixType> DirectSPDSolverType;
//        typedef Eigen::SparseLU<SparseMatrixType> DirectSquareSolverType;
//    #endif
//        typedef Eigen::ConjugateGradient<SparseMatrixType> IterativeSPDSolverType;
//        typedef Eigen::BiCGSTAB<SparseMatrixType> IterativeSquareSolverType;

        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef MicrostructureContainer<dim> MicrostructureContainerType;

        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;

        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;


//        typedef typename ClusterDynamicsFEM<dim>::FiniteElementType FiniteElementType;
//        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
//        static constexpr int dVorder=4;
//        typedef IntegrationDomain<FiniteElementType,0,dVorder,GaussLegendre> VolumeIntegrationDomainType;
        
//        typedef typename ClusterDynamicsFEM<dim>::MobileTrialType MobileTrialType;
//        typedef typename ClusterDynamicsFEM<dim>::MobileGradType MobileGradType;
//        typedef typename ClusterDynamicsFEM<dim>::MobileFluxType MobileFluxType;
//        typedef BilinearForm<MobileGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;
//        typedef BilinearWeakForm<MobileBilinearFormType,VolumeIntegrationDomainType> MobileBilinearWeakFormType;
//
//        typedef TrialFunction<'d',mSize,FiniteElementType> MobileIncrementTrialType;
//        typedef TrialGrad<MobileIncrementTrialType> MobileIncrementGradType;
//        typedef TrialProd<FluxMatrix<dim>,MobileIncrementGradType> MobileIncrementFluxType;
//        typedef BilinearForm<MobileIncrementGradType,TrialProd<Constant<double,1,1>,MobileIncrementFluxType>> MobileIncrementBilinearFormType;
//        typedef BilinearWeakForm<MobileIncrementBilinearFormType,VolumeIntegrationDomainType> MobileIncrementBilinearWeakFormType;

//        typedef typename ClusterDynamicsFEM<dim>::ImmobileTrialType ImmobileTrialType;

        
        typedef UniformController<mSize*dim> UniformControllerType;
        typedef std::map<size_t,const std::unique_ptr<UniformControllerType>> UniformControllerContainerType; // one per each grain
        
        ClusterDynamics(MicrostructureContainerType& mc);
        
        const ClusterDynamicsParameters<dim> cdp;

        const bool useClusterDynamicsFEM;
        const std::unique_ptr<ClusterDynamicsFEM<dim>> clusterDynamicsFEM;
        UniformControllerContainerType uniformControllers;

        
//        const int nodeListInternalExternal;
//        MobileIncrementTrialType mobileClustersIncrement;
//        VolumeIntegrationDomainType dV;
//        MobileBilinearWeakFormType mBWF;
//        MobileIncrementBilinearWeakFormType dmBWF;
//
////        FixedDirichletSolver<MobileBilinearWeakFormType> mSolver;
//        FixedDirichletSolver mSolver;
//
//        const Eigen::VectorXd cascadeGlobalProduction;
        
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels) override;
        void solve() override;
        double getDt() const override;
        void output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const override;
        void updateConfiguration() override;
        MatrixDim averagePlasticDistortion() const override ;
        MatrixDim averagePlasticDistortionRate() const override;
        VectorDim displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
        MatrixDim stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
        MatrixDim averageStress() const override;
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        VectorMSize mobileConcentration(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        void applyBoundaryConditions();

    
        
        
        static UniformControllerContainerType getUniformControllers(const DislocationDynamicsBase<dim>& ddBase,const ClusterDynamicsParameters<dim>& cdp);

    };



    
}
#endif

