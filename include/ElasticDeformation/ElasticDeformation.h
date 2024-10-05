/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformation_H_
#define model_ElasticDeformation_H_

#ifdef MODELIB_CHOLMOD // SuiteSparse Cholmod module
#include <Eigen/CholmodSupport>
#endif

#include <MicrostructureContainer.h>
#include <ElasticDeformationFEM.h>
#include <IntegrationDomain.h>
#include <GaussLegendre.h>
#include <FEMfaceEvaluation.h>
#include <UniformController.h>

namespace model
{
    template<int dim>
    struct ElasticDeformation : public MicrostructureBase<dim>
{
    typedef typename MicrostructureBase<dim>::ElementType ElementType;
    typedef typename MicrostructureBase<dim>::FiniteElementType FiniteElementType;
    typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
    typedef typename MicrostructureBase<dim>::NodeType NodeType;
    typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
    typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
    typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
    typedef MicrostructureContainer<dim> MicrostructureContainerType;
    static constexpr int imageTractionIntegrationOrder=3;
    typedef IntegrationDomain<FiniteElementType,1,imageTractionIntegrationOrder,GaussLegendre> TractionIntegrationDomainType;
    typedef IntegrationList<1,FEMfaceEvaluation<ElementType,dim,dim>> TractionIntegrationListType;
    typedef UniformController<SymmetricVoigtTraits<dim>::voigtSize> UniformControllerType;
    
#ifdef CHOLMOD_H // SuiteSparse Cholmod module
    typedef Eigen::CholmodSupernodalLLT<typename ElasticDeformationFEM<dim>::SparseMatrixType> SPDsolverType;
#else
    typedef Eigen::SimplicialLLT<typename ElasticDeformationFEM<dim>::SparseMatrixType> SPDsolverType;
#endif

    
    const bool useElasticDeformationFEM;
    const std::unique_ptr<ElasticDeformationFEM<dim>> elasticDeformationFEM;
    const std::unique_ptr<UniformControllerType> uniformLoadController;
    const double inertiaReliefPenaltyFactor;
    SPDsolverType directSolver;
    TractionIntegrationDomainType ndA;
    TractionIntegrationListType tractionList;
    Eigen::VectorXd zDot;
    bool solverInitialized;
    
    ElasticDeformation(MicrostructureContainerType& mc);
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
    void initializeSolver();
};
}
#endif

