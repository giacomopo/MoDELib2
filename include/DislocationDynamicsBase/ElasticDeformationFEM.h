/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformationFEM_H_
#define model_ElasticDeformationFEM_H_


#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <VoigtTraits.h>


namespace model
{

    template<int dim>
    struct ElasticDeformationFEM
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef typename DislocationDynamicsBase<dim>::VectorDim VectorDim;
        typedef typename DislocationDynamicsBase<dim>::MatrixDim MatrixDim;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<'u',dim,FiniteElementType> TrialFunctionType;
        typedef TrialGrad<TrialFunctionType> TrialGradType;
        typedef TrialDef<TrialFunctionType> TrialDefType;
        typedef Eigen::Matrix<double,dim*(dim+1)/2,dim*(dim+1)/2> CmatrixType;
        typedef Constant<CmatrixType,dim*(dim+1)/2,dim*(dim+1)/2> CconstantType;
        typedef TrialProd<CconstantType,TrialDefType> TrialStressType;
        typedef BilinearForm<TrialDefType,TrialStressType> BilinearFormType;
        typedef IntegrationDomain<FiniteElementType,0,4,GaussLegendre> IntegrationDomainType;
        typedef BilinearWeakForm<BilinearFormType,IntegrationDomainType> BilinearWeakFormType;
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        
        typedef Quadrature<dim,4,GaussLegendre> MomentumQuadratureType;
//        typedef Eigen::Matrix<double,1,ElementType::nodesPerElement> MomentumKernelMatrixType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        typedef Eigen::Matrix<double,2*dim,dim*ElementType::nodesPerElement> MomentumMatrixType;
        
        typedef TrialFunction<'z',dim,FiniteElementType> DiffusiveTrialType;


        
        MomentumMatrixType elementMomentumKernelMatrix(const VectorDim& abscissa,const ElementType& ele) const;
        MomentumMatrixType elementMomentumMatrix(const ElementType& ele) const;

        DislocationDynamicsBase<dim>& ddBase;
//        const SymmetricVoigtTraits<dim>& voigtTraits;

        CmatrixType C; // matrix of elastic moduli
        TrialFunctionType u; // displacement u=[u1 u2 u3]'
        TrialGradType  b;  // displacement gradient b=[u1,1 u1,2 u1,3 u2,1 u2,2 u2,3 u3,1 u3,2 u3,3]'
        TrialDefType e; // strain s=[e11 e22 e33 e12 e23 e13]'
        TrialStressType s;     // stress s=[s11 s22 s33 s12 s23 s13]'
        IntegrationDomainType dV; // volume element
        BilinearWeakFormType bWF; // bilinear weak form (test(e),s)*dV
        SparseMatrixType A; // stiffness matrix of bWF
                
        DiffusiveTrialType z; // diffusive displacement

        ElasticDeformationFEM(DislocationDynamicsBase<dim>& ddBase_in);
        
        CmatrixType get_C(const double& mu, const double& nu) const;
        
    };
    
}
#endif

