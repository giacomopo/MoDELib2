/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsBase_H_
#define model_ClusterDynamicsBase_H_


#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <EvalFunction.h>

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

    template<int dim>
    struct ClusterDynamicsBase
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;
        typedef TrialFunction<'m',mSize,FiniteElementType> MobileTrialType;
        typedef TrialFunction<'i',iSize,FiniteElementType> ImmobileTrialType;
//        typedef TrialFunction<'z',dim,FiniteElementType> DiffusiveTrialType;
        typedef TrialGrad<MobileTrialType> MobileGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileGradType> MobileFluxType;


        const DislocationDynamicsBase<dim>& ddBase;
        const ClusterDynamicsParameters<dim> cdp;
        
        MobileTrialType mobileClusters;
        MobileGradType mobileGrad;
        MobileFluxType mobileFlux;


        ImmobileTrialType immobileClusters;
                


        ClusterDynamicsBase(const DislocationDynamicsBase<dim>& ddBase_in);
        
    };
    
}
#endif

