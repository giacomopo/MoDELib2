/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsBase_cpp_
#define model_ClusterDynamicsBase_cpp_

#include <cmath>
#include <ClusterDynamicsBase.h>

namespace model
{

    template<int dim>
    FluxMatrix<dim>::FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */cdp(cdp_in)
    {
        
    }
        
    template<int dim>
    const typename FluxMatrix<dim>::MatrixType FluxMatrix<dim>::operator() (const ElementType& elem, const BaryType&) const
    {
        const size_t grainID(elem.simplex.region->regionID);
        MatrixType D(MatrixType::Zero());
        for(size_t k=0;k<mSize;++k)
        {
            D.template block<dim,dim>(k*dim,k*dim)=-cdp.D.at(grainID)[k]/cdp.omega;
        }
        return D;
    }

template struct FluxMatrix<3>;


    template<int dim>
    ClusterDynamicsBase<dim>::ClusterDynamicsBase(const DislocationDynamicsBase<dim>& ddBase_in) :
    /* init */ ddBase(ddBase_in)
    /* init */,cdp(ddBase_in)
    /* init */,mobileClusters(ddBase.fe->template trial<'m',mSize>())
    /* init */,mobileGrad(grad(mobileClusters))
    /* init */,mobileFlux(FluxMatrix<dim>(cdp)*mobileGrad)
    /* init */,immobileClusters(ddBase.fe->template trial<'i',iSize>())
    {
        mobileClusters.setConstant(cdp.equilibriumMobileConcentration(0.0).matrix().transpose());
        immobileClusters.setConstant(Eigen::Matrix<double,iSize,1>::Zero());
    }
    
template struct ClusterDynamicsBase<3>;

}
#endif

