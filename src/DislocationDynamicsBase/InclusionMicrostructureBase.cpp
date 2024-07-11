/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionMicrostructureBase_cpp_
#define model_InclusionMicrostructureBase_cpp_

#include <map>

#include <InclusionMicrostructureBase.h>


namespace model
{
    template<int dim>
    InclusionMicrostructureBase<dim>::InclusionMicrostructureBase(DislocationDynamicsBase<dim>& ddBase_in) :
    /* init */ ddBase(ddBase_in)
    {
    }

    template<int dim>
    typename InclusionMicrostructureBase<dim>::EshelbyInclusionContainerType& InclusionMicrostructureBase<dim>::eshelbyInclusions()
    {
        return *this;
    }

    template<int dim>
    const typename InclusionMicrostructureBase<dim>::EshelbyInclusionContainerType& InclusionMicrostructureBase<dim>::eshelbyInclusions() const
    {
        return *this;
    }

    template<int dim>
    typename InclusionMicrostructureBase<dim>::PolyhedronInclusionNodeContainerType& InclusionMicrostructureBase<dim>::polyhedronInclusionNodes()
    {
        return *this;
    }

    template<int dim>
    const typename InclusionMicrostructureBase<dim>::PolyhedronInclusionNodeContainerType& InclusionMicrostructureBase<dim>::polyhedronInclusionNodes() const
    {
        return *this;
    }

    template struct InclusionMicrostructureBase<3>;
}
#endif

