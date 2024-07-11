/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionMicrostructureBase_H_
#define model_InclusionMicrostructureBase_H_

#include <map>

#include <DislocationDynamicsBase.h>
#include <EshelbyInclusionBase.h>
#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>


namespace model
{

    template<int dim>
    struct InclusionMicrostructureBase : public std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>>
    /*                                */,public std::map<size_t,PolyhedronInclusionNodeIO<dim>>

    {
        typedef typename DislocationDynamicsBase<dim>::VectorDim VectorDim;
        typedef typename DislocationDynamicsBase<dim>::MatrixDim MatrixDim;
        typedef std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>> EshelbyInclusionContainerType;
        typedef std::map<size_t,PolyhedronInclusionNodeIO<dim>> PolyhedronInclusionNodeContainerType;


        DislocationDynamicsBase<dim>& ddBase;
        InclusionMicrostructureBase(DislocationDynamicsBase<dim>& ddBase_in);
        
        const EshelbyInclusionContainerType& eshelbyInclusions() const;
        EshelbyInclusionContainerType& eshelbyInclusions();
        const PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes() const;
        PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes();

    };
    
}
#endif

