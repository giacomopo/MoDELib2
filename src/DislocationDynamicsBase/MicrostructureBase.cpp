/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureBase_cpp_
#define model_MicrostructureBase_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <map>
#include <DDconfigIO.h>
#include <DDauxIO.h>
#include <MicrostructureBase.h>
#include <MicrostructureContainer.h>

namespace model
{
    template <int dim>
    MicrostructureBase<dim>::MicrostructureBase(const std::string& tag_in,MicrostructureContainer<dim>& m_in):
    /* init */ tag(tag_in)
    /* init */,microstructures(m_in)
    {
        
    }

    template <int dim>
    double MicrostructureBase<dim>::timeSinceLastUpdate() const
    {
        return microstructures.ddBase.simulationParameters.totalTime-lastUpdateTime;
    }

template <int dim>
Eigen::Matrix<double,Eigen::Dynamic,dim> MicrostructureBase<dim>::displacement(const Eigen::Matrix<double,Eigen::Dynamic,dim>& points) const
{
    Eigen::Matrix<double,Eigen::Dynamic,dim> temp(Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(points.rows(),dim));
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(long int k=0;k<points.rows();++k)
    {
        temp.row(k)=displacement(points.row(k),nullptr,nullptr,nullptr);
    }
    return temp;
}

template <int dim>
std::set<const Grain<dim>*> MicrostructureBase<dim>::pointGrains(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
{
    std::set<const Grain<dim>*> temp;
    if(node)
    {
        for(const auto& nodeEle : *node)
        {
            temp.emplace(&microstructures.ddBase.poly.grain(nodeEle->simplex.region->regionID));
        }
    }
    else
    {
        if(ele)
        {
            temp.emplace(&microstructures.ddBase.poly.grain(ele->simplex.region->regionID));
        }
        else
        {
            const std::pair<bool,const Simplex<dim,dim>*> found(microstructures.ddBase.mesh.searchWithGuess(x,guess));
            if(found.first)
            {
                temp.emplace(&microstructures.ddBase.poly.grain(found.second->region->regionID));
            }
        }
    }
    return temp;
}

    template struct MicrostructureBase<3>;

}
#endif
