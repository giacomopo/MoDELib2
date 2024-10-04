/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshLine_cpp_
#define model_MeshLine_cpp_

#include <numbers>
#include <cfloat>
#include <map>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

#include <MeshModule.h>
#include <MeshRegion.h>
#include <MeshRegionObserver.h>
#include <PlaneLineIntersection.h>
#include <Polygon2D.h>

namespace model
{
    
template <int dim>
FiniteLineSegment<dim> MeshLine<dim>::getFiniteLineSegment(const SimplicialMesh<dim>& mesh,
                                                   const VectorDim& lineP,
                                                   const VectorDim& lineD)
{
    const double dNorm(lineD.norm());
    if(dNorm<FLT_EPSILON)
    {
        throw std::runtime_error("Line direction has zero norm");
    }
    const VectorDim unitD(lineD/dNorm);
    std::map<float,VectorDim> rootsMap;
    for(const auto& region : mesh.regions())
    {
        for(const auto& face : region.second->faces())
        {
            const auto plane(face.second->asPlane());
            PlaneLineIntersection<dim> pli(plane.P,plane.unitNormal,lineP,unitD);
            if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
            {
                
                const auto& faceHull(face.second->convexHull());
                std::vector<Eigen::Matrix<double,dim-1,1>> localHull;
                for(const auto& vtx : faceHull)
                {
                    localHull.push_back(plane.localPosition(vtx->P0));
                }
                const int wn(Polygon2D::windingNumber(plane.localPosition(pli.P),localHull));
                if(wn!=0)
                {
                    const double u((pli.P-lineP).dot(unitD));
                    rootsMap.emplace(u,pli.P);
                }
            }
        }
    }
    if(rootsMap.size()>=2)
    {
        return FiniteLineSegment<dim>(rootsMap.begin()->second,rootsMap.rbegin()->second);
    }
    else
    {
        throw std::runtime_error("FiniteLineSegment::getFiniteLineSegment failed ("+std::to_string(rootsMap.size())+" intersections)");
        return FiniteLineSegment<dim>(VectorDim::Zero(),VectorDim::Zero());
    }
}

//        const std::pair<int,int> regionIDs;
//        const BoundingMeshSegments<dim> meshIntersections;
template <int dim>
MeshLine<dim>::MeshLine(const SimplicialMesh<dim>& mesh,
                        const VectorDim& lineP,
                        const VectorDim& lineD):
/* init */ FiniteLineSegment<dim>(getFiniteLineSegment(mesh,lineP,lineD))
{
    
}
    
    template struct MeshLine<3>;
//    template struct MeshLine<2>;

}
#endif
