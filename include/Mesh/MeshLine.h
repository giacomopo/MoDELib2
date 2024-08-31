/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshLine_H_
#define model_MeshLine_H_


#include <cfloat>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <SimplicialMesh.h>
#include <PlanarMeshFace.h>
#include <FiniteLineSegment.h>
#include <StaticID.h>
#include <PlaneSegmentIntersection.h>
#include <PlanePlaneIntersection.h>
#include <MeshBoundarySegment.h>
#include <LineLineIntersection.h>
//#include <EmbeddedPolygonTriangulation.h>

namespace model
{
    
    
    template <int dim>
    struct MeshLine : public FiniteLineSegment<dim>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        static FiniteLineSegment<dim> getFiniteLineSegment(const SimplicialMesh<dim>& mesh,
                                                           const VectorDim& p,
                                                           const VectorDim& d);
        
//        const std::pair<int,int> regionIDs;
//        const BoundingMeshSegments<dim> meshIntersections;

        MeshLine(const SimplicialMesh<dim>& mesh,
                  const VectorDim& p,
                  const VectorDim& d);
    };
    
//    template <int dim,class T>
//    T& operator << (T& os, const MeshPlane<dim>& gp)
//    {
//        for (const auto& x : gp.meshIntersections)
//        {
//            os<<gp.sID<<" "<<*x;
//        }
//        return os;
//    }
    
}
#endif
