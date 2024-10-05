/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneSegmentIntersection_H_
#define model_PlaneSegmentIntersection_H_

#include <tuple>
#include <Eigen/Dense>
#include <Plane.h>
#include <FiniteLineSegment.h>

namespace model
{
    
    template <int dim>
    struct PlaneSegmentIntersection
    {
        
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT,OFFSET};
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findIntersection(const VectorDimD& P0,
                                             const VectorDimD& Normal,
                                             const VectorDimD& v0,
                                             const VectorDimD& v1);
        
        
    public:
        
        const IntersectionType& type;
        const VectorDimD x0;
        const VectorDimD x1;
        
        /**********************************************************************/
        PlaneSegmentIntersection(const VectorDimD& P0,
                                 const VectorDimD& N,
                                 const VectorDimD& v0,
                                 const VectorDimD& v1);
        
        /**********************************************************************/
        PlaneSegmentIntersection(const Plane<dim>& plane,
                                 const FiniteLineSegment<dim>& seg);
        
    };

}
#endif
