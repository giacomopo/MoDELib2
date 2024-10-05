/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineSimplexIntersection_H_
#define model_LineSimplexIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    class LineSimplexIntersection
    {
        
        /*\todo MAKE THIS CLASS GENERAL FOR SIMPLICES.
         */
        
        typedef Eigen::Matrix<double,dim,1>     VectorDim;
        typedef Eigen::Matrix<double,dim,dim>   MatrixDim;
        
    public:
        
        /**********************************************************************/
        static int lineTriangleIntersection(const VectorDim& X,const VectorDim& s,
                                            const VectorDim& P1,const VectorDim& P2,const VectorDim& P3);
        
    };
    
}
#endif
