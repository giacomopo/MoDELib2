/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LowerCorner_H_
#define model_LowerCorner_H_

#include <Eigen/Dense>

namespace model
{
    class LowerCorner
    {
        
        const Eigen::MatrixXd xMin;
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        LowerCorner(const FiniteElementType& fe) :
        xMin(fe.xMin())
        {

        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return (node.P0-xMin).norm()<FLT_EPSILON;
            
        }
        
    };
    
    
}	// close namespace
#endif

