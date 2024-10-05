/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Fix_H_
#define model_Fix_H_

#include <Eigen/Dense>


namespace model
{


    struct Fix
    {        

        /**************************************/
        template <typename NodeType,int dofPerNode>
        Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType&,
                                                       Eigen::Matrix<double,dofPerNode,1>& val) const
        {
            val.setZero();
            return val;
        }
        
        
    };
    
    
}	// close namespace
#endif
