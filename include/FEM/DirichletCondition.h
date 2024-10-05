/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DirichletCondition_H_
#define model_DirichletCondition_H_

#include <Eigen/Dense>

namespace model
{

    struct DirichletCondition
    {
        
        const double value;
        
        DirichletCondition(const double& value_in) :
        /* init list */ value(value_in)
        {
        
        }

        /**************************************/
        template <typename NodeType,int dofPerNode>
        Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType&,
                                                       Eigen::Matrix<double,dofPerNode,1>& val) const
        {
            val.setConstant(value);
            return val;
        }
        
//        /**************************************/
//        template <typename NodeType,typename TrialFunctionType>
//        Eigen::Matrix<typename TrialFunctionType::Scalar,TrialFunctionType::dofPerNode,1>
//        operator()(const NodeType&, const TrialFunctionType&) const
//        {
//            return Eigen::Matrix<typename TrialFunctionType::Scalar,TrialFunctionType::dofPerNode,1>::Constant(value);
//        }
        
    };
    
    
}	// close namespace
#endif

