/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FEMbaseEvaluation_H_
#define model_FEMbaseEvaluation_H_

#include <Eigen/Dense>

namespace model
{
    
 
    /******************************************************************************/
    template <typename ElementType,int rows,int cols>
    struct FEMbaseEvaluation : public Eigen::Matrix<double,rows,cols>
    {
        
        typedef Eigen::Matrix<double,ElementType::dim,1> VectorDim;
        
        const VectorDim P;
        
        /**********************************************************************/
        FEMbaseEvaluation(const VectorDim& _P) :
        /* init */ Eigen::Matrix<double,rows,cols>(Eigen::Matrix<double,rows,cols>::Zero()),
        /* init */ P(_P)
        {
        }
        
    };

}
#endif
