/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TensorToVoigt_H_
#define model_TensorToVoigt_H_

#include <Eigen/Dense>

namespace model
{
    template <int rows,int cols>
    struct TensorToVoigt
    {
        
        static int t2v(const int& i,const int& j)
        {
            assert(i>=0 && i<rows);
            assert(j>=0 && j<cols);
            return cols*i+j;
        }
        
    };
}
#endif




