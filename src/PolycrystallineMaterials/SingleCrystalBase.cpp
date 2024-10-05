/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_cpp_
#define model_SingleCrystalBase_cpp_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <iterator>

#include <SingleCrystalBase.h>

namespace model
{
    template<int dim>
    SingleCrystalBase<dim>::SingleCrystalBase(const MatrixDim& A,
                                      const MatrixDim& C2G) :
    /* init */ LatticeType(A,C2G)
    {
    }

    template<int dim>
    SingleCrystalBase<dim>::~SingleCrystalBase(){}

    template struct SingleCrystalBase<3>;

}
#endif
