/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneKey_H
#define model_GlidePlaneKey_H

#include <array>
#include <Eigen/Dense>
#include <LatticeModule.h>


namespace model
{
    
    template <int dim>
    using GlidePlaneKey = LatticePlaneKey<dim>;
    
}
#endif

