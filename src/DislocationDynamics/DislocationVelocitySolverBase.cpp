/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationVelocitySolverBase_cpp_
#define model_DislocationVelocitySolverBase_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationVelocitySolverBase.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    DislocationVelocitySolverBase<DislocationNetworkType>::DislocationVelocitySolverBase(const DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    {
        
    }
    
    template struct DislocationVelocitySolverBase<DislocationNetwork<3,0>>;

}
#endif
