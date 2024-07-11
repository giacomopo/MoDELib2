/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
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
