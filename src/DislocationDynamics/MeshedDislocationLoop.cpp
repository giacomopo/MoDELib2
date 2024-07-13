/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshedDislocationLoop_cpp_
#define model_MeshedDislocationLoop_cpp_


#include <MeshedDislocationLoop.h>

namespace model
{
    
MeshedDislocationLoop::MeshedDislocationLoop(const Plane<3>& plane_in,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize):
/* init */ plane(plane_in)
{
    
    std::deque<Eigen::Matrix<double,2,1>> localBndPts;
    std::deque<Eigen::Matrix<double,2,1>> internalPts;
    for(const auto& gpt : globalBndPts)
    {
        localBndPts.push_back(plane.localPosition(gpt));
    }
    this->reMesh(localBndPts,internalPts,meshSize);

    
    for(const auto& v : this->vertices())
    {
        points.push_back(plane.globalPosition(v));
    }
}

}
#endif
