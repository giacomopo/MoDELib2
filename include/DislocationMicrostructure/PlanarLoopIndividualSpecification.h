/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarLoopIndividualSpecification_H_
#define model_PlanarLoopIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct PlanarLoopIndividualSpecification : public MicrostructureSpecificationBase
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        Eigen::Matrix<double,Eigen::Dynamic,3> loopPoints;
        VectorDim burgers;
        VectorDim normal;
        
        PlanarLoopIndividualSpecification();
        PlanarLoopIndividualSpecification(const std::string& fileName);
    };
}
#endif
