/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopIndividualSpecification_cpp_
#define model_PrismaticLoopIndividualSpecification_cpp_

#include <PrismaticLoopIndividualSpecification.h>

namespace model
{
    PrismaticLoopIndividualSpecification::PrismaticLoopIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("PrismaticLoop","Individual")
    {
        
    }

    PrismaticLoopIndividualSpecification::PrismaticLoopIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("PrismaticLoop","Individual",fileName)
    /* init */,slipSystemIDs(this->parser->readArray<int>("slipSystemIDs",true))
    /* init */,loopRadii(this->parser->readArray<double>("loopRadii_SI",true))
    /* init */,loopCenters(this->parser->readMatrix<double>("loopCenters",slipSystemIDs.size(),3,true))
    /* init */,glideSteps(this->parser->readArray<double>("glideSteps",true))
    {
        
    }
}
#endif
