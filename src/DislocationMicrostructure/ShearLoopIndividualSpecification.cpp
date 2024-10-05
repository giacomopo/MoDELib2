/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopIndividualSpecification_cpp_
#define model_ShearLoopIndividualSpecification_cpp_

#include <ShearLoopIndividualSpecification.h>

namespace model
{
    ShearLoopIndividualSpecification::ShearLoopIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("ShearLoop","Individual")
    {
        
    }

    ShearLoopIndividualSpecification::ShearLoopIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ShearLoop","Individual",fileName)
    {
        slipSystemIDs=this->parser->readArray<int>("slipSystemIDs",true);
        if(slipSystemIDs.size())
        {
            loopRadii=this->parser->readArray<double>("loopRadii_SI",true);
            loopCenters=this->parser->readMatrix<double>("loopCenters",slipSystemIDs.size(),3,true);
            loopSides=this->parser->readArray<int>("loopSides",true);
        }
    }
}
#endif
