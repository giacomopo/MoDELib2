/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsIndividualSpecification_cpp_
#define model_FrankLoopsIndividualSpecification_cpp_

#include <FrankLoopsIndividualSpecification.h>

namespace model
{
    FrankLoopsIndividualSpecification::FrankLoopsIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("FrankLoops","Individual")
    {
        
    }

    FrankLoopsIndividualSpecification::FrankLoopsIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("FrankLoops","Individual",fileName)
    {
        planeIDs=this->parser->readArray<int>("planeIDs",true);
        if(planeIDs.size())
        {
            loopRadii=this->parser->readArray<double>("loopRadii_SI",true);
            loopCenters=this->parser->readMatrix<double>("loopCenters",planeIDs.size(),3,true);
            loopSides=this->parser->readArray<int>("loopSides",true);
            isVacancyLoop=this->parser->readArray<int>("isVacancyLoop",true);
        }
    }
}
#endif
