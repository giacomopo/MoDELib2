/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleIndividualSpecification_cpp_
#define model_PeriodicDipoleIndividualSpecification_cpp_

#include <PeriodicDipoleIndividualSpecification.h>

namespace model
{
    PeriodicDipoleIndividualSpecification::PeriodicDipoleIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("PeriodicDipole","Individual")
    {
        
    }

    PeriodicDipoleIndividualSpecification::PeriodicDipoleIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("PeriodicDipole","Individual",fileName)
    /* init */,slipSystemIDs(this->parser->readArray<int>("slipSystemIDs",true))
    /* init */,exitFaceIDs(this->parser->readArray<int>("exitFaceIDs",true))
    /* init */,dipoleCenters(this->parser->readMatrix<double>("dipoleCenters",slipSystemIDs.size(),3,true))
    /* init */,dipoleHeights(this->parser->readArray<double>("dipoleHeights",true))
    /* init */,nodesPerLine(this->parser->readArray<int>("nodesPerLine",true))
    /* init */,glideSteps(this->parser->readArray<double>("glideSteps",true))
    {
        
    }
}
#endif
