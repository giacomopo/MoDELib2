/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StackingFaultTetrahedraIndividualSpecification_cpp_
#define model_StackingFaultTetrahedraIndividualSpecification_cpp_

#include <StackingFaultTetrahedraIndividualSpecification.h>

namespace model
{
    StackingFaultTetrahedraIndividualSpecification::StackingFaultTetrahedraIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("StackingFaultTetrahedra","Individual")
    {
        
    }

    StackingFaultTetrahedraIndividualSpecification::StackingFaultTetrahedraIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("StackingFaultTetrahedra","Individual",fileName)
    /* init */,planeIDs(this->parser->readArray<int>("planeIDs",true))
    /* init */,areInverted(this->parser->readArray<int>("areInverted",true))
    /* init */,sizes(this->parser->readArray<double>("sizes",true))
    /* init */,basePoints(this->parser->readMatrix<double>("basePoints",planeIDs.size(),3,true))
    {
        if(int(planeIDs.size())!=basePoints.rows())
        {
            std::cout<<"planeIDs.size()="<<planeIDs.size()<<std::endl;
            std::cout<<"basePoints.rows()="<<basePoints.rows()<<std::endl;
            throw std::runtime_error("You must provide one point for each SFT. Each point is a row of the matrix basePoints.");
        }
        
        if(planeIDs.size()!=sizes.size())
        {
            std::cout<<"planeIDs.size()="<<planeIDs.size()<<std::endl;
            std::cout<<"sizes.size()="<<sizes.size()<<std::endl;
            throw std::runtime_error("You must provide one size for each SFT. Each size is an element of the vector sizes.");
        }
        
        if(planeIDs.size()!=areInverted.size())
        {
            std::cout<<"planeIDs.size()="<<planeIDs.size()<<std::endl;
            std::cout<<"areInverted()="<<sizes.size()<<std::endl;
            throw std::runtime_error("You must provide one boolean value of areInverted for each SFT. Each value is an element of the vector areInverted.");
        }
    }
}
#endif
