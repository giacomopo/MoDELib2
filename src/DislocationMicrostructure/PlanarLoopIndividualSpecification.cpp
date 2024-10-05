/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarLoopIndividualSpecification_cpp_
#define model_PlanarLoopIndividualSpecification_cpp_

#include <PlanarLoopIndividualSpecification.h>

namespace model
{
    PlanarLoopIndividualSpecification::PlanarLoopIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("PlanarLoop","Individual")
    /* init */,burgers(VectorDim::Zero())
    /* init */,normal(VectorDim::Zero())
    {
        
    }

    PlanarLoopIndividualSpecification::PlanarLoopIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("PlanarLoop","Individual",fileName)
    /* init */,loopPoints(this->parser->readMatrixCols<double>("loopPoints",3,true))
    /* init */,burgers(loopPoints.rows()>2? this->parser->readMatrix<double>("burgers",1,3,true) : VectorDim::Zero())
    /* init */,normal(loopPoints.rows()>2? this->parser->readMatrix<double>("normal",1,3,true) : VectorDim::Zero())
    {
        
    }
}
#endif
