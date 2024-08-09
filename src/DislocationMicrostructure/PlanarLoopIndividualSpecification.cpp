/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
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
    /* init */,burgers(this->parser->readMatrix<double>("burgers",1,3,true))
    /* init */,normal(this->parser->readMatrix<double>("normal",1,3,true))
    {
        
    }
}
#endif
