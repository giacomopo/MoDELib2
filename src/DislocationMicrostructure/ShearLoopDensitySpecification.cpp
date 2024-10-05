/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopDensitySpecification_cpp_
#define model_ShearLoopDensitySpecification_cpp_

#include <ShearLoopDensitySpecification.h>

namespace model
{
    ShearLoopDensitySpecification::ShearLoopDensitySpecification():
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density")
    /* init */,targetDensity(0.0)
    /* init */,numberOfSides(0)
    /* init */,radiusDistributionMean(0.0)
    /* init */,radiusDistributionStd(0.0)
    {
        
    }

    ShearLoopDensitySpecification::ShearLoopDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity",true))
    /* init */,numberOfSides(targetDensity>0.0? this->parser->readScalar<int>("numberOfSides",true) : 0)
    /* init */,radiusDistributionMean(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionMean",true) : 0.0)
    /* init */,radiusDistributionStd(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionStd",true) : 0.0)
    {
        
    }
}
#endif
