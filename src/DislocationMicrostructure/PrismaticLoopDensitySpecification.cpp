/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopDensitySpecification_cpp_
#define model_PrismaticLoopDensitySpecification_cpp_

#include <PrismaticLoopDensitySpecification.h>

namespace model
{
    PrismaticLoopDensitySpecification::PrismaticLoopDensitySpecification():
    /* init */ MicrostructureSpecificationBase("PrismaticLoop","Density")
    /* init */,targetDensity(0.0)
    /* init */,radiusDistributionMean(0.0)
    /* init */,radiusDistributionStd(0.0)
    {
        
    }

    PrismaticLoopDensitySpecification::PrismaticLoopDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("PrismaticLoop","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity_SI",true))
    /* init */,radiusDistributionMean(this->parser->readScalar<double>("radiusDistributionMean_SI",true))
    /* init */,radiusDistributionStd(this->parser->readScalar<double>("radiusDistributionStd_SI",true))
    {
        
    }
}
#endif
