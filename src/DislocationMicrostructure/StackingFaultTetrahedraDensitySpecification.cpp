/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StackingFaultTetrahedraDensitySpecification_cpp_
#define model_StackingFaultTetrahedraDensitySpecification_cpp_

#include <StackingFaultTetrahedraDensitySpecification.h>

namespace model
{
    StackingFaultTetrahedraDensitySpecification::StackingFaultTetrahedraDensitySpecification():
    /* init */ MicrostructureSpecificationBase("StackingFaultTetrahedra","Density")
    /* init */,targetDensity(0.0)
    /* init */,sizeDistributionMean(0.0)
    /* init */,sizeDistributionStd(0.0)
    {
        
    }

    StackingFaultTetrahedraDensitySpecification::StackingFaultTetrahedraDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("StackingFaultTetrahedra","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity_SI",true))
    /* init */,sizeDistributionMean(targetDensity>0.0? this->parser->readScalar<double>("sizeDistributionMean_SI",true) : 0.0)
    /* init */,sizeDistributionStd(targetDensity>0.0? this->parser->readScalar<double>("sizeDistributionStd_SI",true) : 0.0)
    {
        
    }
}
#endif
