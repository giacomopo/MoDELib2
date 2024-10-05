/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleDensitySpecification_cpp_
#define model_PeriodicDipoleDensitySpecification_cpp_

#include <PeriodicDipoleDensitySpecification.h>

namespace model
{
    PeriodicDipoleDensitySpecification::PeriodicDipoleDensitySpecification():
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density")
    /* init */,targetDensity(0.0)
    {
        
    }

    PeriodicDipoleDensitySpecification::PeriodicDipoleDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity",true))
    {
        
    }
}
#endif
