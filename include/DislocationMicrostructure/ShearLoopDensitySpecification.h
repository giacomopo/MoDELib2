/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopDensitySpecification_H_
#define model_ShearLoopDensitySpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct ShearLoopDensitySpecification : public MicrostructureSpecificationBase
    {
        double targetDensity;
        int numberOfSides;
        double radiusDistributionMean;
        double radiusDistributionStd;
        
        ShearLoopDensitySpecification();
        ShearLoopDensitySpecification(const std::string& fileName);
    };
}
#endif
