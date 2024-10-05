/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StackingFaultTetrahedraDensitySpecification_H_
#define model_StackingFaultTetrahedraDensitySpecification_H_

#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{

struct StackingFaultTetrahedraDensitySpecification : public MicrostructureSpecificationBase
{
    double targetDensity;
    double sizeDistributionMean;
    double sizeDistributionStd;

    StackingFaultTetrahedraDensitySpecification();
    StackingFaultTetrahedraDensitySpecification(const std::string& fileName);
};

}
#endif
