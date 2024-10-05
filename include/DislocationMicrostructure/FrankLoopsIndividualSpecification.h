/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsIndividualSpecification_H_
#define model_FrankLoopsIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct FrankLoopsIndividualSpecification : public MicrostructureSpecificationBase
    {
        std::vector<int> planeIDs;
        std::vector<double> loopRadii;
        Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters;
        std::vector<int> loopSides;
        std::vector<int> isVacancyLoop;
        
        FrankLoopsIndividualSpecification();
        FrankLoopsIndividualSpecification(const std::string& fileName);
    };
}
#endif
