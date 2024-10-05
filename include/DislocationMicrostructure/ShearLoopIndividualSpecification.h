/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopIndividualSpecification_H_
#define model_ShearLoopIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct ShearLoopIndividualSpecification : public MicrostructureSpecificationBase
    {
        std::vector<int> slipSystemIDs;
        std::vector<double> loopRadii;
        Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters;
        std::vector<int> loopSides;
        
        ShearLoopIndividualSpecification();
        ShearLoopIndividualSpecification(const std::string& fileName);
    };
}
#endif
