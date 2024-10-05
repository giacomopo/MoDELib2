/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleIndividualSpecification_H_
#define model_PeriodicDipoleIndividualSpecification_H_

#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{

struct PeriodicDipoleIndividualSpecification : public MicrostructureSpecificationBase
{
    
    
    std::vector<int> slipSystemIDs;
    std::vector<int> exitFaceIDs;
    Eigen::Matrix<double,Eigen::Dynamic,3> dipoleCenters;
    std::vector<double> dipoleHeights;
    std::vector<int> nodesPerLine;
    std::vector<double> glideSteps;

    
    PeriodicDipoleIndividualSpecification();
    PeriodicDipoleIndividualSpecification(const std::string& fileName);
};

}
#endif
