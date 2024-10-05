/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionIndividualSpecification_H_
#define model_SphericalInclusionIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct SphericalInclusionIndividualSpecification : public MicrostructureSpecificationBase
    {
        
        static constexpr int dim=3;

        std::vector<double> radii_SI;
        Eigen::Matrix<double,Eigen::Dynamic,dim> centers;
        Eigen::Matrix<double,Eigen::Dynamic,dim*dim> eigenDistortions;
        std::vector<double> velocityReductionFactors;
        std::vector<int> phaseIDs;
        bool allowOverlap;
        bool allowOutside;


        SphericalInclusionIndividualSpecification();
        SphericalInclusionIndividualSpecification(const std::string& fileName);
    };
}
#endif
