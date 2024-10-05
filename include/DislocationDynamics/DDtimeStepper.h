/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeStepper_H_
#define model_DDtimeStepper_H_

#include <chrono>

#include <Eigen/Dense>
#include <TextFileParser.h>
#include <DislocationNetwork.h>

namespace model
{
	
    template <typename DislocationNetworkType>
    struct DDtimeStepper
    {
//        static constexpr auto tag="vMax integrator";
        const DislocationNetworkType& DN;
        const std::string timeSteppingMethod;
        const double dxMax;
//        const double shearWaveSpeedFraction;
//        const double dtMax;
        
        DDtimeStepper(const DislocationNetworkType& DN_in);
        double getDt(const double& vRef, const double& vFactor) const;
    };

} // end namespace
#endif

