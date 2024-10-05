/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityHEXpyramidal_h_
#define _model_DislocationMobilityHEXpyramidal_h_

#include <DislocationMobilityBase.h>

namespace model
{
    
    struct DislocationMobilityHEXpyramidal : public DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double kB;
        
        /**********************************************************************/
        DislocationMobilityHEXpyramidal(const PolycrystallineMaterialBase& material) ;
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) override;
        
    };
    
}
#endif
