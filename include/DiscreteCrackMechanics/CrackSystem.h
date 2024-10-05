/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackSystem_H_
#define model_CrackSystem_H_

#include <Eigen/Dense>

namespace model
{
    template <int dim>
    class CrackSystem
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
    public:
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim&) const
        {/*!\param[in] P position vector
          * \returns The displacement field in the DefectiveCrystal at P
          */

            return VectorDim::Zero();
        }
        
        /**********************************************************************/
        template<typename ElementType>
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>&) const
        {

        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim&) const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */

            return MatrixDim::Zero();
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */

            return MatrixDim::Zero();
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            return MatrixDim::Zero();
        }
        
    };
}
#endif
