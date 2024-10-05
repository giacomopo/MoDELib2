/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_FieldBase_h
#define _model_FieldBase_h

#include <Eigen/Dense>
//#include <SpatialCellObserver.h>


namespace model
{
    

    template<typename _Scalar, int _rows, int _cols>
    struct FieldBase
    {
        typedef _Scalar Scalar;
        enum{rows=_rows};
        enum{cols=_cols};
        
        typedef Eigen::Matrix<Scalar,rows,cols> MatrixType;

        
        template <typename ParticleType>
        static MatrixType addSourceContribution(const ParticleType&)
        {
            return MatrixType::Zero();
        }
        
        
    };
} // end namespace
#endif
