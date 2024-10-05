/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_JGNselector_H_
#define model_JGNselector_H_

#include <Eigen/Dense>

namespace model
{
    template <int cols>
    struct JGNselector
    {
		template <int dim>
        static Eigen::Matrix<double,dim,1> jGN(const Eigen::Matrix<double,dim,1>& jgn)
        {
            static_assert(cols==dim,"DerivedEval MUST HAVE EXACTLY either one or dim COLUMNS");
            return jgn;
		}
        
    };
    
    template <>
    struct JGNselector<1>
    {
		template <int dim>
        static double jGN(const Eigen::Matrix<double,dim,1>& jgn)
        {
            return jgn.norm();
		}
        
    };
    
}	// close namespace
#endif
