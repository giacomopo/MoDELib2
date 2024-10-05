/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EntireDomain_H_
#define model_EntireDomain_H_

//#include <deque>
//#include <utility>      // std::pair, std::make_pair
//#include <Eigen/Dense>
//#include <Simplex.h>
#include <IntegrationDomain.h>

namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
	struct EntireDomain
    {

        template <typename FiniteElementType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
        static IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> domain(const FiniteElementType& fe)
        {
            IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> temp;

            for (const auto& ele : fe.elements())
            {
                temp.emplace_back(&ele.second);
            }
            
            return temp;
        }
        
    };

}	// close namespace
#endif
