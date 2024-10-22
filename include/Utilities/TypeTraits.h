/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_TYPETRAITS_H_
#define model_TYPETRAITS_H_

#include <assert.h>

namespace model
{
	
	template <typename T>
	struct TypeTraits
    {
		
		/*! \brief A helper struct for predefinition of types in model classes. 
		 */
		
		TypeTraits()
        {
			assert(0 && "YOU MUST DEFINE A TypeTraits SPECIALIZATION FOR YOUR CLASS.");
		}
	};
	
}
#endif


