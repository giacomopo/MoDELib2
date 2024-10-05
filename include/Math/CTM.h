/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CTM_H_
#define model_CTM_H_

//#include <Factorial.h>

namespace model
{
	
    /*\brief A structure which wraps compile-time math functions
     */
    struct CTM
    {
        static constexpr int binomial(int N, int k)
        {/*\brief A compile-time implementation of binomial coefficients
          */
            return factorial(N) / (factorial(N-k) * factorial(k));
        }
    
        static constexpr int factorial(int n)
        {/*!\brief A compile-time implementation of the factorial function         */
            return n <= 1? 1 : (n * factorial(n - 1));
        }
        

        static constexpr size_t gcd(size_t A, size_t B)
        {/*!\bief A compile-time calculation of the Greater Common Divisors (gcd).
                  */
            return B==0? A : gcd(B, A % B );
        }
        
        /*!\brief A compile-time implementation of the power function.
         */
        static constexpr int pow(int x, int n)
        {
            return n==0? 1 : x * pow(x,n-1);
        }
        
    };
    
}

#endif
