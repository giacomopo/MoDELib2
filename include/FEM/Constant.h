/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Constatnt_H_
#define model_Constatnt_H_

#include <utility> // for std::move

#include <Eigen/Dense>
//#include <EvalExpression.h>
#include <EvalFunction.h>


namespace model
{
    template <typename T, int _rows, int _cols>
	struct Constant : public EvalFunction<Constant<T,_rows,_cols> >
    {
        constexpr static int rows=_rows;
        constexpr static int cols=_cols;
   
        const T c;
        Constant(const T& c_in) : c(c_in)
        {
        }
        
        template<typename ElementType, typename BaryType>
        const T& operator() (const ElementType&, const BaryType&) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the constant c.
          */
            return c;
        }
    };
    
    // Operators
    template <int rows, int cols>
    Constant<Eigen::Matrix<double,rows,cols>,rows,cols> make_constant(const Eigen::Matrix<double,rows,cols>& c)
    {
        return Constant<Eigen::Matrix<double,rows,cols>,rows,cols>(c);
    }

    static inline Constant<double,1,1> make_constant(const double& c)
    {
        return Constant<double,1,1>(c);
    }
     
}	// close namespace
#endif
