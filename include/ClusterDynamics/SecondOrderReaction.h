/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SecondOrderReaction_H_
#define model_SecondOrderReaction_H_

#include <Eigen/Dense>
#include <EvalFunction.h>
#include <ClusterDynamicsParameters.h>

namespace model
{

    template <typename TrialFunctionType>
    struct SecondOrderReaction : public EvalFunction<SecondOrderReaction<TrialFunctionType>>
    {
        
        constexpr static int rows=TrialFunctionType::rows;
        constexpr static int cols=rows;
        constexpr static int dim=TrialFunctionType::dim;
        constexpr static int mSize=rows;
        typedef EvalExpression<TrialFunctionType> EvalFunctionType;
        
        const ClusterDynamicsParameters<dim>& cdp;
        const TrialFunctionType& c;       // concentration field
        const EvalFunctionType ce;
        const Eigen::Matrix<double,rows,1> inFactors;
        const Eigen::Matrix<double,rows,1> outFactors;
        
        /**********************************************************************/
        SecondOrderReaction(const TrialFunctionType& c_in, const ClusterDynamicsParameters<dim>& cdp_in) :
        /* init */ cdp(cdp_in),
        /* init */ c(c_in),
        /* init */ ce(c_in),
        /* init */ inFactors((1.0/cdp.msVector.abs()).matrix()),
        /* init */ outFactors(cdp.msVector.abs().matrix())
        {
            
        }
        
        template<typename ElementType, typename BaryType>
        const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the matrix having in row k the product R_2k*c
          *
          * Note that the coefficients in cdp.R2 are meant to multiply cluster concentrations, not point defect concentrations.
          * Therefore, we first convert ce(ele,bary) to cluster concentrations using the inFactors, and then re-convert the result to
          * point defect concentrations using the outFactors.
          */
            const Eigen::Matrix<double,rows,1> cvalue(inFactors.asDiagonal()*ce(ele,bary));
            Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());
            for(size_t r=0;r<rows;++r)
            {
                temp.row(r)=cdp.R2[r]*cvalue;
            }
            return  outFactors.asDiagonal()*(temp)*inFactors.asDiagonal();
        }
        
    };

}
#endif

