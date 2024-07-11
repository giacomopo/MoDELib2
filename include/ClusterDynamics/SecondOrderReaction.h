/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SecondOrderReaction_H_
#define model_SecondOrderReaction_H_

//#include <MicrostructureContainer.h>
//#include <ClusterDynamicsParameters.h>
//#include <ClusterDynamicsBase.h>
//#include <SimpleNullSpaceSolver.h>
//#include <FixedDirichletSolver.h>
//#include <UniformController.h>
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

    const ClusterDynamicsParameters<dim>& cdp;
//    const ODESystem& odeS;
    
//    typedef TrialFunction<'c',4,FiniteElementType> TrialFunctionType;
    typedef EvalExpression<TrialFunctionType> EvalFunctionType;
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
    
    std::vector<Eigen::Matrix<double,mSize,mSize>> buildR2(const Eigen::Matrix<double,mSize,1>& cin) const
    {
        Eigen::Array<double,1,mSize> rn(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            if(fabs(cdp.msVector(k)-1)<FLT_EPSILON || fabs(cdp.msVector(k)+1)<FLT_EPSILON)
            {
                rn(k)=pow(3.0*cdp.omega/4.0/M_PI,1.0/3.0);
            }
            else
            {
                rn(k)=pow(fabs(cdp.msVector(k)*cdp.omega/cdp.b/M_PI),1.0/2.0);
            }
        }
        
        const std::vector<Eigen::Matrix<double,3,3>> Dlocal(cdp.getDlocal());
        
        std::vector<Eigen::Matrix<double,mSize,mSize>> tempR2(cdp.R2);
        
        if(cdp.use0DsinkStrength)
        {
            for(int k=0;k<mSize;++k)
            {
                // Add 1D interaction rate
                for(const auto& pair: cdp.reactionMap)
                {
                    const auto key=pair.first;
                    if(k==key.first || k==key.second)
                    {// (key.first, key.second) consume k_specie
                        if( cdp.msVector(key.second)>1.0 )
                        { // 3D immobile + 1D mobile
                            tempR2[k](key.first ,key.second) -= pair.second*M_PI*M_PI*pow(rn(key.second),4.0)*(Dlocal[key.second](0,0))/cdp.omega/cdp.omega*cin(key.first);
                            tempR2[k](key.second,key.first)  -= pair.second*M_PI*M_PI*pow(rn(key.second),4.0)*(Dlocal[key.second](0,0))/cdp.omega/cdp.omega*cin(key.first);
                        }
                    }
                    else if( fabs(cdp.msVector(key.first)+cdp.msVector(key.second)-cdp.msVector(k))<FLT_EPSILON )
                    {// (key.first, key.second) generate k_specie
                        if( cdp.msVector(key.second)>1.0 )
                        { // 3D mobile + 1D immobile
                            const double same= (key.first==key.second) ? 0.5 : 1.0;
                            tempR2[k](key.first ,key.second) += same*pair.second*M_PI*M_PI*pow(rn(key.second),4.0)*(Dlocal[key.second](0,0))/cdp.omega/cdp.omega*cin(key.first);
                            tempR2[k](key.second,key.first)  += same*pair.second*M_PI*M_PI*pow(rn(key.second),4.0)*(Dlocal[key.second](0,0))/cdp.omega/cdp.omega*cin(key.first);
                        }
                    }
                }
            }
        }
        
        return tempR2;
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the current stiffness C, which in general is a funciton of C0 and grad(u)
      */
        // dot(c_a) = R1_ab*cb ->  t_a*dot(c_a) = t_a*R1_ab* d(cb)
        
        // Assume R2_abg=R2_agb
        // dot(c_a) = 1/2*R2_abg*cb*cg -> t_a*dot(c_a) = 1/2*t_a*(R2_abg*d(cb)*cg+R2_abg*cb*d(cg)) = t_a*(R2_abg*cg)*d(cb)
        const Eigen::Matrix<double,rows,1> cvalue(inFactors.asDiagonal()*ce(ele,bary));
        std::vector<Eigen::Matrix<double,rows,cols>> R2V(buildR2(cvalue));
        
        Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());
        for(size_t r=0;r<rows;++r)
        {
            temp.row(r)=R2V[r]*cvalue;
        }

        return  outFactors.asDiagonal()*(temp)*inFactors.asDiagonal();
    }
    
};
    
}
#endif

