/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialProd_H_
#define model_TrialProd_H_

#include <utility> // for std::move
#include <TypeTraits.h>
#include <Simplex.h>
#include <TrialBase.h>
#include <TrialGrad.h>
#include <TrialExpressionBase.h>
#include <Constant.h>
//#include <EvalExpression.h>
#include <EvalFunction.h>
#include <ExpressionRef.h>


namespace model
{
    
    template<int rows1,int cols1, int rows2>
    struct TrialProdRows
    {
        static_assert(cols1==rows2,"YOU ARE MULTIPLYING TrialExpressionBaseS WITH DIFFERENT SIZE.");
        constexpr static int rows=rows1;
    };
    
    template<int rows2>
    struct TrialProdRows<1,1,rows2>
    {
        constexpr static int rows=rows2;
    };
    
    template<>
    struct TrialProdRows<1,1,1>
    {
        constexpr static int rows=1;
    };
    
    template <typename T1, typename T2>
    struct TrialProd : //public TrialBase<typename T2::TrialFunctionType>,
    /*             */ public TrialExpressionBase<TrialProd<T1,T2> >
    {
        typedef typename T2::TrialFunctionType TrialFunctionType;
        typedef TrialBase<TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=TrialProdRows<T1::rows,T1::cols,T2::rows>::rows;
        constexpr static int cols=T2::cols;
        
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows*dim,dofPerElement> ShapeFunctionGradMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;
        
        
        //        const T1& op1;  // first  operand (the Constant).
        //        const T2& op2;
        
        ExpressionRef<T1> op1;
        ExpressionRef<T2> op2;

        TrialProd(const T1& x,
                  const T2& y) :
        /* init list           */ op1(x),
        /* init list           */ op2(y)
        {/*!
          */
        }
        
        /**********************************************************************/
        TrialProd(T1&& x,
                  const T2& y) :
        //        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(std::move(x)),
        /* init list           */ op2(y)
        {/*!
          */
        }
        
        /**********************************************************************/
        TrialProd(const T1& x,
                  T2&& y) :
        //        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(x),
        /* init list           */ op2(std::move(y))
        {/*!
          */
        }
        
        /**********************************************************************/
        TrialProd(T1&& x,
                  T2&& y) :
        //        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(std::move(x)),
        /* init list           */ op2(std::move(y))
        {/*!
          */
        }
        
        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                    const BaryType& bary) const
        {
            return op1()(ele,bary)*op2().sfm(ele,bary);
        }
        
        ShapeFunctionGradMatrixType sfmGrad(const ElementType& ele,
                                                                       const BaryType& bary) const
        {
            const auto t1(op1()(ele,bary));
            Eigen::Matrix<double,T1::rows*dim,TrialGrad<T2>::rows> t1Ex(Eigen::Matrix<double,T1::rows*dim,TrialGrad<T2>::rows>::Zero());
            for(int I=0;I<T1::rows;++I)
            {
                for(int J=0;J<T1::cols;++J)
                {
                    for(int d=0;d<dim;++d)
                    {
                        const int i(I*dim+d);
                        const int j(J*dim+d);
                        t1Ex(i,j)=t1(I,J);
                    }
                }
            }
            return t1Ex*op2().sfmGrad(ele,bary);
        }
    };
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(const EvalFunction<DerivedEval>& evalFunc,
                                                  const TrialExpressionBase<DerivedTrial>& trialExp)
    {
        return TrialProd<DerivedEval,DerivedTrial>(evalFunc.derived(),trialExp.derived());
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(EvalFunction<DerivedEval>&& evalFunc,
                                                  const TrialExpressionBase<DerivedTrial>& trialExp)
    {
        return TrialProd<DerivedEval,DerivedTrial>(std::move(evalFunc.derived()),trialExp.derived());
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(const EvalFunction<DerivedEval>& evalFunc,
                                                  TrialExpressionBase<DerivedTrial>&& trialExp)
    {
        return TrialProd<DerivedEval,DerivedTrial>(evalFunc.derived(),std::move(trialExp.derived()));
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(EvalFunction<DerivedEval>&& evalFunc,
                                                  TrialExpressionBase<DerivedTrial>&& trialExp)
    {
        return TrialProd<DerivedEval,DerivedTrial>(std::move(evalFunc.derived()),std::move(trialExp.derived()));
    }
    
    template <typename DerivedTrial>
    TrialProd<Constant<double,1,1>,DerivedTrial> operator*(const double& d,
                                                           const TrialExpressionBase<DerivedTrial>& trialExp)
    {
        return TrialProd<Constant<double,1,1>,DerivedTrial>(make_constant(d),trialExp.derived());
    }
    
    template <typename DerivedTrial>
    TrialProd<Constant<double,1,1>,DerivedTrial> operator*(const double& d,
                                                           TrialExpressionBase<DerivedTrial>&& trialExp)
    {
        return TrialProd<Constant<double,1,1>,DerivedTrial>(make_constant(d),std::move(trialExp.derived()));
    }
    
    template <int rows1, int cols1, typename T2>
    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(const Eigen::Matrix<double,rows1,cols1>& c,
                                                                                    const TrialExpressionBase<T2>& op2)
    {
        return TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2>(make_constant(c),op2.derived());
    }
    
    template <int rows1, int cols1, typename T2>
    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(const Eigen::Matrix<double,rows1,cols1>& c,
                                                                                    TrialExpressionBase<T2>&& op2)
    {
        return TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2>(make_constant(c),std::move(op2.derived()));
    }

}	// close namespace
#endif


//    
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(const Constant<T1,rows1,cols1>& c,
//                                                     const TrialExpressionBase<T2>& op2)
//    {
//                std::cout<<"TrialProd operator* 1a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(c,op2);
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(Constant<T1,rows1,cols1>&& c,
//                                                     const TrialExpressionBase<T2>& op2)
//    {
//                std::cout<<"TrialProd operator* 2a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(std::move(c),op2);
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(const Constant<T1,rows1,cols1>& c,
//                                                     TrialExpressionBase<T2>&& op2)
//    {
//                std::cout<<"TrialProd operator* 3a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(c,std::move(op2));
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(Constant<T1,rows1,cols1>&& c,
//                                                     TrialExpressionBase<T2>&& op2)
//    {
//                std::cout<<"TrialProd operator* 4a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(std::move(c),std::move(op2));
//
//    }
