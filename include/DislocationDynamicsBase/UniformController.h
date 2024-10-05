/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniformController_H_
#define model_UniformController_H_

#include <cmath>
#include <Eigen/Dense>

namespace model
{
    
    template <int voigtSize>
    struct UniformController
    {
        typedef Eigen::Matrix<double,voigtSize,1>          VectorVoigt;
        typedef Eigen::Matrix<double,voigtSize,voigtSize>  MatrixVoigt;
        
        const double t0;
        const MatrixVoigt A;
        const VectorVoigt alpha;
        const VectorVoigt g0;
        const VectorVoigt g0Dot;
        const VectorVoigt f0;
        const VectorVoigt f0Dot;
        const MatrixVoigt D;
        const MatrixVoigt MG;
        const MatrixVoigt MF;
        VectorVoigt gs;

        UniformController(const double& t0_in,
                          const MatrixVoigt&  A_in,const VectorVoigt& alpha_in,
                          const VectorVoigt& g0_in,const VectorVoigt& g0Dot_in,
                          const VectorVoigt& f0_in,const VectorVoigt& f0Dot_in) :
        /* init */ t0(t0_in)
        /* init */,A(A_in)
        /* init */,alpha(alpha_in)
        /* init */,g0(g0_in)
        /* init */,g0Dot(g0Dot_in)
        /* init */,f0(f0_in)
        /* init */,f0Dot(f0Dot_in)
        /* init */,D((alpha.array()*A.diagonal().array()).matrix().asDiagonal())
        /* init */,MG((A+D).inverse())
        /* init */,MF(A*MG)
        /* init */,gs(VectorVoigt::Zero())
        {

//            std::cout<<"A=\n"<<A<<std::endl;
//            std::cout<<"D=\n"<<D<<std::endl;
//            std::cout<<"MG=\n"<<MG<<std::endl;
//            std::cout<<"MF=\n"<<MF<<std::endl;
//            std::cout<<"MG.norm()="<<MG.norm()<<std::endl;
//            std::cout<<"MG finite="<<std::isfinite(MG.norm())<<std::endl;
//
//            const double norm2MG(MG.squaredNorm());
//            if(norm2MG!=norm2MG)
//            {
//                throw std::runtime_error("UniformController: MG is not finite.");
//            }

        }
        
        VectorVoigt flux(const double& t) const
        {
            return MF*(f0+f0Dot*(t-t0)
                       +D*(g0+g0Dot*(t-t0)-gs)
                       );
        }
        
        VectorVoigt grad(const double& t) const
        {
            return MG*(f0+f0Dot*(t-t0)
                       +D*(g0+g0Dot*(t-t0))
                       +A*gs
                       );
        }
        
//        VectorComp operator()(const VectorDim& x,const double& t) const
//        {
//            return grad(t).reshaped(dim,comp).transpose()*x;
//        }
        
    };
        
}
#endif
