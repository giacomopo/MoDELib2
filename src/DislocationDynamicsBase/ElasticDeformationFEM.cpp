/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformationFEM_cpp_
#define model_ElasticDeformationFEM_cpp_

#include <cmath>
#include <ElasticDeformationFEM.h>

namespace model
{

    template<int dim>
    ElasticDeformationFEM<dim>::ElasticDeformationFEM(DislocationDynamicsBase<dim>& ddBase_in) :
    /* init */ ddBase(ddBase_in)
//    /* init */,voigtTraits((typename SymmetricVoigtTraits<3>::VoigtSizeMatrixType()<<0,0,1,1,2,2,0,1,1,2,0,2).finished())
//    /* init */,voigtTraits(ddBase.voigtTraits)
    /* init */,C(get_C(ddBase.poly.mu,ddBase.poly.nu))
    /* init */,u(ddBase.fe->template trial<'u',dim>())
    /* init */,b(grad(u))
    /* init */,e(def(u))
    /* init */,s(C*e)
    /* init */,dV(ddBase.fe->template domain<EntireDomain,4,GaussLegendre>())
    /* init */,bWF((test(e),s),dV)
    /* init */,z(ddBase.fe->template trial<'z',dim>())
    {
        u.setConstant(Eigen::Matrix<double,dim,1>::Zero());        
        z.setConstant(Eigen::Matrix<double,dim,1>::Zero());
    }

template<int dim>
typename ElasticDeformationFEM<dim>::CmatrixType ElasticDeformationFEM<dim>::get_C(const double& mu, const double& nu) const
{
    const double lam=2.0*mu*nu/(1.0-2.0*nu);
    const double C11(lam+2.0*mu);
    const double C12(lam);
    const double C44(mu); // C multiplies engineering strain
    
    CmatrixType temp;
    temp<<C11, C12, C12, 0.0, 0.0, 0.0,
    /***/ C12, C11, C12, 0.0, 0.0, 0.0,
    /***/ C12, C12, C11, 0.0, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
    return temp;
}

template<int dim>
typename ElasticDeformationFEM<dim>::MomentumMatrixType ElasticDeformationFEM<dim>::elementMomentumKernelMatrix(const VectorDim& abscissa,const ElementType& ele) const
{
    const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(abscissa));
    const auto NJ(ele.sf(bary)*ele.absJ(bary)/this->ddBase.mesh.volume());
    MomentumMatrixType km(MomentumMatrixType::Zero());
    const VectorDim x(ele.position(bary));
    const MatrixDim X((MatrixDim()<<0.0,-x(2),x(1),x(2),0.0,-x(0),-x(1),x(0),0.0).finished());
    for(int j=0;j<ElementType::nodesPerElement;++j)
    {
        km.template block<dim,dim>(0,dim*j)=MatrixDim::Identity()*NJ(j);
        km.template block<dim,dim>(dim,dim*j)=X*NJ(j);
    }
    return km;
}


template<int dim>
typename ElasticDeformationFEM<dim>::MomentumMatrixType ElasticDeformationFEM<dim>::elementMomentumMatrix(const ElementType& ele) const
{
    MomentumMatrixType km(MomentumMatrixType::Zero());
    MomentumQuadratureType::integrate(this,km,&ElasticDeformationFEM<dim>::elementMomentumKernelMatrix,ele);
    return km;
}
    
    template struct ElasticDeformationFEM<3>;

}
#endif

