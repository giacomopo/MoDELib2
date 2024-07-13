/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StressStraight_cpp_
#define model_StressStraight_cpp_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <cfloat>
#include <StressStraight.h>

namespace model
{

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::MatrixDim StressStraight<dim,Scalar>::nonSymmStress_kernel(const VectorDim& r) const
    {
    #if _MODEL_NON_SINGULAR_DD_ == 0
        /* Devincre, B., & Condat, M. (1992). Model validation of a 3D
         * simulation of dislocation dynamics: Discretization and line tension
         * effects. Acta Metallurgica Et Materialia, 40(10), 2629–2637.
         */
        
        static_assert(0,"THERE IS NO CONSISTENT REGULARIZATION OF THIS EXPRESSION. USE _MODEL_NON_SINGULAR_DD_=1");
        
        //            const Scalar L(r.dot(t));
        //            const Scalar R(r.norm()+DislocationFieldBase<dim>::a);
        //            const VectorDim rho(r-L*t); // distance vector to the line
        //            const VectorDim Y((L+R)*t+rho); // = R*t + r
        //            // Y2=R^2+R^2+2R*(r*t)=2R*(R+L)
        ////            const Scalar Y2(Y.squaredNorm()+DislocationFieldBase<dim>::a2);
        ////            return (material.C1* b.cross(Y)*t.transpose()
        ////            /*                 */ -b.cross(t)*Y.transpose()
        ////            /*                 */ -b.dot(Y.cross(t))*(2.0/Y2*rho*Y.transpose()+0.5*(MatrixDim::Identity()+t*t.transpose()+2.0/Y2*L/R*Y*Y.transpose()))
        ////                    )*2.0/Y2;
        //
        //
        ////            const Scalar f1(2.0/Y.squaredNorm());
        //            const Scalar f1(1.0/(R*(R+L))); // =2/Y2=2/(2R*(R+L))
        //
        ////            const Scalar f1(2.0/(Y.squaredNorm()+DislocationFieldBase<dim>::a2));
        ////            const Scalar f1(2.0/Y2);
        //            const Scalar f2(material.C1*f1);
        //            const Scalar bDotYcrosst(b.dot(Y.cross(t)));
        //            const Scalar f3(bDotYcrosst*f1*f1);
        //            const Scalar f4(0.5*f1*bDotYcrosst);
        //            const Scalar f5(f4*f1*L/R);
        //            //            const Scalar f5(0.5*f3*L/R);
        //
        //            return  f2*b.cross(Y)*t.transpose()
        //            /*  */ -f1*b.cross(t)*Y.transpose()
        //            /*  */ -f3*rho*Y.transpose()
        //            /*  */ -f4*(MatrixDim::Identity()+t*t.transpose())
        //            /*  */ -f5*Y*Y.transpose();
        
        //            const Scalar R(r.norm());
        const Scalar Ra(sqrt(r.squaredNorm()+DislocationFieldBase<dim>::a2));
        const VectorDim Y(r+Ra*t);
        const Scalar Yt(Y.dot(t));
        //            const Scalar Yta(Yt+DislocationFieldBase<dim>::a);
        const Scalar Yta(sqrt(Yt*Yt+DislocationFieldBase<dim>::a2));
        const Scalar Y2(Y.squaredNorm()+DislocationFieldBase<dim>::a2);
        const Scalar bYt(b.cross(Y).dot(t));
        
        
//        const Scalar f1(2.0/Y2);
//        const Scalar f1(EwaldLength>FLT_EPSILON? 2.0*erfc(Ra/EwaldLength)/Y2 : 2.0/Y2);
        const Scalar f1(EwaldLength>FLT_EPSILON? 2.0*erfc(Y2/EwaldLength/EwaldLength)/Y2 : 2.0/Y2);

        
        return  f1*material.C1*t*(b.cross(Y)).transpose()
        /*   */-f1*Y*bCt.transpose()
        /*   */-f1*bYt/Yta*t*r.transpose()
        /*   */-0.5*f1*bYt*MatrixDim::Identity()
        /*   */-f1*bYt*(R+Yt)/Ra/Y2*r*r.transpose()
        /*   */-0.5*f1*bYt*R/Yta*t*t.transpose();
        
        
    #elif _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */
        
        //            const double Ra2=r.squaredNorm()+DislocationFieldBase<dim>::a2;
        //            const double Ra=sqrt(Ra2);
        //            const double Ra3=std::pow(Ra,3);
        //            const double rdt=r.dot(t);
        //            const double rdt2=std::pow(rdt,2);
        //            const double A1=-(rdt*(3.0*Ra2-rdt2))/(std::pow(Ra2-rdt2,2)*Ra3);
        //            const double A2=1.0/Ra3-rdt*A1;
        //            const double A6=-rdt/((Ra2-rdt2)*Ra);
        //            const double A3=-rdt/Ra3+A6+rdt2*A1;
        //            const double A4=A6+DislocationFieldBase<dim>::a2*A1;
        //            const double A5=-material.C1*A6-0.5*DislocationFieldBase<dim>::a2*material.C1*A1;
        //            const double A7=material.nu/Ra-rdt*A6-0.5*DislocationFieldBase<dim>::a2*material.C1*A2;
        //
        //            const double rbt(r.cross(b).dot(t));
        //
        //            return  0.5*rbt*A1*r*r.transpose()
        //            /*  */ +    rbt*A2*t*r.transpose()
        //            /*  */ +0.5*rbt*A3*t*t.transpose()
        //            /*  */ +0.5*rbt*A4*MatrixDim::Identity()
        //            /*  */ +A5*r.cross(b)*t.transpose()
        //            /*  */ +A6*t.cross(b)*r.transpose()
        //            /*  */ +A7*t.cross(b)*t.transpose();
        
        const Scalar Ra2=r.squaredNorm()+DislocationFieldBase<dim>::a2;
        const Scalar Ra(sqrt(Ra2));
        const VectorDim Ya(r+Ra*t);
        const Scalar Yat(Ya.dot(t));
        const Scalar Ya2a2(Ya.squaredNorm()+DislocationFieldBase<dim>::a2);
        const VectorDim bYa(b.cross(Ya));
        const Scalar bYat(bYa.dot(t));
        
        
//        const Scalar f1(2.0/Ya2a2);
//        const Scalar f1(EwaldLength>FLT_EPSILON? 2.0*erfc(Ra/EwaldLength)/Ya2a2 : 2.0/Ya2a2);
        const Scalar f1(EwaldLength>FLT_EPSILON? 2.0*erfc(Ya2a2/EwaldLength/EwaldLength)/Ya2a2 : 2.0/Ya2a2);


        return f1*material.C1*(1.0+DislocationFieldBase<dim>::a2/Ya2a2)*t*bYa.transpose()
        /*  */+f1*material.C1*0.5*DislocationFieldBase<dim>::a2/Ra2*t*b.cross(r).transpose()
        /*  */-f1*Ya*bCt.transpose()
        /*  */-f1*bYat/Yat*t*r.transpose()
        /*  */-0.5*f1*bYat*(1.0+2.0*DislocationFieldBase<dim>::a2/Ya2a2+DislocationFieldBase<dim>::a2/Ra2)*MatrixDim::Identity()
        /*  */-f1*bYat*(Ra+Yat)/Ra/Ya2a2*r*r.transpose()
        /*  */-0.5*f1*bYat*Ra/Yat*t*t.transpose();
        
    #elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
        static_assert(0,"NOT IMPLEMENTED. USE _MODEL_NON_SINGULAR_DD_=1");
    #else
    #error Unsupported choice of field regularization
    #endif
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::VectorDim StressStraight<dim,Scalar>::displacement_kernel(const VectorDim& r) const
    {
        const Scalar Ra(sqrt(r.squaredNorm()+DislocationFieldBase<dim>::a2));
        const VectorDim Ya(r+Ra*t);
        const Scalar Yat(Ya.dot(t));
        return -(2.0-0.5/material.C1+(2.0-1.0/material.C1)*log(Yat)-DislocationFieldBase<dim>::a2/Ra/Yat)/8.0/M_PI*bCt
        /*  */ +bCt.dot(r)/Yat/Ra/material.C1/8.0/M_PI*Ya;
    }

    template <int dim,typename Scalar>
    StressStraight<dim,Scalar>::StressStraight(const PolycrystallineMaterialBase& material_in,const VectorDim& _P0,const VectorDim& _P1, const VectorDim& _b,
                                               const double& EwaldLength_in) :
    /* init list */ material(material_in),
    /* init list */ P0(_P0),
    /* init list */ P1(_P1),
    /* init list */ b(_b),
    /* init list */ chord(P1-P0),
    /* init list */ length(chord.norm()),
    /* init list */ t(chord/length),
    /* init  */ EwaldLength(EwaldLength_in),
    /* init list */ bCt(b.cross(t))
    {/*!\param[in] _P0 starting point of the segment
      * \param[in] _P0 ending point of the segment
      * \param[in] _b Burgers vector of the segment
      */
        
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::MatrixDim StressStraight<dim,Scalar>::nonSymmStress(const VectorDim& x) const
    {
        return length>FLT_EPSILON? nonSymmStress_kernel(P1-x)-nonSymmStress_kernel(P0-x) : MatrixDim::Zero().eval();
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::MatrixDim StressStraight<dim,Scalar>::stress(const VectorDim& x) const
    {
        const MatrixDim temp = nonSymmStress(x);
        return material.C2*(temp+temp.transpose());
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::VectorDim StressStraight<dim,Scalar>::displacement(const VectorDim& x) const
    {/*!\returns the line-integral part of the displacement contribution of this straight segment.
      * Note: the return value  does NOT include the solid-angle contribution to the displacement field.
      */
        return displacement_kernel(P1-x)-displacement_kernel(P0-x);
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::ConcentrationMatrixType StressStraight<dim,Scalar>::concentrationMatrices(const VectorDim& x, const size_t& grainID, const VectorDim& sourceDir, const VectorDim& sinkDir, const ClusterDynamicsParameters<dim>& icp) const
    {
        const auto& Dinv(icp.invD.at(grainID));
        const auto& Ddet(icp.detD.at(grainID));
        ConcentrationMatrixType M(ConcentrationMatrixType::Zero());
        if(length>FLT_EPSILON)
        {
            Eigen::Array<double,mSize,1> a(Eigen::Array<double,mSize,1>::Zero());
            Eigen::Array<double,mSize,1> b(Eigen::Array<double,mSize,1>::Zero());
            Eigen::Array<double,mSize,1> c(Eigen::Array<double,mSize,1>::Zero());
            for(size_t k=0; k<mSize; k++)
            {
                a(k) = chord.dot( Dinv[k]*chord );
                b(k) =-2.0*(x-P0).dot( Dinv[k]*chord  );
                c(k) = (x-P0).dot( Dinv[k]*(x-P0) )+ DislocationFieldBase<dim>::a2/pow(Ddet(k),1.0/3.0);
            }
            const Eigen::Array<double,mSize,1> ba(b/a);
            const Eigen::Array<double,mSize,1> ca(c/a);
            const Eigen::Array<double,mSize,1> sqbca(sqrt(1.0+ba+ca));
            const Eigen::Array<double,mSize,1> sqca(sqrt(ca));
            const Eigen::Array<double,mSize,1> logTerm(log((2.0*sqbca+2.0+ba)/(2.0*sqca+ba)));
            const double bxtSource(bCt.dot(sourceDir));
            const double bxtSink  (bCt.dot(sinkDir));
            const Eigen::Array<double,mSize,1> I0((1.0+0.5*ba)*logTerm-sqbca+sqca);
            const Eigen::Array<double,mSize,1> I1(     -0.5*ba*logTerm+sqbca-sqca);
            M.col(0)=bxtSource*length/(4.0*M_PI)*(I0/sqrt(a*Ddet)).matrix();
            M.col(1)=bxtSink*  length/(4.0*M_PI)*(I1/sqrt(a*Ddet)).matrix();
        }
        return M;
    }

    template <int dim,typename Scalar>
    typename StressStraight<dim,Scalar>::ConcentrationVectorType StressStraight<dim,Scalar>::clusterConcentration(const VectorDim& x, const size_t& grainID, const VectorDim& sourceDir, const Eigen::Array<double,1,mSize>& sourceVScalar, const VectorDim& sinkDir, const Eigen::Array<double,1,mSize>& sinkVScalar, const ClusterDynamicsParameters<dim>& icp) const
    {
        ConcentrationVectorType temp(ConcentrationVectorType::Zero());
        const ConcentrationMatrixType cM(concentrationMatrices(x, grainID, sourceDir, sinkDir, icp));
        for(int k=0; k<mSize; k++)
        {
            temp(k)=cM.row(k)*(Eigen::Matrix<double,2,1>()<<sourceVScalar(k),sinkVScalar(k)).finished();
        }
        return temp;
    }

    template class StressStraight<3,double>;

}
#endif
