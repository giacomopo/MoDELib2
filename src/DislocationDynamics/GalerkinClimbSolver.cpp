/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinClimbSolver_cpp_
#define model_GalerkinClimbSolver_cpp_

#include <ClusterDynamicsParameters.h>
#include <GalerkinClimbSolver.h>
#include <TerminalColors.h>

namespace model
{

    template <typename DislocationNetworkType>
    GalerkinClimbSolver<DislocationNetworkType>::GalerkinClimbSolver(const DislocationNetworkType& DN_in,const ClusterDynamics<dim>* const CD_in) :
    /* init */ DislocationClimbSolverBase<DislocationNetworkType>(DN_in,CD_in)
    {
        std::cout<<greenBoldColor<<"Creating GalerkinClimbSolver"<<defaultColor<<std::endl;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::ForceVectorMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterForceKernel(const int& k,const NetworkLinkType& fieldSegment) const
    {
        const Eigen::Matrix<double,1,mSize> deltaC(fieldSegment.quadraturePoint(k).cDD-fieldSegment.quadraturePoint(k).cCD);
        const double u(fieldSegment.quadraturePoint(k).abscissa);
        return  (ForceVectorMatrixType()<<(1.0-u)*deltaC,u*deltaC).finished();
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::ForceVectorMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterForceVector(const NetworkLinkType& fieldSegment) const
    {
        ForceVectorMatrixType F(ForceVectorMatrixType::Zero());
        const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,F,&GalerkinClimbSolver<DislocationNetworkType>::clusterForceKernel,fieldSegment);
        F.row(0)*=bxt.dot(fieldSegment.source->climbDirection());
        F.row(1)*=bxt.dot(fieldSegment.sink->climbDirection());
        return F;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::StiffnessMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessKernel(const int& k,const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType temp(StiffnessMatrixType::Zero());
        if(sourceSegment.grains().size() == 1)
        {
            const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
            const double u(fieldSegment.quadraturePoint(k).abscissa);
            const Eigen::Matrix<double,2,1> m(( Eigen::Matrix<double,2,1>()<<(1.0-u)*bxt.dot(fieldSegment.source->climbDirection()),
                                               /*                               */ u*bxt.dot(fieldSegment.sink->climbDirection())).finished());
            const auto sourceConcentratoinMatrices(sourceSegment.concentrationMatrices(fieldSegment.quadraturePoint(k).r,this->CD->cdp));
            for(int kc=0; kc < mSize; kc++)
            {
                temp.template block<2,2>(kc*2,0)+=m*sourceConcentratoinMatrices.row(kc);
            }
        }
        return temp;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::StiffnessMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessMatrix(const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType K(StiffnessMatrixType::Zero());
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,K,&GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessKernel,fieldSegment,sourceSegment);
        return K;
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocitiesPipe() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        return nodeVelocities;
    }

    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocities()
    {
        std::cout<<", climbSolver "<<std::flush;
        computeClimbScalarVelocitiesBulk();
    }

    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocitiesBulk()
    {
        std::vector<std::vector<Eigen::Triplet<double>>> lhsT(mSize);
        std::vector<Eigen::VectorXd> Fc(mSize,Eigen::VectorXd::Zero(this->DN.networkNodes().size()));
        bool useLumpedSolver(true);
        
        for(const auto& fieldLink : this->DN.networkLinks())
        {// sum line-integral part of displacement field per segment
            if(   !fieldLink.second.lock()->hasZeroBurgers()
               && fieldLink.second.lock()->isSessile()
               && !fieldLink.second.lock()->isBoundarySegment()
               && !fieldLink.second.lock()->isGrainBoundarySegment()
               &&  fieldLink.second.lock()->chordLength()>FLT_EPSILON
               )
            {
                const size_t i0(fieldLink.second.lock()->source->gID());
                const size_t i1(fieldLink.second.lock()->  sink->gID());
                
                const ForceVectorMatrixType fc(clusterForceVector(*fieldLink.second.lock()));
                for(int kc=0; kc<mSize; ++kc)
                {
                    Fc[kc](i0)+=fc(0,kc);
                    Fc[kc](i1)+=fc(1,kc);
                }
                
                for(const auto& sourceLink : this->DN.networkLinks())
                {// sum line-integral part of displacement field per segment
                    if(   !sourceLink.second.lock()->hasZeroBurgers()
                       && !sourceLink.second.lock()->isBoundarySegment()
                       && !sourceLink.second.lock()->isGrainBoundarySegment()
                       &&  sourceLink.second.lock()->chordLength()>FLT_EPSILON
                       )
                    {
                        const size_t j0(sourceLink.second.lock()->source->gID());
                        const size_t j1(sourceLink.second.lock()->  sink->gID());
                        const StiffnessMatrixType kcc(clusterStiffnessMatrix(*fieldLink.second.lock(),*sourceLink.second.lock()));
                        
                        for(int kc=0; kc<mSize; ++kc)
                        {
                            const Eigen::Matrix<double,2,2> kccs(kcc.template block<2,2>(2*kc,0));
                            
                            if(useLumpedSolver)
                            {
                                lhsT[kc].emplace_back(i0,i0,0.5*kccs(0,0));
                                lhsT[kc].emplace_back(j0,j0,0.5*kccs(0,0));
                                
                                lhsT[kc].emplace_back(i0,i0,0.5*kccs(0,1));
                                lhsT[kc].emplace_back(j1,j1,0.5*kccs(0,1));
                                
                                lhsT[kc].emplace_back(i1,i1,0.5*kccs(1,0));
                                lhsT[kc].emplace_back(j0,j0,0.5*kccs(1,0));
                                
                                lhsT[kc].emplace_back(i1,i1,0.5*kccs(1,1));
                                lhsT[kc].emplace_back(j1,j1,0.5*kccs(1,1));
                            }
                            else
                            {
                                const bool forceSym(true);// at the moment we need to symmtrize, since numerical intgration on fieldLink is not equivalent to analytical integration on sourceLink
                                if(forceSym)
                                {
                                    lhsT[kc].emplace_back(i0,j0,0.5*kccs(0,0));
                                    lhsT[kc].emplace_back(i0,j1,0.5*kccs(0,1));
                                    lhsT[kc].emplace_back(i1,j0,0.5*kccs(1,0));
                                    lhsT[kc].emplace_back(i1,j1,0.5*kccs(1,1));
                                    
                                    lhsT[kc].emplace_back(j0,i0,0.5*kccs(0,0));
                                    lhsT[kc].emplace_back(j1,i0,0.5*kccs(0,1));
                                    lhsT[kc].emplace_back(j0,i1,0.5*kccs(1,0));
                                    lhsT[kc].emplace_back(j1,i1,0.5*kccs(1,1));
                                }
                                else
                                {
                                    lhsT[kc].emplace_back(i0,j0,kccs(0,0));
                                    lhsT[kc].emplace_back(i0,j1,kccs(0,1));
                                    lhsT[kc].emplace_back(i1,j0,kccs(1,0));
                                    lhsT[kc].emplace_back(i1,j1,kccs(1,1));
                                }
                            }
                        }
                    }
                }
            }
        }
        
        Eigen::SparseMatrix<double> Kcc(this->DN.networkNodes().size(),this->DN.networkNodes().size());
        std::vector<Eigen::Array<double,1,mSize>> nodeV(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        
        for(int kc=0; kc<mSize; ++kc)
        {
            Kcc.setFromTriplets(lhsT[kc].begin(),lhsT[kc].end());
            
            if(size_t(Kcc.rows())!=this->DN.networkNodes().size() || size_t(Kcc.cols())!=this->DN.networkNodes().size())
            {
                throw std::runtime_error("the Stiffness Matrix size is not equal to the node size.");
            }
            
            if(useLumpedSolver)
            {
                Eigen::VectorXd Kccd(Kcc.diagonal());
                for (size_t n=0; n<this->DN.networkNodes().size(); n++)
                {
                    if(std::fabs(Kccd(n))>FLT_EPSILON)
                    {
                        nodeV[n](kc)=Fc[kc](n)/Kccd(n);
                    }
                }
            }
        }

        this->scalarVelocities().resize(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        for(size_t k=0; k<this->DN.networkNodes().size();++k)
        {
            this->scalarVelocities()[k]=nodeV[k];
        }
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocitiesBulk() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        size_t k=0;
        for(auto& node: this->DN.networkNodes())
        {
            const double vScalarTotal(-1.0*(this->scalarVelocities()[k]*this->CD->cdp.msVector/this->CD->cdp.msVector.abs()).matrix().sum());
            nodeVelocities.template segment<dim>(k*dim)= vScalarTotal*node.second.lock()->climbDirection();
            k++;
        }
        return nodeVelocities;
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        return getNodeVelocitiesBulk()+getNodeVelocitiesPipe();
    }

    template class GalerkinClimbSolver<DislocationNetwork<3,0>>;

}
#endif
