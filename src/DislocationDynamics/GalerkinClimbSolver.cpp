/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinClimbSolver_cpp_
#define model_GalerkinClimbSolver_cpp_

#include <deque>


#include <ClusterDynamicsParameters.h>
#include <GalerkinClimbSolver.h>
#include <TerminalColors.h>
#include <EqualIteratorRange.h>

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
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocities(const bool& updateStiffness)
    {
        std::cout<<", climbSolver "<<std::flush;
        computeClimbScalarVelocitiesBulk(updateStiffness);
    }

    // Retained temporarily for an easy source-level comparison with the
    // optimized assembly below. The complete original file is also backed up
    // as GalerkingClimbSolver_BK.cpp.
#if 0
    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocitiesBulk(const bool& updateStiffness)
    {
//        std::vector<std::vector<Eigen::Triplet<double>>> lhsT(mSize);
        std::vector<Eigen::VectorXd> Fc(mSize,Eigen::VectorXd::Zero(this->DN.networkNodes().size()));
        std::vector<Eigen::VectorXd> KKc(mSize,Eigen::VectorXd::Zero(this->DN.networkNodes().size()));
//        bool useLumpedSolver(true);
        
#ifdef _OPENMP
        const size_t nThreads(omp_get_max_threads());
//        std::vector<std::vector<std::vector<Eigen::Triplet<double>>>> lhsTV(nThreads,lhsT);
        std::vector<std::vector<Eigen::VectorXd>> KKcT(nThreads,KKc);
        std::vector<std::vector<Eigen::VectorXd>> FcT(nThreads,Fc);
        const EqualIteratorRange<typename DislocationNetworkType::NetworkLinkContainerType::const_iterator> eir(this->DN.networkLinks().begin(),this->DN.networkLinks().end(),nThreads);
#pragma omp parallel for
        for(size_t thread=0;thread<eir.size();++thread)
        {
            auto& Fc_ref(FcT[thread]);
            auto& KKc_ref(KKcT[thread]);
            for(auto fieldLinkIter=eir[thread].first;fieldLinkIter!=eir[thread].second;++fieldLinkIter)
//            for(const auto& fieldLink : this->DN.networkLinks())
            {// sum line-integral part of displacement field per segment
                const auto& fieldLink(*fieldLinkIter);
#else
                auto& Fc_ref(Fc);
                auto& KKc_ref(KKc);
                for(const auto& fieldLink : this->DN.networkLinks())
                {
#endif
////#pragma omp parallel for
//        for(size_t thread=0;thread< eir.size();++thread)
//        {
//            for(auto fieldLinkIter=eir[thread].first;fieldLinkIter!=eir[thread].second;++fieldLinkIter)
////            for(const auto& fieldLink : this->DN.networkLinks())
//            {// sum line-integral part of displacement field per segment
//                const auto& fieldLink(*fieldLinkIter);
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
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                        {
                            Fc_ref[kc](i0)+=fc(0,kc);
                            Fc_ref[kc](i1)+=fc(1,kc);
//                        }
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
                                
//                                if(useLumpedSolver)
//                                {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                                    {
                                        KKc_ref[kc](i0)+=0.5*kccs(0,0)+0.5*kccs(0,1);
                                        KKc_ref[kc](j0)+=0.5*kccs(0,0)+0.5*kccs(1,0);
                                        KKc_ref[kc](i1)+=0.5*kccs(1,0)+0.5*kccs(1,1);
                                        KKc_ref[kc](j1)+=0.5*kccs(0,1)+0.5*kccs(1,1);

                                        
                                        //lhsT[kc].emplace_back(i0,i0,0.5*kccs(0,0));
                                        //lhsT[kc].emplace_back(j0,j0,0.5*kccs(0,0));
                                        
                                        //lhsT[kc].emplace_back(i0,i0,0.5*kccs(0,1));
//                                        lhsT[kc].emplace_back(j1,j1,0.5*kccs(0,1));
                                        
//                                        lhsT[kc].emplace_back(i1,i1,0.5*kccs(1,0));
                                        //lhsT[kc].emplace_back(j0,j0,0.5*kccs(1,0));
                                        
  //                                      lhsT[kc].emplace_back(i1,i1,0.5*kccs(1,1));
//                                        lhsT[kc].emplace_back(j1,j1,0.5*kccs(1,1));
//                                    }
//                                }
//                                else
//                                {
//                                    const bool forceSym(true);// at the moment we need to symmtrize, since numerical intgration on fieldLink is not equivalent to analytical integration on sourceLink
//                                    if(forceSym)
//                                    {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                                        {
//                                            lhsT[kc].emplace_back(i0,j0,0.5*kccs(0,0));
//                                            lhsT[kc].emplace_back(i0,j1,0.5*kccs(0,1));
//                                            lhsT[kc].emplace_back(i1,j0,0.5*kccs(1,0));
//                                            lhsT[kc].emplace_back(i1,j1,0.5*kccs(1,1));
//                                            
//                                            lhsT[kc].emplace_back(j0,i0,0.5*kccs(0,0));
//                                            lhsT[kc].emplace_back(j1,i0,0.5*kccs(0,1));
//                                            lhsT[kc].emplace_back(j0,i1,0.5*kccs(1,0));
//                                            lhsT[kc].emplace_back(j1,i1,0.5*kccs(1,1));
//                                        }
//                                    }
//                                    else
//                                    {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                                        {
//                                            lhsT[kc].emplace_back(i0,j0,kccs(0,0));
//                                            lhsT[kc].emplace_back(i0,j1,kccs(0,1));
//                                            lhsT[kc].emplace_back(i1,j0,kccs(1,0));
//                                            lhsT[kc].emplace_back(i1,j1,kccs(1,1));
//                                        }
//                                    }
//                                }
                            }
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
            for(size_t thread=0;thread< eir.size();++thread)
            {// recombine contributions of different threads
                for(int kc=0; kc<mSize; ++kc)
                {
                    Fc[kc]+=FcT[thread][kc];
                    KKc[kc]+=KKcT[thread][kc];
//                    lhsT[kc].insert(lhsT[kc].end(), lhsTV[thread][kc].begin(), lhsTV[thread][kc].end());
                }
            }
#endif
        
//        Eigen::SparseMatrix<double> Kcc(this->DN.networkNodes().size(),this->DN.networkNodes().size());
        std::vector<Eigen::Array<double,1,mSize>> nodeV(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        
        for(int kc=0; kc<mSize; ++kc)
        {
//            Kcc.setFromTriplets(lhsT[kc].begin(),lhsT[kc].end());
            
//            if(size_t(Kcc.rows())!=this->DN.networkNodes().size() || size_t(Kcc.cols())!=this->DN.networkNodes().size())
//            {
//                throw std::runtime_error("the Stiffness Matrix size is not equal to the node size.");
//            }
            
//            if(useLumpedSolver)
//            {
//                Eigen::VectorXd Kccd(Kcc.diagonal());
                for (size_t n=0; n<this->DN.networkNodes().size(); n++)
                {
                    if(std::fabs(KKc[kc](n))>FLT_EPSILON)
                    {
                        nodeV[n](kc)=Fc[kc](n)/KKc[kc](n);
                    }
                }
//            }
//            else
//            {
//                throw std::runtime_error("GalerkinClimbSolver: only the lumped solver is implemented.");
//            }
        }

        this->scalarVelocities().resize(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        for(size_t k=0; k<this->DN.networkNodes().size();++k)
        {
            this->scalarVelocities()[k]=nodeV[k];
        }
    }

#endif

    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocitiesBulk(const bool& updateStiffness)
    {
        const size_t nodeCount(this->DN.networkNodes().size());
        Fc.assign(mSize,Eigen::VectorXd::Zero(nodeCount));
        if(updateStiffness)
        {
            KKc.assign(mSize,Eigen::VectorXd::Zero(nodeCount));
        }

        using NetworkLinkPtrType=std::shared_ptr<NetworkLinkType>;
        struct SourceData
        {
            NetworkLinkPtrType link;
            size_t j0;
            size_t j1;
        };
        struct FieldQuadratureData
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            VectorDim r;
            double weightedM0;
            double weightedM1;
        };
        using FieldQuadratureContainer=std::vector<FieldQuadratureData,Eigen::aligned_allocator<FieldQuadratureData>>;
        struct FieldData
        {
            NetworkLinkPtrType link;
            size_t i0;
            size_t i1;
            ForceVectorMatrixType force;
            FieldQuadratureContainer quadrature;
        };

        // Lock and filter each network link once, outside the O(N^2) loop.
        std::vector<SourceData> sources;
        std::vector<FieldData> fields;
        sources.reserve(this->DN.networkLinks().size());
        fields.reserve(this->DN.networkLinks().size());
        for(const auto& networkLink : this->DN.networkLinks())
        {
            const auto link(networkLink.second.lock());
            if(!link)
            {
                continue;
            }

            const bool isActive(!link->hasZeroBurgers()
                                && !link->isBoundarySegment()
                                && !link->isGrainBoundarySegment()
                                && link->chordLength()>FLT_EPSILON);
            // clusterStiffnessKernel returned zero for all other sources.
            if(isActive && link->grains().size()==1)
            {
                sources.push_back({link,link->source->gID(),link->sink->gID()});
            }

            if(isActive && link->isSessile())
            {
                FieldData field{link,link->source->gID(),link->sink->gID(),clusterForceVector(*link),{}};
                field.quadrature.reserve(link->quadraturePoints().size());

                // These quantities depend only on the field segment and its
                // quadrature point, so compute them once rather than per source.
                const auto bxt(link->burgers().cross(link->chord()));
                const double sourceProjection(bxt.dot(link->source->climbDirection()));
                const double sinkProjection(bxt.dot(link->sink->climbDirection()));
                for(const auto& qp : link->quadraturePoints())
                {
                    field.quadrature.push_back({qp.r,
                                                qp.weight*(1.0-qp.abscissa)*sourceProjection,
                                                qp.weight*qp.abscissa*sinkProjection});
                }
                fields.emplace_back(std::move(field));
            }
        }

#ifdef _OPENMP
        const size_t nThreads(omp_get_max_threads());
        std::vector<std::vector<Eigen::VectorXd>> KKcT(nThreads,KKc);
        std::vector<std::vector<Eigen::VectorXd>> FcT(nThreads,Fc);
#pragma omp parallel for schedule(guided)
        for(size_t fieldIndex=0;fieldIndex<fields.size();++fieldIndex)
        {
            const size_t thread(omp_get_thread_num());
            auto& Fc_ref(FcT[thread]);
            auto& KKc_ref(KKcT[thread]);
            const auto& field(fields[fieldIndex]);
#else
        auto& Fc_ref(Fc);
        auto& KKc_ref(KKc);
        for(const auto& field : fields)
        {
#endif
            for(int kc=0;kc<mSize;++kc)
            {
                Fc_ref[kc](field.i0)+=field.force(0,kc);
                Fc_ref[kc](field.i1)+=field.force(1,kc);
            }

            if(updateStiffness)
            {
                for(const auto& source : sources)
                {
                    Eigen::Array<double,mSize,1> i0Term(Eigen::Array<double,mSize,1>::Zero());
                    Eigen::Array<double,mSize,1> j0Term(Eigen::Array<double,mSize,1>::Zero());
                    Eigen::Array<double,mSize,1> i1Term(Eigen::Array<double,mSize,1>::Zero());
                    Eigen::Array<double,mSize,1> j1Term(Eigen::Array<double,mSize,1>::Zero());

                    // Assemble exactly the four lumped contributions consumed by
                    // the solver; do not materialize the 2*mSize-by-2 pair matrix.
                    for(const auto& qp : field.quadrature)
                    {
                        const auto concentration(source.link->concentrationMatrices(qp.r,this->CD->cdp));
                        const auto c0(concentration.col(0).array());
                        const auto c1(concentration.col(1).array());
                        i0Term+=0.5*qp.weightedM0*(c0+c1);
                        j0Term+=0.5*(qp.weightedM0+qp.weightedM1)*c0;
                        i1Term+=0.5*qp.weightedM1*(c0+c1);
                        j1Term+=0.5*(qp.weightedM0+qp.weightedM1)*c1;
                    }

                    for(int kc=0;kc<mSize;++kc)
                    {
                        KKc_ref[kc](field.i0)+=i0Term(kc);
                        KKc_ref[kc](source.j0)+=j0Term(kc);
                        KKc_ref[kc](field.i1)+=i1Term(kc);
                        KKc_ref[kc](source.j1)+=j1Term(kc);
                    }
                }
            }
        }

#ifdef _OPENMP
        for(size_t thread=0;thread<nThreads;++thread)
        {
            for(int kc=0;kc<mSize;++kc)
            {
                Fc[kc]+=FcT[thread][kc];
                if(updateStiffness)
                {
                    KKc[kc]+=KKcT[thread][kc];
                }
            }
        }
#endif

        std::vector<Eigen::Array<double,1,mSize>> nodeV(nodeCount,Eigen::Array<double,1,mSize>::Zero());
        for(int kc=0;kc<mSize;++kc)
        {
            for(size_t n=0;n<nodeCount;++n)
            {
                if(std::fabs(KKc[kc](n))>FLT_EPSILON)
                {
                    nodeV[n](kc)=Fc[kc](n)/KKc[kc](n);
                }
            }
        }

        this->scalarVelocities().resize(nodeCount,Eigen::Array<double,1,mSize>::Zero());
        for(size_t k=0;k<nodeCount;++k)
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

    template class GalerkinClimbSolver<DislocationNetwork<3>>;

}
#endif
