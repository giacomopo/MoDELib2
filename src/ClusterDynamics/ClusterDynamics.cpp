/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamics_cpp_
#define model_ClusterDynamics_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ClusterDynamics.h>
#include <ExternalAndInternalBoundary.h>
#include <Fix.h>


namespace model
{

    //template <int dim>
    //struct BoundaryConcentration
    //{
    //    const IrradiationClimbParameters<dim>& icp;
    //    const Eigen::Matrix<double,dim,dim>& stress;
    //
    //    /**********************************************************************/
    //    BoundaryConcentration(const IrradiationClimbParameters<dim>& icp_in,const Eigen::Matrix<double,dim,dim>& stress_in) :
    //    /* */ icp(icp_in)
    //    /* */,stress(stress_in)
    //    {
    //
    //    }
    //
    //    /**********************************************************************/
    //    template <typename NodeType,int dofPerNode>
    //    Eigen::Matrix<double,dofPerNode,1> operator()(const NodeType& node,
    //                                                  Eigen::Matrix<double,dofPerNode,1>& val) const
    //    {
    //
    //        const auto v1(icp.equilibriumClusterConcentration(stress.trace()));
    //
    //        const auto outNormal(node.outNormal()); // used to compute traction
    //        const double tn = outNormal.dot(stress*outNormal);
    //
    //        val=(v1.array()*exp(-tn*icp.omega/icp.kB/icp.T*icp.speciesVector)).matrix().transpose();
    //
    //        return val;
    //
    //    }
    //
    //};

template <int dim>
typename ClusterDynamics<dim>::UniformControllerContainerType ClusterDynamics<dim>::getUniformControllers(const DislocationDynamicsBase<dim>& ddBase,const ClusterDynamicsParameters<dim>& cdp)
{
    typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
    typedef typename UniformControllerType::VectorVoigt VectorVoigt;
    
    UniformControllerContainerType temp;
    for(const auto& pair : cdp.D)
    {
        MatrixVoigt C(MatrixVoigt::Zero());
        const auto& grainID(pair.first);
        for(size_t k=0;k<pair.second.size();++k)
        {
            const auto& Dm(pair.second[k]);
            C.template block<dim,dim>(k*dim,k*dim)=Dm;
        }
        
        const auto f0(VectorVoigt::Zero()); // flux
        const auto f0Dot(VectorVoigt::Zero()); // flux rate
        const auto g0(VectorVoigt::Zero()); // Concentration gradient
        const auto g0Dot(VectorVoigt::Zero()); // Concentration gradient rate
        const auto stiffnessRatio(VectorVoigt::Ones()*1.0e20); // impose concentration gradient
        const double t0(0.0);

        temp.emplace(grainID,new UniformControllerType(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
    }
    return temp;
}


    template<int dim>
    ClusterDynamics<dim>::ClusterDynamics(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("ClusterDynamics",mc)
    /* init */,ClusterDynamicsBase<dim>(this->microstructures.ddBase)
    /* init */,useClusterDynamicsFEM(this->ddBase.fe? bool(TextFileParser(this->ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("useClusterDynamicsFEM",true)) : false )
    /* init */,uniformControllers(useClusterDynamicsFEM? UniformControllerContainerType() : getUniformControllers(this->microstructures.ddBase,this->cdp))
    /* init */,nodeListInternalExternal(this->microstructures.ddBase.isPeriodicDomain ? -1 : this->microstructures.ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
    /* init */,mobileClustersIncrement(this->microstructures.ddBase.fe->template trial<'d',mSize>())
    /* init */,dV(this->microstructures.ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
    /* init */,mBWF((test(this->mobileGrad),-this->microstructures.ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,dmBWF((test(grad(mobileClustersIncrement)),-this->microstructures.ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)
    /* init */,mSolver(true,FLT_EPSILON)
//    /* init */,mSolver(mBWF,true,FLT_EPSILON)
    /* init */,cascadeGlobalProduction(((test(this->mobileClusters),make_constant(this->cdp.G.transpose().eval()))*dV).globalVector())
{
        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        
        
//        std::cout<<cascadeGlobalProduction<<std::endl;
        
    }

    template<int dim>
    void ClusterDynamics<dim>::applyBoundaryConditions()
    {
     
        const auto& nodesInternalExternal(this->mobileClusters.fe().nodeList(nodeListInternalExternal));
        
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t k=0;k<nodesInternalExternal.size();++k)
        {
            const auto& node(nodesInternalExternal[k]);
            const auto outNormal(node->outNormal()); // used to compute traction
            const MatrixDim sigma(this->microstructures.stress(node->P0,node,nullptr,nullptr));
            const double normalTraction(outNormal.dot(sigma*outNormal));
            const auto bndConcentration(this->cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
            
            VectorMSize otherConcentration(VectorMSize::Zero());
            for(const auto& microstructure : this->microstructures)
            {
                if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
                {// not the ClusterDynamics physics
                    otherConcentration += microstructure->mobileConcentration(node->P0,node,nullptr,nullptr);
                }
            }
            for(int k=0;k<mSize;++k)
            {
                this->mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
            }
        }
    }

    template<int dim>
    void ClusterDynamics<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(uniformControllers.size())
        {// nothing to initialize for uniformControllers
            
        }
        else
        {
            if(size_t(configIO.cdMatrix().size())==this->mobileClusters.gSize()+this->immobileClusters.gSize())
            {
                const size_t nNodes(this->mobileClusters.fe().nodes().size());
                this->mobileClusters=configIO.cdMatrix().block(0,0,nNodes,this->mSize).transpose().reshaped(this->mobileClusters.gSize(),1);
                this->immobileClusters=configIO.cdMatrix().block(0,this->mSize,nNodes,this->iSize).transpose().reshaped(this->immobileClusters.gSize(),1);
            }
            else
            {
                if(configIO.cdMatrix().size())
                {
                    throw std::runtime_error("ClusterDynamics: TrialFunctions initializatoin size mismatch");
                }
            }
            
            std::array<bool,mSize> allComps;
            for(int k=0;k<mSize;++k)
            {
                allComps[k]=true;
            }

            this->mobileClusters.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
            mobileClustersIncrement.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
            mSolver.compute(mBWF); // call this after assigning the BCs
        }
    }

    template<int dim>
    void ClusterDynamics<dim>::solve()
    {
        if(uniformControllers.size())
        {// nothing to solve for uniformControllers
            
        }
        else
        {
            solveMobileClusters();
            solveImmobileClusters();
        }
    }

    template<int dim>
    void ClusterDynamics<dim>::updateConfiguration()
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
    }

    template<int dim>
    double ClusterDynamics<dim>::getDt() const
    {
        return this->ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ClusterDynamics<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(uniformControllers.size())
        {// nothing to solve for uniformControllers
            
        }
        else
        {
            const size_t nNodes(this->mobileClusters.fe().nodes().size());
            configIO.cdMatrix().resize(nNodes,this->mSize+this->iSize);
            configIO.cdMatrix().block(0,0,nNodes,this->mSize)=this->mobileClusters.dofVector().reshaped(this->mSize,nNodes).transpose();
            configIO.cdMatrix().block(0,this->mSize,nNodes,this->iSize)=this->immobileClusters.dofVector().reshaped(this->iSize,nNodes).transpose();
        }
    }

    template<int dim>
    void ClusterDynamics<dim>::solveMobileClusters()
    {
        
        // Find new mobileConcentration
        std::cout<<" mobile BCs,"<<std::flush;
        applyBoundaryConditions();
        std::cout<<" mobile solver,"<<std::flush;
        this->mobileClusters=mSolver.solve(cascadeGlobalProduction);
        
        if(this->cdp.computeReactions)
        {
            const double cTol(1.0e-5);
            double cError(1.0);
            while(cError>cTol)
            {
                const auto R1((this->cdp.R1cd).eval());
                auto bWF_R1((test(this->mobileClustersIncrement),R1*(-1.0*mobileClustersIncrement))*dV); // THIS SHOULD BE STORED SINCE IT IS ALWAYS THE SAME
                auto lWF_R1((test(this->mobileClustersIncrement),eval(R1*this->mobileClusters))*dV);

                SecondOrderReaction<MobileTrialType> R2(this->mobileClusters,this->cdp);
                auto bWF_R2((test(this->mobileClustersIncrement),R2*(-1.0*mobileClustersIncrement))*dV);
                auto lWF_R2((test(this->mobileClustersIncrement),eval(R2*(0.5*this->mobileClusters)))*dV);
                
                // Missing immobile sinks

                Eigen::SparseMatrix<double,Eigen::RowMajor> AcIR;
                AcIR.resize(mobileClustersIncrement.gSize(),mobileClustersIncrement.gSize());
                std::vector<Eigen::Triplet<double>> globalTripletsR((bWF_R1+bWF_R2).globalTriplets());
                AcIR.setFromTriplets(globalTripletsR.begin(),globalTripletsR.end());

                FixedDirichletSolver rSolver(false,FLT_EPSILON);
                rSolver.compute(dmBWF+bWF_R1+bWF_R2);
                mobileClustersIncrement=rSolver.solve(cascadeGlobalProduction-mSolver.getA()*this->mobileClusters.dofVector()+(lWF_R1+lWF_R2).globalVector());

                Eigen::MatrixXd cOld(this->mobileClusters.dofVector());
                cOld.resize(mSize,this->mobileClusters.gSize()/mSize);
                this->mobileClusters += mobileClustersIncrement.dofVector();

                Eigen::MatrixXd cNew(this->mobileClusters.dofVector());
                cNew.resize(mSize,this->mobileClusters.gSize()/mSize);
                
                const Eigen::VectorXd absErr((cNew-cOld).rowwise().norm());
                const Eigen::VectorXd cNewNorm((cNew.rowwise().norm().array()+1.e-20).matrix());
                const Eigen::VectorXd relErr((absErr.array()/cNewNorm.array()).matrix());

                cError=relErr.maxCoeff();//aError/cInorm;
                std::cout<<"max values="<<cNew.rowwise().maxCoeff().transpose()<<std::endl;
                std::cout<<"min values="<<cNew.rowwise().minCoeff().transpose()<<std::endl;
                std::cout<<"absolute errors="<<absErr.transpose()<<std::endl;
                std::cout<<"solution norms="<<cNewNorm.transpose()<<std::endl;
                std::cout<<"relative error="<<relErr.transpose()<<std::endl;
                std::cout<<"convergenceError="<<cError<<std::endl;

                
            }
        }
        // Find immobile rate
        
        
        // Find diffusive-displacement rate
        
    }

    template<int dim>
    void ClusterDynamics<dim>::solveImmobileClusters()
    {
        
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        Eigen::Matrix<double,dim*mSize,1> speciesFlux(Eigen::Matrix<double,dim*mSize,1>::Zero());
        if(node)
        {
            speciesFlux=eval(this->mobileFlux)(*node);
        }
        else
        {
            if(ele)
            {
                speciesFlux=eval(this->mobileFlux)(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                speciesFlux=eval(this->mobileFlux)(x,guess);
            }
        }
        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
        for(int i=0; i<mSize; i++)
        {
            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
            netFlux += speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
        }
        return netFlux;
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortion() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorMSize ClusterDynamics<dim>::mobileConcentration(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformControllers.size())
        {
            const auto pointGrains(this->pointGrains(x,node,ele,guess));
            VectorMSize averageVal(VectorMSize::Zero());
            for(const auto& grain : pointGrains)
            {
                const auto grainID(grain->grainID);
                const auto& controller(uniformControllers.at(grainID));
                averageVal += controller->grad(this->microstructures.ddBase.simulationParameters.totalTime).reshaped(dim,mSize).transpose()*x;
            }
            return averageVal/pointGrains.size();
        }
        else
        {
            if(node)
            {
                return eval(this->mobileClusters)(*node);
            }
            else
            {
                if(ele)
                {
                    return eval(this->mobileClusters)(*ele,ele->simplex.pos2bary(x));
                }
                else
                {
                    return eval(this->mobileClusters)(x,guess);
                }
            }
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averageStress() const
    {
        return MatrixDim::Zero();
    }

    template struct ClusterDynamics<3>;
}
#endif


//template<int dim>
//void ClusterDynamics<dim>::solveDiffusiveDisplacement()
//{
//    std::cout<<"Solving diffusiveDisplacementRate..."<<std::flush;
//    const auto t0= std::chrono::system_clock::now();
//    diffusiveDisplacementRate=Eigen::VectorXd::Zero(this->diffusiveDisplacement.gSize());
//    for(const auto& node: this->microstructures.ddBase.fe->nodes())
//    {
//        const Eigen::Matrix<double,dim*mSize,1> speciesFlux(eval(this->mobileFlux)(node));
//        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
//        for(int i=0; i<mSize; i++)
//        {
//            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
//            netFlux+= speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
//        }
//        diffusiveDisplacementRate.template segment<dim>(dim*node.gID)=netFlux;
//    }
//    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//}
