/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamics_cpp_
#define model_ClusterDynamics_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ClusterDynamics.h>
//#include <ExternalAndInternalBoundary.h>
//#include <Fix.h>


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
    /* init */,cdp(this->microstructures.ddBase)
//    /* init */,ClusterDynamicsBase<dim>(this->microstructures.ddBase)
    /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.fe? bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("useClusterDynamicsFEM",true)) : false )
    /* init */,clusterDynamicsFEM(useClusterDynamicsFEM?  new ClusterDynamicsFEM<dim>(this->microstructures.ddBase,cdp) : nullptr)
    /* init */,uniformControllers(useClusterDynamicsFEM? UniformControllerContainerType() : getUniformControllers(this->microstructures.ddBase,this->cdp))
//    /* init */,nodeListInternalExternal(this->microstructures.ddBase.isPeriodicDomain ? -1 : this->microstructures.ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
//    /* init */,mobileClustersIncrement(this->microstructures.ddBase.fe->template trial<'d',mSize>())
//    /* init */,dV(this->microstructures.ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-this->microstructures.ddBase.poly.Omega*this->mobileFlux)*dV)
//    /* init */,dmBWF((test(grad(mobileClustersIncrement)),-this->microstructures.ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)
//    /* init */,mSolver(true,FLT_EPSILON)
////    /* init */,mSolver(mBWF,true,FLT_EPSILON)
//    /* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    ///* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G.transpose().eval()))*dV).globalVector())
{
//        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        
        
//        std::cout<<cascadeGlobalProduction<<std::endl;
        
    }

    

    template<int dim>
    void ClusterDynamics<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(clusterDynamicsFEM)
        {
            clusterDynamicsFEM->initializeConfiguration(configIO,f_file,F_labels);
        }
        else
        {// nothing to initialize for uniformControllers
            
        }
    }

template<int dim>
void ClusterDynamics<dim>::applyBoundaryConditions()
{
 
    const auto& nodesInternalExternal(clusterDynamicsFEM->mobileClusters.fe().nodeList(clusterDynamicsFEM->nodeListInternalExternal));
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t k=0;k<nodesInternalExternal.size();++k)
    {
        const auto& node(nodesInternalExternal[k]);
        const auto outNormal(node->outNormal()); // used to compute traction
        const MatrixDim sigma(this->microstructures.stress(node->P0,node,nullptr,nullptr));
        const double normalTraction(outNormal.dot(sigma*outNormal));
        const auto bndConcentration(cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
        
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
            clusterDynamicsFEM->mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
        }
    }
}

    template<int dim>
    void ClusterDynamics<dim>::solve()
    {
        if(clusterDynamicsFEM)
        {
            if(!clusterDynamicsFEM->solverInitialized)
            {
                clusterDynamicsFEM->initializeSolver();
            }
            std::cout<<", mobile BCs"<<std::flush;
            applyBoundaryConditions();
            clusterDynamicsFEM->solve();
        }
        else
        {

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
        return this->microstructures.ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ClusterDynamics<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(clusterDynamicsFEM)
        {
            const size_t nNodes(clusterDynamicsFEM->mobileClusters.fe().nodes().size());
            configIO.cdMatrix().resize(nNodes,mSize+iSize);
            configIO.cdMatrix().block(0,0,nNodes,mSize)=clusterDynamicsFEM->mobileClusters.dofVector().reshaped(mSize,nNodes).transpose();
            configIO.cdMatrix().block(0,mSize,nNodes,iSize)=clusterDynamicsFEM->immobileClusters.dofVector().reshaped(iSize,nNodes).transpose();
        }
        else
        {
            
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(clusterDynamicsFEM)
        {
            return clusterDynamicsFEM->inelasticDisplacementRate(x,node,ele,guess);
        }
        else
        {
            return VectorDim::Zero();
        }
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
                return eval(clusterDynamicsFEM->mobileClusters)(*node);
            }
            else
            {
                if(ele)
                {
                    return eval(clusterDynamicsFEM->mobileClusters)(*ele,ele->simplex.pos2bary(x));
                }
                else
                {
                    return eval(clusterDynamicsFEM->mobileClusters)(x,guess);
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
