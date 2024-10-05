/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_cpp_
#define model_ClusterDynamicsFEM_cpp_

#include <cmath>
#include <ClusterDynamicsFEM.h>
#include <ExternalAndInternalBoundary.h>
#include <Fix.h>

namespace model
{

    template<int dim>
    FluxMatrix<dim>::FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */cdp(cdp_in)
    {
        
    }

    template<int dim>
    const typename FluxMatrix<dim>::MatrixType FluxMatrix<dim>::operator() (const ElementType& elem, const BaryType&) const
    {
        const size_t grainID(elem.simplex.region->regionID);
        MatrixType D(MatrixType::Zero());
        for(size_t k=0;k<mSize;++k)
        {
            D.template block<dim,dim>(k*dim,k*dim)=-cdp.D.at(grainID)[k]/cdp.omega;
        }
        return D;
    }

    template struct FluxMatrix<3>;

template<int dim>
InvDscaling<dim>::InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in) :
/* init */cdp(cdp_in)
{
    
}

template<int dim>
const typename InvDscaling<dim>::MatrixType InvDscaling<dim>::operator() (const ElementType& elem, const BaryType&) const
{
    const size_t grainID(elem.simplex.region->regionID);
    MatrixType InvDscaling(MatrixType::Zero());
    for(size_t k=0;k<mSize;++k)
    {
        InvDscaling(k,k)=3.0/(cdp.D.at(grainID)[k].trace()/cdp.omega);
//        InvDscaling(k,k)=(cdp.D.at(grainID)[k].trace()/cdp.omega)/3.0;
//        InvDscaling(k,k)=1.0;

    }
    return InvDscaling;
}

template struct InvDscaling<3>;


    template<int dim>
    ClusterDynamicsFEM<dim>::ClusterDynamicsFEM(const DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */ ddBase(ddBase_in)
    /* init */,cdp(cdp_in)
    /* init */,iDs(cdp)
    /* init */,mobileClusters(ddBase.fe->template trial<'m',mSize>())
    /* init */,mobileGrad(grad(mobileClusters))
    /* init */,mobileFlux(FluxMatrix<dim>(cdp)*mobileGrad)
    /* init */,immobileClusters(ddBase.fe->template trial<'i',iSize>())
    /* init */,nodeListInternalExternal(ddBase.isPeriodicDomain ? -1 : ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
    /* init */,mobileClustersIncrement(ddBase.fe->template trial<'d',mSize>())
    /* init */,dV(ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,mBWF((test(grad(iDs*mobileClusters)),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,dmBWF((test(grad(iDs*mobileClustersIncrement)),-ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)
    /* init */,mSolver(true,FLT_EPSILON)
    /* init */,solverInitialized(false)
//    /* init */,cascadeGlobalProduction(((test(this->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    /* init */,cascadeGlobalProduction(((test(iDs*this->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    {
        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        mobileClusters.setConstant(cdp.equilibriumMobileConcentration(0.0).matrix().transpose());
        immobileClusters.setConstant(Eigen::Matrix<double,iSize,1>::Zero());
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveMobileClusters()
    {
        std::cout<<", mobile solver "<<std::flush;
        mobileClusters=mSolver.solve(cascadeGlobalProduction);
        
        if(this->cdp.computeReactions)
        {
            const double cTol(1.0e-5);
            double cError(1.0);
            while(cError>cTol)
            {
                const auto R1((this->cdp.R1cd).eval());
                auto bWF_R1((test(iDs*mobileClustersIncrement),R1*(-1.0*mobileClustersIncrement))*dV); // THIS SHOULD BE STORED SINCE IT IS ALWAYS THE SAME
                auto lWF_R1((test(iDs*mobileClustersIncrement),eval(R1*mobileClusters))*dV);
                
                SecondOrderReaction<MobileTrialType> R2(mobileClusters,this->cdp);
                auto bWF_R2((test(iDs*mobileClustersIncrement),R2*(-1.0*mobileClustersIncrement))*dV);
                auto lWF_R2((test(iDs*mobileClustersIncrement),eval(R2*(0.5*mobileClusters)))*dV);
                
                // Missing immobile sinks
                
                Eigen::SparseMatrix<double,Eigen::RowMajor> AcIR;
                AcIR.resize(mobileClustersIncrement.gSize(),mobileClustersIncrement.gSize());
                std::vector<Eigen::Triplet<double>> globalTripletsR((bWF_R1+bWF_R2).globalTriplets());
                AcIR.setFromTriplets(globalTripletsR.begin(),globalTripletsR.end());
                
                MobileReactionSolverType rSolver(false,FLT_EPSILON);
                rSolver.compute(dmBWF+bWF_R1+bWF_R2);
                mobileClustersIncrement=rSolver.solve(cascadeGlobalProduction-mSolver.getA()*mobileClusters.dofVector()+(lWF_R1+lWF_R2).globalVector());
                
                Eigen::MatrixXd cOld(mobileClusters.dofVector());
                cOld.resize(mSize,mobileClusters.gSize()/mSize);
                mobileClusters += mobileClustersIncrement.dofVector();
                
                Eigen::MatrixXd cNew(mobileClusters.dofVector());
                cNew.resize(mSize,mobileClusters.gSize()/mSize);
                
                const Eigen::VectorXd absErr((cNew-cOld).rowwise().norm());
                const Eigen::VectorXd cNewNorm((cNew.rowwise().norm().array()+1.e-50).matrix());
                const Eigen::VectorXd relErr((absErr.array()/cNewNorm.array()).matrix());
                
                cError=relErr.maxCoeff();//aError/cInorm;
                if(false)
                {
                    std::cout<<"max values="<<cNew.rowwise().maxCoeff().transpose()<<std::endl;
                    std::cout<<"min values="<<cNew.rowwise().minCoeff().transpose()<<std::endl;
                    std::cout<<"absolute errors="<<absErr.transpose()<<std::endl;
                    std::cout<<"solution norms="<<cNewNorm.transpose()<<std::endl;
                    std::cout<<"relative error="<<relErr.transpose()<<std::endl;
                }
                std::cout<<"convergenceError="<<cError<<std::endl;
            }
        }
        // Find immobile rate
        
        
        // Find diffusive-displacement rate
        
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveImmobileClusters()
    {
        
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solve()
    {
        solveMobileClusters();
        solveImmobileClusters();
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeSolver()
    {
        std::cout<<" decomposing"<<std::flush;
        std::array<bool,mSize> allComps;
        for(int k=0;k<mSize;++k)
        {
            allComps[k]=true;
        }
        
        mobileClusters.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mobileClustersIncrement.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mSolver.compute(mBWF); // call this after assigning the BCs
        solverInitialized=true;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        
        if(size_t(configIO.cdMatrix().size())==mobileClusters.gSize()+immobileClusters.gSize())
        {
            const size_t nNodes(mobileClusters.fe().nodes().size());
            mobileClusters=configIO.cdMatrix().block(0,0,nNodes,mSize).transpose().reshaped(mobileClusters.gSize(),1);
            immobileClusters=configIO.cdMatrix().block(0,mSize,nNodes,iSize).transpose().reshaped(immobileClusters.gSize(),1);
        }
        else
        {
            if(configIO.cdMatrix().size())
            {
                throw std::runtime_error("ClusterDynamics: TrialFunctions initializatoin size mismatch");
            }
        }
        
        
        
    }

    template<int dim>
    typename ClusterDynamicsFEM<dim>::VectorDim ClusterDynamicsFEM<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
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

    template struct ClusterDynamicsFEM<3>;

}
#endif



//template<int dim>
//void ClusterDynamicsFEM<dim>::applyBoundaryConditions()
//{
//
//    const auto& nodesInternalExternal(mobileClusters.fe().nodeList(nodeListInternalExternal));
//
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for(size_t k=0;k<nodesInternalExternal.size();++k)
//    {
//        const auto& node(nodesInternalExternal[k]);
//        const auto outNormal(node->outNormal()); // used to compute traction
//        const MatrixDim sigma(microstructures.stress(node->P0,node,nullptr,nullptr));
//        const double normalTraction(outNormal.dot(sigma*outNormal));
//        const auto bndConcentration(this->cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
//
//        VectorMSize otherConcentration(VectorMSize::Zero());
//        for(const auto& microstructure : this->microstructures)
//        {
//            if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
//            {// not the ClusterDynamics physics
//                otherConcentration += microstructure->mobileConcentration(node->P0,node,nullptr,nullptr);
//            }
//        }
//        for(int k=0;k<mSize;++k)
//        {
//            mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
//        }
//    }
//}
