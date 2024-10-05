/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformation_cpp_
#define model_ElasticDeformation_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ElasticDeformation.h>
#include <FEMfaceEvaluation.h>
#include <LinearWeakList.h>

namespace model
{

    bool kd(const int& i,const int& j)
    {
        return i==j;
    }

    template <int dim>
    std::unique_ptr<UniformController<SymmetricVoigtTraits<dim>::voigtSize>> getUniformEDcontroller(const DislocationDynamicsBase<dim>& ddBase)
    {
        
        typedef UniformController<SymmetricVoigtTraits<dim>::voigtSize> UniformControllerType;
        typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
        typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
        
        // Build stiffness matrix
        const double mu(ddBase.poly.mu);
        const double nu(ddBase.poly.nu);
        const double lam(2.0*nu*mu/(1.0-2.0*nu));
        MatrixVoigt C(MatrixVoigt::Zero()); // C_ijkl=lam*d_ij*d_kl+mu*(d_ik*d_jl+d_il*d_jk)
        for(int vI=0;vI<SymmetricVoigtTraits<dim>::voigtSize;++vI)
        {
            const auto& i(ddBase.voigtTraits.tensorIndex(vI,0));
            const auto& j(ddBase.voigtTraits.tensorIndex(vI,1));
            for(int vJ=0;vJ<SymmetricVoigtTraits<dim>::voigtSize;++vJ)
            {
                const auto& k(ddBase.voigtTraits.tensorIndex(vJ,0));
                const auto& l(ddBase.voigtTraits.tensorIndex(vJ,1));
                C(vI,vJ)= lam*kd(i,j)*kd(k,l)+mu*(kd(i,k)*kd(j,l)+kd(i,l)*kd(j,k));
            }
        }
        
        TextFileParser parser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt");
        
        const auto f0(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStress0",true));
        const auto f0Dot(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStressRate",true));
        const auto g0(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrain0",true));
        const auto g0Dot(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrainRate",true));
        const auto stiffnessRatio(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("stiffnessRatio",true));
        const double t0(0.0);
        return std::unique_ptr<UniformControllerType>(new UniformControllerType(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
    }

    template<int dim>
    ElasticDeformation<dim>::ElasticDeformation(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("ElasticDeformation",mc)
    /* init */,useElasticDeformationFEM(this->microstructures.ddBase.fe? bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("useElasticDeformationFEM",true)) : false )
    /* init */,elasticDeformationFEM(useElasticDeformationFEM?  new ElasticDeformationFEM<dim>(this->microstructures.ddBase) : nullptr)
    /* init */,uniformLoadController(useElasticDeformationFEM? std::unique_ptr<UniformControllerType>(nullptr) : getUniformEDcontroller(this->microstructures.ddBase))
    /* init */,inertiaReliefPenaltyFactor(uniformLoadController? 0.0 : TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("inertiaReliefPenaltyFactor",true))
    /* init */,ndA(this->microstructures.ddBase.fe? this->microstructures.ddBase.fe->template boundary<ExternalBoundary,imageTractionIntegrationOrder,GaussLegendre>() : TractionIntegrationDomainType())
    /* init */,tractionList(ndA.template integrationList<FEMfaceEvaluation<ElementType,dim,dim>>())
    /* init */,solverInitialized(false)
{
        
    }

    template<int dim>
    void ElasticDeformation<dim>::initializeSolver()
    {
        std::cout<<" gSize="<<elasticDeformationFEM->u.gSize()<<std::flush;
        elasticDeformationFEM->A.resize(elasticDeformationFEM->u.gSize(),elasticDeformationFEM->u.gSize());
        std::vector<Eigen::Triplet<double> > globalTriplets(elasticDeformationFEM->bWF.globalTriplets());
        if(inertiaReliefPenaltyFactor>0.0)
        {
            //                std::cout<<"Adding inertia relief triplets"<<std::endl;
            for(const auto& elePair : this->microstructures.ddBase.fe->elements())
            {
                const auto& ele(elePair.second);
                const auto C(elasticDeformationFEM->elementMomentumMatrix(ele));
                
                const auto Ct(C.block(  0,0,dim,dim*ElementType::nodesPerElement)); // rigid translation contraint
                const auto Cr(C.block(dim,0,dim,dim*ElementType::nodesPerElement)); // rigid rotation contraint
                
                //                    const auto CTC((C.transpose()*C).eval());
                const auto CTC((Ct.transpose()*Ct+Cr.transpose()*Cr).eval());
                
                for(int i=0;i<dim*ElementType::nodesPerElement;++i)
                {
                    for(int j=0;j<dim*ElementType::nodesPerElement;++j)
                    {
                        if(CTC(i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dim);
                            const size_t nodeDof_I(i%dim);
                            const size_t gI= ele.node(nodeID_I).gID*dim+nodeDof_I;
                            
                            const size_t  nodeID_J(j/dim);
                            const size_t nodeDof_J(j%dim);
                            const size_t gJ=ele.node(nodeID_J).gID*dim+nodeDof_J;
                            
                            globalTriplets.emplace_back(gI,gJ,inertiaReliefPenaltyFactor*CTC(i,j));
                        }
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error("inertiaReliefPenaltyFactor<=0.0");
        }
        std::cout<<", assembling"<<std::flush;
        elasticDeformationFEM->A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
        //            elasticDeformationFEM->A.prune(elasticDeformationFEM->A.norm()/elasticDeformationFEM->A.nonZeros(),FLT_EPSILON);
        std::cout<<", decomposing"<<std::flush;
        directSolver.compute(elasticDeformationFEM->A);
        if(directSolver.info()!=Eigen::Success)
        {
            throw std::runtime_error("Decomposing global stiffness matrix FAILED.");
        }
        solverInitialized=true;
    }

    template<int dim>
    void ElasticDeformation<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(uniformLoadController)
        {
            const MatrixDim apd(this->microstructures.averagePlasticDistortion());
            const MatrixDim aps(0.5*(apd+apd.transpose()));
            uniformLoadController->gs=this->microstructures.ddBase.voigtTraits.m2v(aps,true);
        }
        else
        {
            if(size_t(configIO.displacementMatrix().size())==elasticDeformationFEM->u.gSize()+elasticDeformationFEM->z.gSize())
            {
                std::cout<<" from configIO,"<<std::flush;
                const size_t nNodes(elasticDeformationFEM->u.fe().nodes().size());
                elasticDeformationFEM->u=configIO.displacementMatrix().block(0,0,nNodes,dim).transpose().reshaped(elasticDeformationFEM->u.gSize(),1);
                elasticDeformationFEM->z=configIO.displacementMatrix().block(0,dim,nNodes,dim).transpose().reshaped(elasticDeformationFEM->z.gSize(),1);
            }
            else
            {
                if(configIO.displacementMatrix().size())
                {
                    throw std::runtime_error("ElasticDeformation: TrialFunction u initializatoin size mismatch");
                }
            }
        }
    }

    template<int dim>
    void ElasticDeformation<dim>::solve()
    {
        if(uniformLoadController)
        {
            const MatrixDim apd(this->microstructures.averagePlasticDistortion());
            const MatrixDim aps(0.5*(apd+apd.transpose()));
            uniformLoadController->gs=this->microstructures.ddBase.voigtTraits.m2v(aps,true);
        }
        else
        {
            if(!solverInitialized)
            {
                initializeSolver();
                
            }
            //            elasticDeformationFEM->u.clearDirichletConditions();
            
            //
            //            lc->update(DN);
            //            lc->addDirichletConditions(DN);
            
            // Modify Dirichlet conditions
            //            std::vector<FEMnodeEvaluation<ElementType,dim,1>> fieldPoints;
            //            fieldPoints.reserve(displacement().dirichletNodeMap().size());
            //            for (const auto& pair : displacement().dirichletNodeMap())
            //            {
            //                const auto& gID(pair.first);
            //                fieldPoints.emplace_back(gID,fe.node(gID).P0);
            //            }
            //            for(const auto& microstructure : this->microstructures)
            //            {
            //                if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
            //                {// not the ElasticDeformation physics
            //                    qPoint.stress += microstructure->stress(qPoint.r);
            //                }
            //            }
            //            const auto t0= std::chrono::system_clock::now();
            std::cout<<", image traction"<<std::flush;
            //            auto ndA=this->microstructures.ddBase.fe->template boundary<ExternalBoundary,imageTractionIntegrationOrder,GaussLegendre>();
            //            auto eb_list = ndA.template integrationList<FEMfaceEvaluation<ElementType,dim,dim>>(); // TO DO: make this a member data to be able to output
            
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
            for(size_t k=0;k<tractionList.size();++k)
            {
                auto& pt(tractionList[k]);
                pt.setZero();
                for(const auto& microstructure : this->microstructures)
                {
                    if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
                    {// not the ElasticDeformation physics
                        pt += microstructure->stress(pt.P,nullptr,&pt.ele,nullptr);
                    }
                }
            }
            auto tractionWF=(test(elasticDeformationFEM->u),tractionList);
            std::cout<<", directSolver"<<std::flush;
            elasticDeformationFEM->u=directSolver.solve(-tractionWF.globalVector());
            
            // Compute zDot
            std::cout<<", zDot "<<std::flush;
            zDot.setZero(elasticDeformationFEM->z.gSize());
            for(const auto& node : elasticDeformationFEM->z.fe().nodes())
            {
                zDot.template segment<dim>(dim*node.gID)=this->microstructures.inelasticDisplacementRate(node.P0,&node,nullptr,nullptr);
            }
        }
    }

    template<int dim>
    void ElasticDeformation<dim>::updateConfiguration()
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(uniformLoadController)
        {// already updated in solve()
            
        }
        else
        {
            // Update z
            elasticDeformationFEM->z += zDot*this->timeSinceLastUpdate();
        }
    }

    template<int dim>
    double ElasticDeformation<dim>::getDt() const
    {
        return this->microstructures.ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ElasticDeformation<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(uniformLoadController)
        {
            const auto epsil(uniformLoadController->grad(this->microstructures.ddBase.simulationParameters.totalTime));
            const auto sigma(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime));
            
            f_file<<epsil.transpose()<<" "<<sigma.transpose()<<" ";
                        
            if(this->microstructures.ddBase.simulationParameters.runID==0)
            {
                for(int k=0;k<SymmetricVoigtTraits<dim>::voigtSize;++k)
                {
                    const auto& i(this->microstructures.ddBase.voigtTraits.tensorIndex(k,0)+1); // 1-based notation
                    const auto& j(this->microstructures.ddBase.voigtTraits.tensorIndex(k,1)+1); // 1-based notation
                    const std::string lab("e_"+std::to_string(i)+std::to_string(j));
                    F_labels<<lab<<"\n";
                }
                for(int k=0;k<SymmetricVoigtTraits<dim>::voigtSize;++k)
                {
                    const auto& i(this->microstructures.ddBase.voigtTraits.tensorIndex(k,0)+1); // 1-based notation
                    const auto& j(this->microstructures.ddBase.voigtTraits.tensorIndex(k,1)+1); // 1-based notation
                    const std::string lab("s_"+std::to_string(i)+std::to_string(j));
                    F_labels<<lab<<"\n";
                }
                F_labels<<std::endl;
            }
        }
        else
        {
            const size_t nNodes(elasticDeformationFEM->u.fe().nodes().size());
            configIO.displacementMatrix().resize(nNodes,dim+dim);
            configIO.displacementMatrix().block(0,0,nNodes,dim)=elasticDeformationFEM->u.dofVector().reshaped(dim,nNodes).transpose();
            configIO.displacementMatrix().block(0,dim,nNodes,dim)=elasticDeformationFEM->z.dofVector().reshaped(dim,nNodes).transpose();
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorDim ElasticDeformation<dim>::inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averagePlasticDistortion() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorDim ElasticDeformation<dim>::displacement(const VectorDim& x,const NodeType* const  node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->grad(this->microstructures.ddBase.simulationParameters.totalTime),true)*x;
        }
        else
        {
            if(node)
            {
                return eval(elasticDeformationFEM->u)(*node);
            }
            else
            {
                if(ele)
                {
                    return eval(elasticDeformationFEM->u)(*ele,ele->simplex.pos2bary(x));
                }
                else
                {
                    return eval(elasticDeformationFEM->u)(x,guess);
                }
            }
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::stress(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime),false);
        }
        else
        {
            if(node)
            {
                return this->microstructures.ddBase.voigtTraits.v2m(eval(elasticDeformationFEM->s)(*node),false);
            }
            else
            {
                if(ele)
                {
                    return this->microstructures.ddBase.voigtTraits.v2m(eval(elasticDeformationFEM->s)(*ele,ele->simplex.pos2bary(x)),false);
                }
                else
                {
                    return this->microstructures.ddBase.voigtTraits.v2m(eval(elasticDeformationFEM->s)(x,guess),false);
                }
            }
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorMSize ElasticDeformation<dim>::mobileConcentration(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorMSize::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averageStress() const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime),false);
        }
        else
        {
            return MatrixDim::Zero();
        }
    }

    template struct ElasticDeformation<3>;
}
#endif


//    template <int dim>
//    std::unique_ptr<ExternalLoadControllerBase<dim>> selectExternalLoadController(const DislocationDynamicsBase<dim>& ddBase)
//    {
//        if(ddBase.simulationParameters.externalLoadControllerName=="UniformExternalLoadController")
//        {
//            return std::unique_ptr<ExternalLoadControllerBase<dim>>(new UniformExternalLoadController<dim>(ddBase));
//        }
//        else
//        {
//            std::cout<<"Unknown externalLoadController name "<<ddBase.simulationParameters.externalLoadControllerName<<"No controller applied."<<std::endl;
//            return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
//        }
//    }


//template <int dim>
//typename VoigtTraits<dim,dim>::VoigtStorageType getVoigtOrder()
//{
//    typename VoigtTraits<dim,dim>::VoigtStorageType voigtOrder;
//    for(int i=0;i<dim;++i)
//    {
//        for(int j=0;j<dim;++j)
//        {
//            voigtOrder(dim*i+j,0)=i;
//            voigtOrder(dim*i+j,1)=j;
//        }
//    }
//    return voigtOrder;
//}

//template <int dim>
//std::unique_ptr<UniformController<dim,dim>> getUniformController(const DislocationDynamicsBase<dim>& ddBase, const VoigtTraits<dim,dim>& vt)
//{
//    typedef typename UniformController<dim,dim>::MatrixGrad MatrixGrad;
//    typedef typename UniformController<dim,dim>::MatrixGrad MatrixGrad;
//
//    const double mu(ddBase.poly.mu);
//    const double nu(ddBase.poly.nu);
//    const double lam(2.0*nu*mu/(1.0-2.0*nu));
//
//    MatrixGrad C(MatrixGrad::Zero()); // C_ijkl=lam*d_ij*d_kl+mu*(d_ik*d_jl+d_il*d_jk)
//    for(int vI=0;vI<VoigtTraits<dim,dim>::voigtSize;++vI)
//    {
//        const auto& i(vt.tensorIndex(vI,0));
//        const auto& j(vt.tensorIndex(vI,1));
//        for(int vJ=0;vJ<VoigtTraits<dim,dim>::voigtSize;++vJ)
//        {
//            const auto& k(vt.tensorIndex(vJ,0));
//            const auto& l(vt.tensorIndex(vJ,1));
//            C(vI,vJ)= lam*kd(i,j)*kd(k,l)+mu*(kd(i,k)*kd(j,l)+kd(i,l)*kd(j,k));
//        }
//    }
//
//    TextFileParser parser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/uniformExternalLoadControllerNew.txt");
//
//    SymmetricVoigtTraits<dim> svt(parser.readMatrix<size_t,SymmetricVoigtTraits<dim>::voigtSize,2>("voigtOrder",true));
//    const auto f0(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStress0",true),false)));
//    const auto f0Dot(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStressRate",true),false)));
//    const auto g0(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrain0",true),true)));
//    const auto g0Dot(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrainRate",true),true)));
//    const auto stiffnessRatio(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("stiffnessRatio",true),false)));
//    const double t0(0.0);
//    std::cout<<f0<<std::endl<<std::endl;
//    std::cout<<f0Dot<<std::endl<<std::endl;
//    std::cout<<g0<<std::endl<<std::endl;
//    std::cout<<g0Dot<<std::endl<<std::endl;
//
//    return std::unique_ptr<UniformController<dim,dim>>(new UniformController<dim,dim>(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
//}
