/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_CPP_
#define model_DislocationNetwork_CPP_

#include <numbers>

#include <DislocationNetwork.h>
#include <algorithm>


namespace model
{
    template <int dim, short unsigned int corder>
    DislocationNetwork<dim,corder>::DislocationNetwork(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("DislocationDynamics",mc)
    /* init */,glideStepsSinceLastClimb(0)
    /* init */,ddBase(this->microstructures.ddBase)
    /* init */,networkRemesher(*this)
    /* init */,junctionsMaker(*this)
    /* init */,crossSlipModel(DislocationCrossSlip<DislocationNetwork<dim,corder>>::getModel(ddBase.poly,ddBase.simulationParameters.traitsIO))
    /* init */,crossSlipMaker(*this)
    /* init */,nodeContractor(*this)
    /* init */,timeStepper(*this)
    /* init */,stochasticForceGenerator(ddBase.simulationParameters.use_stochasticForce? new StochasticForceGenerator(ddBase.simulationParameters.traitsIO) : nullptr)
    /* init */,computeDDinteractions(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("computeDDinteractions",true))
    /* init */,outputQuadraturePoints(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputQuadraturePoints",true))
    /* init */,outputLinkingNumbers(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputLinkingNumbers",true))
    /* init */,outputLoopLength(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputLoopLength",true))
    /* init */,outputSegmentPairDistances(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputSegmentPairDistances",true))
    /* init */,outputPlasticDistortionPerSlipSystem(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputPlasticDistortionPerSlipSystem",true))
    /* init */,computeElasticEnergyPerLength(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("computeElasticEnergyPerLength",true))
    /* init */,alphaLineTension(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("alphaLineTension",true))
    /* init */,use_velocityFilter(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("use_velocityFilter",true))
    /* init */,velocityReductionFactor(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("velocityReductionFactor",true))
    /* init */,verboseDislocationNode(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseDislocationNode",true))
    {
        assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        
        LoopNetworkType::verboseLevel=TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseLoopNetwork",false);
        verboseDislocationNetwork=TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseDislocationNetwork",false);
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim, corder>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream&,const std::ofstream&)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        LoopType::initFromFile(ddBase.simulationParameters.traitsIO.ddFile);
        LoopNodeType::initFromFile(ddBase.simulationParameters.traitsIO.ddFile);
        LoopLinkType::initFromFile(ddBase.simulationParameters.traitsIO.ddFile);
        NetworkLinkType::initFromFile(ddBase.simulationParameters.traitsIO.ddFile);
        DislocationFieldBase<dim>::initFromFile(ddBase.simulationParameters.traitsIO.ddFile);
        
        glideSolver=DislocationGlideSolverFactory<DislocationNetwork<dim,corder>>::getGlideSolver(*this,TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readString("glideSolverType",true));
        climbSolver=DislocationClimbSolverFactory<DislocationNetwork<dim,corder>>::getClimbSolver(*this,TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readString("climbSolverType",true));
        _inclusions=this->microstructures.template getUniqueTypedMicrostructure<InclusionMicrostructure<dim>>();
                
        setConfiguration(configIO);
        updateGeometry();
    }

    //New Version
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::setConfiguration(const DDconfigIO<dim>& evl)
    {
        this->loopLinks().clear(); // erase base network to clear current config
        
        // Create Loops
        std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
        size_t loopNumber=1;
        for(const auto& loop : evl.loops())
        {
            const bool faulted(ddBase.poly.grain(loop.grainID).singleCrystal->rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false);
            VerboseDislocationNetwork(1,"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;);

            //std::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
            const size_t loopIDinFile(loop.sID);
            LoopType::set_count(loopIDinFile);
            
            GlidePlaneKey<dim> loopPlaneKey(loop.P, ddBase.poly.grain(loop.grainID).singleCrystal->reciprocalLatticeDirection(loop.N));
            tempLoops.push_back(this->loops().create(loop.B, ddBase.glidePlaneFactory.getFromKey(loopPlaneKey)));
            assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
            loopNumber++;
            
            
            //            switch (loop.loopType)
            //            {
            //                case DislocationLoopIO<dim>::GLISSILELOOP:
            //                {
            //                    GlidePlaneKey<dim> loopPlaneKey(loop.P, ddBase.poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
            //                    tempLoops.push_back(this->loops().create(loop.B, ddBase.glidePlaneFactory.getFromKey(loopPlaneKey)));
            //                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
            //                    loopNumber++;
            //                    break;
            //                }
            //                case DislocationLoopIO<dim>::SESSILELOOP:
            //                {
            //                    tempLoops.push_back(this->loops().create(loop.B,loop.grainID,loop.loopType ));
            //                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
            //                    loopNumber++;
            //                    break;
            //                }
            //                default:
            //                    assert(false && "Unknown DislocationLoop type");
            //                    break;
            //            }
        }
        
        // Create NetworkNodes
        std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
        size_t netNodeNumber=1;
        for(const auto& node : evl.nodes())
        {
            VerboseDislocationNetwork(1,"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;);
            const size_t nodeIDinFile(node.sID);
            NetworkNodeType::set_count(nodeIDinFile);
            tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.climbVelocityScalar,node.velocityReduction));
            assert(this->networkNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            netNodeNumber++;
        }
        
        // Create LoopNodes
        std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
        size_t loopNodeNumber=1;
        for(const auto& node : evl.loopNodes())
        {
            VerboseDislocationNetwork(1,"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::endl;);
            const size_t nodeIDinFile(node.sID);
            LoopNodeType::set_count(nodeIDinFile);
            const auto loop(this->loops().get(node.loopID));
            const auto netNode(this->networkNodes().get(node.networkNodeID));
            const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
            const auto periodicPatchEdge((periodicPatch && node.edgeIDs.first>=0)? (node.edgeIDs.second>=0 ? std::make_pair(periodicPatch->edges()[node.edgeIDs.first],
                                                                                                                            periodicPatch->edges()[node.edgeIDs.second]):
                                                                                    std::make_pair(periodicPatch->edges()[node.edgeIDs.first],nullptr)):std::make_pair(nullptr,nullptr));
            tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
            assert(this->loopNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            loopNodeNumber++;
        }
        
        // Insert Loops
        std::map<size_t,std::map<size_t,size_t>> loopMap;
        for(const auto& looplink : evl.loopLinks())
        {// Collect LoopLinks by loop IDs
            loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
        }
        assert(loopMap.size()==evl.loops().size());
        
        for(const auto& loop : evl.loops())
        {// for each loop in the DDconfigIO<dim> object
            
            const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
            assert(loopFound!=loopMap.end());
            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
            loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first));
            for(size_t k=0;k<loopFound->second.size();++k)
            {
                const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                if(k<loopFound->second.size()-1)
                {
                    loopNodes.push_back(this->loopNodes().get(nodeFound->second));
                }
                else
                {
                    assert(nodeFound->second==loopNodes[0]->sID);
                }
            }
            //        std::cout<<" Inserting loop "<<loop.sID<<std::endl;
            this->insertLoop(this->loops().get(loop.sID),loopNodes);
        }
        updateGeometry();
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::updateGeometry()
    {
        VerboseDislocationNetwork(2,"DislocationNetwork::updateGeometry"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateGeometry();
        }
        // updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateGeometry DONE"<<std::endl;);
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::averagePlasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->averagePlasticDistortion();
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::averagePlasticStrain() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        const MatrixDim apd(averagePlasticDistortion());
        return 0.5*(apd+apd.transpose());
    }

    template <int dim, short unsigned int corder>
    std::map<std::pair<int,int>,double> DislocationNetwork<dim,corder>::slipSystemAveragePlasticDistortion() const
    {
        std::map<std::pair<int,int>,double> temp; // <grainID,slipSystemID>
        for(const auto& weakloop : this->loops())
        {
            const auto loop(weakloop.second.lock());
            if(loop->slipSystem())
            {
                const double val = loop->slipSystem()->unitSlip.transpose()*loop->averagePlasticDistortion()*loop->slipSystem()->unitNormal;
                const std::pair<int,int> key(std::make_pair(loop->grain.grainID,loop->slipSystem()->sID));
                auto iter(temp.find(key));
                if(iter!=temp.end())
                {// iter found
                    iter->second += val ;
                }
                else
                {// iter not found
                    temp.emplace(key,val);
                }
            }
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::averagePlasticDistortionRate() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->averagePlasticDistortionRate();
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::averagePlasticStrainRate() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        const MatrixDim apdr(averagePlasticDistortionRate());
        return 0.5*(apdr+apdr.transpose());
    }

    template <int dim, short unsigned int corder>
    std::tuple<double,double,double,double> DislocationNetwork<dim,corder>::networkLength() const
    {/*!\returns the total line length of the DislocationNetwork. The return
      * value is a tuple, where the first value is the length of bulk glissile
      * dislocations, the second value is the length of bulk sessile
      * dislocations, and the third value is the length accumulated on
      * the mesh boundary.
      */
        double bulkGlissileLength(0.0);
        double bulkSessileLength(0.0);
        double boundaryLength(0.0);
        double grainBoundaryLength(0.0);
        
        for(auto& loop : this->loops())
        {
            for(const auto& loopLink : loop.second.lock()->loopLinks())
            {
                if(loopLink->networkLink())
                {
                    if(!loopLink->networkLink()->hasZeroBurgers())
                    {
                        if(loopLink->networkLink()->isBoundarySegment())
                        {
                            boundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else if(loopLink->networkLink()->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else
                        {
                            if(loopLink->networkLink()->isSessile())
                            {
                                bulkSessileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                        }
                    }
                }
            }
        }
        return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
    }


    template <int dim, short unsigned int corder>
    const std::shared_ptr<InclusionMicrostructure<dim>>& DislocationNetwork<dim,corder>::inclusions() const
    {
        return _inclusions;
    }

    template <int dim, short unsigned int corder>
    bool DislocationNetwork<dim,corder>::contract(std::shared_ptr<NetworkNodeType> nA,
                                                  std::shared_ptr<NetworkNodeType> nB)
    {
        return nodeContractor.contract(nA,nB);
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::VectorDim DislocationNetwork<dim,corder>::displacement(const VectorDim& x,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        VectorDim temp(VectorDim::Zero());
        
        for(const auto& loop : this->loops())
        {// sum solid angle of each loop
            for(const auto& shift : ddBase.periodicShifts)
            {
                temp-=loop.second.lock()->solidAngle(x+shift)/4.0/std::numbers::pi*loop.second.lock()->burgers();
            }
        }
        
        for(const auto& link : this->networkLinks())
        {// sum line-integral part of displacement field per segment
            if(   !link.second.lock()->hasZeroBurgers())
            {
                for(const auto& shift : ddBase.periodicShifts)
                {
                    temp+=link.second.lock()->straight.displacement(x+shift);
                }
            }
        }
        
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::averageStress() const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        return MatrixDim::Zero();
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::VectorMSize DislocationNetwork<dim,corder>::mobileConcentration(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        VectorMSize temp(VectorMSize::Zero());
        if(climbSolver)
        {
            const auto pointGrains(this->pointGrains(x,node,ele,guess));
//            std::set<const Grain<dim>*> pointGrains;
//            if(node)
//            {
//                for(const auto& nodeEle : *node)
//                {
//                    pointGrains.emplace(&ddBase.poly.grain(nodeEle->simplex.region->regionID));
//                }
//            }
//            else
//            {
//                if(ele)
//                {
//                    pointGrains.emplace(&ddBase.poly.grain(ele->simplex.region->regionID));
//                }
//                else
//                {
//                    const std::pair<bool,const Simplex<dim,dim>*> found(ddBase.mesh.searchWithGuess(x,guess));
//                    if(found.first)
//                    {
//                        pointGrains.emplace(&ddBase.poly.grain(found.second->region->regionID));
//                    }
//                }
//            }
            
            if(pointGrains.size())
            {
                for(const auto& link : this->networkLinks())
                {// sum stress field per segment
                    if(   !link.second.lock()->hasZeroBurgers()
                       )
                    {
                        
                        const auto segment(link.second.lock());
                        const auto segGrains(segment->grains());
                        std::set<const Grain<dim>*> intersect;
                        std::set_intersection(pointGrains.begin(), pointGrains.end(), segGrains.begin(), segGrains.end(),
                                              std::inserter(intersect, intersect.begin()));
                        if(intersect.size()==1)
                        {
//                            const int grainID((*intersect.begin())->region.regionID);
//                            StressStraight<3> ss(ddBase.poly,segment->source->get_P(),segment->sink->get_P(),segment->burgers(),ddBase.EwaldLength);
//                            for(const auto& shift : ddBase.periodicShifts)
//                            {
                                temp+=segment->clusterConcentration(x,climbSolver->CD->cdp);

//                                temp+=ss.clusterConcentration(x+shift,grainID, segment->source->climbDirection(), segment->source->climbVelocityScalar , segment->sink->climbDirection(), segment->sink->climbVelocityScalar, climbSolver->CD->cdp);
//                            }
                        }
                    }
                }
            }
        }
        return temp;
    }


    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::stress(const VectorDim& x,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& link : this->networkLinks())
        {// sum stress field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               )
            {
                for(const auto& shift : ddBase.periodicShifts)
                {
                    temp+=link.second.lock()->straight.stress(x+shift);
                }
            }
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim, corder>::VectorDim DislocationNetwork<dim, corder>::inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim, corder>::updateConfiguration()
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        moveNodes(ddBase.simulationParameters.dt);
        executeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
        if(isClimbStep())
        {
            glideStepsSinceLastClimb=0;
        }
        else
        {
            glideStepsSinceLastClimb++;
        }
        updateGeometry();
    }

    template <int dim, short unsigned int corder>
    double DislocationNetwork<dim, corder>::getDt() const
    {
        if(isClimbStep())
        {
            return timeStepper.getDt(climbSolver->vClimbRef,1.0e-7);
            
        }
        else
        {
            //        climbSolver->CD->cdp
            return timeStepper.getDt(ddBase.poly.cs,1.0e-7);
        }
    }

    template <int dim, short unsigned int corder>
    bool DislocationNetwork<dim, corder>::isClimbStep() const
    {
        if(glideSolver)
        {
            return climbSolver?  (glideStepsSinceLastClimb && averagePlasticDistortionRate().norm()<climbSolver->glideEquilibriumRate): false;
        }
        else
        {
            return true;
        }
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim, corder>::solve()
    {
        
        const bool isClimbingStep(isClimbStep());
if(isClimbingStep)
{
    std::cout<<" climbStep"<<std::flush;
}
else
{
    std::cout<<" glideStep"<<std::flush;
}
        double maxVelocity = 0.0;
        for (const auto &nodeIter : this->networkNodes())
        {
            const double vNorm(nodeIter.second.lock()->get_V().norm());
            if (vNorm > maxVelocity)
            {
                maxVelocity = vNorm;
            }
        }
        
    #ifdef _OPENMP
        const size_t nThreads = omp_get_max_threads();
    #else
        const size_t nThreads = 1;
    #endif
        
        //! -1 Compute the interaction StressField between dislocation particles
        std::map<int, int> velocityBinMap;
        for (const auto &binVal : ddBase.simulationParameters.subcyclingBins)
        {
            velocityBinMap.emplace(binVal, 0);
        }
        
        if (corder == 0)
        { // For straight segments use analytical expression of stress field
            std::cout <<" creating qPoints "<< std::flush;
            for (const auto &links : this->networkLinks())
            {
                
                const int velGroup((ddBase.simulationParameters.useSubCycling && !isClimbingStep) ? links.second.lock()->velocityGroup(maxVelocity, ddBase.simulationParameters.subcyclingBins) : 1);
                auto velocityBinIter(velocityBinMap.find(velGroup));
                assert(velocityBinIter != velocityBinMap.end());
                velocityBinIter->second++;
                
                if ((ddBase.simulationParameters.runID % velGroup) == 0)
                {
                    links.second.lock()->createQuadraturePoints(isClimbingStep);
                }
            }
            
            std::cout <<" ,updating qPoints (" << nThreads << " threads) " << std::flush;
    #ifdef _OPENMP
    #pragma omp parallel for
            for (size_t k = 0; k < this->networkLinks().size(); ++k)
            {
                auto linkIter(this->networkLinks().begin());
                std::advance(linkIter, k);
                const int velGroup((ddBase.simulationParameters.useSubCycling && !isClimbingStep) ? linkIter->second.lock()->velocityGroup(maxVelocity, ddBase.simulationParameters.subcyclingBins) : 1);
                
                if ((ddBase.simulationParameters.runID % velGroup) == 0)
                {
                    linkIter->second.lock()->updateQuadraturePoints(isClimbingStep);
                }
                //            else
                //            {
                //                linkIter->second.lock()->assembleGlide(false);
                //            }
            }
    #else
            for (auto &linkIter : this->networkLinks())
            {
                const int velGroup(ddBase.simulationParameters.useSubCycling ? linkIter.second.lock()->velocityGroup(maxVelocity, ddBase.simulationParameters.subcyclingBins) : 1);
                
                if ((ddBase.simulationParameters.runID % velGroup) == 0)
                {
                    linkIter.second.lock()->updateQuadraturePoints(isClimbingStep);
                }
                //            else
                //            {
                //                linkIter.second.lock()->assembleGlide(false);
                //            }
            }
    #endif
        }
        else
        { // For curved segments use quandrature integration of stress field
            //        assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
            throw std::runtime_error("DislocationNetwork::SolveNodalVelocities not implemented for corder>0.");
        }
        
        //    const DislocationVelocitySolverBase<DislocationNetwork<dim,corder>>* const vSolver(isClimbingStep? static_cast<const DislocationVelocitySolverBase<DislocationNetwork<dim,corder>>*>(climbSolver.get())
        //                                                                                                  : static_cast<const DislocationVelocitySolverBase<DislocationNetwork<dim,corder>>*>(glideSolver.get()));
        //
        //    if(vSolver)
        //    {
        //        const Eigen::VectorXd X(vSolver->getNodeVelocities());
        //        if(int(NdofXnode*this->networkNodes().size())==X.size())
        //        {
        //            size_t k=0;
        //            for (auto& networkNode : this->networkNodes())
        //            {
        ////                std::cout<<X.segment(NdofXnode*k,NdofXnode)<<std::endl;
        //                networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode)); // double cast to remove some numerical noise
        //                ++k;
        //            }
        //        }
        //        else
        //        {
        //            std::cout<<"NdofXnode*this->networkNodes().size()="<<NdofXnode*this->networkNodes().size()<<std::endl;
        //            std::cout<<"vSolver->getNodeVelocities().size()="<<X.size()<<std::endl;
        //            throw std::runtime_error("vSolver returned wrong velocity vector size.");
        //        }
        //
        //        VerboseDislocationNetwork(2,"DislocationNetwork::updateRates"<<std::endl;);
        //        for(auto& loop : this->loops())
        //        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
        //            loop.second.lock()->updateRates();
        //        }
        //        // updatePlasticDistortionRateFromAreas();
        //        VerboseDislocationNetwork(3,"DislocationNetwork::updateRates DONE"<<std::endl;);
        //
        //    }
        //    else
        //    {
        ////        std::cout<<"No vSolver"<<std::endl;
        //    }
        
        Eigen::VectorXd X(Eigen::VectorXd::Zero(0));
        if(isClimbingStep)
        {
            climbSolver->computeClimbScalarVelocities();
            size_t k=0;
            for (auto& networkNode : this->networkNodes())
            {
                networkNode.second.lock()->climbVelocityScalar=climbSolver->scalarVelocities()[k];
                ++k;
            }
            X=climbSolver->getNodeVelocities();
        }
        else
        {
            if(glideSolver)
            {
                size_t k=0;
                for (auto& networkNode : this->networkNodes())
                {
                    networkNode.second.lock()->climbVelocityScalar.setZero();
                    ++k;
                }
                X=glideSolver->getNodeVelocities();
            }
            
        }
        if(int(NdofXnode*this->networkNodes().size())==X.size())
        {
            size_t k=0;
            for (auto& networkNode : this->networkNodes())
            {
                networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode)); // double cast to remove some numerical noise
                ++k;
            }
        }
        else
        {
            std::cout<<"NdofXnode*this->networkNodes().size()="<<NdofXnode*this->networkNodes().size()<<std::endl;
            std::cout<<"vSolver->getNodeVelocities().size()="<<X.size()<<std::endl;
            throw std::runtime_error("vSolver returned wrong velocity vector size.");
        }
        
        VerboseDislocationNetwork(2,"DislocationNetwork::updateRates"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateRates();
        }
        // updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateRates DONE"<<std::endl;);
        
        storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
        
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::moveNodes(const double & dt_in)
    {/*! Moves all nodes in the DislocationNetwork using the stored glide velocity and current dt
      */
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
        danglingBoundaryLoopNodes.clear();
        for(auto& node : this->networkNodes())
        {
            node.second.lock()->trySet_P(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in);
        }
        updateBoundaryNodes();
        //    if(ddBase.isPeriodicDomain)
        //    {
        //        danglingBoundaryLoopNodes.clear();
        //        for(auto& node : this->networkNodes())
        //        {
        //            node.second.lock()->trySet_P(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in);
        //        }
        //        updateBoundaryNodes();
        //    }
        //    else
        //    {
        //        for (auto& nodeIter : this->networkNodes())
        //        {
        //            nodeIter.second.lock()->set_P(nodeIter.second.lock()->get_P()+nodeIter.second.lock()->get_V()*dt_in);
        //        }
        //    }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::storeSingleGlideStepDiscreteEvents(const long int&)
    {
        crossSlipMaker.findCrossSlipSegments();
        crossSlipMaker.execute(); // this is now performed before moving, since internally it computes forces and velocities on the new segments
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::executeSingleGlideStepDiscreteEvents(const long int& runID)
    {
        
        //    crossSlipMaker.execute();
        //        //! 13- Node redistribution
        networkRemesher.remesh(runID);
        //        //! 12- Form Junctions
        junctionsMaker.formJunctions(3.0*networkRemesher.Lmin);
        //Calling remesh again so that any other topological changes created by junctions whihc are otherwise removable can be removed
        //        //! 13- Node redistribution
        //        this->io().output(runID);
        
        networkRemesher.remesh(runID);
        //        updateVirtualBoundaryLoops();
        
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        
        for(const auto& loop : this->loops())
        {
            configIO.loops().emplace_back(*loop.second.lock());
        }
        
        // Store LoopNodes
        for(const auto& node : this->loopNodes())
        {
            configIO.loopNodes().emplace_back(*node.second.lock());
        }
        
        // Store LoopLinks
        for(const auto& link : this->loopLinks())
        {
            configIO.loopLinks().emplace_back(link.second);
        }
        
        // Store NetworkNodes
        for(const auto& node : this->networkNodes())
        {
            configIO.nodes().emplace_back(*node.second.lock());
        }
        
//        for(const auto& node : this->polyhedronInclusionNodes())
//        {
//            configIO.polyhedronInclusionNodes().emplace_back(node.second);
//        }
//        
//        // Store Eshelby Inclusions
//        for(const auto& ei : this->eshelbyInclusions())
//        {
//            
//            auto* sphericalDerived = dynamic_cast<SphericalInclusion<dim>*>(ei.second.get());
//            if (sphericalDerived)
//            {
//                configIO.sphericalInclusions().emplace_back(*sphericalDerived);
//            }
//            
//            auto* polyhedronDerived = dynamic_cast<PolyhedronInclusion<dim>*>(ei.second.get());
//            if (polyhedronDerived)
//            {
//                configIO.polyhedronInclusions().emplace_back(*polyhedronDerived);
//                for(const auto& face : polyhedronDerived->faces)
//                {
//                    for(size_t k=0;k<face.second.size();++k)
//                    {
//                        const size_t k1(k<face.second.size()-1? k+1 : 0);
//                        configIO.polyhedronInclusionEdges().emplace_back(polyhedronDerived->sID,face.first,face.second[k].first,face.second[k1].first);
//                    }
//                }
//            }
//        }
        
        
        // AuxIO
        if (this->outputQuadraturePoints)
        {
            for (const auto& link : this->networkLinks())
            {
                for(const auto& qPoint : link.second.lock()->quadraturePoints())
                {
                    auxIO.quadraturePoints().push_back(qPoint);
                }
            }
        }
        
        const std::tuple<double,double,double,double> length(this->networkLength());
        const double densityFactor(1.0/this->ddBase.mesh.volume()/std::pow(this->ddBase.poly.b_SI,2));
        f_file<<std::get<0>(length)*densityFactor<<" "<<std::get<1>(length)*densityFactor<<" "<<std::get<2>(length)*densityFactor<<" "<<std::get<3>(length)*densityFactor<<" ";
        if(ddBase.simulationParameters.runID==0)
        {
            F_labels<<"glissile density [m^-2]\n";
            F_labels<<"sessile density [m^-2]\n";
            F_labels<<"boundary density [m^-2]\n";
            F_labels<<"grain boundary density [m^-2]\n";
        }
        
        if(this->outputPlasticDistortionPerSlipSystem)
        {
            const auto ssapd(slipSystemAveragePlasticDistortion());
            for(const auto& pair : ssapd)
            {
                f_file<<pair.second<<" ";
                if(ddBase.simulationParameters.runID==0)
                {
                    F_labels<<"betaP_g"+std::to_string(pair.first.first)+"_ss"+std::to_string(pair.first.second)+" [-]\n";
                }
            }
        }
        
        if(this->computeElasticEnergyPerLength)
        {
            double eE(0.0);
            double eC(0.0);
            
            for(const auto& linkIter : this->networkLinks())
            {// Collect LoopLinks by loop IDs
                const auto link(linkIter.second.lock());
                for(const auto& qPoint : link->quadraturePoints())
                {
                    eE+=qPoint.elasticEnergyPerLength*qPoint.dL;
                    eC+=qPoint.coreEnergyPerLength*qPoint.dL;
                }
            }
            f_file<<eE<<" "<<eC<<" ";
            if(ddBase.simulationParameters.runID==0)
            {
                F_labels<<"dislocation elastic energy [mu b^3]\n";
                F_labels<<"dislocation core energy [mu b^3]\n";
            }
        }
    }

    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::updateBoundaryNodes()
    {
        
        /*!Step 1. Before removing populate the junction information
         *       This populates the network node where the junctions are needed to be preserved
         *       To pupulate junction information we create a map with key=pair(sourceNetNode,sinkNetNode) and val = set(loopIDs passing through nodes)
         */
        std::map<std::pair<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>>,std::set<size_t>> networkNodeLoopMap; //Size_t corresponds to the loopID
        for (const auto& ln : this->loopNodes())
        {
            const auto sharedLNptr(ln.second.lock());
            if (sharedLNptr->periodicPlaneEdge.first)
            {//Node on a boundary: possible junction is between sharedLNptr->periodicPrev() and sharedLNptr->periodicNext()
                if (sharedLNptr->networkNode->loopNodes().size()>1)
                {//junction node
                    const auto loopsThis (sharedLNptr->networkNode->loopIDs());
                    
                    const LoopNodeType *pPrev(sharedLNptr->periodicPrev());
                    const LoopNodeType *pNext(sharedLNptr->periodicNext());
                    
                    
                    const auto pPrevNetwork (pPrev->networkNode);
                    const auto pNextNetwork (pNext->networkNode);
                    
                    assert(pPrevNetwork!=nullptr);
                    assert(pNextNetwork!=nullptr);
                    
                    const auto loopspPrev(pPrevNetwork->loopIDs());
                    const auto loopspNext(pNextNetwork->loopIDs());
                    
                    std::set<size_t> tempPrev;
                    std::set<size_t> tempNext;
                    std::set_intersection(loopspPrev.begin(), loopspPrev.end(), loopsThis.begin(), loopsThis.end(), std::inserter(tempPrev, tempPrev.begin()));
                    std::set_intersection(loopsThis.begin(), loopsThis.end(), loopspNext.begin(), loopspNext.end(), std::inserter(tempNext, tempNext.begin()));
                    
                    if (tempPrev!=tempNext)
                    {
                        std::cout<<"For bnd network node"<<sharedLNptr->networkNode->sID<<" loops are "<<std::flush;
                        for (const auto& loop : loopsThis)
                        {
                            std::cout<<loop<<", ";
                        }
                        std::cout<<std::endl;
                        
                        std::cout<<"For prev network node"<<pPrevNetwork->sID<<" loops are "<<std::flush;
                        for (const auto& loop : loopspPrev)
                        {
                            std::cout<<loop<<", ";
                        }
                        std::cout<<std::endl;
                        
                        std::cout<<"For next network node"<<pNextNetwork->sID<<" loops are "<<std::flush;
                        for (const auto& loop : loopspNext)
                        {
                            std::cout<<loop<<", ";
                        }
                        std::cout<<std::endl;
                        throw std::runtime_error("BND node must have the same loops as the common loops between the internal nodes");
                    }
                    else
                    {
                        if (pPrevNetwork->sID < pNextNetwork->sID)
                        {
                            networkNodeLoopMap.emplace(std::make_pair(pPrevNetwork, pNextNetwork), tempPrev);
                        }
                        else
                        {
                            networkNodeLoopMap.emplace(std::make_pair(pNextNetwork, pPrevNetwork), tempPrev);
                        }
                    }
                }
            }
            else
            {//Node not on a boundary: possible junction is between sharedLNptr and sharedLNptr->periodicNext()
                if (sharedLNptr->networkNode->loopNodes().size()>1)
                {// a junction node
                    if (sharedLNptr->boundaryNext().size()==0 && sharedLNptr->periodicNext()->networkNode->loopNodes().size()>1)
                    {
                        if (sharedLNptr->periodicPlanePatch()!=sharedLNptr->periodicNext()->periodicPlanePatch())
                        {// Junction is across boundary
                            const auto netLink (sharedLNptr->next.second->networkLink());
                            if(netLink)
                            {
                                std::set<size_t> netLinkLoopIDs (netLink->loopIDs());
                                if (netLink->loopLinks().size()>=2)
                                {
                                    //a junction node moving out
                                    if (sharedLNptr->networkNode->sID < sharedLNptr->periodicNext()->networkNode->sID)
                                    {
                                        networkNodeLoopMap.emplace(std::make_pair(sharedLNptr->networkNode, sharedLNptr->periodicNext()->networkNode), netLinkLoopIDs);
                                    }
                                    else
                                    {
                                        networkNodeLoopMap.emplace(std::make_pair(sharedLNptr->periodicNext()->networkNode,sharedLNptr->networkNode), netLinkLoopIDs);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*!Step 2. Remove nodes in danglingBoundaryLoopNodes
         */
        if (danglingBoundaryLoopNodes.size())
        {
            //       std::cout << "Removing bnd Nodes" << std::endl;
            VerboseDislocationNetwork(1, "Removing bnd Nodes"<< std::endl;);
            
            for (const auto &node : danglingBoundaryLoopNodes)
            {
                this->removeLoopNode(node->sID);
            }
        }
        
        //    std::cout << "Inserting new boundary nodes" << std::endl;
        VerboseDislocationNetwork(1, "Inserting new boundary nodes"<< std::endl;);
        
        /*!Step 3. Create a map of new boundary NetworkNode to be inserted
         *       key = tuple<souceNetNode,
         *                sinkNetNode,
         *                set<all mesh faces that current link is crossing>,
         *                set<current mesh faces contining new boundary NetworkNode>,
         *                number of mesh intersection crossed by link beyond  new boundary NetworkNode>
         *       value= new  boundary NetworkNode
         */
        std::map<std::tuple<const std::shared_ptr<NetworkNodeType>, const std::shared_ptr<NetworkNodeType>, std::set<const PlanarMeshFace<dim> *>, std::set<const PlanarMeshFace<dim> *>,size_t>, std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
        for (const auto &weakLoop : this->loops())
        {
            
            const auto loop(weakLoop.second.lock());
            if (loop->periodicGlidePlane)
            {
                VerboseDislocationNetwork(1, " DisloationNetwork::DislocationLoop " << loop->sID << " Updating boundary nodes " << std::endl;);
                
                std::vector<std::pair<VectorLowerDim, const LoopNodeType *const>> loopNodesPos;
                std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *> bndNodesMap;
                for (const auto &loopLink : loop->linkSequence())
                {
                    if (!loopLink->source->periodicPlaneEdge.first)
                    {
                        loopNodesPos.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(loopLink->source->get_P()), loopLink->source.get());
                    }
                    else
                    {
                        
                        if (loopLink->source->periodicPlaneEdge.second)
                        {// loopLink->source is a corner boundary node. For nodes intersecting the corner, storage should be from smaller edgeId and larger edgeID
                            if (loopLink->source->periodicPlaneEdge.first->edgeID<loopLink->source->periodicPlaneEdge.second->edgeID)
                            {
                                bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.first.get(), loopLink->source->periodicPlaneEdge.second.get())), loopLink->source.get());
                            }
                            else
                            {
                                bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.second.get(),loopLink->source->periodicPlaneEdge.first.get())), loopLink->source.get());
                            }
                        }
                        else
                        {
                            bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.first.get(), loopLink->source->periodicPlaneEdge.second.get())), loopLink->source.get());
                        }
                    }
                }
                
                if (loopNodesPos.size())
                {
                    const auto polyInt(loop->periodicGlidePlane->polygonPatchIntersection(loopNodesPos)); // Note: last element of polyInt is always an internal node
                    // 2dPos,shift,edgeIDs ,set of edgeiD(all edges crossed),size_t(number of boundary nodes on same loop link after current node),LoopNodeType* (from loopNodesPos if internal, nullptr otherwise)
                    //edges still needed to be traversed will be reverse of the value if junction formation includes link in the opposite direction
                    std::map<const LoopNodeType *, size_t> polyIntMap; // ID of internal loopNodes into polyInt
                    for (size_t p = 0; p < polyInt.size(); ++p)
                    {
                        if (std::get<5>(polyInt[p]))
                        {// an internal node
                            polyIntMap.emplace(std::get<5>(polyInt[p]), p);
                        }
                    }
                    
                    for (size_t k = 0; k < loopNodesPos.size(); ++k)
                    {
                        const auto periodicPrev(loopNodesPos[k].second);
                        const auto periodicPrevNetwork(periodicPrev->networkNode);
                        const auto polyIter(polyIntMap.find(periodicPrev));
                        assert(polyIter != polyIntMap.end());
                        const size_t p(polyIter->second);
                        
                        const size_t k1(k < loopNodesPos.size() - 1 ? k + 1 : 0);
                        const auto periodicNext(loopNodesPos[k1].second);
                        const auto periodicNextNetwork(periodicNext->networkNode);
                        const auto polyIter1(polyIntMap.find(periodicNext));
                        assert(polyIter1 != polyIntMap.end());
                        const size_t p1(polyIter1->second);
                        
                        const auto periodicNetworkSource(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicPrevNetwork : periodicNextNetwork);
                        const auto periodicNetworkSink(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicNextNetwork : periodicPrevNetwork);
                        
                        VerboseDislocationNetwork(2, " PeriodicNetworkSource " << periodicNetworkSource->sID << std::endl;);
                        VerboseDislocationNetwork(2, " PeriodicNetworkSink " << periodicNetworkSink->sID << std::endl;);
                        // std::cout<<"periodicPrev->periodicNext"<<periodicPrev->tag()<<"==>"<<periodicNext->tag()<<std::endl;
                        
                        const LoopNodeType *currentSource(periodicPrev);
                        size_t p2 = (p + 1) % polyInt.size();
                        while (p2 < p1)
                        {
                            const auto periodicPatch(loop->periodicGlidePlane->getPatch(std::get<1>(polyInt[p2])));
                            const auto periodicPatchEdge(std::get<2>(polyInt[p2]).second < 0 ? std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], nullptr) : std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], periodicPatch->edges()[std::get<2>(polyInt[p2]).second]));
                            
                            VerboseDislocationNetwork(2, " First EdgeID on periodicPatchEdge " << periodicPatchEdge.first->edgeID << std::endl;);
                            if (periodicPatchEdge.second)
                            {
                                VerboseDislocationNetwork(2, " Second EdgeID on periodicPatchEdge " << periodicPatchEdge.second->edgeID << std::endl;);
                            }
                            
                            //                            typename std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *>::iterator bndIter(bndNodesMap.end());
                            auto bndIter(bndNodesMap.end());
                            
                            //For the case where the second edge exists, the pair of edge should be from minimum to maximum.
                            if (periodicPatchEdge.second)
                            {
                                if (periodicPatchEdge.first->edgeID < periodicPatchEdge.second->edgeID)
                                {
                                    bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.first.get(), periodicPatchEdge.second.get())));
                                }
                                else
                                {
                                    bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.second.get(), periodicPatchEdge.first.get())));
                                }
                            }
                            else
                            {
                                bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.first.get(), periodicPatchEdge.second.get())));
                            }
                            
                            if (bndIter != bndNodesMap.end())
                            { // exising bnd node found
                                // std::cout<<" Using boundary node "<<bndIter->second->tag()<<std::endl;
                                currentSource = bndIter->second;
                            }
                            else
                            {
                                // Original
                                
                                // const VectorDim loopNodePos(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                                // const VectorDim networkNodePos(loopNodePos + std::get<1>(polyInt[p2]));
                                
                                const auto periodicPatchEdgesAll(std::get<3>(polyInt[p2])); // Get all the edges which the nodes are intesecting
                                // std::cout<<"periodicPatchEdgesAll size is "<<periodicPatchEdgesAll.size()<<std::endl;
                                //get the mesh faces corresponding to the edges from periodicPatchEdgesAll
                                std::set<const PlanarMeshFace<dim> *> allMeshFaces;
                                size_t totalNumEdgesCrossed(0);
                                for (const auto& allPatchMap : periodicPatchEdgesAll)
                                {
                                    // std::cout<<"periodicPatchEdgesAll second size is "<<allPatchMap.second.size()<<std::endl;
                                    totalNumEdgesCrossed+=allPatchMap.second.size();
                                    const auto periodicPatchTemp(loop->periodicGlidePlane->getPatch(allPatchMap.first)); //From shift grab the patch
                                    
                                    for (const auto& allPatchEdges : allPatchMap.second)
                                    {
                                        const auto periodicPatchEdgeTemp(allPatchEdges.second < 0 ? std::make_pair(periodicPatchTemp->edges()[allPatchEdges.first], nullptr)
                                                                         : std::make_pair(periodicPatchTemp->edges()[allPatchEdges.first], periodicPatchTemp->edges()[allPatchEdges.second]));
                                        
                                        for (const auto &pmface : periodicPatchEdgeTemp.first->meshIntersection->faces)
                                        {
                                            allMeshFaces.emplace(pmface);
                                        }
                                        if (periodicPatchEdgeTemp.second)
                                        {
                                            for (const auto &pmface : periodicPatchEdgeTemp.second->meshIntersection->faces)
                                            {
                                                allMeshFaces.emplace(pmface);
                                            }
                                        }
                                    }
                                }
                                
                                
                                
                                
                                const auto networkLoopMapIter(networkNodeLoopMap.find(std::make_pair(periodicNetworkSource, periodicNetworkSink)));
                                // std::set<LoopType *> commonLoops;
                                bool loopBelongtoCommonLoop(false);
                                //                                bool firstLoopInJunction(false); //only for the first loop we need to insert a new network node,otherwise we should be able to find an already existent network node
                                if (networkLoopMapIter != networkNodeLoopMap.end())
                                {
                                    // Insert the loops as determined from the map to preserve the junction information
                                    const auto networkLoopSetIter (networkLoopMapIter->second.find(loop->sID));
                                    loopBelongtoCommonLoop = (networkLoopSetIter != networkLoopMapIter->second.end());
                                    //                                    firstLoopInJunction = (networkLoopSetIter == networkLoopMapIter->second.begin()); //If this is true then we need to create a new network node
                                }
                                else
                                {
                                    // Insert the common loop
                                    loopBelongtoCommonLoop = false;
                                }
                                
                                const VectorDim loopNodePostemp(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                                
                                VectorDim networkNodePos(periodicPatchEdge.first->meshIntersection->snap(loopNodePostemp + std::get<1>(polyInt[p2])));
                                //                                std::set<std::shared_ptr<PeriodicPlanePatch<dim>>> auxiliaryPatches; //Aux patches with only be populated if the second patch edge exists
                                // i.e. The interseection is taking place diagonally
                                if (periodicPatchEdge.second)
                                {
                                    SegmentSegmentDistance<dim> ssd(periodicPatchEdge.first->meshIntersection->P0, periodicPatchEdge.first->meshIntersection->P1,
                                                                    periodicPatchEdge.second->meshIntersection->P0, periodicPatchEdge.second->meshIntersection->P1);
                                    assert(ssd.dMin < FLT_EPSILON && "Two edges must intersect");
                                    networkNodePos = 0.5 * (ssd.x0 + ssd.x1);
                                    
                                    //                                    auxiliaryPatches.insert(loop->periodicGlidePlane->getPatch(periodicPatch->shift+periodicPatchEdge.first->deltaShift));
                                    //                                    auxiliaryPatches.insert(loop->periodicGlidePlane->getPatch(periodicPatch->shift+periodicPatchEdge.second->deltaShift));
                                    
                                    //                                    if ((networkNodePos - (loopNodePostemp + std::get<1>(polyInt[p2]))).norm() > FLT_EPSILON)
                                    //                                    {
                                    //                                        std::cout<<" PeriodicNetworkSource "<<periodicNetworkSource->sID<<std::endl;
                                    //                                        std::cout<<" PeriodicNetworkSink "<<periodicNetworkSink->sID<<std::endl;
                                    //                                        std::cout<<" First edge mesh intersection "<<std::endl;
                                    //                                        for (const auto& mf : periodicPatchEdge.first->meshIntersection->faces)
                                    //                                        {
                                    //                                            std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                    //                                        }
                                    //                                        std::cout<<std::endl;
                                    //                                        std::cout<<" Second edge mesh intersection "<<std::endl;
                                    //                                        for (const auto& mf : periodicPatchEdge.second->meshIntersection->faces)
                                    //                                        {
                                    //                                            std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                    //                                        }
                                    //                                        std::cout<<std::endl;
                                    //                                        std::cout<<" Total edges crossed "<<totalNumEdgesCrossed<<" numEdges yet to be crossed "<<(std::get<4>(polyInt[p2]))<<std::endl;
                                    //                                        std::cout<<" Trying to set network Node position to "<<networkNodePos.transpose()<<std::endl;
                                    //                                        std::cout<<" Actual network node position  "<<(loopNodePostemp + std::get<1>(polyInt[p2])).transpose()<<std::endl;
                                    //                                    }
                                    assert((networkNodePos - (loopNodePostemp + std::get<1>(polyInt[p2]))).norm() < FLT_EPSILON && "Position mismatch");
                                }
                                //                                else
                                //                                {
                                //                                    auxiliaryPatches.insert(std::shared_ptr<PeriodicPlanePatch<dim>>{nullptr});
                                //                                    networkNodePos = periodicPatchEdge.first->meshIntersection->snap(loopNodePostemp + std::get<1>(polyInt[p2]));
                                //                                }
                                const VectorDim loopNodePos(networkNodePos - std::get<1>(polyInt[p2]));
                                
                                const auto currentLoopLink(currentSource->next.second);
                                
                                const auto currentNetworkLink(currentLoopLink->networkLink());
                                
                                std::set<const PlanarMeshFace<dim> *> tmpMeshFaces;
                                for (const auto &pmface : periodicPatchEdge.first->meshIntersection->faces)
                                {
                                    tmpMeshFaces.emplace(pmface);
                                }
                                if (periodicPatchEdge.second)
                                {
                                    for (const auto &pmface : periodicPatchEdge.second->meshIntersection->faces)
                                    {
                                        tmpMeshFaces.emplace(pmface);
                                    }
                                }
                                
                                const size_t edgesStillRemainingtoCross(std::get<4>(polyInt[p2]));
                                VerboseDislocationNetwork(2, " edgesStillRemainingtoCross " << edgesStillRemainingtoCross << std::endl;);
                                VerboseDislocationNetwork(2, " total edges crossed " << totalNumEdgesCrossed << std::endl;);
                                
                                const int u(std::min(edgesStillRemainingtoCross,totalNumEdgesCrossed - edgesStillRemainingtoCross - 1));
                                //                                size_t u(0);
                                //
                                //                                if (loopBelongtoCommonLoop)
                                //                                {
                                //                                    if (periodicNetworkSource == periodicNetworkSink)
                                //                                    {
                                //                                        if (firstLoopInJunction)
                                //                                        {
                                //                                            //We can go as it is...i.e. from periodicPrev to periodicNext
                                //                                            u = edgesStillRemainingtoCross;
                                //                                        }
                                //                                        else
                                //                                        {
                                //                                            //need to check the alignment with the first loop junction nodes
                                //                                            const VectorDim currentLoopLinkChord(periodicNext->get_P()-periodicPrev->get_P());
                                //                                            const double currentLoopLinkChordLength(currentLoopLinkChord.norm());
                                //                                            VectorDim firstLoopLinkChord (VectorDim::Zero());
                                //                                            for (const auto& ln : periodicNetworkSource->loopNodes())
                                //                                            {
                                //                                                if (ln->loop()->sID==*(networkLoopMapIter->second.begin()))
                                //                                                {
                                //                                                    //Can grab any oriented direction
                                //                                                    if (ln->periodicNext()->networkNode==periodicNetworkSource)
                                //                                                    {
                                //                                                        //We are at the source of the loop link
                                //                                                        // std::cout<<" Case A "<<ln->tag()<<ln->get_P().transpose()<<std::endl;
                                //                                                        // std::cout<<" Case A "<<ln->periodicNext()->tag()<<ln->periodicNext()->get_P().transpose()<<std::endl;
                                //                                                        firstLoopLinkChord=ln->periodicNext()->get_P()-ln->get_P();
                                //                                                    }
                                //                                                    else if (ln->periodicPrev()->networkNode==periodicNetworkSource)
                                //                                                    {
                                //                                                        //We are at the sink of the loop link
                                //                                                        // std::cout<<" Case B "<<ln->tag()<<ln->get_P().transpose()<<std::endl;
                                //                                                        // std::cout<<" Case B "<<ln->periodicPrev()->tag()<<ln->periodicPrev()->get_P().transpose()<<std::endl;
                                //                                                        firstLoopLinkChord=ln->get_P()-ln->periodicPrev()->get_P();
                                //                                                    }
                                //                                                    // if (ln->periodicNext()->networkNode==periodicNetworkSource)
                                //                                                    // assert(false && "FINISH HERE");
                                //                                                }
                                //                                            }
                                //                                            const double firstLoopLinkChordLength(firstLoopLinkChord.norm());
                                //                                            assert(firstLoopLinkChordLength>FLT_EPSILON && "First looplink chord must be finite length");
                                //                                            assert(fabs(firstLoopLinkChordLength - currentLoopLinkChordLength)<FLT_EPSILON && "Chord length must be same");
                                //                                            // Determine the alignment
                                //                                            const VectorDim currentLoopLinkDir(currentLoopLinkChord / currentLoopLinkChordLength);
                                //                                            const VectorDim firstLoopLinkDir(firstLoopLinkChord / firstLoopLinkChordLength);
                                //                                            const double dirRef(currentLoopLinkDir.dot(firstLoopLinkDir));
                                //                                            if (fabs(dirRef + 1.0) < FLT_EPSILON)
                                //                                            {
                                //                                                //anti aligned
                                //                                                u = totalNumEdgesCrossed - edgesStillRemainingtoCross - 1; // 1 is to make it to go to 0 (at the end 0 edges should be crossed)
                                //                                            }
                                //                                            else
                                //                                            {
                                //                                                //aligned
                                //                                                u = edgesStillRemainingtoCross;
                                //                                            }
                                //                                        }
                                ////                                        assert(false && "First case occuring... check how the implementation works here");
                                //                                    }
                                //                                    else
                                //                                    {
                                //                                        if (periodicPrevNetwork == periodicNetworkSource)
                                //                                        {
                                //                                            // This condition will give problem if perioicNetworkSource and Sink are same
                                //                                            u = edgesStillRemainingtoCross;
                                //                                        }
                                //                                        else if (periodicNextNetwork == periodicNetworkSource)
                                //                                        {
                                //                                            // This condition will give problem if perioicNetworkSource and Sink are same
                                //                                            u = totalNumEdgesCrossed - edgesStillRemainingtoCross - 1; // 1 is to make it to go to 0 (at the end 0 edges should be crossed)
                                //                                        }
                                //                                    }
                                //                                }
                                
                                // std::cout<<"Total number of edges crossed "<<totalNumEdgesCrossed<<std::endl;
                                const auto key(std::make_tuple(periodicNetworkSource, periodicNetworkSink, allMeshFaces, tmpMeshFaces,u));
                                // std::cout<<" Key is "<<periodicNetworkSource->tag()<<" "<<periodicNetworkSink->tag()<<" "<<u<<" "<<std::endl;
                                // std::cout<<"All Meshfaces "<<std::endl;
                                // for (const auto& mf : allMeshFaces)
                                // {
                                //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                // }
                                // std::cout<<std::endl<< " Temp mesh faces "<<std::endl;
                                // for (const auto& mf : tmpMeshFaces)
                                // {
                                //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                // }
                                // std::cout<<std::endl;
                                if (loopBelongtoCommonLoop)
                                {
                                    //Either we need to insert a new node or we grab an already existent node
                                    const auto networkNodeIter(newNetworkNodesMap.find(key));
                                    if (networkNodeIter == newNetworkNodesMap.end())
                                    {
                                        // Insert a new node
                                        const auto newNetNode(this->networkNodes().create(networkNodePos)); // TODO compute velocity and velocityReduction by interpolation
                                        // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                        VerboseDislocationNetwork(2, " Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                        
                                        newNetworkNodesMap.emplace(key, newNetNode);
                                        const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    }
                                    else
                                    {
                                        // A node has already been inserted... use that node
                                        const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                        const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    }
                                    
                                    
                                    //                                    if (firstLoopInJunction)
                                    //                                    {
                                    //                                        //May be a self annihilation case...For the self annihilation case use the already existent key
                                    //                                        const auto networkNodeIter(newNetworkNodesMap.find(key));
                                    //                                        if (networkNodeIter == newNetworkNodesMap.end())
                                    //                                        {
                                    //                                            // Insert a new node
                                    //                                            const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                    //                                            // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                    //                                            VerboseDislocationNetwork(2, " Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                    //
                                    //                                            newNetworkNodesMap.emplace(key, newNetNode);
                                    //                                            const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                    //                                            currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    //                                        }
                                    //                                        else
                                    //                                        {
                                    //                                            // A node has already been inserted... use that node
                                    //                                            const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                    //                                            const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                    //                                            currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    //                                        }
                                    //                                    }
                                    //                                    else
                                    //                                    {
                                    //                                        //grab an already existent node
                                    //                                        const auto networkNodeIter(newNetworkNodesMap.find(key));
                                    //                                        if (networkNodeIter == newNetworkNodesMap.end())
                                    //                                        {
                                    //                                            std::cout<<"Printing common loops "<<std::endl;
                                    //                                            for (const auto& loop : networkLoopMapIter->second)
                                    //                                            {
                                    //                                                std::cout<<loop<<std::endl;
                                    //                                            }
                                    //                                            std::cout<<"Printing already existent keys "<< newNetworkNodesMap.size()<<" new network node map size"<<std::endl;
                                    //                                            for (const auto& nnMap : newNetworkNodesMap)
                                    //                                            {
                                    //                                                std::cout<<std::get<0>(nnMap.first)->sID<<" "<<std::get<1>(nnMap.first)->sID<<" All Mesh face "<<std::endl;
                                    //                                                std::cout << "All Meshfaces " << std::endl;
                                    //                                                for (const auto &mf : std::get<2>(nnMap.first))
                                    //                                                {
                                    //                                                    std::cout << mf->sID << " with outnormal" << mf->outNormal().transpose() << "\t " << std::flush;
                                    //                                                }
                                    //                                                std::cout<<std::endl<< " Temp mesh faces "<<std::endl;
                                    //                                                for (const auto &mf : std::get<3>(nnMap.first))
                                    //                                                {
                                    //                                                    std::cout << mf->sID << " with outnormal" << mf->outNormal().transpose() << "\t " << std::flush;
                                    //                                                }
                                    //                                                std::cout << std::endl;
                                    //                                                std::cout<<" NetworkNode is "<<nnMap.second->sID<<std::endl;
                                    //                                            }
                                    //
                                    //                                            std::cout<<"Printing bnd Nodes Map "<< bndNodesMap.size()<<" bndNodemap size"<<std::endl;
                                    //                                            // std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *> bndNodesMap;
                                    //                                            for (const auto& bndNode : bndNodesMap )
                                    //                                            {
                                    //                                                std::cout<<"Prev Node "<<(std::get<0>(bndNode.first))->tag()<<std::endl;
                                    //                                                std::cout<<"Next Node "<<(std::get<1>(bndNode.first))->tag()<<std::endl;
                                    //                                                std::cout<<"First Edge ID "<<(std::get<2>(bndNode.first)).first->edgeID<<std::endl;
                                    //                                                if ((std::get<2>(bndNode.first)).second)
                                    //                                                    std::cout<<"Second Edge ID "<<(std::get<2>(bndNode.first)).second->edgeID<<std::endl;
                                    //                                                std::cout<<" Stored Node "<<bndNode.second->tag()<<std::endl;
                                    //
                                    //                                                // std::cout<<(std::get<0>(bndNode.first))->tag()<<" "<<(std::get<1>(bndNode.first))->tag()<<" "<<(std::get<2>(bndNode.first)).first->edgeID<<" "<<(std::get<2>(bndNode.first)).second->edgeID<<" "<<
                                    //                                                // bndNode.second->tag()<<std::endl;
                                    //                                            }
                                    //
                                    //                                        }
                                    //                                        assert(networkNodeIter != newNetworkNodesMap.end() && "Inserting bnd node corresponding to a junction... bnd node should be present already");
                                    //                                        /* The network node position must be commensurate */
                                    //                                        if ((networkNodeIter->second->get_P() - networkNodePos).norm() > 10000 * FLT_EPSILON) // THis condition is just to check for widely different positions
                                    //                                        {
                                    //                                            std::cout << std::scientific << std::setprecision(15) << " PeriodicNetwork Source " << periodicNetworkSource->sID << "P= " << periodicNetworkSource->get_P().transpose() << "\n PeriodicNetwork Sink " << periodicNetworkSink->sID << " P = " << periodicNetworkSink->get_P().transpose() << std::endl;
                                    //                                            std::cout << " Position of the network node " << std::scientific << std::setprecision(15) << networkNodeIter->second->get_P().transpose() << std::endl;
                                    //                                            std::cout << " Actual Position of the network node should be " << std::scientific << std::setprecision(15) << networkNodePos.transpose() << std::endl;
                                    //                                            std::cout << " Position difference " << std::setprecision(15) << (networkNodeIter->second->get_P() - networkNodePos).transpose() << std::endl;
                                    //                                            std::cout << " Position difference norm " << std::setprecision(15) << (networkNodeIter->second->get_P() - networkNodePos).norm() << std::endl;
                                    //                                            assert(false && "Loop node position and network node position mismatch");
                                    //                                        }
                                    //                                        // Here we want the loop node position to be commensurate with the network node (This is important for minimizing accumulating error)
                                    //                                        VerboseDislocationNetwork(2, " Junction case....Grabbing an already existent node " << networkNodeIter->second->tag() << std::endl;);
                                    //                                        const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                    //                                        const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                    //                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    //                                    }
                                }
                                else
                                {
                                    //create a new node
                                    const auto newNetNode(this->networkNodes().create(networkNodePos)); // TODO compute velocity and velocityReduction by interpolation
                                    VerboseDislocationNetwork(2, " non-Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                    
                                    // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                    newNetworkNodesMap.emplace(key, newNetNode);
                                    const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                    currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                }
                            }
                            p2 = (p2 + 1) % polyInt.size();
                        }
                    }
                }
            }
        }
        
        danglingBoundaryLoopNodes.clear();
        
    }

    template <int dim, short unsigned int corder>
    int DislocationNetwork<dim,corder>::verboseDislocationNetwork=0;

    template class DislocationNetwork<3,0>;

}
#endif
