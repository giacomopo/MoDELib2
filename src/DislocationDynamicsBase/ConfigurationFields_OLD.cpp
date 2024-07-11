/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfigurationFields_cpp_
#define model_ConfigurationFields_cpp_

#include <numbers>
#include <ConfigurationFields.h>
#include <StressStraight.h>


namespace model
{

    template <typename T>
    int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    template <int dim>
    ConfigurationFields<dim>::ConfigurationFields(DislocationDynamicsBase<dim>& ddBase_in,const DDconfigIO<dim>& configIO_in):
//    /* init */ ClusterDynamicsBase<dim>(ddBase_in)
//    /* init */ ElasticDeformationBase<dim>(ddBase_in)
    /* init */,ddBase(ddBase_in)
    /* init */,configIO(configIO_in)
    /* init */,MicrostructureContainerType(ddBase_in)
    /* init */,f_file(this->ddBase.simulationParameters.traitsIO.fFile,std::ios_base::app)
    /* init */,F_labels(this->ddBase.simulationParameters.traitsIO.flabFile,std::ios_base::app)
    {
        if (!f_file.is_open())
        {
            throw std::runtime_error("Cannot open file "+this->ddBase.simulationParameters.traitsIO.fFile);
        }
        if (!F_labels.is_open())
        {
            throw std::runtime_error("Cannot open file "+this->ddBase.simulationParameters.traitsIO.flabFile);
        }
        
        if(ddBase.simulationParameters.useClusterDynamics)
        {
            this->emplace_back(new ClusterDynamicsType(*this));
        }
        if(ddBase.simulationParameters.useDislocations)
        {
            this->emplace_back(new DislocationNetworkType(*this));
        }
        if(ddBase.simulationParameters.useElasticDeformation)
        {
            this->emplace_back(new ElasticDeformationType(*this));
        }
        
//        this->initializeConfiguration(configIO,f_file,F_labels);

    }


    template <int dim>
    void ConfigurationFields<dim>::updateConfiguration()
    {
        // clean up
        loopPatches().clear();
        segments()=configIO.segments();
        polyhedronInclusionNodes().clear();
        eshelbyInclusions().clear();
        EshelbyInclusionBase<dim>::force_count(0);
        
        // Update loop patches
        for(const auto& pair : configIO.loopNodeSequence())
        {
            const auto& loopID(pair.first);
            const auto& loopIO(configIO.loops()[configIO.loopMap().at(loopID)]);
            const auto& grain(ddBase.poly.grain(loopIO.grainID));
            GlidePlaneKey<3> glidePlaneKey(loopIO.P, grain.singleCrystal->reciprocalLatticeDirection(loopIO.N));
            std::shared_ptr<PeriodicGlidePlane<3>> periodicGlidePlane(ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            
            std::vector<Eigen::Matrix<double,3,1>> nodeShifts;
            std::vector<Eigen::Matrix<double,3,1>> nodePos;
            
            for(const auto& loopNodeID : pair.second)
            {
                const auto& loopNodeIO(configIO.loopNodes()[configIO.loopNodeMap().at(loopNodeID)]);
                nodeShifts.push_back(loopNodeIO.periodicShift);
                if(loopNodeIO.edgeIDs.first<0 && loopNodeIO.edgeIDs.second<0)
                {
                    nodePos.push_back(loopNodeIO.P);
                }
            }
            
            DislocationLoopPatches<3> currentPatches(periodicGlidePlane);
            currentPatches.update(nodeShifts,nodePos);
            loopPatches().emplace(loopID,currentPatches);
        }
        
        // update inclusions
        for(const auto& inclusion : configIO.sphericalInclusions())
        {
            //        std::cout<<"Creating spherical inclusion "<<inclusion.inclusionID<<std::endl;
            const std::pair<bool,const Simplex<dim,dim>*> searchPair(ddBase.mesh.search(inclusion.C));
            if(searchPair.first)
            {
                
                const auto& grain(ddBase.poly.grain(searchPair.second->region->regionID));
                if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                {
                    const auto secondPhase(grain.singleCrystal->secondPhases().at(inclusion.phaseID));
                    EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                    
                    
                    std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new SphericalInclusion<dim>(inclusion.C,inclusion.a,inclusion.eT,ddBase.poly.nu,ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                    
                    eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                }
                else
                {
                    throw std::runtime_error("phaseID does not exist in grain.");
                }
            }
        }
        
        for(const auto& piNode : configIO.polyhedronInclusionNodes())
        {
            polyhedronInclusionNodes().emplace(piNode.nodeID,piNode);
        }
        
        std::map<size_t,std::map<size_t,std::vector<size_t>>> faceMap;
        for(const auto& edge : configIO.polyhedronInclusionEdges())
        {
            const size_t& iID(edge.inclusionID);
            const size_t& fID(edge.faceID);
            const size_t& sourceID(edge.sourceID);
            faceMap[iID][fID].push_back(sourceID);
        }
        
        
        for(const auto& inclusion : configIO.polyhedronInclusions())
        {
            //        std::cout<<"Creating polyhedron inclusion "<<inclusion.inclusionID<<std::endl;
            
            const auto faceIter(faceMap.find(inclusion.inclusionID));
            if(faceIter!=faceMap.end())
            {
                const auto& faces(faceIter->second);
                //                std::cout<<"    #faces= "<<faces.size()<<std::endl;
                std::set<const PolyhedronInclusionNodeIO<dim>*> uniquePolyNodes;
                for(const auto& pair : faces)
                {
                    for(const auto& nodeID : pair.second)
                    {
                        uniquePolyNodes.emplace(&polyhedronInclusionNodes().at(nodeID));
                    }
                }
                //                std::cout<<"    #nodes= "<<uniquePolyNodes.size()<<std::endl;
                if(uniquePolyNodes.size()>=dim+1)
                {
                    // Find grain
                    std::set<size_t> grainIDs;
                    for(const auto& nodePtr : uniquePolyNodes)
                    {
                        const std::pair<bool,const Simplex<dim,dim>*> searchPair(ddBase.mesh.search(nodePtr->P));
                        if(searchPair.first)
                        {
                            grainIDs.insert(searchPair.second->region->regionID);
                        }
                        else
                        {
                            throw std::runtime_error("inclusion node outside mesh");
                        }
                    }
                    
                    // Add inclusion
                    if(grainIDs.size()==1)
                    {
                        const auto& grain(ddBase.poly.grain(*grainIDs.begin()));
                        if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                        {
                            const auto secondPhase(grain.singleCrystal->secondPhases().at(inclusion.phaseID));
                            EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                            std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new PolyhedronInclusion<dim>( polyhedronInclusionNodes(),faces,inclusion.eT,ddBase.poly.nu,ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                            eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                        }
                        else
                        {
                            throw std::runtime_error("phaseID does not exist in grain.");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("inclusion across grain boundaries");
                    }
                }
                else
                {
                    throw std::runtime_error("inclusion does not have enough nodes");
                }
            }
            else
            {
                throw std::runtime_error("inclusionID not found in faceMap");
            }
        }
        
        
//        const int edNodes(configIO.displacementMatrix().rows());
//        if(edNodes*dim == elasticDeformation().u.dofVector().size())
//        {
//            elasticDeformation().u=configIO.displacementMatrix().block(0,0,edNodes,dim).transpose().reshaped(elasticDeformation().u.dofVector().size(),1);
//        }
//        else
//        {   
//            std::cout<<"elasticDeformation DOF mismatch."<<std::endl;
////            throw std::runtime_error("elasticDeformation DOF mismatch.");
//        }


//        const int cdNodes(configIO.cdMatrix().rows());
//        if(cdNodes*this->mSize == clusterDynamics().mobileClusters.dofVector().size())
//        {
//            clusterDynamics().mobileClusters=configIO.cdMatrix().block(0,0,cdNodes,this->mSize).transpose().reshaped(clusterDynamics().mobileClusters.dofVector().size(),1);
//        }
//        else
//        {
//            std::cout<<"mobileSpecies DOF mismatch."<<std::endl;
////            throw std::runtime_error("mobileSpecies DOF mismatch.");
//        }
//        if(cdNodes*this->iSize == clusterDynamics().immobileClusters.dofVector().size())
//        {
//            clusterDynamics().immobileClusters=configIO.cdMatrix().block(0,this->mSize,cdNodes,this->iSize).transpose().reshaped(clusterDynamics().immobileClusters.dofVector().size(),1);
//        }
//        else
//        {
//            std::cout<<"immobileSpecies DOF mismatch."<<std::endl;
////            throw std::runtime_error("immobileSpecies DOF mismatch.");
//        }
    }

    template <int dim>
    const std::map<size_t,DislocationLoopPatches<dim>>& ConfigurationFields<dim>::loopPatches() const
    {
        return *this;
    }

    template <int dim>
    std::map<size_t,DislocationLoopPatches<dim>>& ConfigurationFields<dim>::loopPatches()
    {
        return *this;
    }

    template <int dim>
    const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& ConfigurationFields<dim>::segments() const
    {
        return *this;
    }

    template <int dim>
    std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& ConfigurationFields<dim>::segments()
    {
        return *this;
    }

    template <int dim>
    const typename ConfigurationFields<dim>::PolyhedronInclusionNodeContainerType& ConfigurationFields<dim>::polyhedronInclusionNodes() const
    {
        return *this;
    }

    template <int dim>
    typename ConfigurationFields<dim>::PolyhedronInclusionNodeContainerType& ConfigurationFields<dim>::polyhedronInclusionNodes()
    {
        return *this;
    }

    template <int dim>
    const typename ConfigurationFields<dim>::EshelbyInclusionContainerType& ConfigurationFields<dim>::eshelbyInclusions() const
    {
        return *this;
    }

    template <int dim>
    typename ConfigurationFields<dim>::EshelbyInclusionContainerType& ConfigurationFields<dim>::eshelbyInclusions()
    {
        return *this;
    }

//template <int dim>
//const ClusterDynamicsBase<dim>& ConfigurationFields<dim>::clusterDynamics() const
//{
//    return *this;
//}
//
//template <int dim>
//ClusterDynamicsBase<dim>& ConfigurationFields<dim>::clusterDynamics()
//{
//    return *this;
//}

//template <int dim>
//const ElasticDeformationBase<dim>& ConfigurationFields<dim>::elasticDeformation() const
//{
//    return *this;
//}
//
//template <int dim>
//ElasticDeformationBase<dim>& ConfigurationFields<dim>::elasticDeformation()
//{
//    return *this;
//}

    template <int dim>
    double ConfigurationFields<dim>::solidAngle(const VectorDim& x) const
    {
        double temp(0.0);
        for(const auto& patch : loopPatches())
        {
            for(const auto& shift : ddBase.periodicShifts)
            {
                temp+=patch.second.solidAngle(x+shift);
            }
        }
        return temp;
    }

    template <int dim>
    typename ConfigurationFields<dim>::VectorDim ConfigurationFields<dim>::dislocationPlasticDisplacement(const VectorDim& x) const
    {
        VectorDim temp(VectorDim::Zero());
        for(const auto& patch : loopPatches())
        {
            const auto& loop(configIO.loop(patch.first));
            for(const auto& shift : ddBase.periodicShifts)
            {
                temp-=patch.second.solidAngle(x+shift)/4.0/std::numbers::pi*loop.B;
            }
        }
        return temp;
    }

    template <int dim>
    typename  ConfigurationFields<dim>::ConcentrationVectorType ConfigurationFields<dim>::dislocationMobileConentrations(const VectorDim& x,
                                                                                                                         const Simplex<dim,dim>* guess) const
    {
        ConcentrationVectorType temp(ConcentrationVectorType::Zero());
        const std::pair<bool,const Simplex<dim,dim>*> found(ddBase.mesh.searchWithGuess(x,guess));
        if(found.first)
        {
            guess = found.second; // return (reassign) guess for possible next search
            const int grainID(found.second->region->regionID);
            
            for (const auto& segment : segments())
            {
                if(segment.second.grainIDs.find(grainID)!=segment.second.grainIDs.end())
                {
                    auto itSource(configIO.nodeMap().find(segment.second.sourceID)); //source
                    auto   itSink(configIO.nodeMap().find(segment.second.sinkID)); //sink
                    if(itSource!=configIO.nodeMap().end() && itSink!=configIO.nodeMap().end())
                    {
                        const auto& sourceNode(configIO.nodes()[itSource->second]);
                        const auto&   sinkNode(configIO.nodes()[itSink->second]);
                        StressStraight<3> ss(ddBase.poly,sourceNode.P,sinkNode.P,segment.second.b,ddBase.EwaldLength);
                        
                        const double sourceVnorm(sourceNode.V.norm());
                        const double sourceVsum(sourceNode.climbVelocityScalar.sum());
                        if(std::abs(sourceVnorm-std::abs(sourceVsum))>FLT_EPSILON)
                        {
                            throw std::runtime_error("sourceVnorm != abs(sourceVsum)");
                        }
                        const VectorDim sourceVDir(sourceVnorm==0.0 ? VectorDim::Zero() : (sourceNode.V/sourceVsum).eval());
                        const double sinkVnorm(sinkNode.V.norm());
                        const double sinkVsum(sinkNode.climbVelocityScalar.sum());
                        if(std::abs(sinkVnorm-std::abs(sinkVsum))>FLT_EPSILON)
                        {
                            throw std::runtime_error("sinkVnorm != abs(sinkVsum)");
                        }
                        const VectorDim sinkVDir(sinkVnorm==0.0 ? VectorDim::Zero() : (sinkNode.V/sinkVsum).eval());

                        for(const auto& shift : ddBase.periodicShifts)
                        {
                            temp+=ss.clusterConcentration(x+shift,grainID, sourceVDir, sourceNode.climbVelocityScalar , sinkVDir, sinkNode.climbVelocityScalar, this->cdp);
                        }
                    }
                }
            }
        }
        return temp;
    }

template <int dim>
typename  ConfigurationFields<dim>::MatrixDim ConfigurationFields<dim>::dislocationStress(const VectorDim& x) const
{
    MatrixDim temp(MatrixDim::Zero());
    for (const auto& segment : segments())
    {
        if(segment.second.meshLocation==0 || segment.second.meshLocation==2)
        {// segment inside mesh or on grain boundary
            auto itSource(configIO.nodeMap().find(segment.second.sourceID)); //source
            auto   itSink(configIO.nodeMap().find(segment.second.sinkID)); //sink
            if(itSource!=configIO.nodeMap().end() && itSink!=configIO.nodeMap().end())
            {
                const auto& sourceNode(configIO.nodes()[itSource->second]);
                const auto&   sinkNode(configIO.nodes()[itSink->second]);
                StressStraight<3> ss(ddBase.poly,sourceNode.P,sinkNode.P,segment.second.b,ddBase.EwaldLength);
                for(const auto& shift : ddBase.periodicShifts)
                {
                    temp+=ss.stress(x+shift);
                }
            }
            else
            {
                
            }
        }
    }
    return temp;
}

    template <int dim>
    typename  ConfigurationFields<dim>::MatrixDim ConfigurationFields<dim>::inclusionStress(const VectorDim& x) const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& inclusion : eshelbyInclusions() )
        {
            for(const auto& shift : ddBase.periodicShifts)
            {
                temp+=inclusion.second->stress(x+shift);
            }
        }
        return temp;
    }

    template class ConfigurationFields<3>;

}
#endif
