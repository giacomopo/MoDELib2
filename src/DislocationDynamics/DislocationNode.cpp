/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_cpp_
#define model_DislocationNode_cpp_

#include <DislocationNode.h>
#include <PlanesIntersection.h>

namespace model
{
    
    template <int dim>
    DislocationNode<dim>::DislocationNode(LoopNetworkType* const net,
                                                                   const VectorDim& P,
                                                                   const VectorDim& V,
                                                                   const ClimbVelocityScalarType& cvs,
                                                                   const double& vrc) :
    /* init */ NetworkNode<DislocationNode>(net)
    /* init */,SplineNodeType(P)
    /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
    /* init */,climbVelocityScalar(cvs)
    /* init */,velocity(V)
    /* init */,vOld(velocity)
    /* init */,velocityReductionCoeff(vrc)
    {
        VerboseDislocationNode(1, "  Creating Network Node " << this->tag() <<" @ "<<this->get_P().transpose() << std::endl;);        
    }

template <int dim>
DislocationNode<dim>::DislocationNode(LoopNetworkType* const net,
                                                               const VectorDim& P) :
/* init */ NetworkNode<DislocationNode>(net)
/* init */,SplineNodeType(P)
/* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
/* init */,climbVelocityScalar(ClimbVelocityScalarType::Zero())
/* init */,velocity(VectorDim::Zero())
/* init */,vOld(velocity)
/* init */,velocityReductionCoeff(1.0)
{
    VerboseDislocationNode(1, "  Creating Network Node " << this->tag() <<" @ "<<this->get_P().transpose() << std::endl;);
}

template <int dim>
typename DislocationNode<dim>::VectorDim DislocationNode<dim>::climbDirection() const
{
    VectorDim temp(VectorDim::Zero());
    for(const auto& link : this->neighbors())
    {// weighted average of segments climb directions
        if(std::get<1>(link.second)->isSessile())
        {
            temp+=std::get<1>(link.second)->climbDirection()*std::get<1>(link.second)->chordLength();
        }
    }
    const double tempNorm(temp.norm());
    return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
}

    template <int dim>
    DislocationNode<dim>::~DislocationNode()
    {
        VerboseDislocationNode(1, "  Destroying Network Node " << this->tag() << std::endl;);

    }
    
    template <int dim>
    std::shared_ptr<DislocationNode<dim>> DislocationNode<dim>::clone() const
    {
        return std::shared_ptr<DislocationNode<dim>>(new DislocationNode(this->p_network(),this->get_P(),get_V(),climbVelocityScalar,velocityReduction()));
    }
    
    template <int dim>
    const Simplex<dim,dim>* DislocationNode<dim>::get_includingSimplex(const VectorDim& X,const Simplex<dim,dim>* const guess) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
        if (guess==NULL)
        {
            temp=this->network().ddBase.mesh.search(X);
        }
        else
        {
            const auto grains(this->grains());
            if(grains.size()==1)
            {// node only in one region
                if((*grains.begin())->grainID!=guess->region->regionID)
                {
                    temp=this->network().ddBase.mesh.searchRegion((*grains.begin())->grainID,X);
                }
                else
                {
                    temp=this->network().ddBase.mesh.searchRegionWithGuess(X,guess);
                }
            }
            else
            {
                temp=this->network().ddBase.mesh.searchWithGuess(X,guess);
            }
        }
        if(!temp.first) // PlanarDislocationNode not found inside mesh
        {// Detect if the PlanarDislocationNode is sligtly outside the boundary
            int faceID;
            const double baryMin(temp.second->pos2bary(X).minCoeff(&faceID));
            const bool isApproxOnBoundary(std::fabs(baryMin)<1.0e3*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
            if(!isApproxOnBoundary)
            {
                std::cout<<"PlanarDislocationNode "<<this->sID<<" @ "<<X.transpose()<<std::endl;
                std::cout<<"Simplex "<<temp.second->xID<<std::endl;
                std::cout<<"bary "<<temp.second->pos2bary(X)<<std::endl;
                std::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                std::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                std::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                throw std::runtime_error("DISLOCATION NODE OUTSIDE MESH.");
//                assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
            }
        }
        return temp.second;
    }

    
    template <int dim>
    void DislocationNode<dim>::projectVelocity(const bool& isClimbingStep)
    {
        
        VectorOfNormalsType temp;
//        temp.push_back(VectorDim::UnitZ());

        
        for(const auto& loopNode : this->loopNodes())
        {
            const auto& loop(loopNode->loop());
            if(loop)
            {
                if(loop->glidePlane)
                {
                    temp.push_back(loop->glidePlane->unitNormal);
                }
                if(!isClimbingStep && loop->loopType==DislocationLoopIO<dim>::SESSILELOOP)
                {
                    velocity.setZero();
                    break;
                }
//                if(loop->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
//                {
//                    if(loop->glidePlane)
//                    {
//                        temp.push_back(loop->glidePlane->unitNormal);
//                    }
//                }
//                else if(loop->loopType==DislocationLoopIO<dim>::SESSILELOOP)
//                {
//                    velocity.setZero();
//                    break;
//                }
            }
//            else
//            {
//                assert(false && "no loop in node");
//            }

        }
        
//        if(velocity.squaredNorm()>FLT_EPSILON)
//        {
            
            
            if(!this->network().ddBase.isPeriodicDomain)
            {// Use boundary planes to confine velocity in case of non-periodic simulation
                
                for(const auto& face : this->meshFaces())
                {
                    temp.push_back(face->asPlane().unitNormal);
                }
            }
            
            for(int r=0;r<this->network().nodalVelocityConstraints.rows();++r)
            {
                temp.push_back(this->network().nodalVelocityConstraints.row(r));
            }
            

            GramSchmidt::orthoNormalize(temp);
            
            for(const auto& vec : temp)
            {
                velocity-=velocity.dot(vec)*vec;
            }
            
//        }
    }
    
    template <int dim>
    void DislocationNode<dim>::set_V(const VectorDim& vNew,const bool& isClimbingStep)
    {
        
        double emaStrength=0.0;
                
        vOld=velocity; // store current value of velocity before updating
        
        //velocity=vNew;
        velocity=(1.0-emaStrength)*vNew+emaStrength*velocity;
        projectVelocity(isClimbingStep);

//        if (this->network().capMaxVelocity && !this->network().timeIntegrator.isTimeStepControlling(*this))
//        {
//            //Need to cap the max velocity based on the time step
//            const double vmax=this->network().timeIntegrator.dxMax/this->network().timeIntegrator.dtMax;
//            assert(fabs(vmax)>FLT_EPSILON && "Max velocity should be greater than 0");
//            if (velocity.norm()>=FLT_EPSILON && velocity.norm()>fabs(vmax))
//            {
//                const double velredFac(velocity.norm()/fabs(vmax));
//                assert(velredFac>=1.0 && "Velocity reduction factor must be greater than 1");
//                velocity/=velredFac;
//                totalCappedNodes+=1;
//            }
//        }
        
//        if(this->network().use_velocityFilter)
//        {
//            if(!isBoundaryNode())
//            {
//                const double filterThreshold=0.05*velocity.norm()*vOld.norm()+FLT_EPSILON;
//                
//                if(velocity.dot(vOld)<-filterThreshold)
//                {
//                    velocityReductionCoeff*=this->network().velocityReductionFactor;
//                }
//                else if(velocity.dot(vOld)>filterThreshold)
//                {
//                    velocityReductionCoeff/=this->network().velocityReductionFactor;
//                }
//                else
//                {
//                    // don't change velocityReductionCoeff
//                }
//                if(velocityReductionCoeff>1.0)
//                {
//                    velocityReductionCoeff=1.0;
//                }
//                if(velocityReductionCoeff<0.005)
//                {
//                    velocityReductionCoeff=0.005;
//                }
//                velocity*=velocityReductionCoeff;
//            }
//            else
//            {
//                // TO DO BOUNDARY NODES NEED TO INTERPOLATE ALSO THE VELOCTY FILTER
//            }
//            
//        }
    }

    template <int dim>
    void DislocationNode<dim>::updateGeometry()
    {
        for(const auto& neighboor : this->neighbors())
        {
            std::get<1>(neighboor.second)->updateGeometry();
        }
    }

    template <int dim>
    bool DislocationNode<dim>::trySet_P(const typename DislocationNode<dim>::VectorDim& newP)
    {
        VerboseDislocationNode(1, " Try  Setting P for Network Node " << this->tag()<<" to"<<newP.transpose()<< std::endl;);
        const VectorDim snappedPosition(snapToGlidePlanesinPeriodic(newP));
        VerboseDislocationNode(2, " snappedPosition= " << snappedPosition.transpose()<< std::endl;);
        std::pair<bool, const Simplex<dim,dim>*> temp(this->network().ddBase.mesh.searchRegionWithGuess(newP,p_Simplex));
        if((this->get_P()-snappedPosition).norm()>FLT_EPSILON && (this->network().ddBase.isPeriodicDomain || temp.first))
        {
            for (auto &loopNode : this->loopNodes())
            {
                loopNode->set_P(loopNode->periodicPlanePatch() ? snappedPosition - loopNode->periodicPlanePatch()->shift : snappedPosition);
            }
        }

        return true;
    }
    
    template <int dim>
    bool DislocationNode<dim>::set_P(const typename DislocationNode<dim>::VectorDim& newP)
    {
        VerboseDislocationNode(1, "  Setting P for Network Node " << this->tag() << std::endl;);

        if((SplineNodeType::get_P()-newP).norm()>FLT_EPSILON)   //Check with the base classs
        {
            SplineNodeType::set_P(newP);
            // this->updateConfinement(this->get_P());
            p_Simplex=get_includingSimplex(this->get_P(),p_Simplex); // update including simplex
            
        }
        
        return (this->get_P()-newP).norm()<FLT_EPSILON;
    }


    
    template <int dim>
    const typename DislocationNode<dim>::VectorDim& DislocationNode<dim>::get_V() const
    {/*! The nodal velocity vector
      */
        return velocity;
    }

    template <int dim>
    std::set<typename DislocationNode<dim>::LoopType*> DislocationNode<dim>::sessileLoops() const
    {
        std::set<LoopType*> temp;
        for (const auto& loop : this->loops())
        {
            if (loop->loopType == DislocationLoopIO<dim>::SESSILELOOP)
            {
                temp.insert(loop);
            }
        }
        return temp;
    }


    template <int dim>
    const double& DislocationNode<dim>::velocityReduction() const
    {
        return velocityReductionCoeff;
    }
    
    template <int dim>
    typename DislocationNode<dim>::MeshLocation DislocationNode<dim>::meshLocation() const
    {/*!\returns the position of *this relative to the bonudary:
      * 1 = inside mesh
      * 2 = on mesh boundary
      */
        MeshLocation temp = MeshLocation::outsideMesh;

//        if(!isVirtualBoundaryNode())
//        {
            if(isBoundaryNode())
            {
                temp=MeshLocation::onMeshBoundary;
            }
            else
            {
                if(isGrainBoundaryNode())
                {
                    temp=MeshLocation::onRegionBoundary;
                }
                else
                {
                    temp=MeshLocation::insideMesh;
                }
            }
//        }
        return temp;
    }
    
//    template <int dim>
//    const std::shared_ptr<typename DislocationNode<dim>::NetworkNodeType>& DislocationNode<dim>::virtualBoundaryNode() const
//    {
//        return virtualNode;
//    }
    
    template <int dim>
    typename DislocationNode<dim>::VectorDim DislocationNode<dim>::invariantDirectionOfMotion() const
    {/*!\returns the direction of alignment if all links connected to this node are geometrically aligned.
      * Otherwise it returns the zero vector.
      */
        VectorDim temp(this->neighbors().size()? std::get<1>(this->neighbors().begin()->second)->chord().normalized() : VectorDim::Zero());
        for (const auto& neighborIter : this->neighbors())
        {
            if(std::get<1>(neighborIter.second)->chord().normalized().cross(temp).norm()>FLT_EPSILON)
            {
                temp.setZero();
            }
        }
        return temp;
    }
    
    template <int dim>
    bool DislocationNode<dim>::isMovableTo(const VectorDim& X) const
    {
        for(const auto& loopNode : this->loopNodes())
        {
            const VectorDim shift(loopNode->periodicPlanePatch()? loopNode->periodicPlanePatch()->shift : VectorDim::Zero());
            
            if(!loopNode->isMovableTo(X-shift))
            {
                VerboseDislocationNode(2,"  DislocationNode "<<this->tag()<< " NOT movable "<<std::endl;);
                return false;
            }
        }
        return true;
        
    }
    
    template <int dim>
    const Simplex<dim,dim>* DislocationNode<dim>::includingSimplex() const
    {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
      */
        return p_Simplex;
    }
    
//    template <int dim>
//    bool DislocationNode<dim>::isVirtualBoundaryNode() const
//    {
//        return masterNode && this->meshFaces().size()==0;
////        return false;
//    }
    


    template <int dim>
    bool DislocationNode<dim>::isBoundaryNode() const
    {
        return isOnExternalBoundary();
    }
    
    template <int dim>
    bool DislocationNode<dim>::isGrainBoundaryNode() const
    {
        return this->isOnInternalBoundary();
    }

    template <int dim>
    typename DislocationNode<dim>::GlidePlaneContainerType DislocationNode<dim>::glidePlanes() const
    {
        GlidePlaneContainerType temp;
        for (const auto &ln : this->loopNodes())
        {
            if (ln->periodicPlanePatch())
            {
                temp.emplace(ln->periodicPlanePatch()->glidePlane.get());
            }
        }
        return temp;
    }

    template <int dim>
    typename DislocationNode<dim>::PlanarMeshFaceContainerType DislocationNode<dim>::meshFaces() const
    {
        PlanarMeshFaceContainerType temp;
        for (const auto &ln : this->loopNodes())
        {
            if (ln->periodicPlaneEdge.first)
            {
                temp = ln->periodicPlaneEdge.first->meshIntersection->faces;
            }
            if (ln->periodicPlaneEdge.second)
            {
                for (const auto &face : ln->periodicPlaneEdge.second->meshIntersection->faces)
                {
                    temp.emplace(face);
                }
            }
        }
        return temp;
    }

    template <int dim>
    bool DislocationNode<dim>::isOnExternalBoundary() const
    { /*!\returns _isOnExternalBoundarySegment.
       */
        bool _isOnExternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first == face->regionIDs.second)
            {
                _isOnExternalBoundary = true;
                break;
            }
        }

        return _isOnExternalBoundary;
    }

    template <int dim>
    bool DislocationNode<dim>::isOnInternalBoundary() const
    {
        bool _isOnInternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first != face->regionIDs.second)
            {
                _isOnInternalBoundary = true;
                break;
            }
        }
        return _isOnInternalBoundary;
    }

    template <int dim>
    bool DislocationNode<dim>::isOnBoundary() const
    {
        return isOnExternalBoundary() || isOnInternalBoundary();
    }

    template <int dim>
    typename DislocationNode<dim>::VectorDim DislocationNode<dim>::bndNormal() const
    {
        VectorDim _outNormal(VectorDim::Zero());
        for (const auto &face : meshFaces())
        {
            _outNormal += face->outNormal();
        }
        const double _outNormalNorm(_outNormal.norm());
        if (_outNormalNorm > FLT_EPSILON)
        {
            _outNormal /= _outNormalNorm;
        }
        else
        {
            _outNormal.setZero();
        }
        return _outNormal;
    }

    template <int dim>
    typename DislocationNode<dim>::VectorDim DislocationNode<dim>::snapToGlidePlanesinPeriodic(const VectorDim &x) const
    {
        GlidePlaneContainerType gps(glidePlanes());
        if(gps.size())
        {
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));

            int k=0;
            for(const auto& plane : gps)
            {
                N.col(k)=plane->unitNormal;
                P.col(k)=plane->P;
                ++k;
            }
            PlanesIntersection<dim> pInt(N,P,FLT_EPSILON);

            const std::pair<bool,VectorDim> snapped(pInt.snap(x));
            if(snapped.first)
            {
                return snapped.second;
            }
            else
            {
                std::cout<<"glidePlanes.size()="<<gps.size()<<std::endl;
                std::cout<<"N="<<N<<std::endl;
                std::cout<<"P="<<P<<std::endl;
                throw std::runtime_error("DislocationNode, cannot snap, glidePlanes dont intersect.");
                return snapped.second;
            }
        }
        else
        {
            if(sessileLoops().size()==this->loops().size())
            {
                return x;
            }
            else
            {
                throw std::runtime_error("All loops must be sessile if there are no glide planes.");
                return x;
            }
        }
    }

    template <int dim>
    typename DislocationNode<dim>::GrainContainerType DislocationNode<dim>::grains() const
    {
        GrainContainerType temp;
        for (const auto &glidePlane : glidePlanes())
        {
            temp.insert(&glidePlane->grain);
        }
        return temp;
    }
    
    template <int dim>
    int DislocationNode<dim>::totalCappedNodes=0;
    
    template class DislocationNode<3>;
    
}
#endif
