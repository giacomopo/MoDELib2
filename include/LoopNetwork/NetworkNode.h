/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkNode_H_
#define model_NetworkNode_H_

#include <iostream>
#include <set>
#include <memory>
#include <assert.h>
#include <tuple>
#include <limits.h>

#include <StaticID.h>
#include <CRTP.h>
#include <NetworkBase.h>


#ifndef NDEBUG
#define VerboseNetworkNode(N,x) if(verboseLevel>=N){std::cout<<redColor<<x<<defaultColor;}
#else
#define VerboseNetworkNode(N,x)
#endif

namespace model
{
    
    
    template<typename Derived>
    class NetworkNode : public StaticID<Derived>
    /*            */,public CRTP<Derived>
    /*            */,public NetworkBase<Derived,size_t>
    /*            */,public std::set<typename TypeTraits<Derived>::LoopNodeType*>
    /*            */,private std::map<size_t,std::tuple<Derived* const ,typename TypeTraits<Derived>::NetworkLinkType* const>>
    {
        
    public:
        
        typedef Derived NetworkNodeType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef typename std::set<LoopLinkType*> LoopLinkContainerType;
                typedef typename TypeTraits<Derived>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef std::set<LoopNodeType*> LoopNodeContainerType;
        typedef NetworkBase<Derived,size_t> NetworkBaseType;
        typedef std::tuple<Derived* const ,NetworkLinkType* const>                NeighborType;
        typedef std::map<size_t,NeighborType>                            NeighborContainerType;
        
    public:

        static int verboseLevel;

        NetworkNode(LoopNetworkType* const loopNetwork_in) :
        /* init */ NetworkBaseType(loopNetwork_in,&loopNetwork_in->networkNodes(),this->sID)
        {
            VerboseNetworkNode(1,"Constructing NetworkNode "<<tag()<<std::endl);
        }
        
        /**********************************************************************/
        ~NetworkNode()
        {
            VerboseNetworkNode(1,"Destroying NetworkNode "<<tag()<<std::endl);
            assert(neighbors().empty());
        }
        
         LoopLinkContainerType outLoopLinks() const
        {
            LoopLinkContainerType temp;
            for (const auto &ln : loopNodes())
            {
                if (ln->next.second)
                {
                    temp.insert(ln->next.second);
                }
            }
            return temp;
        }

        LoopLinkContainerType inLoopLinks() const
        {
            LoopLinkContainerType temp;
            for (const auto &ln : loopNodes())
            {
                if (ln->prev.second)
                {
                    temp.insert(ln->prev.second);
                }
            }
            return temp;
        }

        std::set<size_t> loopIDs() const
        {
            std::set<size_t> temp;
            for (const auto &loopIter : loops())
            {
                temp.insert(loopIter->sID);
            }
            return temp;
        }

        std::set<LoopType*> loops() const
        {
            std::set<LoopType*> temp;
            for (const auto& ln : loopNodes())
            {
                temp.insert(ln->loop().get());
            }
            return temp;
        }

        const LoopNodeContainerType& loopNodes() const
        {
            return *this;
        }
        
        LoopNodeContainerType& loopNodes()
        {
            return *this;
        }
        
        void addLoopNode(LoopNodeType* const pN)
        {
            const bool success(loopNodes().insert(pN).second);
            if(!success)
            {
                throw std::runtime_error("Duplicate LoopNode could not be added");
            }
//            assert(success && "Duplicate LoopNode could not be added");
        }
        
        void removeLoopNode(LoopNodeType* const pN)
        {
            const size_t erased(loopNodes().erase(pN));
            if(erased!=1)
            {
                throw std::runtime_error("LoopNode could not be erased");
            }
        }
        
        bool isContractableTo(const std::shared_ptr<Derived>&) const
        {
            return true;
        }
        
        size_t gID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            const auto nodeIter(this->network().networkNodes().find(this->sID));
            if(nodeIter!=this->network().networkNodes().end())
            {
                return std::distance(this->network().networkNodes().begin(),nodeIter);
            }
            else
            {
                throw std::runtime_error("NetworkNode "+std::to_string(this->sID)+" not found in calling NetworkNode::gID().");
                return 0;
            }
//            return ;
//            size_t globalNodeID(const size_t& sID) const
//            {
//                
//            }
//            
//            return this->network().globalNodeID(this->sID);
        }
        
        const NeighborContainerType& neighbors() const
        {
            return *this;
        }

        NeighborContainerType& neighbors()
        {
            return *this;
        }
        
        void addToNeighborhood(NetworkLinkType* const pL)
        {/*!@param[in] pL a pointer to a LinkType edge
          */
            if (pL->source->sID==this->sID)
            {// this vertex is the source of edge *pL, so sink of *pL is the neighbor
                const NeighborType temp(pL->sink.get(),pL);
                // std::cout<<" Addng to neighborhood "<<pL->tag()<<std::endl;
                // std::cout<<" Current Neighborhood "<<std::endl;
                // for (const auto& neigh : neighbors())
                // {
                //     std::cout<<std::get<0>(neigh.second)->tag()<<"=>"<<std::get<1>(neigh.second)->tag()<<std::endl;
                // }
                const bool success=neighbors().emplace(pL->sink->sID,temp).second;
                if(!success)
                {
                    throw std::runtime_error("CANNOT INSERT IN NEIGHBORHOOD.");
                }
//                assert(success && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else if (pL->sink->sID==this->sID)
            {// this vertex is the sink of edge *pL, so source of *pL is the neighbor
                const NeighborType temp(pL->source.get(),pL);
                const bool success=neighbors().emplace(pL->source->sID,temp).second;
                if(!success)
                {
                    throw std::runtime_error("CANNOT INSERT IN NEIGHBORHOOD.");
                }
//                assert(success  && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else
            {
                throw std::runtime_error("CANNOT INSERT NON-INCIDENT EDGE");
//
//                assert(0 && "CANNOT INSERT NON-INCIDENT EDGE");
            }
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(NetworkLinkType* const pL)
        {
            if (pL->source->sID==this->sID)
            {
                const size_t key=pL->sink->sID;
                const int success=neighbors().erase(key);
                if(success!=1)
                {
                    throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
                }
//                assert(success==1);
            }
            else if (pL->sink->sID==this->sID)
            {
                const size_t key=pL->source->sID;
                const int success=neighbors().erase(key);
                if(success!=1)
                {
                    throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
                }
            }
            else
            {
                throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
            }
        }
                
        std::string tag() const
        {/*!\returns the string "i" where i is this->sID
          */
            return std::to_string(this->sID) ;
        }
    };
    
    template<typename Derived>
    int NetworkNode<Derived>::verboseLevel=0;
}
#endif
