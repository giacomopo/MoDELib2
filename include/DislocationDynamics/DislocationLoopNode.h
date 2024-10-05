/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNode_h_
#define model_DislocationLoopNode_h_

#include <PeriodicGlidePlane.h>
// #include <PeriodicGlidePlane.hpp>

#include <DislocationDynamicsModule.h>

#ifndef NDEBUG
#define VerboseDislocationLoopNode(N,x) if(verboseDislocationLoopNode>=N){std::cout<<x;}
#else
#define VerboseDislocationLoopNode(N,x)
#endif

namespace model
{
    
    template <int dim, short unsigned int corder>
    class DislocationLoopNode : public LoopNode<DislocationLoopNode<dim,corder>>
    /*                       */,public SplineNode<DislocationLoopNode<dim,corder>,dim,corder,Hermite>
    {
        
        
        std::shared_ptr<PeriodicPlanePatch<dim>> _periodicPlanePatch;
        
        public:
        
        typedef DislocationLoopNode<dim,corder> DislocationLoopNodeType;
        typedef TypeTraits<DislocationLoopNodeType> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;
        typedef SplineNode<DislocationLoopNode<dim,corder>,dim,corder,Hermite> SplineNodeType;

        
        static int verboseDislocationLoopNode;
        const std::pair<const std::shared_ptr<PeriodicPlaneEdge<dim>>,const std::shared_ptr<PeriodicPlaneEdge<dim>>> periodicPlaneEdge; 
        //First will always be populated if the ndoe is on an edge..Second will only be populated if the node belongs to multiple edges        

    public:
        DislocationLoopNode(LoopNetworkType* const,
                            const std::shared_ptr<LoopType>&,
                            const std::shared_ptr<NetworkNodeType>&,
                            const VectorDim&,
                            const std::shared_ptr<PeriodicPlanePatch<dim>>&,
                            const std::pair<const std::shared_ptr<PeriodicPlaneEdge<dim>>,const std::shared_ptr<PeriodicPlaneEdge<dim>>>&
                            );

        DislocationLoopNode(LoopNetworkType* const,
                            const std::shared_ptr<LoopType>&,
                            const std::shared_ptr<NetworkNodeType>&,
                            const LoopLinkType* const);

        
        
        std::shared_ptr<DislocationLoopNode> clone(const std::shared_ptr<LoopType>&,
                                             const std::shared_ptr<NetworkNodeType>&) const;
        
        const DislocationLoopNodeType* periodicPrev() const;
        const DislocationLoopNodeType* periodicNext() const;


        std::vector<DislocationLoopNodeType*> boundaryPrev() const;
        std::vector<DislocationLoopNodeType*> boundaryNext() const;
        
        void updateGeometry();

        void updateConfinedGeometry();

        void set_P(const VectorDim&);
        void set_P(const VectorLowerDim&);
        std::shared_ptr<PeriodicPlanePatch<dim>> periodicPlanePatch() const;
        std::pair<bool,size_t> isRemovable(const double& Lmin,const double& Lmax, const double& relAreaTh);
        // std::pair<bool,std::set<size_t>> isRemovable(const double& Lmin, const double& relAreaTh);
        bool isMovableTo(const VectorDim&) const;
        bool isContractableTo(const LoopNodeType* const other) const;
        bool isGeometricallyRemovable(const double& Lmin,const double& Lmax, const double& relAreaTh);
        static void initFromFile(const std::string&);

    };
    

    
}
#endif
