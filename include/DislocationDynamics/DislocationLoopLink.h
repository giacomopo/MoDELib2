/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopLink_h_
#define model_DislocationLoopLink_h_

#include <DislocationDynamicsModule.h>
#include <CatmullRomSplineSegment.h>

//#include <LatticeVector.h>
//#include <LatticePlaneBase.h>
//#include <LatticePlane.h>
//#include <Grain.h>
//#include <GlidePlane.h>
////#include <DislocationLoopIO.h>
////#include <PlanarDislocationLoop.h>
//#include <SlipSystem.h>

//#include <PlanarPolygon.h>
//#include <DislocationNetwork.h>

#ifndef NDEBUG
#define VerboseDislocationLoopLink(N,x) if(verboseDislocationLoopLink>=N){std::cout<<x;}
#else
#define VerboseDislocationLoopLink(N,x)
#endif


namespace model
{
    template <int dim, short unsigned int corder>
    class DislocationLoopLink : public LoopLink<DislocationLoopLink<dim,corder>>
    {

    public:
        
        typedef TypeTraits<DislocationLoopLink<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;

        static int verboseDislocationLoopLink;

        
        DislocationLoopLink(const std::shared_ptr<LoopNodeType>&,
                      const std::shared_ptr<LoopNodeType>&,
                      const std::shared_ptr<LoopType>&);
        
        bool hasNetworkLink() const;
        std::shared_ptr<PeriodicPlanePatch<dim>> periodicPlanePatch() const;
        CatmullRomSplineSegment<dim> spline() const;
        
        const LoopLinkType * twinnedLink() const ;
        static void initFromFile(const std::string&);
        
    };
    
}
#endif
