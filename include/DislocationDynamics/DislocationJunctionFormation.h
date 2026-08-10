/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONJUNCTIONFORMATION_H_
#define model_DISLOCATIONJUNCTIONFORMATION_H_

#include <utility> // for std::pair
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <SegmentSegmentDistance.h>
#include <SweepPlane.h>
#include <TypeTraits.h>
#include <StressStraight.h>
#include <cassert>

#include <DislocationNetworkRemesh.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation
    {
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef typename TypeTraits<DislocationNetworkType>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<DislocationNetworkType>::NetworkNodeType NetworkNodeType;
        typedef typename TypeTraits<DislocationNetworkType>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef typename NetworkLinkType::KeyType KeyType;
        typedef std::tuple<KeyType,KeyType,SegmentSegmentDistance<dim>> IntersectionType;
        typedef std::deque<IntersectionType> IntersectionTypeContainerType;
        
        void insertIntersection(std::deque<IntersectionTypeContainerType> &intersectionContainer,
                                const NetworkLinkType *const linkA,
                                const NetworkLinkType *const linkB,
                                const SegmentSegmentDistance<dim> &ssd,
                                const double &currentcCollisionTOL,
                                const bool &bndJunction, const bool &gbndJunction);
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               const size_t& nThreads);
        std::pair<std::shared_ptr<NetworkNodeType>, bool> junctionNode(const double &t,
                                                                       const VectorDim &x,
                                                                       const std::shared_ptr<NetworkLinkType> &L);
        size_t contractJunctions(const std::deque<IntersectionTypeContainerType>& intersectionContainer);
        void glissileJunctions(const double &dx);
        
        void breakJunctions();
        
        DislocationNetworkType& DN;
        
    public:
                
        static double collisionTol;     //! The tolerance (in units of distance) used for collision detection
        const size_t maxJunctionIterations;
        const int verboseJunctions;
        const double infiniteLineLength;
        
        DislocationJunctionFormation(DislocationNetworkType& DN_in);
        static void initFromFile(const std::string& fileName);
        void formJunctions(const double& dx);
    };
    
}
#endif
