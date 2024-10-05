/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalAndInternalBoundary_H_
#define model_ExternalAndInternalBoundary_H_

#include <deque>
#include <utility>      // std::pair, std::make_pair
#include <Eigen/Dense>
#include <Simplex.h>
#include <IntegrationDomain.h>

namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
	struct ExternalAndInternalBoundary //: public IntegrationDomain<dim,1,qOrder,QuadratureRule>
    {

        /**************************************/
        template <typename FiniteElementType>
        ExternalAndInternalBoundary(const FiniteElementType& fe)
        {
            
        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return node.normalTensor().squaredNorm()>FLT_EPSILON;
        }
//
//        template <typename FiniteElementType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
//        static IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary(const FiniteElementType& fe)
//        {
//            IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> temp;
//
//            for (const auto& eIter : fe.elements())
//            {
//                if(eIter.second.isBoundaryElement())
//                {
//                    const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
//                    for (size_t f=0;f<boundaryFaces.size();++f)
//                    {
//                        temp.emplace_back(&eIter.second,boundaryFaces[f]);
//                    }
//                }
//            }
//            
//            return temp;
//        }
        
    };

}	// close namespace
#endif
