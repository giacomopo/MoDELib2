/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_H_
#define model_DislocationNodeContraction_H_

#include <memory>
#include <Eigen/Dense>
#include <TypeTraits.h>

namespace model
{
    template <typename DislocationNetworkType>
    class DislocationNodeContraction
    {
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef typename TypeTraits<DislocationNetworkType>::NetworkNodeType NetworkNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        std::vector<std::set<std::shared_ptr<NetworkNodeType>>> contractBoundary(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        bool contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        bool contractYoungest(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        bool contractToPosition(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB,const VectorDim& X, const double& maxRange);

        DislocationNetworkType& DN;

    public:
        
        const int verboseNodeContraction;

        DislocationNodeContraction(DislocationNetworkType& DN_in);
        bool contract(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
    };
}
#endif

