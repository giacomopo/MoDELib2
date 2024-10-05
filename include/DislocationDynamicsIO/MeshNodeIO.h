/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshNodeIO_H_
#define model_MeshNodeIO_H_

#include <tuple>
#include <FEMnodeEvaluation.h>

namespace model
{
    
    template<short unsigned int dim>
    struct MeshNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        
        size_t nodeID;          // sID
        VectorDim displacement;
        
        /**********************************************************************/
        template<typename ElementType>
        MeshNodeIO(const FEMnodeEvaluation<ElementType,dim,1>& node) :
        /* init */ nodeID(node.pointID),
        /* init */ displacement(node)
        {
         
//            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        MeshNodeIO(const size_t& nodeID_in,          // sID
                   const VectorDim& disp) :
        /* init */ nodeID(nodeID_in),
        /* init */ displacement(disp)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        MeshNodeIO() :
        /* init */ nodeID(0),
        /* init */ displacement(VectorDim::Zero())
        {
            
            
        }
        
        /**********************************************************************/
        MeshNodeIO(std::stringstream& ss) :
        /* init */ nodeID(0),
        /* init */ displacement(VectorDim::Zero())
        {
            ss>>nodeID;
            for(int d=0;d<dim;++d)
            {
                ss>>displacement(d);
            }
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshNodeIO<dim>& ds)
        {
            os  << ds.nodeID<<"\t"
            /**/<< ds.displacement.transpose();
            return os;
        }
        
	};
	
}
#endif

