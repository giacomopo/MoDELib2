/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationGlideSolverFactory_cpp_
#define model_DislocationGlideSolverFactory_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationGlideSolverFactory.h>
#include <GalerkinGlideSolver.h>
#include <PyGlideSolver.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    DislocationGlideSolverBase<DislocationNetworkType>::DislocationGlideSolverBase(const DislocationNetworkType& DN) :
//    /* init */ DN(DN_in)
    /* init */ DislocationVelocitySolverBase<DislocationNetworkType>(DN)
    {
        
    }

template <typename DislocationNetworkType>
std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>> DislocationGlideSolverFactory<DislocationNetworkType>::getGlideSolver(const DislocationNetworkType& DN,const std::string& solverType)
{
        if(solverType=="Galerkin" || solverType=="galerkin")
        {
            return std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>>(new GalerkinGlideSolver<DislocationNetworkType>(DN));
        }
        else if(solverType=="Pybind11" || solverType=="pybind11")
        {
            return std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>>(new PyGlideSolver<DislocationNetworkType>(DN));
        }
    else
    {
        std::cout<<redBoldColor<<"Unknown glide solver type "<<solverType<<". Glide disabled."<<defaultColor<<std::endl;
        return nullptr;
    }
    
}
    
    template struct DislocationGlideSolverBase<DislocationNetwork<3,0>>;
    template struct DislocationGlideSolverFactory<DislocationNetwork<3,0>>;

}
#endif
