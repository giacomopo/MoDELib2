/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationClimbSolverFactory_cpp_
#define model_DislocationClimbSolverFactory_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationClimbSolverFactory.h>
#include <GalerkinClimbSolver.h>
#include <TextFileParser.h>

namespace model
{

    template <typename DislocationNetworkType>
    DislocationClimbSolverBase<DislocationNetworkType>::DislocationClimbSolverBase(const DislocationNetworkType& DN,const ClusterDynamics<dim>* const CD_in) :
    /* init */ DislocationVelocitySolverBase<DislocationNetworkType>(DN)
    /* init */,glideEquilibriumRate(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("glideEquilibriumRate",true))
    /* init */,CD(CD_in)
    /* init */,vClimbRef(getVclimbRef())
    {
    }

    template <typename DislocationNetworkType>
    const typename DislocationClimbSolverBase<DislocationNetworkType>::ScalarVelocitiesContainerType& DislocationClimbSolverBase<DislocationNetworkType>::scalarVelocities() const
    {
        return *this;
    }

    template <typename DislocationNetworkType>
    typename DislocationClimbSolverBase<DislocationNetworkType>::ScalarVelocitiesContainerType& DislocationClimbSolverBase<DislocationNetworkType>::scalarVelocities()
    {
        return *this;
    }

    template <typename DislocationNetworkType>
    double DislocationClimbSolverBase<DislocationNetworkType>::getVclimbRef() const
    {
        double vRef(0.0);
            const auto eqC(CD->cdp.equilibriumMobileConcentration(0.0));
            for(const auto& pair : CD->cdp.D)
            {
                const auto& dVec(pair.second);
                for(int m=0;m<CD->cdp.mSize;++m)
                {
                    vRef=std::max(vRef,dVec[m].trace()*eqC[m]/3.0/this->DN.ddBase.poly.b);
                }
            }
        return  vRef;
    }

//    template <typename DislocationNetworkType>
//    const ClusterDynamics<DislocationClimbSolverFactory<DislocationNetworkType>::dim>*  DislocationClimbSolverFactory<DislocationNetworkType>::getCD(const DislocationNetworkType& DN)
//    {
//        const auto cdMicrostructures(DN.microstructures.template getTypedMicrostructures<ClusterDynamics<dim>>());
//        if(cdMicrostructures.size()==1)
//        {
//            return *cdMicrostructures.begin();
//        }
//        else
//        {
//            std::cout<<"Found "<<cdMicrostructures.size()<<" microstructures of type ClusterDynamics<"<<dim<<">"<<std::endl;
//            return nullptr;
//        }
//    }

    template <typename DislocationNetworkType>
    std::shared_ptr<DislocationClimbSolverBase<DislocationNetworkType>> DislocationClimbSolverFactory<DislocationNetworkType>::getClimbSolver(const DislocationNetworkType& DN,const std::string& solverType)
{
//        const ClusterDynamics<DislocationClimbSolverFactory<DislocationNetworkType>::dim>* CD(getCD(DN));
        const auto CD(DN.microstructures.template getUniqueTypedMicrostructure<ClusterDynamics<dim>>());
        if(CD)
        {
            if(solverType=="Galerkin" || solverType=="galerkin")
            {
                
                return std::shared_ptr<DislocationClimbSolverBase<DislocationNetworkType>>(new GalerkinClimbSolver<DislocationNetworkType>(DN,CD.get()));
            }
            else
            {
                std::cout<<redBoldColor<<"Unknown climb solver type "<<solverType<<". Climb disabled."<<defaultColor<<std::endl;
                return nullptr;
            }
        }
        else
        {
            return nullptr;
        }
    }

    template struct DislocationClimbSolverBase<DislocationNetwork<3,0>>;
    template struct DislocationClimbSolverFactory<DislocationNetwork<3,0>>;

}
#endif
