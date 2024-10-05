/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeStepper_cpp_
#define model_DDtimeStepper_cpp_

#include <cfloat>
#include <DDtimeStepper.h>

namespace model
{

    template <typename DislocationNetworkType>
    DDtimeStepper<DislocationNetworkType>::DDtimeStepper(const DislocationNetworkType& DN_in):
    /* init */ DN(DN_in)
    /* init */,timeSteppingMethod(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readString("timeSteppingMethod", true))
    /* init */,dxMax(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("dxMax", true))
//    /* init */,shearWaveSpeedFraction(1.0e-7)
//    /* init */,dtMax(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("dtMax", true))
    {
        if (dxMax < FLT_EPSILON)
        {
            throw std::runtime_error("dxMax must be > FLT_EPSILON.");
        }
        

    }

    template <typename DislocationNetworkType>
    double DDtimeStepper<DislocationNetworkType>::getDt(const double& vRef, const double& vFactor) const
    {
        if(timeSteppingMethod=="fixed" || timeSteppingMethod=="Fixed" || timeSteppingMethod=="FIXED")
        {
            return DN.ddBase.simulationParameters.dtMax;
        }
        else if(timeSteppingMethod=="adaptive" || timeSteppingMethod=="Adaptive" || timeSteppingMethod=="ADAPTIVE")
        {
            double vmax = 0.0;
            long int vmaxID=-1;
            for (const auto &nodeIter : DN.networkNodes())
            {
                if(  !nodeIter.second.lock()->isBoundaryNode()
                   && nodeIter.second.lock()->glidePlanes().size()==1)
                {
                    const double vNorm(nodeIter.second.lock()->get_V().norm());
                    if (vNorm > vmax)
                    {
                        vmax = vNorm;
                        vmaxID=nodeIter.second.lock()->sID;
                    }
                }
            }
//            const double vRef(1.0);
            const double vEff(vmax+vRef*std::exp(-vmax/(vRef*vFactor)));
            std::cout<<" (vMax="<<vmax<<", vMaxID="<<vmaxID<<", vRef="<<vRef<<", vEff="<<vEff<<")"<<std::flush;
            return std::min(std::ceil(dxMax / vEff),DN.ddBase.simulationParameters.dtMax);
        }
        else
        {
            throw std::runtime_error("Unkonwn timeSteppingMethod "+timeSteppingMethod);
            return 0.0;
        }
    }

template struct DDtimeStepper<DislocationNetwork<3,0>>;

} // end namespace
#endif

