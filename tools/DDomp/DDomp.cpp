/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
// Define the non-singluar method used for calculations
#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <DefectiveCrystal.h>

using namespace model;

int main (int argc, char* argv[])
{

#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    pybind11::scoped_interpreter guard{};
#endif

    const std::string folderName(argc>1? std::string(argv[1]) : "./");
    
    try 
    {
        DislocationDynamicsBase<3> ddBase(folderName);
        DDconfigIO<3> configIO(ddBase.simulationParameters.traitsIO.evlFolder);
        configIO.read(ddBase.simulationParameters.runID);
        DefectiveCrystal<3> DC(ddBase);
        DC.initializeConfiguration(configIO);
        DC.runSteps();
    } 
    catch (const std::exception& e)
    {
        std::cout<<e.what()<<std::endl;
    }
    
    return 0;
}
