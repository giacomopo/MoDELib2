/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <DislocationDynamicsBase.h>
#include <MicrostructureGenerator.h>

using namespace model;

int main(int argc, char** argv)
{
#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    pybind11::scoped_interpreter guard{};
#endif
    
    const std::string folderName(argc>1? std::string(argv[1]) : "./");
    
    DislocationDynamicsBase<3> ddBase(folderName);
    MicrostructureGenerator mg(ddBase);
    mg.readMicrostructureFile();
    mg.writeConfigFiles(0);
    return 0;
}
