/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_pyMoDELib_cpp_
#define model_pyMoDELib_cpp_

#ifdef _MODEL_PYBIND11_
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#endif

#include <string>
#include <vector>
#include <memory>

#include <MicrostructureBase.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <DefectiveCrystal.h>


namespace model
{

#ifdef _MODEL_PYBIND11_
PYBIND11_MODULE(pyMoDELib,m)
{
    namespace py=pybind11;
    
    py::class_<model::DislocationDynamicsBase<3>>(m,"DislocationDynamicsBase")
        .def(py::init<const std::string&>())
    ;
    
    // https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
    py::class_<model::MicrostructureBase<3>>(m,"MicrostructureBase")
        .def("displacement", static_cast<Eigen::Matrix<double,Eigen::Dynamic,3> (model::MicrostructureBase<3>::*)(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const>(&model::MicrostructureBase<3>::displacement))
    ;
    
    py::class_<model::MicrostructureContainer<3>,model::MicrostructureBase<3>>(m,"MicrostructureContainer")
        .def(py::init<model::DislocationDynamicsBase<3>&>())
    ;
 
    py::class_<model::DefectiveCrystal<3>,model::MicrostructureContainer<3>>(m,"DefectiveCrystal")
        .def(py::init<model::DislocationDynamicsBase<3>&>())
    ;
}
#endif

}

int main(int argc, char** argv)
{
    
    return 0;
}


#endif
