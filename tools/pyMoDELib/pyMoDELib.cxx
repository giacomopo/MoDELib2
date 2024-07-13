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
    
    py::class_<model::DefectiveCrystal<3>>(m,"DefectiveCrystal")
        .def(py::init<model::DislocationDynamicsBase<3>&>())
    //        .def("displacement", static_cast<Eigen::Matrix<double,Eigen::Dynamic,3> (model::MicrostructureBase<3>::*)(const Eigen::Matrix<double,Eigen::Dynamic,3>&) const>(&model::DefectiveCrystal<3>::displacement))
        .def("displacement", &model::DefectiveCrystal<3>::displacement)
    ;
}
#endif

}

int main(int argc, char** argv)
{
    
    return 0;
}


#endif
