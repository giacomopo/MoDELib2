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

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <set>
#include <memory>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <PolycrystallineMaterialBase.h>
#include <Polycrystal.h>
#include <MicrostructureBase.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <DefectiveCrystal.h>


namespace model
{

// https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html

#ifdef _MODEL_PYBIND11_
PYBIND11_MODULE(pyMoDELib,m)
{
    namespace py=pybind11;

    py::class_<model::DefectiveCrystalParameters>(m,"DefectiveCrystalParameters")
        .def(py::init<const std::string&>())
        .def_readwrite("runID", &model::DefectiveCrystalParameters::runID)
        .def_readonly("useFEM", &model::DefectiveCrystalParameters::useFEM)
    ;
    
    py::class_<model::MeshRegionObserver<model::MeshRegion<3>>>(m,"MeshRegionObserver")
        .def(py::init<>())
    ;
    
    py::class_<model::SimplexReader<3>>(m,"SimplexReader")
        .def(py::init<>())
    ;
    
    py::class_<std::map<typename model::SimplexTraits<3,3>::SimplexIDType,const model::Simplex<3,3>>>(m,"SimplexIDMap")
        .def(py::init<>())
    ;
    
    py::class_<std::map<std::pair<size_t,size_t>,model::MeshRegionBoundary<3>>>(m,"MeshRegionBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<model::SimplicialMesh<3>,
    /*      */ model::MeshRegionObserver<model::MeshRegion<3>>,
    /*      */ model::SimplexReader<3>,
    /*      */ std::map<typename model::SimplexTraits<3,3>::SimplexIDType,const model::Simplex<3,3>>,
    /*      */ std::map<std::pair<size_t,size_t>,model::MeshRegionBoundary<3>>>(m,"SimplicialMesh")
        .def(py::init<>())
//        .def(py::init<const std::string&,const Eigen::Matrix<double,3,3>&,const Eigen::Matrix<double,3,1>&,const std::set<int>&>())
        .def("xMin",static_cast<const Eigen::Matrix<double,3,1>& (model::SimplicialMesh<3>::*)() const>(&model::SimplicialMesh<3>::xMin))
        .def("xMax",static_cast<const Eigen::Matrix<double,3,1>& (model::SimplicialMesh<3>::*)() const>(&model::SimplicialMesh<3>::xMax))
        .def("volume",&model::SimplicialMesh<3>::volume)
    ;

    py::class_<model::PolycrystallineMaterialBase>(m,"PolycrystallineMaterialBase")
        .def(py::init<const std::string&,const double&>())
    ;

    py::class_<std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"GrainBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<model::Grain<3>,std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"Grain")
        .def(py::init<const MeshRegion<3>&,const PolycrystallineMaterialBase&,const std::string& >())
//            .def_readonly("grainID", &model::Grain<3>::grainID)

    ;
    
    py::class_<model::Polycrystal<3>,PolycrystallineMaterialBase>(m,"Polycrystal")
        .def(py::init<const std::string&,const SimplicialMesh<3>&>())
//        .def("mesh", &model::Polycrystal<3>::mesh)
        .def("randomPoint", &model::Polycrystal<3>::randomPoint)
        .def_readonly("grains", &model::Polycrystal<3>::grains)
        .def("grain", &model::Polycrystal<3>::grain)
    ;
    
    py::class_<model::DislocationDynamicsBase<3>>(m,"DislocationDynamicsBase")
        .def(py::init<const std::string&>())
//        .def("getMesh", &model::DislocationDynamicsBase<3>::getMesh)
        .def_readonly("poly", &model::DislocationDynamicsBase<3>::poly)
//        .def_readonly("mesh", &model::DislocationDynamicsBase<3>::mesh)
    ;
    
    py::class_<model::MicrostructureBase<3>>(m,"MicrostructureBase")
        .def("displacement", static_cast<Eigen::Matrix<double,Eigen::Dynamic,3> (model::MicrostructureBase<3>::*)(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const>(&model::MicrostructureBase<3>::displacement))
    ;
    
    py::class_<model::MicrostructureContainer<3>,model::MicrostructureBase<3>>(m,"MicrostructureContainer",py::multiple_inheritance())
        .def(py::init<model::DislocationDynamicsBase<3>&>())
    ;
 
    py::class_<model::DefectiveCrystal<3>,model::MicrostructureContainer<3>>(m,"DefectiveCrystal")
        .def(py::init<model::DislocationDynamicsBase<3>&>())
        .def("dislocationNetwork", &model::DefectiveCrystal<3>::dislocationNetwork)
    ;
    
//    py::class_<model::DislocationNetwork<3,0>,model::MicrostructureBase<3>,model::LoopNetwork<DislocationNetwork<3,0>>>(m,"DislocationNetwork")
//        .def(py::init<model::MicrostructureContainer<3>&>())
//    ;
//    
//    py::class_<model::LoopNetwork<model::DislocationNetwork<3,0>>,model::WeakPtrFactory<model::DislocationNetwork<3,0>,typename TypeTraits<model::DislocationNetwork<3,0>>::LoopType>>(m,"LoopNetwork",py::multiple_inheritance())
//        .def(py::init<>())
//        .def("loops", static_cast<const model::WeakPtrFactory<model::DislocationNetwork<3,0>,typename TypeTraits<model::DislocationNetwork<3,0>>::LoopType>& (model::LoopNetwork<model::DislocationNetwork<3,0>>::*)()const>(&model::LoopNetwork<model::DislocationNetwork<3,0>>::loops))
//    ;
}
#endif
}

int main(int argc, char** argv)
{
    
    return 0;
}


#endif
