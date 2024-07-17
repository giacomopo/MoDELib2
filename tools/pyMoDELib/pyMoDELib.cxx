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
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/stl_bind.h>
#endif

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <set>
#include <memory>
#include <map>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <PolycrystallineMaterialBase.h>
#include <Polycrystal.h>
#include <MicrostructureBase.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <DefectiveCrystal.h>


using namespace model;
// https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html

typedef Eigen::Matrix<double,3,1> VectorDim;
typedef Eigen::Matrix<double,3,3> MatrixDim;
typedef GlidePlane<3> GlidePlaneType;
typedef DislocationNetwork<3,0> DislocationNetworkType;
typedef typename TypeTraits<DislocationNetworkType>::LoopNodeType LoopNodeType;
typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
typedef typename TypeTraits<DislocationNetworkType>::NetworkNodeType NetworkNodeType;

#ifdef _MODEL_PYBIND11_

PYBIND11_MAKE_OPAQUE(std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>);
PYBIND11_MAKE_OPAQUE(std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MeshedDislocationLoop>>);

PYBIND11_MODULE(pyMoDELib,m)
{
    namespace py=pybind11;

    py::class_<DefectiveCrystalParameters>(m,"DefectiveCrystalParameters")
        .def(py::init<const std::string&>())
        .def_readwrite("runID", &DefectiveCrystalParameters::runID)
        .def_readonly("useFEM", &DefectiveCrystalParameters::useFEM)
    ;
    
    py::class_<MeshRegionObserver<MeshRegion<3>>>(m,"MeshRegionObserver")
        .def(py::init<>())
    ;
    
    py::class_<SimplexReader<3>>(m,"SimplexReader")
        .def(py::init<>())
    ;
    
    py::class_<std::map<typename SimplexTraits<3,3>::SimplexIDType,const Simplex<3,3>>>(m,"SimplexIDMap")
        .def(py::init<>())
    ;
    
    py::class_<std::map<std::pair<size_t,size_t>,MeshRegionBoundary<3>>>(m,"MeshRegionBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<SimplicialMesh<3>,
    /*      */ MeshRegionObserver<MeshRegion<3>>,
    /*      */ SimplexReader<3>,
    /*      */ std::map<typename SimplexTraits<3,3>::SimplexIDType,const Simplex<3,3>>,
    /*      */ std::map<std::pair<size_t,size_t>,MeshRegionBoundary<3>>>(m,"SimplicialMesh")
        .def(py::init<>())
//        .def(py::init<const std::string&,const Eigen::Matrix<double,3,3>&,const Eigen::Matrix<double,3,1>&,const std::set<int>&>())
        .def("xMin",static_cast<const Eigen::Matrix<double,3,1>& (SimplicialMesh<3>::*)() const>(&SimplicialMesh<3>::xMin))
        .def("xMax",static_cast<const Eigen::Matrix<double,3,1>& (SimplicialMesh<3>::*)() const>(&SimplicialMesh<3>::xMax))
        .def("volume",&SimplicialMesh<3>::volume)
    ;

    py::class_<PolycrystallineMaterialBase>(m,"PolycrystallineMaterialBase")
        .def(py::init<const std::string&,const double&>())
    ;

    py::class_<std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"GrainBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<Grain<3>,std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"Grain")
        .def(py::init<const MeshRegion<3>&,const PolycrystallineMaterialBase&,const std::string& >())
//            .def_readonly("grainID", &Grain<3>::grainID)
    ;
    
    py::class_<Polycrystal<3>,PolycrystallineMaterialBase>(m,"Polycrystal")
        .def(py::init<const std::string&,const SimplicialMesh<3>&>())
//        .def("mesh", &Polycrystal<3>::mesh)
        .def("randomPoint", &Polycrystal<3>::randomPoint)
        .def_readonly("grains", &Polycrystal<3>::grains)
        .def("grain", &Polycrystal<3>::grain)
    ;
    
    py::class_<DislocationDynamicsBase<3>>(m,"DislocationDynamicsBase")
        .def(py::init<const std::string&>())
//        .def("getMesh", &DislocationDynamicsBase<3>::getMesh)
        .def_readonly("poly", &DislocationDynamicsBase<3>::poly)
//        .def_readonly("mesh", &DislocationDynamicsBase<3>::mesh)
    ;
    
    py::class_<MicrostructureBase<3>>(m,"MicrostructureBase")
        .def("displacement", static_cast<Eigen::Matrix<double,Eigen::Dynamic,3> (MicrostructureBase<3>::*)(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const>(&MicrostructureBase<3>::displacement))
    ;
    
    py::class_<MicrostructureContainer<3>,MicrostructureBase<3>>(m,"MicrostructureContainer",py::multiple_inheritance())
        .def(py::init<DislocationDynamicsBase<3>&>())
    ;
    
    // LoopNode
    py::bind_map<std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>>(m, "LoopNodeWeakPtrMap");

    py::class_<WeakPtrFactory<DislocationNetworkType,LoopNodeType>
    /*      */,std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>
    /*      */>(m,"LoopNodeWeakPtrFactory")
        .def(py::init<>())
        .def("getRef",&WeakPtrFactory<DislocationNetworkType,LoopNodeType>::getRef,pybind11::return_value_policy::reference)
    ;
    
    py::class_<LoopNodeType
    /*      */ >(m,"LoopNode")
        .def(py::init<typename TypeTraits<DislocationNetworkType>::LoopNetworkType* const,
             const std::shared_ptr<LoopType>&,
             const std::shared_ptr<NetworkNodeType>&,
             const typename TypeTraits<DislocationNetworkType>::VectorDim&,
             const std::shared_ptr<PeriodicPlanePatch<3>>&,
             const std::pair<const std::shared_ptr<PeriodicPlaneEdge<3>>,const std::shared_ptr<PeriodicPlaneEdge<3>>>&>())
    ;
    
    //Loop
    py::bind_vector<std::vector<std::shared_ptr<MeshedDislocationLoop>>>(m, "MeshedDislocationLoopVector");

    py::bind_map<std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>>(m, "LoopWeakPtrMap");

    py::class_<WeakPtrFactory<DislocationNetworkType,LoopType>
    /*      */,std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>
    /*      */>(m,"LoopWeakPtrFactory")
        .def(py::init<>())
        .def("getRef",&WeakPtrFactory<DislocationNetworkType,LoopType>::getRef,pybind11::return_value_policy::reference)
    ;
    
    py::class_<LoopType
    /*      */ >(m,"Loop")
        .def(py::init<DislocationNetworkType* const,
             const VectorDim&,
             const std::shared_ptr<GlidePlaneType>&>())
        .def("solidAngle",&LoopType::solidAngle)
        .def("meshed",&LoopType::meshed)
    ;
    
    py::class_<MeshedDislocationLoop
    /*      */ >(m,"MeshedDislocationLoop")
        .def(py::init<const VectorDim&,
             const std::vector<VectorDim>&,
             const Plane<3>&,const std::vector<Eigen::Matrix<double,3,1>>&,
             const double&>())
//        .def("solidAngle",&LoopType::solidAngle)
//        .def("meshed",&LoopType::meshed)
    ;

    py::class_<LoopNetwork<DislocationNetworkType>
    /*      */,WeakPtrFactory<DislocationNetworkType,LoopType>
    /*      */,WeakPtrFactory<DislocationNetworkType,LoopNodeType>
    /*      */ >(m,"LoopNetwork",py::multiple_inheritance())
        .def(py::init<>())
        .def("loops", static_cast<const WeakPtrFactory<DislocationNetworkType,LoopType>& (LoopNetwork<DislocationNetworkType>::*)()const>(&LoopNetwork<DislocationNetworkType>::loops),pybind11::return_value_policy::reference)
        .def("loopNodes", static_cast<const WeakPtrFactory<DislocationNetworkType,LoopNodeType>& (LoopNetwork<DislocationNetworkType>::*)()const>(&LoopNetwork<DislocationNetworkType>::loopNodes),pybind11::return_value_policy::reference)
    ;
    
    py::class_<DislocationNetworkType
    /*      */,MicrostructureBase<3>
    /*      */,LoopNetwork<DislocationNetworkType>>(m,"DislocationNetwork")
        .def(py::init<MicrostructureContainer<3>&>())
    ;
    
    py::class_<DefectiveCrystal<3>
    /*      */,MicrostructureContainer<3>
    /*      */>(m,"DefectiveCrystal")
        .def(py::init<DislocationDynamicsBase<3>&>())
        .def("dislocationNetwork", &DefectiveCrystal<3>::dislocationNetwork,pybind11::return_value_policy::reference)
    ;

}
#endif


int main(int argc, char** argv)
{
    return 0;
}

#endif
