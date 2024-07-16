/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityPy_h_
#define _model_DislocationMobilityPy_h_

#include <DislocationMobilityBase.h>

#ifdef _MODEL_PYBIND11_
    #undef slots
    #include <pybind11/embed.h>
    #include <pybind11/eigen.h>
    #include <pybind11/stl.h>
    #define slots Q_SLOTS
#endif

//#ifdef _MODEL_PYBIND11_
//    #ifdef slots
//        #define SAVED_SLOTS slots
//        #undef slots
//    #endif
//    #include <pybind11/embed.h>
//    #include <pybind11/eigen.h>
//    #include <pybind11/stl.h>
//    #ifdef SAVED_SLOTS
//        #define slots SAVED_SLOTS
//        #undef SAVED_SLOTS
//    #else
//        #define slots Q_SLOTS
//    #endif
//#endif


namespace model
{
    
    struct DislocationMobilityPy : public DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        const double kB;
        const std::string pyModuleName;

#ifdef _MODEL_PYBIND11_
        pybind11::module pyModule;
#endif

        DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in);

        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) override;
        
    };
    
}
#endif
