/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionDensitySpecification_cpp_
#define model_SphericalInclusionDensitySpecification_cpp_

#include <SphericalInclusionDensitySpecification.h>

namespace model
{
    SphericalInclusionDensitySpecification::SphericalInclusionDensitySpecification():
    /* init */ MicrostructureSpecificationBase("SphericalInclusion","Density")
    /* init */,targetDensity(0.0)
    /* init */,diameterLognormalDistribution_M(0.0)
    /* init */,diameterLognormalDistribution_S(0.0)
    /* init */,diameterLognormalDistribution_A(0.0)
    /* init */,transformationEigenDistortion(Eigen::Matrix<double,1,dim*dim>::Zero())
    /* init */,patternVector_SI(Eigen::Matrix<double,1,dim>::Zero())
    /* init */,allowOverlap(false)
    /* init */,allowOutside(false)
    /* init */,velocityReductionFactor(1.0)
    /* init */,phaseIDs(-1)
    {
        
    }

    SphericalInclusionDensitySpecification::SphericalInclusionDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("SphericalInclusion","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity",true))
    /* init */,diameterLognormalDistribution_M(targetDensity>0.0? this->parser->readScalar<double>("diameterLognormalDistribution_M",true) : 0.0 )
    /* init */,diameterLognormalDistribution_S(targetDensity>0.0? this->parser->readScalar<double>("diameterLognormalDistribution_S",true) : 0.0)
    /* init */,diameterLognormalDistribution_A(targetDensity>0.0? this->parser->readScalar<double>("diameterLognormalDistribution_A",true) : 0.0)
    /* init */,transformationEigenDistortion(targetDensity>0.0? this->parser->readMatrix<double>("transformationEigenDistortion",1,dim*dim,true) : Eigen::Matrix<double,1,dim*dim>::Zero())
    /* init */,patternVector_SI(targetDensity>0.0? this->parser->readMatrix<double>("patternVector_SI",1,dim,true) : Eigen::Matrix<double,1,dim>::Zero())
    /* init */,allowOverlap(targetDensity>0.0? this->parser->readScalar<int>("allowOverlap",true) : false)
    /* init */,allowOutside(targetDensity>0.0? this->parser->readScalar<int>("allowOutside",true) : false)
    /* init */,velocityReductionFactor(targetDensity>0.0? this->parser->readScalar<double>("velocityReductionFactor",true) : 1.0)
    /* init */,phaseIDs(targetDensity>0.0? this->parser->readScalar<int>("phaseID",true) : -1)
    {
        
    }
}
#endif
