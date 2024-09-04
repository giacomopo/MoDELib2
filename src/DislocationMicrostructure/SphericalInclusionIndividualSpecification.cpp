/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionIndividualSpecification_cpp_
#define model_SphericalInclusionIndividualSpecification_cpp_

#include <SphericalInclusionIndividualSpecification.h>

namespace model
{
    SphericalInclusionIndividualSpecification::SphericalInclusionIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("SphericalInclusion","Individual") 
    /* init */,allowOverlap(false)
    /* init */,allowOutside(false)
    {
        
    }

    SphericalInclusionIndividualSpecification::SphericalInclusionIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("SphericalInclusion","Individual",fileName)
    /* init */,radii_SI(this->parser->readArray<double>("radii_SI",true))
    /* init */,centers(radii_SI.size()? this->parser->readMatrix<double>("centers",radii_SI.size(),dim,true) : Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(0,dim))
    /* init */,eigenDistortions(radii_SI.size()? this->parser->readMatrix<double>("eigenDistortions",radii_SI.size(),dim*dim,true) : Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(0,dim))
    /* init */,velocityReductionFactors(radii_SI.size()? this->parser->readArray<double>("velocityReductionFactors",true) : std::vector<double>())
    /* init */,phaseIDs(radii_SI.size()? this->parser->readArray<int>("phaseIDs",true) : std::vector<int>())
    /* init */,allowOverlap(radii_SI.size()? this->parser->readScalar<int>("allowOverlap",true) : false)
    /* init */,allowOutside(radii_SI.size()? this->parser->readScalar<int>("allowOutside",true) : false)
    {
        if(radii_SI.size())
        {
//            const std::vector<int> periodicDipoleExitFaceIDs(this->parser->readArray<int>("periodicDipoleExitFaceIDs",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,dim> centers(this->parser->readMatrix<double>("centers",radii_SI.size(),dim,true));
//            const Eigen::Matrix<double,Eigen::Dynamic,dim*dim> eigenDistortions(this->parser->readMatrix<double>("eigenDistortions",radii_SI.size(),dim*dim,true));
//            const std::vector<double> velocityReductionFactors(this->parser->readArray<double>("inclusionVelocityReductionFactors",true));
//            const std::vector<int> phaseIDs(this->parser->readArray<int>("phaseIDs",true));

            if(int(radii_SI.size())!=centers.rows())
            {
                throw std::runtime_error("radii_SI.size()="+std::to_string(radii_SI.size())+" NOT EQUAL TO centers.rows()="+std::to_string(centers.rows()));
            }
            if(int(radii_SI.size())!=eigenDistortions.rows())
            {
                throw std::runtime_error("radii_SI.size()="+std::to_string(radii_SI.size())+" NOT EQUAL TO eigenDistortions.rows()="+std::to_string(eigenDistortions.rows()));
            }
            if(radii_SI.size()!=velocityReductionFactors.size())
            {
                throw std::runtime_error("radii_SI.size()="+std::to_string(radii_SI.size())+" NOT EQUAL TO velocityReductionFactors.size()="+std::to_string(velocityReductionFactors.size()));
            }
            if(radii_SI.size()!=phaseIDs.size())
            {
                throw std::runtime_error("radii_SI.size()="+std::to_string(radii_SI.size())+" NOT EQUAL TO phaseIDs.size()="+std::to_string(phaseIDs.size()));
            }
        }
    }
}
#endif
