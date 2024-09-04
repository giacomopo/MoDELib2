/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionDensitySpecification_H_
#define model_SphericalInclusionDensitySpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct SphericalInclusionDensitySpecification : public MicrostructureSpecificationBase
    {
        
        static constexpr int dim=3;
        
        double targetDensity;
        double diameterLognormalDistribution_M;
        double diameterLognormalDistribution_S;
        double diameterLognormalDistribution_A;
        Eigen::Matrix<double,1,dim*dim> transformationEigenDistortion;
        Eigen::Matrix<double,1,dim> patternVector_SI;
        bool allowOverlap;
        bool allowOutside;
        double velocityReductionFactor;
        int phaseIDs;

        SphericalInclusionDensitySpecification();
        SphericalInclusionDensitySpecification(const std::string& fileName);
    };
}
#endif
