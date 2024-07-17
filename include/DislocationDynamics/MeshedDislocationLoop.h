/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshedDislocationLoop_h_
#define model_MeshedDislocationLoop_h_




#ifndef NDEBUG
#define VerboseMeshedDislocationLoop(N,x) if(verboseMeshedDislocationLoop>=N){std::cout<<x;}
#else
#define VerboseMeshedDislocationLoop(N,x)
#endif

#include <vector>
#include <TriangularMesh.h>
#include <Plane.h>

namespace model
{
    struct MeshedDislocationLoop : public TriangularMesh
    {
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        const VectorDim burgers;
        const std::vector<VectorDim> periodicShifts;
        const Plane<3> plane;
        std::deque<VectorDim> points;
        
        MeshedDislocationLoop(const VectorDim& burgers_in,const std::vector<VectorDim>& periodicShifts_in,const Plane<3>& plane_in,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize);
        VectorDim plasticDisplacement(const VectorDim&) const;
        Eigen::Matrix<double,Eigen::Dynamic,3> plasticDisplacement(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const;
        VectorDim triangleAreaVector(const Eigen::Vector3i&) const;
        double solidAngle(const VectorDim& x) const;
        
        template <typename T>
        static int sgn(const T& val)
        {
            return (val > T(0)) - (val < T(0));
        }

    };
    
}
#endif
