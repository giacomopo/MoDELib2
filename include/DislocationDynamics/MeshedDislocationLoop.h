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

#include <deque>
#include <vector>
#include <Plane.h>
#include <DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DislocationLoopPatches.h>

namespace model
{
    struct MeshedDislocationLoop //: public TriangularMesh
    {
        typedef Eigen::Matrix<double,3,1> VectorDim;
        typedef CompareVectorsByComponent<double,2,float> CompareType;
        std::map<Eigen::Matrix<double,2,1>,int,CompareType> uniquePointsIDs;

        VectorDim burgers;
        std::vector<VectorDim> periodicShifts;
        //Plane<3> plane;
        std::deque<VectorDim> points;
        std::deque<VectorDim> displacements;

        std::deque<Eigen::Vector3i> triangles;
        
    //    MeshedDislocationLoop(const VectorDim& burgers_in,const DislocationDynamicsBase<3>& ddBase,const GlidePlane<3>& plane_in,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize,const double& localMeshSize);
//        MeshedDislocationLoop(const std::shared_ptr<PeriodicPlanePatch<3>>& patch,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const DislocationNetwork<3,0>& DN,const double& meshSize,const double& localMeshSize);
        MeshedDislocationLoop(const VectorDim& burgers_in,const GlidePlane<3>& plane,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const DislocationNetwork<3,0>& DN,const double& meshSize,const double& localMeshSize);


        
        VectorDim plasticDisplacementKernel(const Eigen::Ref<const VectorDim>& x) const;
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
