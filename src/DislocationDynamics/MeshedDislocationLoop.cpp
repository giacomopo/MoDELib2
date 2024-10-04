/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshedDislocationLoop_cpp_
#define model_MeshedDislocationLoop_cpp_

#include <numbers>
#include <MeshedDislocationLoop.h>
#include <TriangularMesh.h>


namespace model
{



//MeshedDislocationLoop::MeshedDislocationLoop(const VectorDim& burgers_in,const DislocationDynamicsBase<3>& ddBase,const Plane<3>& plane,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize):
///* init */ burgers(burgers_in)
///* init */,periodicShifts(DN.ddBase.periodicShifts)
/////* init */,plane(plane_in)
//{
//
//    std::deque<Eigen::Matrix<double,2,1>> localBndPts;
//    std::deque<Eigen::Matrix<double,2,1>> internalPts;
//    for(const auto& gpt : globalBndPts)
//    {
//        localBndPts.push_back(plane.localPosition(gpt));
//    }
//    TriangularMesh triMesh;
//    triMesh.reMesh(localBndPts,internalPts,meshSize);
//    triangles=triMesh.triangles();
//
//    for(const auto& v : triMesh.vertices())
//    {
//        points.push_back(plane.globalPosition(v));
//    }
//}

    MeshedDislocationLoop::MeshedDislocationLoop(const VectorDim& burgers_in,const GlidePlane<3>& plane,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const DislocationNetwork<3,0>& DN,const double& meshSize,const double& localMeshSize):
    //MeshedDislocationLoop::MeshedDislocationLoop(const VectorDim& burgers_in,const DislocationDynamicsBase<3>& ddBase,const GlidePlane<3>& plane,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize,const double& localMeshSize):
    /* init */ burgers(burgers_in)
    /* init */,periodicShifts(DN.ddBase.periodicShifts)
    {
        std::vector<Eigen::Matrix<double,2,1>> localBndPts;
        localBndPts.reserve(globalBndPts.size());
        for(const auto& gpt : globalBndPts)
        {
            const Eigen::Matrix<double,2,1> localPos(plane.localPosition(gpt));
            const auto pointIter(uniquePointsIDs.find(localPos));
            if(pointIter==uniquePointsIDs.end())
            {// localPos does not exist
                uniquePointsIDs.emplace(localPos,localBndPts.size());
                localBndPts.emplace_back(localPos);
            }
        }
        const auto originalBndPts(localBndPts);
        
        std::vector<Eigen::Matrix<int,2,1>> segments;
        segments.reserve(localBndPts.size());
        for(size_t k=0;k<localBndPts.size();++k)
        {
            const size_t k1(k<localBndPts.size()-1? k+1 : 0); // cycle
            const auto& pointK(localBndPts[k]);
            const auto& pointK1(localBndPts[k1]);
            const auto iterK(uniquePointsIDs.find(pointK));
            const auto iterK1(uniquePointsIDs.find(pointK1));
            
            if(iterK!=uniquePointsIDs.end() && iterK1!=uniquePointsIDs.end())
            {// localPos exists
                segments.emplace_back((Eigen::Matrix<int,2,1>()<<iterK->second,iterK1->second).finished());
            }
            else
            {
                throw std::runtime_error("iterK and iterK1 must exist");
            }
        }
        
        if(localMeshSize>0.0)
        {
            const double L(std::pow(DN.ddBase.mesh.volume(),1.0/3.0));
            
            for(const auto& otherWeakLoop : DN.loops())
            {
                for(const auto& otherPatchPair : otherWeakLoop.second.lock()->patches().localPatches())
                {
                    const auto otherGlidePlane(otherPatchPair.first->glidePlane);
                    if(otherGlidePlane.get()!=&plane)
                    {// intersection of current plane with another plane
                        const PlanePlaneIntersection ppi(plane,*otherGlidePlane);
                        if(ppi.type==PlanePlaneIntersection<3>::INCIDENT)
                        {
                            const auto localLineOrigin(plane.localPosition(ppi.P));
                            const auto localLineDir(plane.localDirection(ppi.d));
                            
                            std::map<float,Eigen::Matrix<double,2,1>> patchLineIntersections;
                            for(size_t k=0;k<originalBndPts.size();++k)
                            {
                                const size_t k1(k<originalBndPts.size()-1? k+1 : 0); // cycle
                                const auto& pointK(originalBndPts[k]);
                                const auto& pointK1(originalBndPts[k1]);
                                SegmentSegmentDistance<2> ssd(localLineOrigin-2.0*L*localLineDir,localLineOrigin+2.0*L*localLineDir,pointK,pointK1);
                                if(ssd.dMin<FLT_EPSILON)
                                {
                                    patchLineIntersections.emplace(ssd.t,0.5*(ssd.x0+ssd.x1));
                                }
                            }
                            
                            if(patchLineIntersections.size()>=2)
                            {
                                const Eigen::Matrix<double,2,1>& p0(patchLineIntersections.begin()->second);
                                const Eigen::Matrix<double,2,1>& p1(patchLineIntersections.rbegin()->second);
                                const Eigen::Matrix<double,2,1> chord(p1-p0);
                                const int np=chord.norm()/localMeshSize;
                                
                                bool oldPointIncluded=false;
                                size_t oldPointID=0;
                                for(int n=0;n<np+1;++n)
                                {
                                    const double u(double(n)/np);
                                    const Eigen::Matrix<double,2,1> P(p0+u*chord);
                                    const Eigen::Matrix<double,2,1> otherP(otherGlidePlane->localPosition(plane.globalPosition(P)));
                                    const int wn(Polygon2D::windingNumber(otherP,otherPatchPair.second));
                                    const bool newPointIncluded(wn!=0);
                                    const auto pliIter(uniquePointsIDs.find(P));
                                    const bool newPointFound(pliIter!=uniquePointsIDs.end());
                                    const size_t newPointID(newPointFound? pliIter->second : localBndPts.size() );
                                    
                                    if(newPointIncluded)
                                    {
                                        if(!newPointFound)
                                        {// a new point
                                            uniquePointsIDs.emplace(P,localBndPts.size());
                                            localBndPts.emplace_back(P);
                                        }
                                        if(oldPointIncluded)
                                        {// both new and old points are inlcuded, add segment
                                            segments.emplace_back((Eigen::Matrix<int,2,1>()<<oldPointID,newPointID).finished());
                                        }
                                    }
                                    oldPointIncluded=newPointIncluded;
                                    oldPointID=newPointID;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        TriangularMesh triMesh;
        triMesh.reMesh(localBndPts,segments,meshSize,"pazq");
        triangles=triMesh.triangles();
        
        for(const auto& v : triMesh.vertices())
        {
            points.push_back(plane.globalPosition(v)); // initialize points
            displacements.push_back(VectorDim::Zero()); // initialize displacements
        }
        
        for(const auto& tri : triMesh.triangles())
        {
            const auto& v0(points[tri(0)]);
            const auto& v1(points[tri(1)]);
            const auto& v2(points[tri(2)]);
            const auto c(1.0/3.0*(v0+v1+v2));
            const auto h(c+(v1-v0).cross(v2-v1));
            
            MatrixDim vertexMatrix(MatrixDim::Zero());
            vertexMatrix.col(0)=v0-h;
            vertexMatrix.col(1)=v1-h;
            vertexMatrix.col(2)=v2-h;
            
            heightPoints.push_back(h); // initialize points
            invVertexMatrix.push_back(vertexMatrix.inverse());
            heightDisplacements.push_back(VectorDim::Zero()); // initialize displacements
        }
        
        defGradients.resize(triangles.size(),MatrixDim::Identity());
    }

void MeshedDislocationLoop::update()
{
    for(size_t k=0;k<points.size();++k)
    {
        points[k]+=displacements[k];
    }
    
    for(size_t k=0;k<heightPoints.size();++k)
    {
        heightPoints[k]+=heightDisplacements[k];
    }
    
    MatrixDim vertexMatrix(MatrixDim::Zero());
    for(size_t k=0;k<triangles.size();++k)
    {
        const auto& tri(triangles[k]);
                
        const auto& v0(points[tri(0)]);
        const auto& v1(points[tri(1)]);
        const auto& v2(points[tri(2)]);
        const auto h(heightPoints[k]);
        vertexMatrix.col(0)=v0-h;
        vertexMatrix.col(1)=v1-h;
        vertexMatrix.col(2)=v2-h;
        defGradients[k]=vertexMatrix*invVertexMatrix[k];
    }
}


    typename MeshedDislocationLoop::VectorDim MeshedDislocationLoop::triangleAreaVector(const Eigen::Vector3i& tri) const
    {
        VectorDim nA(VectorDim::Zero());
        const VectorDim& Ps(points[tri(0)]);
        for(int k0=0;k0<tri.size();++k0)
        {
            const size_t k1(k0+1<tri.size()? k0+1 : 0);
            const VectorDim& P0(points[tri(k0)]);
            const VectorDim& P1(points[tri(k1)]);
            nA+= 0.5*(P0-Ps).cross(P1-P0);
        }
        return nA;
    }

    double MeshedDislocationLoop::solidAngle(const VectorDim& x,const size_t& triID) const
    {
        const double a2(1.0e-2);
        const double oneA2=sqrt(1.0+a2);
        double temp(0.0);
        //for(const auto& tri : triangles)
//        for(size_t triID=0;triID<triangles.size();++triID)
//        {
            const auto& tri(triangles[triID]);
            const auto nA(triangleAreaVector(tri));
            const double triangleArea(nA.norm());
            if(triangleArea>FLT_EPSILON)
            {
                const VectorDim triangleNormal(nA/triangleArea);
                const VectorDim& planePoint(points[tri(0)]);
                const double posNorm((x-planePoint).norm());
                const double dotProd((x-planePoint).dot(triangleNormal));
                if(std::fabs(dotProd)>FLT_EPSILON*posNorm)
                {
                    const VectorDim s(sgn(dotProd)*triangleNormal); // s points along +n for points above, and along -n for points below
                    for(int k=0;k<tri.size();++k)
                    {
                        const size_t k1(k+1<tri.size()? k+1 : 0);
                        VectorDim e1(points[tri(k)]-x);
                        const double e1Norm(e1.norm());
                        if(e1Norm>FLT_EPSILON)
                        {
                            e1/=e1Norm;
                            VectorDim Y1(points[tri(k1)]-x);
                            const double Y1norm(Y1.norm());
                            if(Y1norm>FLT_EPSILON)
                            {
                                Y1/=Y1norm;
                                VectorDim e3(e1.cross(Y1));
                                const double e3Norm(e3.norm());
                                if(e3Norm>FLT_EPSILON)
                                {// e1 and Y1 are not align. If they are the projection on the unit sphere is a point and therefore there is no contribution to solid angle
                                    e3/=e3Norm; // normalize e3
                                    const VectorDim e2(e3.cross(e1));
                                    const double ydy(e1.dot(Y1));
                                    const double w=sqrt((1.0-ydy)/(1.0+ydy));
                                    const double s3(s.dot(e3));
                                    const double s3A2=sqrt(std::pow(s3,2)+a2);
                                    temp+=2.0*s3/oneA2/s3A2*atan(s3A2*w/(oneA2-s.dot(e1)-s.dot(e2)*w));
                                }
                                
                            }
                        }
                    }
                }
            }
//        }
        return temp;
    }


    typename MeshedDislocationLoop::VectorDim MeshedDislocationLoop::plasticDisplacementKernel(const Eigen::Ref<const VectorDim>& x) const
    {
        VectorDim temp(VectorDim::Zero());
        const bool applyDefGradient(true);
        if(applyDefGradient)
        {
            for(const auto& shift : periodicShifts)
            {
                for(size_t triID=0;triID<triangles.size();++triID)
                {
                    temp-=solidAngle(x+shift,triID)/4.0/std::numbers::pi*defGradients[applyDefGradient]*burgers;
                }
            }
        }
        else
        {
            for(const auto& shift : periodicShifts)
            {
                for(size_t triID=0;triID<triangles.size();++triID)
                {
                    temp-=solidAngle(x+shift,triID)/4.0/std::numbers::pi*burgers;
                }
            }
        }
        return temp;
    }

    Eigen::Matrix<double,Eigen::Dynamic,3> MeshedDislocationLoop::plasticDisplacement(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>> points) const
    {
        Eigen::Matrix<double,Eigen::Dynamic,3> temp(Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(points.rows(),3));
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
        for(long int k=0;k<points.rows();++k)
        {
            temp.row(k)=plasticDisplacementKernel(points.row(k));
        }
        return temp;
    }


}
#endif
