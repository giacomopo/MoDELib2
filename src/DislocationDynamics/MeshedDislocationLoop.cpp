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

namespace model
{
    
MeshedDislocationLoop::MeshedDislocationLoop(const VectorDim& burgers_in,const std::vector<VectorDim>& periodicShifts_in,const Plane<3>& plane_in,const std::vector<Eigen::Matrix<double,3,1>>& globalBndPts,const double& meshSize):
/* init */ burgers(burgers_in)
/* init */,periodicShifts(periodicShifts_in)
/* init */,plane(plane_in)
{
    
    std::deque<Eigen::Matrix<double,2,1>> localBndPts;
    std::deque<Eigen::Matrix<double,2,1>> internalPts;
    for(const auto& gpt : globalBndPts)
    {
        localBndPts.push_back(plane.localPosition(gpt));
    }
    this->reMesh(localBndPts,internalPts,meshSize);

    
    for(const auto& v : this->vertices())
    {
        points.push_back(plane.globalPosition(v));
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

double MeshedDislocationLoop::solidAngle(const VectorDim& x) const
{
    const double a2(1.0e-2);
    const double oneA2=sqrt(1.0+a2);
    double temp(0.0);
    for(const auto& tri : this->triangles())
    {
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
    }
    return temp;
}


typename MeshedDislocationLoop::VectorDim MeshedDislocationLoop::plasticDisplacement(const VectorDim& x) const
{
    VectorDim temp(VectorDim::Zero());
    for(const auto& shift : periodicShifts)
    {
        temp-=solidAngle(x+shift)/4.0/std::numbers::pi*burgers;
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
        temp.row(k)=plasticDisplacement(VectorDim(points.row(k)));
    }
    return temp;
}


}
#endif
