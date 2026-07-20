/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopGenerator_cpp_
#define model_PrismaticLoopGenerator_cpp_

#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <limits>
#include <cmath>
#include <numbers> // std::numbers

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PrismaticLoopGenerator.h>
#include <PlanesIntersection.h>
#include <LatticePlane.h>

namespace model
{

//PrismaticLoopGenerator::PrismaticLoopGenerator(const std::string& fileName) :
///* init */ MicrostructureGeneratorBase(fileName)
//{
//    
//}

PrismaticLoopGenerator::PrismaticLoopGenerator(const PrismaticLoopDensitySpecification& spec,MicrostructureGenerator& mg)
{
    
    std::cout<<magentaBoldColor<<"Generating prismatic loop density"<<defaultColor<<std::endl;
    if(spec.targetDensity>0.0)
    {
        std::normal_distribution<double> radiusDistribution(spec.radiusDistributionMean/mg.ddBase.poly.b_SI,spec.radiusDistributionStd/mg.ddBase.poly.b_SI);
        std::mt19937 generator;
        
        std::uniform_int_distribution<> allowedSlipSystemDist(0,spec.allowedSlipSystemIDs.size()-1);

        std::set<int> allowedGrainIDsSet;
        for(const auto& gID : spec.allowedGrainIDs)
        {
            allowedGrainIDsSet.insert(gID);
        }
        const double allowAllGrains(allowedGrainIDsSet.size() ? (*allowedGrainIDsSet.begin())<0 : true);
        
        double density=0.0;
        while(density<spec.targetDensity)
        {
            const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<3> L0=rp.first;
            const size_t grainID=rp.second;
                        
            if(allowAllGrains || allowedGrainIDsSet.find(grainID)!=allowedGrainIDsSet.end())
            {
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID)->slipSystems().size()-1);
                
                const int allowedIndex(allowedSlipSystemDist(generator));
                const int allowedSlipID(spec.allowedSlipSystemIDs[allowedIndex]);
                const int rSS(allowedSlipID<0 ? ssDist(generator) : allowedSlipID); // a random SlipSystem
                const double radius(radiusDistribution(generator));
                try
                {
                    density+=generateSingle(mg,rSS,L0.cartesian(),radius,50.0)/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"prismatic loop density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    std::cout<<e.what()<<std::endl;
                }
            }

        }
    }
}

PrismaticLoopGenerator::PrismaticLoopGenerator(const PrismaticLoopIndividualSpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual prismatic loops"<<defaultColor<<std::endl;
    if(spec.slipSystemIDs.size())
    {
        if(spec.slipSystemIDs.size()!=spec.loopRadii.size())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.loopRadii.size()="+std::to_string(spec.loopRadii.size()));
        }
        if(int(spec.slipSystemIDs.size())!=spec.loopCenters.rows())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
        }
        if(spec.slipSystemIDs.size()!=spec.glideSteps.size())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.glideSteps.size()="+std::to_string(spec.glideSteps.size()));
        }
        
        for(size_t k=0;k<spec.slipSystemIDs.size();++k)
        {
            generateSingle(mg,spec.slipSystemIDs[k],spec.loopCenters.row(k),spec.loopRadii[k]/mg.ddBase.poly.b_SI,spec.glideSteps[k]);
        }
    }
}

double PrismaticLoopGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& guessCenter,const double& radius,const double& L)
{
    std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(guessCenter));
    if(found.first)
    {
        const int grainID(found.second->region->regionID);
        assert(mg.ddBase.poly.grains.size()==1 && "Prismatic dislocations only supported for single crystals");
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        if(rSS>=0 && rSS<int(grain->slipSystems().size()))
        {
            // sort slip systems in the zone axis with incrasing angle from a reference
            const RationalLatticeDirection<3>  b(grain->slipSystems().at(rSS)->s);
            const VectorDimD unitAxis(b.cartesian().normalized());
            const VectorDimD refNormal(Plane<3>::getL2G(unitAxis).col(0));
            std::map<double,ReciprocalLatticeDirection<3>> ssMap;
            for(const auto& ss : grain->slipSystems())
            {
                if((ss.second->s-b).squaredNorm()<FLT_EPSILON)
                {
                    const double cosAngle(ss.second->unitNormal.dot(refNormal));
                    const double sinAngle(ss.second->unitNormal.dot(unitAxis.cross(refNormal)));
                    const double refAngle1(std::atan2(sinAngle,cosAngle));
                    const double refAngle2(std::atan2(-sinAngle,-cosAngle));
                    
                    if(std::fabs(refAngle1)<std::fabs(refAngle2))
                    {
                        ssMap.emplace(refAngle1,ss.second->n*(+1));
                    }
                    else
                    {
                        ssMap.emplace(refAngle2,ss.second->n*(-1));
                    }
                }
            }
            
            // Recenter and define sessile plane

            const VectorDimD center(grain->snapToLattice(guessCenter).cartesian());
            const ReciprocalLatticeDirection<3> axis(grain->reciprocalLatticeDirection(b.cartesian()));
            const int sessilePlaneHeight(axis.planeIndexOfPoint(center));
//            const LatticePlaneKey sessilePlaneKey(axis,sessilePlaneHeight,grain->grainID);
            const LatticePlaneKey sessilePlaneKey(sessilePlaneHeight,axis);
            const auto sessilePlane(mg.ddBase.periodicGlidePlaneFactory.getFromKey(sessilePlaneKey));
            
            // Define prismatic planes.
            std::vector<LatticePlane> planes;
            
            for(const auto& pair : ssMap)
            {// first half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center+pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            for(const auto& pair : ssMap)
            {// second half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center-pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            
            // Intersect two consecutive prismatic planes and sessile plane to find polygon points
            std::vector<VectorDimD> sessileLoopNodePos;
            for(size_t k1=0;k1<planes.size();++k1)
            {
                size_t k2(k1==planes.size()-1? 0 : k1+1);
                const auto& plane1(planes[k1]);
                const auto& plane2(planes[k2]);
                
                Eigen::Matrix<double,3,3> N;
                N.col(0)=sessilePlane->referencePlane->unitNormal;
                N.col(1)=plane1.n.cartesian().normalized();
                N.col(2)=plane2.n.cartesian().normalized();
                
                Eigen::Matrix<double,3,3> P;
                P.col(0)=sessilePlane->referencePlane->P;
                P.col(1)=plane1.planeOrigin();
                P.col(2)=plane2.planeOrigin();
                
                PlanesIntersection<3> pi(N,P,FLT_EPSILON);
                const std::pair<bool,VectorDimD> intPair(pi.snap(VectorDimD::Zero()));
                if(intPair.first)
                {
                    sessileLoopNodePos.push_back(intPair.second);
                }
                else
                {
                    throw std::runtime_error("Cannot determine prismatic polygon.");
                }
                
            }
            
            double loopLength=0.0;
            
            
            const bool sessileLoopInserted(mg.insertJunctionLoop(sessileLoopNodePos,sessilePlane,
                                  b.cartesian(),sessilePlane->referencePlane->unitNormal,
                                  center,grain->grainID,DislocationLoopIO<3>::SESSILELOOP));
            
            if(sessileLoopInserted)
            {
                for(size_t k1=0;k1<sessileLoopNodePos.size();++k1)
                {
                    size_t k2(k1==sessileLoopNodePos.size()-1? 0 : k1+1);
                    loopLength+=(sessileLoopNodePos[k2]-sessileLoopNodePos[k1]).norm();
                }

                const VectorDimD step(L*b.cartesian());
                for(size_t k1=0;k1<planes.size();++k1)
                {
                    size_t k2(k1==planes.size()-1? 0 : k1+1);
                    
                    Eigen::Matrix<double,3,1> shift(Eigen::Matrix<double,3,1>::Zero());
                    for(int c=0;c<mg.ddBase.periodicLatticeBasis.cols();++c)
                    {
                        shift+=mg.ddBase.periodicLatticeBasis.col(c)*std::floor(mg.ddBase.periodicLatticeReciprocalBasis.col(c).dot(sessileLoopNodePos[k2]-mg.ddBase.mesh.xMin()));
                    }
                                        
                    if(mg.ddBase.mesh.searchRegion(grainID,sessileLoopNodePos[k2]-shift).first)
                    {
                        std::vector<VectorDimD> loopNodePos;
                        loopNodePos.push_back(sessileLoopNodePos[k2]-shift);
                        loopNodePos.push_back(sessileLoopNodePos[k1]-shift);
                        loopNodePos.push_back(sessileLoopNodePos[k1]+step-shift);
                        loopNodePos.push_back(sessileLoopNodePos[k2]+step-shift);
                        
                        const int planeHeight(planes[k2].n.planeIndexOfPoint(loopNodePos.front()));
                        const LatticePlaneKey planeKey(planeHeight,planes[k2].n);
                        const auto periodicGlidePlane(mg.ddBase.periodicGlidePlaneFactory.getFromKey(planeKey));
                        
                        mg.insertJunctionLoop(loopNodePos,periodicGlidePlane,
                                              b.cartesian(),periodicGlidePlane->referencePlane->unitNormal,
                                              periodicGlidePlane->referencePlane->P,grain->grainID,DislocationLoopIO<3>::GLISSILELOOP);
                    }
                    else
                    {
                        throw std::runtime_error("Mesh does not contain shifted point.");
                    }
                }
            }
            return loopLength;
        }
        else
        {
            if(rSS<0)
            {
                std::cout<<"Skipping slip system "<<rSS<<std::endl;
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
    }
    else
    {
        std::cout<<"Center outside mesh, skipping prismatic loop"<<std::endl;
    }
    return 0.0;
}

}
#endif
