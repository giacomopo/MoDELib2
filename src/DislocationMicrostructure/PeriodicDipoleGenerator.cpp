/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleGenerator_cpp_
#define model_PeriodicDipoleGenerator_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>
#include <random>
#include <iomanip>

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
#include <PeriodicDipoleGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

//    PeriodicDipoleGenerator::PeriodicDipoleGenerator(const std::string& fileName) :
//    /* init */ MicrostructureGeneratorBase(fileName)
//    {
//        
//    }

PeriodicDipoleGenerator::PeriodicDipoleGenerator(const PeriodicDipoleIndividualSpecification& spec,MicrostructureGenerator& mg)
{
//    const std::vector<int> slipSystemIDs(this->parser.readArray<int>("slipSystemIDs",true));
    std::cout<<magentaBoldColor<<"Generating individual periodic dipoles"<<defaultColor<<std::endl;
    if(spec.slipSystemIDs.size())
    {
        if(spec.slipSystemIDs.size()!=spec.exitFaceIDs.size())
        {
            throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO exitFaceIDs.size()="+std::to_string(spec.exitFaceIDs.size()));
        }
        if(int(spec.slipSystemIDs.size())!=spec.dipoleCenters.rows())
        {
            throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO dipoleCenters.rows()="+std::to_string(spec.dipoleCenters.rows()));
        }
        if(spec.slipSystemIDs.size()!=spec.dipoleHeights.size())
        {
            throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO dipoleHeights.size()="+std::to_string(spec.dipoleHeights.size()));
        }
        if(spec.slipSystemIDs.size()!=spec.nodesPerLine.size())
        {
            throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO nodesPerLine.size()="+std::to_string(spec.nodesPerLine.size()));
        }
        if(spec.slipSystemIDs.size()!=spec.glideSteps.size())
        {
            throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO periodicDipoleGlideSteps.size()="+std::to_string(spec.glideSteps.size()));
        }
        
        for(size_t k=0;k<spec.slipSystemIDs.size();++k)
        {
            generateSingle(mg,spec.slipSystemIDs[k],spec.dipoleCenters.row(k),spec.exitFaceIDs[k],spec.dipoleHeights[k],spec.nodesPerLine[k],spec.glideSteps[k]);
        }
    }
    
}

PeriodicDipoleGenerator::PeriodicDipoleGenerator(const PeriodicDipoleDensitySpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating periodic dipole density"<<defaultColor<<std::endl;
//    const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
    if(spec.targetDensity>0.0)
    {
        std::mt19937 generator;
        double density=0.0;
        while(density<spec.targetDensity)
        {
            const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<3> L0=rp.first;
            const size_t grainID=rp.second;
            std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
            const int rSS(ssDist(generator)); // a random SlipSystem
//                const auto& slipSystem(*poly.grain(grainID).singleCrystal->slipSystems()[rSS]);
            std::uniform_int_distribution<> fDist(0,mg.ddBase.poly.grain(grainID).region.faces().size()-1);
            const int rF(fDist(generator)); // a random face
            auto faceIter(mg.ddBase.poly.grain(grainID).region.faces().begin());
            std::advance(faceIter,rF);

            try
            {
                generateSingle(mg,rSS,L0.cartesian(),faceIter->first,100,0,20.0);
                density+=2.0*faceIter->second->periodicFacePair.first.norm()/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                std::cout<<"periodic dipole density="<<density<<std::endl;
            }
            catch(const std::exception& e)
            {
                
            }
        }
    }
}

//    void PeriodicDipoleGenerator::generateDensity(MicrostructureGenerator& mg)
//    {
//        std::cout<<magentaBoldColor<<"Generating periodic dipole density"<<defaultColor<<std::endl;
//        const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
//        if(targetDensity>0.0)
//        {
//            std::mt19937 generator;
//            double density=0.0;
//            while(density<targetDensity)
//            {
//                const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
//                const LatticeVector<3> L0=rp.first;
//                const size_t grainID=rp.second;
//                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
//                const int rSS(ssDist(generator)); // a random SlipSystem
////                const auto& slipSystem(*poly.grain(grainID).singleCrystal->slipSystems()[rSS]);
//                std::uniform_int_distribution<> fDist(0,mg.ddBase.poly.grain(grainID).region.faces().size()-1);
//                const int rF(fDist(generator)); // a random face
//                auto faceIter(mg.ddBase.poly.grain(grainID).region.faces().begin());
//                std::advance(faceIter,rF);
//
//                try
//                {
//                    generateSingle(mg,rSS,L0.cartesian(),faceIter->first,100,0,20.0);
//                    density+=2.0*faceIter->second->periodicFacePair.first.norm()/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
//                    std::cout<<"periodic dipole density="<<density<<std::endl;
//                }
//                catch(const std::exception& e)
//                {
//                    
//                }
//            }
//        }
//    }

//    void PeriodicDipoleGenerator::generateIndividual(MicrostructureGenerator& mg)
//    {
//        const std::vector<int> slipSystemIDs(this->parser.readArray<int>("slipSystemIDs",true));
//        
//        if(slipSystemIDs.size())
//        {
//            std::cout<<magentaBoldColor<<"Generating individual periodic dipole"<<defaultColor<<std::endl;
//            const std::vector<int> spec.exitFaceIDs(this->parser.readArray<int>("spec.exitFaceIDs",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,3> periodicDipolePoints(this->parser.readMatrix<double>("periodicDipolePoints",slipSystemIDs.size(),dim,true));
//            const std::vector<double> periodicDipoleHeights(this->parser.readArray<double>("periodicDipoleHeights",true));
//            const std::vector<int> periodicDipoleNodes(this->parser.readArray<int>("periodicDipoleNodes",true));
//            const std::vector<double> periodicDipoleGlideSteps(this->parser.readArray<double>("periodicDipoleGlideSteps",true));
//            
//            
//            if(slipSystemIDs.size()!=spec.exitFaceIDs.size())
//            {
//                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(slipSystemIDs.size())+" NOT EQUAL TO spec.exitFaceIDs.size()="+std::to_string(spec.exitFaceIDs.size()));
//            }
//            if(int(slipSystemIDs.size())!=periodicDipolePoints.rows())
//            {
//                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(slipSystemIDs.size())+" NOT EQUAL TO periodicDipolePoints.size()="+std::to_string(periodicDipolePoints.size()));
//            }
//            if(slipSystemIDs.size()!=periodicDipoledipoleHeights.size())
//            {
//                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(slipSystemIDs.size())+" NOT EQUAL TO periodicDipoledipoleHeights.size()="+std::to_string(periodicDipoledipoleHeights.size()));
//            }
//            if(slipSystemIDs.size()!=periodicDipoleNodes.size())
//            {
//                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(slipSystemIDs.size())+" NOT EQUAL TO periodicDipoleNodes.size()="+std::to_string(periodicDipoleNodes.size()));
//            }
//            if(slipSystemIDs.size()!=periodicDipoleGlideSteps.size())
//            {
//                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(slipSystemIDs.size())+" NOT EQUAL TO periodicDipoleGlideSteps.size()="+std::to_string(periodicDipoleGlideSteps.size()));
//            }
//            
//            for(size_t k=0;k<slipSystemIDs.size();++k)
//            {
//                generateSingle(mg,slipSystemIDs[k],periodicDipolePoints.row(k),spec.exitFaceIDs[k],periodicDipoleHeights[k],periodicDipoleNodes[k],periodicDipoleGlideSteps[k]);
//            }
//        }
//        
//    }

    void PeriodicDipoleGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& dipolePoint,const int& exitFaceID,const int& dipoleHeight,const int& dipoleNodes, double glideStep)
    {
        
        
        if(rSS>=0)
        {
            std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(dipolePoint));
            if(!found.first)
            {
                std::cout<<"Point "<<dipolePoint.transpose()<<" is outside mesh. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            const int grainID(found.second->region->regionID);
            assert(mg.ddBase.poly.grains().size()==1 && "Periodic dislocations only supported for single crystals");
            const auto& grain(mg.ddBase.poly.grain(grainID));

            if(rSS<int(grain.singleCrystal->slipSystems().size()))
            {
                const auto periodicFaceIter(grain.region.faces().find(exitFaceID));
                if(periodicFaceIter!=grain.region.faces().end())
                {
                    const auto periodicFaceA(periodicFaceIter->second);
                    const auto periodicFaceB(periodicFaceA->periodicFacePair.second);

                    if(periodicFaceB!=nullptr)
                    {

                        const auto& faceAshift(periodicFaceA->periodicFacePair.first);
                        const auto faceAlatticeShift(grain.singleCrystal->latticeVector(faceAshift));
                        
                        const auto& slipSystem(*grain.singleCrystal->slipSystems()[rSS]);
                        
                        if(slipSystem.n.dot(faceAlatticeShift)==0)
                        {

                            //                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,dipolePoint);
                            
                            const VectorDimD planePoint(dipolePoint-0.5*dipoleHeight*slipSystem.n.interplaneVector());
                            const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(planePoint));
                            GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
                            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
                            //                        const VectorDimD P0(glidePlane->snapToPlane(dipolePoint));
                            const VectorDimD P0(grain.singleCrystal->snapToLattice(planePoint).cartesian()); // WARNING: this may shift the point compared to the input.
                            
                            PlaneLineIntersection<3> pliA(periodicFaceA->center(),periodicFaceA->outNormal(),P0,faceAshift);
                            PlaneLineIntersection<3> pliB(periodicFaceB->center(),periodicFaceB->outNormal(),P0,faceAshift);
                            
                            if(pliA.type==PlaneLineIntersection<3>::INCIDENT && pliB.type==PlaneLineIntersection<3>::INCIDENT)
                            {
                                const VectorDimD AB(pliB.P-pliA.P);
                                if((AB-faceAshift).norm()<FLT_EPSILON)
                                {

                                    GlidePlaneKey<3> parallelGlidePlaneKey(planeIndex+dipoleHeight, slipSystem.n);
                                    std::shared_ptr<PeriodicGlidePlane<3>> parallelglidePlane(mg.ddBase.periodicGlidePlaneFactory.get(parallelGlidePlaneKey));

                                    GlidePlaneKey<3> prismaticPlaneKey(P0, grain.singleCrystal->reciprocalLatticeDirection(glidePlane->referencePlane->unitNormal.cross(AB)));
                                    std::shared_ptr<PeriodicGlidePlane<3>> prismaticGlidePlane(mg.ddBase.periodicGlidePlaneFactory.get(prismaticPlaneKey));

                                    if(parallelglidePlane && prismaticGlidePlane)
                                    {

                                        
                                        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,3,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
                                        // The prismatic loop
                                        const int nShift(2);
                                        const VectorDimD lShift(nShift*AB);
                                        const VectorDimD rShift(lShift+AB);
                                        const VectorDimD startPoint(0.5*(pliB.P+pliA.P));
                                        std::vector<VectorDimD> prismaticNodePos;
                                        prismaticNodePos.push_back(startPoint+lShift); //HERE ADD N*AB
                                        prismaticNodePos.push_back(startPoint+rShift);
                                        prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(startPoint+rShift));
                                        prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(startPoint+lShift));
                                        
                                        const VectorDimD lineP0 = 0.5*(startPoint+lShift+parallelglidePlane->referencePlane->snapToPlane(startPoint+lShift));
                                        const VectorDimD lineP1 = 0.5*(startPoint+rShift+parallelglidePlane->referencePlane->snapToPlane(startPoint+rShift));
                                        
                                        mg.insertJunctionLoop(prismaticNodePos,prismaticGlidePlane,
                                                           slipSystem.s.cartesian(),prismaticGlidePlane->referencePlane->unitNormal,
                                                           P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
                                        
                                        
                                        
                                        const int halfNodeNumber(dipoleNodes/2);
                                        if(halfNodeNumber < 0)
                                        {
                                            throw std::runtime_error("dipoleNodes="+std::to_string(halfNodeNumber)+"is smaller than 0");
                                        }
                                        const double halfDipoleLength(AB.norm()*0.5);
                                        const VectorDimD dipoleDir(AB/AB.norm());
                                        const VectorDimD shiftNodeLength(halfDipoleLength/(1.0*(halfNodeNumber+1))*dipoleDir);
                                        
                                        // First glide loop
    //                                    const double glideStep=50.0;
                                        std::vector<VectorDimD> firstNodePos;
                                        firstNodePos.push_back(startPoint+lShift);
                                        firstNodePos.push_back(startPoint+rShift);
                                        if((rShift-lShift).cross(glideStep*prismaticGlidePlane->referencePlane->unitNormal).dot(glidePlane->referencePlane->unitNormal)>0.0)
                                        {// loop is right-handed to glidePlane->referencePlane->unitNormal
                                            glideStep*=-1;
                                        }
                                        firstNodePos.push_back(startPoint+rShift+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                        for(int k=1; k<=halfNodeNumber; k++)
                                        { // Add nodes on rigth arm
                                            firstNodePos.push_back(startPoint+rShift-k*shiftNodeLength+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                        }
                                        for(int k=halfNodeNumber; k>0; k--)
                                        { // Add nodes on left arm
                                            firstNodePos.push_back(startPoint+lShift+k*shiftNodeLength+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                        }
                                        firstNodePos.push_back(startPoint+lShift+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                                                                
                                        mg.insertJunctionLoop(firstNodePos,glidePlane,
                                                           -slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                                                           P0,grainID,DislocationLoopIO<3>::GLISSILELOOP);


                                        FiniteLineSegment<3> mirrowLine(lineP0,lineP1);
                                        
                                        // Second glide loop
                                        std::vector<VectorDimD> secondNodePos;
                                        for(const auto& pos : firstNodePos)
                                        {
                                            VectorDimD c=mirrowLine.snapToInfiniteLine(pos);
                                            secondNodePos.push_back(2.0*c-pos);
    //                                        secondNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(pos));
                                        }
                                        
                                        
                                        mg.insertJunctionLoop(secondNodePos,parallelglidePlane,
                                                           slipSystem.s.cartesian(),parallelglidePlane->referencePlane->unitNormal,
                                                           parallelglidePlane->referencePlane->snapToPlane(P0),grainID,DislocationLoopIO<3>::GLISSILELOOP);
                                    }
                                    else
                                    {
                                        std::cout<<"Cannot create glide planes"<<std::endl;
                                    }
                                }
                                else
                                {
                                    throw std::runtime_error("periodic line intersection with faces is not the face shift vector");
                                }
                            }
                            else
                            {
                                throw std::runtime_error("periodic line direction does not form an incident intersecitn with periodic faces");
                            }
                        }
                        else
                        {
                            throw std::runtime_error("planeNormal of slipSystem "+std::to_string(rSS)+" is not othogonal to facelatticeShift "+std::to_string(exitFaceID));
                        }
                    }
                    else
                    {
                        throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" is not a periodic face.");
                    }
                }
                else
                {
                    throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" not found.");
                }
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
        else
        {
            std::cout<<"Skipping slip system "<<rSS<<std::endl;
        }
        
    }

}
#endif
