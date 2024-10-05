/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsGenerator_cpp_
#define model_FrankLoopsGenerator_cpp_

#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

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
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <FrankLoopsGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

//FrankLoopsGenerator::FrankLoopsGenerator(const std::string& fileName) :
//    /* init */ MicrostructureGeneratorBase(fileName)
//    {
//        
//    }

    FrankLoopsGenerator::FrankLoopsGenerator(const FrankLoopsDensitySpecification& spec,MicrostructureGenerator& mg)
    {
        std::cout<<magentaBoldColor<<"Generating Frank loop density"<<defaultColor<<std::endl;
//        const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
        if(spec.targetDensity>0.0)
        {
//            const int numberOfSides(this->parser.readScalar<int>("numberOfSides",true));
//            const double radiusDistributionMean(this->parser.readScalar<double>("radiusDistributionMean",true));
//            const double radiusDistributionStd(this->parser.readScalar<double>("radiusDistributionStd",true));
//            const int areVacancyLoops(this->parser.readScalar<int>("areVacancyLoops",true));

            std::normal_distribution<double> radiusDistribution(spec.radiusDistributionMean/mg.ddBase.poly.b_SI,spec.radiusDistributionStd/mg.ddBase.poly.b_SI);
            std::mt19937 generator;
            double density=0.0;
            while(density<spec.targetDensity)
            {
                const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
                const LatticeVector<3> L0=rp.first;
                const size_t grainID=rp.second;
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->planeNormals().size()-1);
                const int pID(ssDist(generator)); // a random SlipSystem
                const double radius(radiusDistribution(generator));
                try
                {
                    const bool isVL(spec.areVacancyLoops>0);
                    generateSingle(mg,pID,L0.cartesian(),radius,spec.numberOfSides,isVL);
                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"Frank loop density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    
                }
            }
        }
    }

FrankLoopsGenerator::FrankLoopsGenerator(const FrankLoopsIndividualSpecification& spec,MicrostructureGenerator& mg)
    {
//        const std::vector<int> planeIDs(this->parser.readArray<int>("planeIDs",true));
        std::cout<<magentaBoldColor<<"Generating individual Frank loops"<<defaultColor<<std::endl;
        if(spec.planeIDs.size())
        {
//            const std::vector<double> loopRadii(this->parser.readArray<double>("loopRadii_SI",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters(this->parser.readMatrix<double>("loopCenters",spec.planeIDs.size(),dim,true));
//            const std::vector<int> numberOfSides(this->parser.readArray<int>("numberOfSides",true));
//            const std::vector<int> isVacancyLoop(this->parser.readArray<int>("isVacancyLoop",true));

            if(spec.planeIDs.size()!=spec.loopRadii.size())
            {
                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.spec.spec.loopRadii.size()="+std::to_string(spec.loopRadii.size()));
            }
            if(int(spec.planeIDs.size())!=spec.loopCenters.rows())
            {
                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
            }
            if(spec.planeIDs.size()!=spec.loopSides.size())
            {
                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopSides.size()="+std::to_string(spec.loopSides.size()));
            }
            if(spec.planeIDs.size()!=spec.isVacancyLoop.size())
            {
                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.isVacancyLoop.size()="+std::to_string(spec.isVacancyLoop.size()));
            }
            for(size_t k=0;k<spec.planeIDs.size();++k)
            {
                
                const bool isVL(spec.isVacancyLoop[k]>0);
                generateSingle(mg,spec.planeIDs[k],spec.loopCenters.row(k),spec.loopRadii[k]/mg.ddBase.poly.b_SI,spec.loopSides[k],isVL);
            }
        }
        
    }

    void FrankLoopsGenerator::generateSingle(MicrostructureGenerator& mg,const int& pID,const VectorDimD& center,const double& radius,const size_t& sides,const bool& isVacancyLoop)
    {
        std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(center));
        if(!found.first)
        {
            std::cout<<"Point "<<center.transpose()<<" is outside mesh. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        
        const int grainID(found.second->region->regionID);
        assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        if(pID>=0 && pID<int(grain.singleCrystal->planeNormals().size()))
        {
            const auto& planeBase(*grain.singleCrystal->planeNormals()[pID]);
            
            
            const long int planeIndex(planeBase.closestPlaneIndexOfPoint(center));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, planeBase);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            const VectorDimD P0(glidePlane->referencePlane->snapToPlane(center));
            
            std::vector<VectorDimD> loopNodePos;
            for(size_t k=0;k< sides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*std::numbers::pi/sides,planeBase.G2L.row(2).transpose())*planeBase.G2L.row(0).transpose()*radius);
            }
            
            
//            std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
            ReciprocalLatticeDirection<3> rp(planeBase);
            const VectorDimD b(isVacancyLoop? (rp.cartesian().normalized()*rp.planeSpacing()).eval() : (-1.0*rp.cartesian().normalized()*rp.planeSpacing()).eval());
            mg.insertJunctionLoop(loopNodePos,glidePlane,
                               b,glidePlane->referencePlane->unitNormal,
                                  P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
        }
        else
        {
            if(pID<0)
            {
                std::cout<<"Skipping planeID "<<pID<<std::endl;
            }
            else
            {
                throw std::runtime_error("planeID "+std::to_string(pID)+" not found, skipping.");
            }
        }
    }

//    void FrankLoopsGenerator::generateDensity(MicrostructureGenerator& mg)
//    {
//        std::cout<<magentaBoldColor<<"Generating Frank loop density"<<defaultColor<<std::endl;
//        const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
//        if(targetDensity>0.0)
//        {
//            const int numberOfSides(this->parser.readScalar<int>("numberOfSides",true));
//            const double radiusDistributionMean(this->parser.readScalar<double>("radiusDistributionMean",true));
//            const double radiusDistributionStd(this->parser.readScalar<double>("radiusDistributionStd",true));
//            const int areVacancyLoops(this->parser.readScalar<int>("areVacancyLoops",true));
//
//            std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
//            std::mt19937 generator;
//            double density=0.0;
//            while(density<targetDensity)
//            {
//                const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
//                const LatticeVector<3> L0=rp.first;
//                const size_t grainID=rp.second;
//                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->planeNormals().size()-1);
//                const int pID(ssDist(generator)); // a random SlipSystem
//                const double radius(radiusDistribution(generator));
//                try
//                {
//                    const bool isVL(areVacancyLoops>0);
//                    generateSingle(mg,pID,L0.cartesian(),radius,numberOfSides,isVL);
//                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
//                    std::cout<<"Frank loop density="<<density<<std::endl;
//                }
//                catch(const std::exception& e)
//                {
//                    
//                }
//            }
//        }
//    }

//    void FrankLoopsGenerator::generateIndividual(MicrostructureGenerator& mg)
//    {
//        const std::vector<int> planeIDs(this->parser.readArray<int>("planeIDs",true));
//        
//        if(spec.planeIDs.size())
//        {
//            std::cout<<magentaBoldColor<<"Generating individual Frank loops"<<defaultColor<<std::endl;
//            const std::vector<double> loopRadii(this->parser.readArray<double>("loopRadii_SI",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters(this->parser.readMatrix<double>("loopCenters",spec.planeIDs.size(),dim,true));
//            const std::vector<int> numberOfSides(this->parser.readArray<int>("numberOfSides",true));
//            const std::vector<int> isVacancyLoop(this->parser.readArray<int>("isVacancyLoop",true));
//
//            if(spec.planeIDs.size()!=spec.spec.spec.loopRadii.size())
//            {
//                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.spec.spec.loopRadii.size()="+std::to_string(spec.spec.spec.loopRadii.size()));
//            }
//            if(int(spec.planeIDs.size())!=spec.loopCenters.rows())
//            {
//                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
//            }
//            if(spec.planeIDs.size()!=spec.loopSides.size())
//            {
//                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopSides.size()="+std::to_string(spec.loopSides.size()));
//            }
//            if(spec.planeIDs.size()!=spec.isVacancyLoop.size())
//            {
//                throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.isVacancyLoop.size()="+std::to_string(spec.isVacancyLoop.size()));
//            }
//            for(size_t k=0;k<spec.planeIDs.size();++k)
//            {
//                
//                const bool isVL(isVacancyLoop[k]>0);
//                generateSingle(mg,planeIDs[k],spec.loopCenters.row(k),loopRadii[k]/mg.ddBase.poly.b_SI,numberOfSides[k],isVL);
//            }
//        }
//        
//    }

}
#endif
