/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopGenerator_cpp_
#define model_ShearLoopGenerator_cpp_

#include <iostream>
#include <ShearLoopGenerator.h>

namespace model
{

//    ShearLoopGenerator::ShearLoopGenerator(const std::string& fileName) :
//    /* init */ MicrostructureGeneratorBase(fileName)
//    {
//        
//    }

ShearLoopGenerator::ShearLoopGenerator(const ShearLoopDensitySpecification& spec,MicrostructureGenerator& mg)
{
    //    const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
    std::cout<<magentaBoldColor<<"Generating shear loop density"<<defaultColor<<std::endl;
    double density=0.0;
    std::cout<<"shear loop density="<<density<<std::endl;
    if(spec.targetDensity>0.0)
    {
        //        const int numberOfSides(this->parser.readScalar<int>("numberOfSides",true));
        //        const double radiusDistributionMean(this->parser.readScalar<double>("radiusDistributionMean",true));
        //        const double radiusDistributionStd(this->parser.readScalar<double>("radiusDistributionStd",true));
        std::normal_distribution<double> radiusDistribution(spec.radiusDistributionMean/mg.ddBase.poly.b_SI,spec.radiusDistributionStd/mg.ddBase.poly.b_SI);
        std::mt19937 generator;
        while(density<spec.targetDensity)
        {
            const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<3> L0=rp.first;
            const size_t grainID=rp.second;
            std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
            const int rSS(ssDist(generator)); // a random SlipSystem
            const double radius(radiusDistribution(generator));
            try
            {
                const bool success=generateSingle(mg,rSS,L0.cartesian(),radius,spec.numberOfSides);
                if(success)
                {
                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"shear loop density="<<density<<std::endl;
                }
            }
            catch(const std::exception& e)
            {
                
            }
        }
    }
}

ShearLoopGenerator::ShearLoopGenerator(const ShearLoopIndividualSpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual shear loops"<<defaultColor<<std::endl;
    
    if(spec.slipSystemIDs.size()!=spec.loopRadii.size())
    {
        throw std::runtime_error("loopSlipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO loopRadii.size()="+std::to_string(spec.loopRadii.size()));
    }
    if(int(spec.slipSystemIDs.size())!=spec.loopCenters.rows())
    {
        throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
    }
    if(spec.slipSystemIDs.size()!=spec.loopSides.size())
    {
        throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO loopSides.size()="+std::to_string(spec.loopSides.size()));
    }
    
    double density=0.0;
    for(size_t k=0;k<spec.slipSystemIDs.size();++k)
    {
        try
        {
            const bool success=generateSingle(mg,spec.slipSystemIDs[k],spec.loopCenters.row(k),spec.loopRadii[k]/mg.ddBase.poly.b_SI,spec.loopSides[k]);
            if(success)
            {
                density+=2.0*std::numbers::pi*spec.loopRadii[k]/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                std::cout<<"shear loop density="<<density<<std::endl;
            }
        }
        catch(const std::exception& e)
        {
            
        }
    }
}

bool ShearLoopGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& center,const double& radius,const size_t& sides)
{
    bool success =false;
    std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(center));
    if(!found.first)
    {
        throw std::runtime_error("ShearLoopGenerator::generateSingle:: point is outside mesh.");
    }
    
    const int grainID(found.second->region->regionID);
    //        std::cout<<"grainID="<<grainID<<std::endl;
    assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
    const auto& grain(mg.ddBase.poly.grain(grainID));
    
    if(rSS>=0 && rSS<int(grain.singleCrystal->slipSystems().size()))
    {
        const auto& slipSystem(*grain.singleCrystal->slipSystems()[rSS]);
        
        
        const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(center));
        GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
        const VectorDimD P0(glidePlane->referencePlane->snapToPlane(center));
        //            std::cout<<"P0="<<P0.transpose()<<std::endl;
        
        std::vector<VectorDimD> loopNodePos;
        for(size_t k=0;k< sides;++k)
        {
            loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*std::numbers::pi/sides,slipSystem.unitNormal)*slipSystem.s.cartesian().normalized()*radius);
        }
        
        
        //            std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
        
        mg.insertJunctionLoop(loopNodePos,glidePlane,
                              slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                              P0,grainID,DislocationLoopIO<3>::GLISSILELOOP);
        success=true;
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
    return success;
}

//    void ShearLoopGenerator::generateDensity(MicrostructureGenerator& mg)
//    {
//        std::cout<<magentaBoldColor<<"Generating shear loop density"<<defaultColor<<std::endl;
//        const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
//        if(targetDensity>0.0)
//        {
//            const int numberOfSides(this->parser.readScalar<int>("numberOfSides",true));
//            const double radiusDistributionMean(this->parser.readScalar<double>("radiusDistributionMean",true));
//            const double radiusDistributionStd(this->parser.readScalar<double>("radiusDistributionStd",true));
//            std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
//            std::mt19937 generator;
//            double density=0.0;
//            while(density<targetDensity)
//            {
//                const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
//                const LatticeVector<3> L0=rp.first;
//                const size_t grainID=rp.second;
//                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
//                const int rSS(ssDist(generator)); // a random SlipSystem
//                const double radius(radiusDistribution(generator));
//                try
//                {
//                    generateSingle(mg,rSS,L0.cartesian(),radius,numberOfSides);
//                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
//                    std::cout<<"shear loop density="<<density<<std::endl;
//                }
//                catch(const std::exception& e)
//                {
//                    
//                }
//            }
//        }
//    }

//    void ShearLoopGenerator::generateIndividual(MicrostructureGenerator& mg)
//    {
//        const std::vector<int> periodicLoopSlipSystemIDs(this->parser.readArray<int>("periodicLoopSlipSystemIDs",true));
//        
//        if(periodicLoopSlipSystemIDs.size())
//        {
//            std::cout<<magentaBoldColor<<"Generating individual periodic loops"<<defaultColor<<std::endl;
//            const std::vector<double> periodicLoopRadii(this->parser.readArray<double>("periodicLoopRadii_SI",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,dim> periodicLoopCenters(this->parser.readMatrix<double>("periodicLoopCenters",periodicLoopSlipSystemIDs.size(),dim,true));
//            const std::vector<int> periodicLoopSides(this->parser.readArray<int>("periodicLoopSides",true));
//            
//            if(periodicLoopSlipSystemIDs.size()!=periodicLoopRadii.size())
//            {
//                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopRadii.size()="+std::to_string(periodicLoopRadii.size()));
//            }
//            if(int(periodicLoopSlipSystemIDs.size())!=periodicLoopCenters.rows())
//            {
//                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopCenters.rows()="+std::to_string(periodicLoopCenters.rows()));
//            }
//            if(periodicLoopSlipSystemIDs.size()!=periodicLoopSides.size())
//            {
//                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopSides.size()="+std::to_string(periodicLoopSides.size()));
//            }
//            
//            for(size_t k=0;k<periodicLoopSlipSystemIDs.size();++k)
//            {
//                generateSingle(mg,periodicLoopSlipSystemIDs[k],periodicLoopCenters.row(k),periodicLoopRadii[k]/mg.ddBase.poly.b_SI,periodicLoopSides[k]);
//            }
//        }
//        
//    }

}
#endif
