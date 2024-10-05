/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalParameters_cpp_
#define model_DefectiveCrystalParameters_cpp_

#include <cfloat>
#include <DefectiveCrystalParameters.h>

namespace model
{
    

        
                /**********************************************************************/
        void DefectiveCrystalParameters::manageRestart()
        {
            // Menage restart
            IDreader<1,200,double> vReader(traitsIO.fFolder+"/F");
            vReader.readLabelsFile(traitsIO.fFolder+"/F_labels.txt");
            if (vReader.isGood(0,true))
            {// F/F_0.txt exists
                vReader.read(0,true);
                if(runID<0)
                {// Restart from last available step
                    if(vReader.size())
                    {// at least a line is available
//                        vReader.readLabelsFile("F/F_labels.txt");
                        runID=vReader.rbegin()->first;
                        totalTime=vReader.last("time [b/cs]");
                        dt=vReader.last("dt [b/cs]");
                    }
                    else
                    {// file is empty, keep default initialization
                        runID=0;
                    }
                }
                else
                {// Restart from specific runID
//                    const auto iter=vReader.find(runID);
                    totalTime=vReader(runID,"time [b/cs]");
                    dt=vReader(runID,"dt [b/cs]");
                }
            }
            else
            {// F/F_0.txt is not there, keep default initialization
                std::cout<<"Unable to read F/F_0.txt"<<std::endl;
                runID=0;
            }
            
            std::cout<<"starting at time step "<<runID<<std::endl;
            std::cout<<"totalTime= "<<totalTime<<std::endl;
            std::cout<<"dt= "<<dt<<std::endl;
        }

        /**********************************************************************/
        std::set<int> DefectiveCrystalParameters::getSubCyclingSet(const std::vector<int> &inpVector)
        {
            std::set<int> temp;
            for (const auto &iv : inpVector)
            {
                temp.insert(iv);
            }
            return temp;
        }

        /**********************************************************************/
        DefectiveCrystalParameters::DefectiveCrystalParameters(const std::string& folderName) :
        /* init */ traitsIO(folderName)
        /* init */,useFEM(TextFileParser(traitsIO.ddFile).readScalar<int>("useFEM",true))
        /* init */,useElasticDeformation(TextFileParser(traitsIO.ddFile).readScalar<int>("useElasticDeformation",true))
        /* init */,useDislocations(TextFileParser(traitsIO.ddFile).readScalar<int>("useDislocations",true))
        /* init */,useClusterDynamics(TextFileParser(traitsIO.ddFile).readScalar<int>("useClusterDynamics",true))
        /* init */,useCracks(TextFileParser(traitsIO.ddFile).readScalar<int>("useCracks",true))
        /* init */,useInclusions(TextFileParser(traitsIO.ddFile).readScalar<int>("useInclusions",true))
        /* init */,Nsteps(TextFileParser(traitsIO.ddFile).readScalar<size_t>("Nsteps",true))
        /* init */,useSubCycling(TextFileParser(traitsIO.ddFile).readScalar<int>("useSubCycling",true))
        /* init */,subcyclingBins(getSubCyclingSet(TextFileParser(traitsIO.ddFile).readArray<int>("subcyclingBins",true)))
        /* init */,use_stochasticForce(TextFileParser(traitsIO.ddFile).readScalar<int>("use_stochasticForce",true))
        /* init */,periodicFaceIDs(TextFileParser(traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true))
        /* init */,dtMax(TextFileParser(traitsIO.ddFile).template readScalar<double>("dtMax",true))
        /* init */,outputFrequency(TextFileParser(traitsIO.ddFile).readScalar<int>("outputFrequency",true))
        /* init */,outputBinary(TextFileParser(traitsIO.ddFile).readScalar<int>("outputBinary",true))
        /* init */,runID(TextFileParser(traitsIO.ddFile).readScalar<long int>("startAtTimeStep",true))
        /* init */,totalTime(0.0)
        /* init */,dt(10.0)
        {
            
            if (Nsteps < 0)
            {
                throw std::runtime_error("Nsteps MUST BE >= 0.");
            }

            if (dtMax < FLT_EPSILON)
            {
                throw std::runtime_error("dtMax must be > FLT_EPSILON.");
            }
            
            manageRestart();

            
        }
        
        
  }
#endif
