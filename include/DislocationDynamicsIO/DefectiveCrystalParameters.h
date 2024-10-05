/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalParameters_H_
#define model_DefectiveCrystalParameters_H_

#include <string>
#include <set>

#include <Eigen/Dense>

#include <TextFileParser.h>
#include <IDreader.h>
#include <DDtraitsIO.h>

namespace model
{
    
    struct DefectiveCrystalParameters
    {
        
        const DDtraitsIO traitsIO;
        const int useFEM;
        const bool useElasticDeformation;
        const bool useDislocations;
        const bool useClusterDynamics;
        const bool useCracks;
        const bool useInclusions;
        const long int Nsteps;
        const int useSubCycling;
        const std::set<int> subcyclingBins; 
        const bool use_stochasticForce;
        const std::set<int> periodicFaceIDs;
        const double dtMax;
        const int outputFrequency;
        const bool outputBinary;
        long int runID;
        double totalTime;
        double dt;
        
        void manageRestart();

    private:

        static std::set<int> getSubCyclingSet(const std::vector<int> &inpVector);

    public:
        
        DefectiveCrystalParameters(const std::string& folderName) ;
        
    };
}
#endif
