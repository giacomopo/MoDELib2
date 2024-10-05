/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilitySelector_h_
#define _model_DislocationMobilitySelector_h_

#include <DislocationMobilityBase.h>
#include <DislocationMobilityFCC.h>
#include <DislocationMobilityBCC.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXpyramidal.h>
#include <DislocationMobilityPy.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilitySelector
    {
        
        const std::string defaultStr;
        
        DislocationMobilitySelector(const std::string& defaultStr_in);
                
        std::shared_ptr<DislocationMobilityBase> getMobility(const std::string& dislocationMobilityType,
                                                             const PolycrystallineMaterialBase& material) const;
        
    };
    
    
}
#endif
