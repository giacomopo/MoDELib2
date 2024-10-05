/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureSpecificationBase_H_
#define model_MicrostructureSpecificationBase_H_

#include <memory>
#include <string>
#include <StaticID.h>
#include <TextFileParser.h>



namespace model
{
    struct MicrostructureSpecificationBase : public StaticID<MicrostructureSpecificationBase>
    {
        const std::string type;
        const std::string style;
        const std::string tag;
        const std::unique_ptr<TextFileParser> parser;
        
        MicrostructureSpecificationBase(const std::string& type_in,const std::string& style_in);
        MicrostructureSpecificationBase(const std::string& type_in,const std::string& style_in,const std::string& fileName);
        
    };
}
#endif
