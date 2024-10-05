/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureSpecificationBase_cpp_
#define model_MicrostructureSpecificationBase_cpp_


#include <string>
#include <MicrostructureSpecificationBase.h>

namespace model
{
    MicrostructureSpecificationBase::MicrostructureSpecificationBase(const std::string& type_in,const std::string& style_in):
    /* init */ type(type_in)
    /* init */,style(style_in)
    /* init */,tag(type+style+std::to_string(this->sID))
    /* init */,parser(nullptr)
    {
        
    }

    MicrostructureSpecificationBase::MicrostructureSpecificationBase(const std::string& type_in,const std::string& style_in,const std::string& fileName):
    /* init */ type(type_in)
    /* init */,style(style_in)
    /* init */,tag(type+style+std::to_string(this->sID))
    /* init */,parser(new TextFileParser(fileName))
    {
        
    }
}
#endif
