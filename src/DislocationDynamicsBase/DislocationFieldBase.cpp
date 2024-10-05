/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationFieldBase_cpp_
#define _model_DislocationFieldBase_cpp_

#include <DislocationFieldBase.h>
#include <TextFileParser.h>

namespace model
{
    
    
        
        template<int dim>
        void DislocationFieldBase<dim>::initFromFile(const std::string& fileName)
        {
            
            a=TextFileParser(fileName).readScalar<double>("coreSize",false);
            assert(a>0.0 && "coreSize MUST BE > 0.");
            a2=a*a;
        }

    
template struct DislocationFieldBase<3>;
    //! Dislocation core size
//    template<int dim>
//    double DislocationFieldBase<dim>::a=1.0;
//    
//    //! Dislocation core size squared
//    template<int dim>
//    double DislocationFieldBase<dim>::a2=1.0;
//    
//    //! Identity matrix
//    template<int dim>
//    const Eigen::Matrix<double,dim,dim> DislocationFieldBase<dim>::I=Eigen::Matrix<double,dim,dim>::Identity();
    
}	// close namespace
#endif



