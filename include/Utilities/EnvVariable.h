/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ENVVARIABLE_H_
#define model_ENVVARIABLE_H_

#include <cstdlib>
#include <string.h>
#include <assert.h>

namespace model {

	/************************************************************/	
	/************************************************************/
	/*! \brief A simple class used to extract environment variables
	 *
	 * Sample usage:
	 * \code
	 * #include <iostream>
	 * #include <EnvVariable.h>
	 * 
	 * int main(){
	 *  model::EnvVariable ev("PATH");
	 *  std::cout<<"PATH is "<<ev<<std::endl;
	 *  return 0;
	 * }
	 * \endcode
	 */

	struct EnvVariable : public std::string {
	
	/******************************/
	static std::string getVariable(const std::string& varname){
  		const char *tmpvar = std::getenv (varname.c_str());
  		assert(tmpvar!=0 && *tmpvar!='\0' && "FAILED TO GET ENVIRONMENT VARIABLE.");
		return tmpvar;
	}

	
	/******************************/
	EnvVariable(const std::string& varname) : std::string(getVariable(varname)) {}

	};

	/************************************************************/	
	/************************************************************/
} // namespace model
#endif
