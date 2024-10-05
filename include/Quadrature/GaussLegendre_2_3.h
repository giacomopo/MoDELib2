 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_2_3_H_
#define model_GAUSSLEGENDRE_2_3_H_

namespace model{

	template <>
	struct GaussLegendre<2,3>
    {
		static Eigen::Matrix<double,3,3> abcsissasAndWeights()
        {
			Eigen::Matrix<double,3,3> U;
			U(0,0)=  1.666666666666667e-01;		U(2,0)= 1.666666666666667e-01;
			U(1,0)=  1.666666666666667e-01;
			U(0,1)=	 6.666666666666667e-01;		U(2,1)= 1.666666666666667e-01;
			U(1,1)=  1.666666666666667e-01;
			U(0,2)=  1.666666666666667e-01;		U(2,2)= 1.666666666666667e-01;
			U(1,2)=  6.666666666666667e-01;
			return U;
		}
	};
}
#endif 

