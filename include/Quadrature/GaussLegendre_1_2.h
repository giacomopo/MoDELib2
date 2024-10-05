 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_2_H_ 
#define model_GAUSSLEGENDRE_1_2_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,2>
   {
       static Eigen::Matrix<double,2,2> abcsissasAndWeights()
       {
           Eigen::Matrix<double,2,2> aw;
           aw<<2.113248654051872e-01, 4.999999999999999e-01, 
               7.886751345948129e-01, 4.999999999999999e-01; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

