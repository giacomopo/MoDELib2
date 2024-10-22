 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_5_H_ 
#define model_GAUSSLEGENDRE_1_5_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,5>
   {
       static Eigen::Matrix<double,2,5> abcsissasAndWeights()
       {
           Eigen::Matrix<double,5,2> aw;
           aw<<4.691007703066818e-02, 1.184634425280944e-01, 
               2.307653449471585e-01, 2.393143352496831e-01, 
               5.000000000000000e-01, 2.844444444444446e-01, 
               7.692346550528416e-01, 2.393143352496835e-01, 
               9.530899229693320e-01, 1.184634425280945e-01; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

