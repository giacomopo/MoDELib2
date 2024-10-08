 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_14_H_ 
#define model_GAUSSLEGENDRE_1_14_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,14>
   {
       static Eigen::Matrix<double,2,14> abcsissasAndWeights()
       {
           Eigen::Matrix<double,14,2> aw;
           aw<<6.858095651593343e-03, 1.755973016587635e-02, 
               3.578255816821341e-02, 4.007904357988003e-02, 
               8.639934246511755e-02, 6.075928534395118e-02, 
               1.563535475941573e-01, 7.860158357909684e-02, 
               2.423756818209231e-01, 9.276919873896886e-02, 
               3.404438155360551e-01, 1.025992318606479e-01, 
               4.459725256463282e-01, 1.076319267315790e-01, 
               5.540274743536719e-01, 1.076319267315792e-01, 
               6.595561844639449e-01, 1.025992318606478e-01, 
               7.576243181790772e-01, 9.276919873896902e-02, 
               8.436464524058430e-01, 7.860158357909743e-02, 
               9.136006575348821e-01, 6.075928534395080e-02, 
               9.642174418317864e-01, 4.007904357988041e-02, 
               9.931419043484058e-01, 1.755973016587551e-02; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

