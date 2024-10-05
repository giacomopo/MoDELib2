/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Plane_H_
#define model_Plane_H_

#include <cfloat>
#include <tuple>
//#include <map>
#include <Eigen/Dense>
#include <iostream>


namespace model
{
    
    template <int dim>
    struct Plane
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        const VectorDim P;
        const VectorDim unitNormal;
        const MatrixDim L2G;

        /**********************************************************************/
        Plane(const VectorDim& p,const VectorDim& n);
        
        /**********************************************************************/
        bool contains(const VectorDim& P0) const;

        /**********************************************************************/
        bool isAbove(const VectorDim& P0) const;

        /**********************************************************************/
        bool isBelow(const VectorDim& P0) const;

        /**********************************************************************/
        VectorDim snapToPlane(const VectorDim& P0) const;
        
        /**********************************************************************/
        double distanceTo(const VectorDim& P0) const;
        
        /**********************************************************************/
        VectorLowerDim localPosition(const VectorDim& point) const;
        
        /**********************************************************************/
        VectorDim globalPosition(const VectorLowerDim& point) const;
        
        VectorLowerDim localDirection(const VectorDim& dir) const;
        
        /**********************************************************************/
        VectorDim globalDirection(const VectorLowerDim& dir) const;
        
        /**********************************************************************/
        static MatrixDim getL2G(VectorDim z);

    };
    
}
#endif
