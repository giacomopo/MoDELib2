/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StraightDislocationSegment_H_
#define model_StraightDislocationSegment_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
#include <PolycrystallineMaterialBase.h>
#include <DislocationFieldBase.h>

namespace model
{
    
    
    template <int dim>
    class StraightDislocationSegment
    {
        typedef double Scalar;
        typedef Eigen::Matrix<Scalar,dim,dim> MatrixDim;
        typedef Eigen::Matrix<Scalar,dim,1>   VectorDim;
        
        MatrixDim nonSymmStress_kernel(const VectorDim& r) const;
        VectorDim displacement_kernel(const VectorDim& r) const;
        double elasticInteractionEnergy_kernel(const VectorDim& z,const VectorDim& tA,const VectorDim& bA) const;
        
    public:
        
        const PolycrystallineMaterialBase& material;
        const VectorDim& P0;
        const VectorDim& P1;
        const VectorDim& b;
        const double& length;
        const VectorDim& t;
        const double& EwaldLength;
        
    private:

        VectorDim bCt;
        
    public:

        StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
                                   const VectorDim& _P0,
                                   const VectorDim& _P1,
                                   const VectorDim& _b,
                                   const double& _length,
                                   const VectorDim& _t,
                                   const double& EwaldLength_in);

        void updateGeometry();
        MatrixDim nonSymmStress(const VectorDim& x) const;
        MatrixDim stress(const VectorDim& x) const;
        VectorDim displacement(const VectorDim& x) const;        
        double elasticInteractionEnergy(const VectorDim& xA,const VectorDim& tA,const VectorDim& bA) const;
        
	};	
	
}
#endif

