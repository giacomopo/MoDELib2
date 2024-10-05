/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StressStraight_H_
#define model_StressStraight_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
#include <PolycrystallineMaterialBase.h>
#include <DislocationFieldBase.h>
#include <ClusterDynamicsParameters.h>

namespace model
{

    template <int dim,typename Scalar=double>
    class StressStraight
    {
    public:
        typedef Eigen::Matrix<Scalar,dim,dim> MatrixDim;
        typedef Eigen::Matrix<Scalar,dim,1>   VectorDim;
//        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
//        typedef Eigen::Matrix<double,mSize,2>   ConcentrationMatrixType;
//        typedef Eigen::Matrix<double,mSize,1>   ConcentrationVectorType;

    private:
        MatrixDim nonSymmStress_kernel(const VectorDim& r) const;
        VectorDim displacement_kernel(const VectorDim& r) const;

    public:
        
        const PolycrystallineMaterialBase& material;
        const VectorDim P0;
        const VectorDim P1;
        const VectorDim b;
        const VectorDim chord;
        const double length;
        const VectorDim t;
        const double EwaldLength;
        const VectorDim bCt;
        
        StressStraight(const PolycrystallineMaterialBase& material_in,const VectorDim& _P0,const VectorDim& _P1, const VectorDim& _b,
                       const double& EwaldLength_in);
        MatrixDim nonSymmStress(const VectorDim& x) const;
        MatrixDim stress(const VectorDim& x) const;
        VectorDim displacement(const VectorDim& x) const;
//        ConcentrationMatrixType concentrationMatrices(const VectorDim& x, const size_t& grainID, const VectorDim& sourceDir, const VectorDim& sinkDir, const ClusterDynamicsParameters<dim>& icp) const;
//        ConcentrationVectorType clusterConcentration(const VectorDim& x, const size_t& grainID, const VectorDim& sourceDir, const Eigen::Array<double,1,mSize>& sourceVScalar, const VectorDim& sinkDir, const Eigen::Array<double,1,mSize>& sinkVScalar, const ClusterDynamicsParameters<dim>& icp) const;
        
	};
	
}
#endif
