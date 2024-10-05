/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VoigtTraits_H_
#define model_VoigtTraits_H_


#include <Eigen/Dense>

namespace model
{
    
    
    template <int dim>
    struct SymmetricVoigtTraits
    {
        static constexpr int voigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<size_t,voigtSize,2> VoigtSizeMatrixType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double,voigtSize,1> VectorVoigt;

        
        const VoigtSizeMatrixType tensorIndex; // i=tensorIndex(k,0), j=tensorIndex(k,1)
//        const Eigen::Matrix<size_t,dim,dim> voigtIndex; // k=voigtIndex(i,j)

        SymmetricVoigtTraits(const VoigtSizeMatrixType& voigtOrder_in);
        
        MatrixDim v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector, const bool& is_strain) const ;
        VectorVoigt m2v(const MatrixDim& input_matrix, const bool& is_strain) const;
        
        
//        static Eigen::Matrix<size_t,dim,dim> getVoigtIndex(const VoigtSizeMatrixType& ti);
        
    };

template <int rows,int cols>
struct VoigtTraits
{
    static constexpr int voigtSize=rows*cols;
    typedef Eigen::Matrix<size_t,voigtSize,2> VoigtStorageType;
    typedef Eigen::Matrix<double,rows,cols>   TensorType;
    typedef Eigen::Matrix<double,voigtSize,1> VoigtVectorType;

    
    const VoigtStorageType tensorIndex; // i=tensorIndex(k,0), j=tensorIndex(k,1)

    VoigtTraits(const VoigtStorageType& voigtOrder_in) :
    /* init */ tensorIndex(voigtOrder_in)
    {
        
    }
    
    TensorType v2m(const VoigtVectorType& voigtvector) const
    {
        //from voigt to matrix format
        TensorType temp(TensorType::Zero());
        for (size_t k=0;k<voigtSize;k++)
        {
            temp(tensorIndex(k,0),tensorIndex(k,1))=voigtvector(k);
        }
        return temp;
    }
    
    VoigtVectorType m2v(const TensorType& input_matrix) const
    {
        //from matrix to voigt format
        VoigtVectorType temp(VoigtVectorType::Zero());
        for (size_t k=0;k<voigtSize;k++)
        {
            temp(k)=input_matrix(tensorIndex(k,0),tensorIndex(k,1));
        }
        return temp;
    }
        
};
        
}
#endif
