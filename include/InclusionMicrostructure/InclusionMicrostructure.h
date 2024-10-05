/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionMicrostructure_H_
#define model_InclusionMicrostructure_H_

#include <MicrostructureContainer.h>
#include <InclusionMicrostructureBase.h>

namespace model
{

    template<int dim>
    struct InclusionMicrostructure : public MicrostructureBase<dim>
    /*                           */, public InclusionMicrostructureBase<dim>
    {
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;

        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef MicrostructureContainer<dim> MicrostructureContainerType;

        InclusionMicrostructure(MicrostructureContainerType& mc);
                
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels) override;
        void solve() override;
        double getDt() const override;
        void output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const override;
        void updateConfiguration() override;
        MatrixDim averagePlasticDistortion() const override ;
        MatrixDim averagePlasticDistortionRate() const override;
        VectorDim displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
        MatrixDim stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
        MatrixDim averageStress() const override;
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        VectorMSize mobileConcentration(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;

    
    };



    
}
#endif

