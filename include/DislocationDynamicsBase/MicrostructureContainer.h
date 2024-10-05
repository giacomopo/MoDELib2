/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureContainer_H_
#define model_MicrostructureContainer_H_

#include <fstream>
#include <memory>
#include <vector>
#include <DislocationDynamicsBase.h>
#include <MicrostructureBase.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

namespace model
{
    template <int dim>
    struct MicrostructureContainer : public MicrostructureBase<dim>
    /*                           */, public std::vector<std::shared_ptr<MicrostructureBase<dim>>>
    {
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;

        typedef std::vector<std::shared_ptr<MicrostructureBase<dim>>> BaseContainerType;
        
        DislocationDynamicsBase<dim>& ddBase;

        MicrostructureContainer(DislocationDynamicsBase<dim>& ddBase_in);
        const BaseContainerType& microstructures() const;
        BaseContainerType& microstructures();
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

//        template<typename MicrostructureType>
//        std::set<MicrostructureType*> getTypedMicrostructures() const
//        {
//            std::set<MicrostructureType*> temp;
//            for(const auto& microstructure : microstructures())
//            {
//                if(dynamic_cast<MicrostructureType*>(microstructure.get()))
//                {
//                    temp.insert(dynamic_cast<MicrostructureType*>(microstructure.get()));
//                }
//            }
//            return temp;
//        }
                
        template<typename MicrostructureType>
        std::set<std::shared_ptr<MicrostructureType>> getTypedMicrostructures() const
        {
            std::set<std::shared_ptr<MicrostructureType>> temp;
            for(const auto& microstructure : microstructures())
            {
                std::shared_ptr<MicrostructureType> derived(std::dynamic_pointer_cast<MicrostructureType>(microstructure));
                if(derived)
                {
                    temp.insert(derived);
                }
            }
            return temp;
        }
        
//        template<typename MicrostructureType>
//        MicrostructureType* getUniqueTypedMicrostructure() const
//        {
//            const auto cdMicrostructures(getTypedMicrostructures<MicrostructureType>());
//            if(cdMicrostructures.size()==1)
//            {
//                return *cdMicrostructures.begin();
//            }
//            else
//            {
//                std::cout<<"Found "<<cdMicrostructures.size()<<" microstructures of requested type"<<std::endl;
//                return nullptr;
//            }
//        }
        
        template<typename MicrostructureType>
        std::shared_ptr<MicrostructureType> getUniqueTypedMicrostructure() const
        {
            const auto cdMicrostructures(getTypedMicrostructures<MicrostructureType>());
            if(cdMicrostructures.size()==1)
            {
                return *cdMicrostructures.begin();
            }
            else
            {
                std::cout<<"Found "<<cdMicrostructures.size()<<" microstructures of requested type"<<std::endl;
                return nullptr;
            }
        }
        
    };
}
#endif
