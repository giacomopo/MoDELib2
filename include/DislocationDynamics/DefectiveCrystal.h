/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <iostream>
#include <MicrostructureContainer.h>
#include <ElasticDeformation.h>
#include <DislocationDynamicsBase.h>
#include <DefectiveCrystalParameters.h>
#include <DislocationDynamicsModule.h>
#include <CrackSystem.h>
#include <ClusterDynamics.h>


namespace model
{
    
    template <int _dim>
    class DefectiveCrystal : public MicrostructureContainer<_dim>
    {
        
        std::ofstream f_file;
        std::ofstream F_labels;

    public:
        static constexpr int dim=_dim; // make dim available outside class
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef MicrostructureContainer<_dim> MicrostructureContainerType;
        
        typedef InclusionMicrostructure<dim> InclusionMicrostructureType;
        typedef DislocationNetwork<dim,0> DislocationNetworkType;
        typedef ClusterDynamics<dim> ClusterDynamicsType;
        typedef ElasticDeformation<dim> ElasticDeformationType;
        
        DefectiveCrystal(DislocationDynamicsBase<_dim>& ddBase_in) ;
        using MicrostructureContainer<_dim>::initializeConfiguration;
        void initializeConfiguration(const DDconfigIO<dim>& configIO);
        void runSteps();
        void runSingleStep();
        
        const DislocationNetwork<_dim,0>& dislocationNetwork() const;
        
    };
}
#endif
