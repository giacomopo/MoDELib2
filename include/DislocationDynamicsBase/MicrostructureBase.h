/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureBase_H_
#define model_MicrostructureBase_H_

#include <vector>
#include <fstream>
#include <DDconfigIO.h>
#include <DDauxIO.h>
#include <DislocationDynamicsBase.h>
#include <ClusterDynamicsParameters.h>

namespace model
{

    template <int dim>
    struct MicrostructureContainer;

    template <int dim>
    struct MicrostructureBase
    {
        
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;

        typedef MicrostructureContainer<dim> MicrostructureContainerType;
        typedef Simplex<dim,dim> SimplexDim;
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;

        typedef typename ElementType::NodeType NodeType;

        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double,mSize,1> VectorMSize;

        MicrostructureBase(const std::string& tag_in,MicrostructureContainer<dim>& m_in);
        virtual ~MicrostructureBase() = default;
        
        const std::string tag;
        MicrostructureContainerType& microstructures;
        
        double lastUpdateTime;
        
        double timeSinceLastUpdate() const;
        std::set<const Grain<dim>*> pointGrains(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const;


        virtual void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels) = 0;
        virtual void solve() = 0;
        virtual double getDt() const = 0;
        virtual void output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const = 0;
        virtual void updateConfiguration() = 0;
        virtual MatrixDim averagePlasticDistortion() const = 0;
        virtual MatrixDim averagePlasticDistortionRate() const = 0;
        virtual VectorDim displacement(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const = 0;
        virtual MatrixDim stress(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const = 0;
        virtual MatrixDim averageStress() const = 0;
        virtual VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const = 0;
        virtual VectorMSize mobileConcentration(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const = 0;

        Eigen::Matrix<double,Eigen::Dynamic,dim> displacement(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,dim>>) const;

    };

}
#endif
