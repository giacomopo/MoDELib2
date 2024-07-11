/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfigurationFields_H_
#define model_ConfigurationFields_H_

#include <map>

#include <Eigen/Dense>

#include <DislocationDynamicsBase.h>
#include <DDconfigIO.h>
#include <DislocationLoopPatches.h>
#include <PeriodicGlidePlaneFactory.h>
#include <EshelbyInclusionBase.h>
#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>
#include <ClusterDynamicsBase.h>
#include <ElasticDeformationBase.h>
#include <StressStraight.h>
#include <MicrostructureContainer.h>
#include <DislocationNetwork.h>
#include <ClusterDynamics.h>
#include <ElasticDeformation.h>

namespace model
{

    template <int dim>
    class ConfigurationFields : private std::map<size_t,DislocationLoopPatches<dim>>
    /*                      */, private std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>
    /*                      */, public std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>>
    /*                      */, public std::map<size_t,PolyhedronInclusionNodeIO<dim>>
//    /*                      */, public ClusterDynamicsBase<dim>
//    /*                      */, public ElasticDeformationBase<dim>
    /*                      */, public MicrostructureContainer<dim>
    {
        
    public:
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef typename StressStraight<dim>::ConcentrationVectorType ConcentrationVectorType;
        typedef std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>> EshelbyInclusionContainerType;
        typedef std::map<size_t,PolyhedronInclusionNodeIO<dim>> PolyhedronInclusionNodeContainerType;
        typedef MicrostructureContainer<dim> MicrostructureContainerType;
        typedef DislocationNetwork<dim,0> DislocationNetworkType;
        typedef ClusterDynamics<dim> ClusterDynamicsType;
        typedef ElasticDeformation<dim> ElasticDeformationType;

        
        DislocationDynamicsBase<dim>& ddBase;
        const DDconfigIO<dim>& configIO;
        
        ConfigurationFields(DislocationDynamicsBase<dim>& ddBase_in,const DDconfigIO<dim>& configIO_in);
        void updateConfiguration();
        const std::map<size_t,DislocationLoopPatches<dim>>& loopPatches() const;
        std::map<size_t,DislocationLoopPatches<dim>>& loopPatches();
        const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments() const;
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments();
        const PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes() const;
        PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes();
        const EshelbyInclusionContainerType& eshelbyInclusions() const;
        EshelbyInclusionContainerType& eshelbyInclusions();
//        const ClusterDynamicsBase<dim>& clusterDynamics() const;
//        ClusterDynamicsBase<dim>& clusterDynamics();
//        const ElasticDeformationBase<dim>& elasticDeformation() const;
//        ElasticDeformationBase<dim>& elasticDeformation();

        double solidAngle(const VectorDim& x) const;
        VectorDim dislocationPlasticDisplacement(const VectorDim& x) const;
        MatrixDim dislocationStress(const VectorDim& x) const;
        MatrixDim inclusionStress(const VectorDim& x) const;
        ConcentrationVectorType dislocationMobileConentrations(const VectorDim& x,const Simplex<dim,dim>* guess) const;

    };

}
#endif
