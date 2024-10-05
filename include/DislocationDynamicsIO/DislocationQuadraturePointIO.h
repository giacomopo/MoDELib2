/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationQuadraturePointIO_H_
#define model_DislocationQuadraturePointIO_H_

#include <Eigen/Dense>
#include <ClusterDynamicsParameters.h>


namespace model
{
    template<int dim>
    struct DislocationQuadraturePointIO
    {

        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        size_t sourceID;
        size_t sinkID;
        int qID;
        VectorDim r;                          // position
        double j;                             // jacobian dl/du
        VectorDim rl;                         // unit tangent dr/dl
        
        MatrixDim stress;
        VectorDim pkForce;
        VectorDim stackingFaultForce;
        VectorDim lineTensionForce;
        VectorDim velocity;
        double elasticEnergyPerLength;
        double coreEnergyPerLength;
        Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize> cCD;
        Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize> cDD;

        /**********************************************************************/
        DislocationQuadraturePointIO() :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,stackingFaultForce(VectorDim::Zero())
        /* init */,lineTensionForce(VectorDim::Zero())
        /* init */,velocity(VectorDim::Zero())
        /* init */,elasticEnergyPerLength(0.0)
        /* init */,coreEnergyPerLength(0.0)
        /* init */,cCD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        /* init */,cDD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        {
        }
        
        /**********************************************************************/
        template<typename DislocationQuadraturePointType>
        DislocationQuadraturePointIO(const DislocationQuadraturePointType& qPoint) :
        /* init */ sourceID(qPoint.sourceID)
        /* init */,sinkID(qPoint.sinkID)
        /* init */,qID(qPoint.qID)
        /* init */,r(qPoint.r)
        /* init */,j(qPoint.j)
        /* init */,rl(qPoint.rl)
        /* init */,stress(qPoint.stress)
        /* init */,pkForce(qPoint.pkForce)
        /* init */,stackingFaultForce(qPoint.stackingFaultForce)
        /* init */,lineTensionForce(qPoint.lineTensionForce)
        /* init */,velocity(qPoint.velocity)
        /* init */,elasticEnergyPerLength(qPoint.elasticEnergyPerLength)
        /* init */,coreEnergyPerLength(qPoint.coreEnergyPerLength)
        /* init */,cCD(qPoint.cCD)
        /* init */,cDD(qPoint.cDD)
        {
            
        }
        
        /**********************************************************************/
        DislocationQuadraturePointIO(std::stringstream& ss) :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,stackingFaultForce(VectorDim::Zero())
        /* init */,lineTensionForce(VectorDim::Zero())
        /* init */,velocity(VectorDim::Zero())
        /* init */,elasticEnergyPerLength(0.0)
        /* init */,coreEnergyPerLength(0.0)
        /* init */,cCD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        /* init */,cDD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        {
            
            ss>>sourceID;
            ss>>sinkID;
            ss>>qID;
            for(int d=0;d<dim;++d)
            {
                ss>>r(d);
            }
            ss>>j;
            for(int d=0;d<dim;++d)
            {
                ss>>rl(d);
            }
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    ss>>stress(i,j);
                }
            }
            for(int d=0;d<dim;++d)
            {
                ss>>pkForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>stackingFaultForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>lineTensionForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>velocity(d);
            }
            ss>>elasticEnergyPerLength;
            ss>>coreEnergyPerLength;
            for(int d=0;d<ClusterDynamicsParameters<dim>::mSize;++d)
            {
                ss>>cCD(d);
            }
            for(int d=0;d<ClusterDynamicsParameters<dim>::mSize;++d)
            {
                ss>>cDD(d);
            }
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationQuadraturePointIO<dim>& ds)
        {
            os  << ds.sourceID<<"\t"
            /**/<< ds.sinkID<<"\t"
            /**/<< ds.qID<<"\t"
            /**/<<std::scientific
            /**/<<std::setprecision(15)
            /**/<< ds.r.transpose()<<"\t"
            /**/<< ds.j<<"\t"
            /**/<< ds.rl.transpose()<<"\t";
            for(int d=0;d<dim;++d)
            {
                os  << ds.stress.row(d)<<" ";
            }
            os  << ds.pkForce.transpose()<<"\t"
            /**/<< ds.stackingFaultForce.transpose()<<"\t"
            /**/<< ds.lineTensionForce.transpose()<<"\t"
            /**/<< ds.velocity.transpose()<<"\t"
            /**/<< ds.elasticEnergyPerLength<<"\t"
            /**/<< ds.coreEnergyPerLength<<"\t"
            /**/<< ds.cCD<<"\t"
            /**/<< ds.cDD;
            return os;
        }
        
    };
    
}
#endif
