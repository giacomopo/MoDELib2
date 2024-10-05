/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureContainer_cpp_
#define model_MicrostructureContainer_cpp_

#include <fstream>
#include <memory>
#include <vector>
#include <chrono>
#include <MicrostructureContainer.h>
#include <TerminalColors.h>


namespace model
{
    template <int dim>
    MicrostructureContainer<dim>::MicrostructureContainer(DislocationDynamicsBase<dim>& ddBase_in) :
    /* init */ MicrostructureBase<dim>("MicrostructureContainer",*this)
    /* init */,ddBase(ddBase_in)
    {
        
    }

    template <int dim>
    const typename MicrostructureContainer<dim>::BaseContainerType& MicrostructureContainer<dim>::microstructures() const
    {
        return *this;
    }

    template <int dim>
    typename MicrostructureContainer<dim>::BaseContainerType& MicrostructureContainer<dim>::microstructures()
    {
        return *this;
    }

    template <int dim>
    void MicrostructureContainer<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        for(auto& pair : microstructures())
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Initializing "<<pair->tag<<":"<<std::flush;
            pair->initializeConfiguration(configIO,f_file,F_labels);
            std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template <int dim>
    void MicrostructureContainer<dim>::solve()
    {
        for(auto& pair : microstructures())
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Solving "<<pair->tag<<":"<<std::flush;
            pair->solve();
            std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template <int dim>
    double MicrostructureContainer<dim>::getDt() const
    {
        const auto t0= std::chrono::system_clock::now();
        double dt(ddBase.simulationParameters.dtMax);
        for(auto& pair : microstructures())
        {
            std::cout<<pair->tag<<":"<<std::flush;
            const double pairDt(pair->getDt());
            std::cout<<" dt="<<pairDt<<std::endl;
            dt=std::min(dt,pairDt);
        }
        std::cout<<"Used dt="<<dt;
        std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        return dt;
    }

    template <int dim>
    void MicrostructureContainer<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        for(const auto& pair : microstructures())
        {
            pair->output(configIO,auxIO,f_file,F_labels);
        }
    }

    template <int dim>
    void MicrostructureContainer<dim>::updateConfiguration()
    {
        for(auto& pair : microstructures())
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Updating "<<pair->tag<<":"<<std::flush;
            pair->updateConfiguration();
            std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template <int dim>
    typename  MicrostructureContainer<dim>::MatrixDim MicrostructureContainer<dim>::averagePlasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->averagePlasticDistortion();
        }
        return temp;
    }

    template <int dim>
    typename  MicrostructureContainer<dim>::MatrixDim MicrostructureContainer<dim>::averagePlasticDistortionRate() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->averagePlasticDistortionRate();
        }
        return temp;
    }

    template <int dim>
    typename  MicrostructureContainer<dim>::VectorDim MicrostructureContainer<dim>::displacement(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        VectorDim temp(VectorDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->displacement(x,node,ele,guess);
        }
        return temp;
    }

    template <int dim>
    typename  MicrostructureContainer<dim>::MatrixDim MicrostructureContainer<dim>::stress(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->stress(x,node,ele,guess);
        }
        return temp;
    }

    template <int dim>
    typename  MicrostructureContainer<dim>::MatrixDim MicrostructureContainer<dim>::averageStress() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->averageStress();
        }
        return temp;
    }

    template<int dim>
    typename MicrostructureContainer<dim>::VectorDim MicrostructureContainer<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        VectorDim temp(VectorDim::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->inelasticDisplacementRate(x,node,ele,guess);
        }
        return temp;
    }

    template<int dim>
    typename MicrostructureContainer<dim>::VectorMSize MicrostructureContainer<dim>::mobileConcentration(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        VectorMSize temp(VectorMSize::Zero());
        for(const auto& pair : microstructures())
        {
            temp += pair->mobileConcentration(x,node,ele,guess);
        }
        return temp;
    }

    template struct MicrostructureContainer<3>;
}
#endif
