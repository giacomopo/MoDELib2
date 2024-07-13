/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_cpp_
#define model_DefectiveCrystal_cpp_

#include <DefectiveCrystal.h>

namespace model
{
    template <int _dim>
    DefectiveCrystal<_dim>::DefectiveCrystal(DislocationDynamicsBase<_dim>& ddBase_in) :
    /* init */ MicrostructureContainerType(ddBase_in)
    /* init */,f_file(this->ddBase.simulationParameters.traitsIO.fFile,std::ios_base::app)
    /* init */,F_labels(this->ddBase.simulationParameters.traitsIO.flabFile,std::ios_base::app)
    {
        if (!f_file.is_open())
        {
            throw std::runtime_error("Cannot open file "+this->ddBase.simulationParameters.traitsIO.fFile);
        }
        if (!F_labels.is_open())
        {
            throw std::runtime_error("Cannot open file "+this->ddBase.simulationParameters.traitsIO.flabFile);
        }
        
        if(this->ddBase.simulationParameters.useInclusions)
        {
            this->emplace_back(new InclusionMicrostructureType(*this));
        }
        if(this->ddBase.simulationParameters.useClusterDynamics)
        {
            this->emplace_back(new ClusterDynamicsType(*this));
        }
        if(this->ddBase.simulationParameters.useDislocations)
        {
            this->emplace_back(new DislocationNetworkType(*this));
        }
        if(this->ddBase.simulationParameters.useElasticDeformation)
        {
            this->emplace_back(new ElasticDeformationType(*this));
        }

        DDconfigIO<dim> configIO(this->ddBase.simulationParameters.traitsIO.evlFolder);
        configIO.read(this->ddBase.simulationParameters.runID);
        this->initializeConfiguration(configIO);
    }

    template <int _dim>
    void DefectiveCrystal<_dim>::initializeConfiguration(const DDconfigIO<dim>& configIO)
    {
        this->initializeConfiguration(configIO,f_file,F_labels);

    }

    template <int _dim>
    void DefectiveCrystal<_dim>::runSingleStep()
    {
        std::cout<<"\n"<<blueBoldColor<< "runID="<<this->ddBase.simulationParameters.runID<<" (of "<<this->ddBase.simulationParameters.Nsteps<<")"
        /*                    */<< ", time="<<this->ddBase.simulationParameters.totalTime<<defaultColor<<std::endl;
        
        this->solve();
        this->ddBase.simulationParameters.dt=this->getDt();
                
        if (!(this->ddBase.simulationParameters.runID%this->ddBase.simulationParameters.outputFrequency))
        {
            DDconfigIO<dim> configIO(this->ddBase.simulationParameters.traitsIO.evlFolder);
            DDauxIO<dim> auxIO(this->ddBase.simulationParameters.traitsIO.auxFolder);
            
            f_file<< this->ddBase.simulationParameters.runID<<" "<<std::setprecision(15)<<std::scientific<<this->ddBase.simulationParameters.totalTime<<" "<<this->ddBase.simulationParameters.dt<<" ";

            const Eigen::Matrix<double,dim,dim>& pD(this->averagePlasticDistortion());
            f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" "<<pD.trace()<<" "<<pD.norm()<<" ";

            const Eigen::Matrix<double,dim,dim>& pDR(this->averagePlasticDistortionRate());
            f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" "<<pDR.trace()<<" "<<pDR.norm()<<" ";
            
            if(this->ddBase.simulationParameters.runID==0)
            {
                F_labels<<"runID\n";
                F_labels<<"time [b/cs]\n";
                F_labels<<"dt [b/cs]\n";
                
                F_labels<<"betaP_11\n";
                F_labels<<"betaP_12\n";
                F_labels<<"betaP_13\n";
                F_labels<<"betaP_21\n";
                F_labels<<"betaP_22\n";
                F_labels<<"betaP_23\n";
                F_labels<<"betaP_31\n";
                F_labels<<"betaP_32\n";
                F_labels<<"betaP_33\n";
                F_labels<<"tr(betaP)\n";
                F_labels<<"norm(betaP)\n";
                
                F_labels<<"dotBetaP_11 [cs/b]\n";
                F_labels<<"dotBetaP_12 [cs/b]\n";
                F_labels<<"dotBetaP_13 [cs/b]\n";
                F_labels<<"dotBetaP_21 [cs/b]\n";
                F_labels<<"dotBetaP_22 [cs/b]\n";
                F_labels<<"dotBetaP_23 [cs/b]\n";
                F_labels<<"dotBetaP_31 [cs/b]\n";
                F_labels<<"dotBetaP_32 [cs/b]\n";
                F_labels<<"dotBetaP_33 [cs/b]\n";
                F_labels<<"tr(dotBetaP) [cs/b]\n";
                F_labels<<"norm(dotBetaP) [cs/b]\n";
            }
            
            this->output(configIO,auxIO,f_file,F_labels);
            f_file<<std::endl;
            
            if(this->ddBase.simulationParameters.runID==0)
            {
                F_labels<<std::flush;
            }
            configIO.write(this->ddBase.simulationParameters.runID,this->ddBase.simulationParameters.outputBinary);
            auxIO.write(this->ddBase.simulationParameters.runID,this->ddBase.simulationParameters.outputBinary);
        }
        
        // updateConfiguration
        this->updateConfiguration();
        
        // Increment runID and time
        this->ddBase.simulationParameters.totalTime+=this->ddBase.simulationParameters.dt;
        ++this->ddBase.simulationParameters.runID;
    }

    template <int _dim>
    void DefectiveCrystal<_dim>::runSteps()
    {/*! Runs a number of simulation time steps defined by simulationParameters.Nsteps
      */
        const auto t0= std::chrono::system_clock::now();
        while (this->ddBase.simulationParameters.runID<this->ddBase.simulationParameters.Nsteps)
        {
            runSingleStep();
        }
        std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<this->ddBase.simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
    }

    template class DefectiveCrystal<3>;
}
#endif
