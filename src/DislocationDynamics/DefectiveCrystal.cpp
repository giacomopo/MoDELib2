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

//    template <int _dim>
//    typename DefectiveCrystal<_dim>::MatrixDim DefectiveCrystal<_dim>::averagePlasticDistortion() const
//    {/*!\param[in] P position vector
//      * \returns The stress field in the DefectiveCrystal at P
//      * Note:
//      */
//        MatrixDim temp(MatrixDim::Zero());
//        for(const auto& pair : microstructures())
//        {
//            temp+pair.second->averagePlasticDistortion();
//        }
//        return temp;
//    }

//    template <int _dim>
//    typename DefectiveCrystal<_dim>::MatrixDim DefectiveCrystal<_dim>::averagePlasticStrain() const
//    {/*!\param[in] P position vector
//      * \returns The stress field in the DefectiveCrystal at P
//      * Note:
//      */
//        MatrixDim temp(averagePlasticDistortion());
//        return 0.5*(temp+temp.transpose());
//    }

//    template <int _dim>
//    typename DefectiveCrystal<_dim>::MatrixDim DefectiveCrystal<_dim>::averagePlasticDistortionRate() const
//    {/*!\param[in] P position vector
//      * \returns The stress field in the DefectiveCrystal at P
//      * Note:
//      */
//
//        MatrixDim temp(MatrixDim::Zero());
//        for(const auto& pair : microstructures())
//        {
//            temp+pair.second->averagePlasticDistortionRate();
//        }
//        return temp;
//    }

//    template <int _dim>
//    typename DefectiveCrystal<_dim>::MatrixDim DefectiveCrystal<_dim>::averagePlasticStrainRate() const
//    {/*!\param[in] P position vector
//      * \returns The stress field in the DefectiveCrystal at P
//      * Note:
//      */
//        MatrixDim temp(averagePlasticDistortionRate());
//        return 0.5*(temp+temp.transpose());
//    }

//    template <int _dim>
//    const typename DefectiveCrystal<_dim>::MicrostructureContainerType& DefectiveCrystal<_dim>::microstructures() const
//    {
//        return *this;
//    }
//
//    template <int _dim>
//    typename DefectiveCrystal<_dim>::MicrostructureContainerType& DefectiveCrystal<_dim>::microstructures()
//    {
//        return *this;
//    }


//    template <int _dim>
//    std::unique_ptr<ExternalLoadControllerBase<_dim>> DefectiveCrystal<_dim>::getExternalLoadController(const DislocationDynamicsBase<dim>& ddBase,const MatrixDim& plasticStrain_in)
//    {
//
//        std::cout<<"gettingExternalLoadController"<<std::endl;
//        if(!ddBase.isPeriodicDomain)
//        {
//            std::cout<<"getExternalLoadController a"<<std::endl;
//            return std::unique_ptr<ExternalLoadControllerBase<_dim>>(nullptr);
//        }
//        else
//        {
//            if(ddBase.simulationParameters.externalLoadControllerName=="UniformExternalLoadController")
//            {
//                std::cout<<"getExternalLoadController b"<<std::endl;
//                return std::unique_ptr<ExternalLoadControllerBase<_dim>>(new UniformExternalLoadController<_dim>(ddBase,plasticStrain_in));
//            }
//            else
//            {
//                std::cout<<"Unknown externalLoadController name "<<ddBase.simulationParameters.externalLoadControllerName<<"No controller applied."<<std::endl;
//                return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
//            }
//        }
//    }

//    template <int _dim>
//    void DefectiveCrystal<_dim>::updateLoadControllers(const long int& runID, const bool& isClimbStep)
//    {/*! Updates bvpSolver using the stress and displacement fields of the
//      *  current DD configuration.
//      */
//        if(bvpSolver)
//        {
//            if (!(runID%bvpSolver->stepsBetweenBVPupdates))
//            {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
//                std::cout<<"Updating bvpSolver ... "<<std::endl;
//                throw std::runtime_error("re-implement BVP update");
//    //                const int quadraturePerTriangle=37;
//    //                bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*DN, isClimbStep);
//            }
//        }
//        if (externalLoadController)
//        {
//            std::cout<<"Updating externalLoadController... "<<std::endl;
//            externalLoadController->update(averagePlasticStrain());
//        }
//    }

//        template <int _dim>
//        double DefectiveCrystal<_dim>::getMaxVelocity() const
//        {
//            double vmax = 0.0;
//
//            for (const auto &nodeIter : DN->networkNodes())
//            {
//                    const double vNorm(nodeIter.second.lock()->get_V().norm());
//                    if (vNorm > vmax)
//                    {
//                        vmax = vNorm;
//                    }
//            }
//            return vmax;
//        }

//        template <int _dim>
//        void DefectiveCrystal<_dim>::singleGlideStep()
//        {
//            if(DN)
//            {
//                std::cout<<"\n"<<blueBoldColor<< "runID="<<ddBase.simulationParameters.runID<<" (of "<<ddBase.simulationParameters.Nsteps<<")"
//                /*                    */<< ", time="<<ddBase.simulationParameters.totalTime<<std::endl;
//                std::cout<< "Glide step: networkNodes="<<DN->networkNodes().size()
//                /*                    */<< ", networkSegments="<<DN->networkLinks().size()
//                /*                    */<< ", loopNodes="<<DN->loopNodes().size()
//                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
//                /*                    */<< ", loops="<<DN->loops().size();
//
//            std::cout<< defaultColor<<std::endl;
//
//                DislocationNode<dim,corder>::totalCappedNodes=0;
//                DN->updateGeometry();
//                updateLoadControllers(ddBase.simulationParameters.runID, false);
////                const double maxVelocity(getMaxVelocity());
//                DN->assembleGlide(ddBase.simulationParameters.runID, maxVelocity);
//                DN->storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
////                DN->solveGlide();
//                DN->solveNodalVelocities(DN->glideSolver.get());
//                ddBase.simulationParameters.dt=DN->timeIntegrator.getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
//                DN->updateRates();
//                DN->io().output(ddBase.simulationParameters.runID);
//                DN->moveGlide(ddBase.simulationParameters.dt);
//                DN->executeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
//                if (DN->capMaxVelocity)
//                {
//                    std::cout<<redBoldColor<<"( "<<(DislocationNode<dim,corder>::totalCappedNodes)<<" total nodes capped "<<defaultColor<<std::endl;
//                    std::cout<<redBoldColor<<", "<<(double(DislocationNode<dim,corder>::totalCappedNodes)/double(DN->networkNodes().size()))<<" fraction of nodes capped "
//                    <<defaultColor<<" )"<<std::endl;
//                }
//                ddBase.simulationParameters.totalTime+=ddBase.simulationParameters.dt;
//                ++ddBase.simulationParameters.runID;
//            }
//        }

//template <int _dim>
//void DefectiveCrystal<_dim>::singleClimbStep()
//{
//    if(DN)
//    {
//        if(DN->climbSolver)
//        {
//            if(DN->plasticDistortionRate().norm()<DN->climbSolver->glideEquilibriumRate)
//            {
//                std::cout<<"\n"<<blueBoldColor<< "runID="<<ddBase.simulationParameters.runID<<" (of "<<ddBase.simulationParameters.Nsteps<<")"
//                /*                    */<< ", time="<<ddBase.simulationParameters.totalTime<<std::endl;
//                std::cout<< "Climb step: networkNodes="<<DN->networkNodes().size()
//                /*                    */<< ", networkSegments="<<DN->networkLinks().size()
//                /*                    */<< ", loopNodes="<<DN->loopNodes().size()
//                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
//                /*                    */<< ", loops="<<DN->loops().size();
//
//            std::cout<< defaultColor<<std::endl;
//
//        //        DislocationNode<dim,corder>::totalCappedNodes=0;
//                DN->updateGeometry();
//                updateLoadControllers(ddBase.simulationParameters.runID, false);
//        //        const double maxVelocity(getMaxVelocity());
//        //        DN->assembleGlide(ddBase.simulationParameters.runID, maxVelocity);
//        //        DN->storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
////                DN->solveClimb();
//                DN->solveNodalVelocities(DN->climbSolver.get());
//                ddBase.simulationParameters.dt=DN->timeIntegrator.getClimbTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
//                DN->updateRates();
//                DN->io().output(ddBase.simulationParameters.runID);
////                DN->moveClimb(ddBase.simulationParameters.dt);
////                DN->executeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
////                if (DN->capMaxVelocity)
////                {
////                    std::cout<<redBoldColor<<"( "<<(DislocationNode<dim,corder>::totalCappedNodes)<<" total nodes capped "<<defaultColor<<std::endl;
////                    std::cout<<redBoldColor<<", "<<(double(DislocationNode<dim,corder>::totalCappedNodes)/double(DN->networkNodes().size()))<<" fraction of nodes capped "
////                    <<defaultColor<<" )"<<std::endl;
////                }
//                ddBase.simulationParameters.totalTime+=ddBase.simulationParameters.dt;
//                ++ddBase.simulationParameters.runID;
//
//            }
//        }
//
//    }
//}
