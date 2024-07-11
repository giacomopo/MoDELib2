/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimpleNullSpaceSolver_H_
#define model_SimpleNullSpaceSolver_H_

#include <iostream>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>


namespace model
{

template<typename TrialFunctionType,typename DirectSolverType, typename IterativeSolverType>
struct SimpleNullSpaceSolver
{// solves Ax=b given x=T*x1+g. Only b, g, and T can change
    
    typedef Eigen::SparseMatrix<double> SparseMatrixType;
    
    const TrialFunctionType& trial;
    const SparseMatrixType& A;
    const double tolerance;
    const bool use_directSolver;
    
    SparseMatrixType T;
    SparseMatrixType A1;
    
    //#ifdef _MODEL_PARDISO_SOLVER_
    //    Eigen::PardisoLLT<SparseMatrixType> directSolver;
    //#else
    //    Eigen::SimplicialLLT<SparseMatrixType> directSolver;
    //#endif
    
    DirectSolverType directSolver;
    IterativeSolverType iterativeSolver;
    
    
    
    SimpleNullSpaceSolver(const TrialFunctionType& trial_in,const SparseMatrixType& A_in,
                          const double& tol, const bool& useDirect) :
    /* init */ trial(trial_in)
    /* init */,A(A_in)
    /* init */,tolerance(tol)
    /* init */,use_directSolver(useDirect)
    {
        
    }
        
    Eigen::VectorXd solve(const Eigen::VectorXd& b, const bool& useGuess)
    {/*!@param[in] b the rhs of A*x=b
      * @param[in] y the guess for x
      */
        
        std::cout<<"Setting up null-space solver..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        const size_t gSize=trial.gSize();
        const size_t cSize(trial.dirichletConditions().size());
        const size_t tSize = gSize-cSize;
        
        std::vector<Eigen::Triplet<double> > tTriplets;
        tTriplets.reserve(tSize);
        
        Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
        Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
        
        size_t startRow=0;
        size_t col=0;
        
        for (const auto& cIter : trial.dirichletConditions())
        {
            const size_t& endRow = cIter.first;
            g(endRow)= cIter.second;
            for (size_t row=startRow;row!=endRow;++row)
            {
                tTriplets.emplace_back(row,col,1.0);
                guess(col)=trial.dofVector()(row);
                col++;
            }
            startRow=endRow+1;
        }
        for (size_t row=startRow;row!=gSize;++row)
        {
            tTriplets.emplace_back(row,col,1.0);
            guess(col)=trial.dofVector()(row);
            col++;
        }
        
        SparseMatrixType T1(gSize,tSize);
        T1.setFromTriplets(tTriplets.begin(),tTriplets.end());
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
        // Check if null-space has changed
        const auto t4= std::chrono::system_clock::now();
        std::cout<<"Checking if null-space has changed... "<<std::flush;
        bool sameT=false;
        if(T.cols()==T1.cols() && T.rows()==T1.rows())
        {
            const double normT2=(T1-T).squaredNorm();
            if(normT2==0.0)
            {
                sameT=true;
            }
            
        }
        std::cout<<!sameT<<std::flush;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t4)).count()<<" sec]"<<defaultColor<<std::endl;
        
        // Update A1 if null-space has changed
        if(!sameT)
        { // need to re-factorized the LLT decomposition
            const auto t1= std::chrono::system_clock::now();
            if(use_directSolver)
            {
                //#ifdef _MODEL_PARDISO_SOLVER_
                //                std::cout<<"PardisoLLT: factorizing..."<<std::flush;
                //#else
                //                std::cout<<"SimplicialLLT: factorizing..."<<std::flush;
                //#endif
                std::cout<<"DirectSolver: factorizing..."<<std::flush;
                
            }
            else
            {
                std::cout<<"IterativeSolver: factorizing..."<<std::flush;
            }
            
            
            
            
            T=T1; // store new T
            A1=T.transpose()*A*T; // store new A1
            
            //                if(gSize==4913)
            //                {
            //                    std::ofstream outA ("A.txt", std::ofstream::out);
            //                    outA<<A.toDense()<<std::endl;
            //                    std::ofstream outT ("T.txt", std::ofstream::out);
            //                    outT<<T.toDense()<<std::endl;
            //                    std::ofstream outA1 ("A1.txt", std::ofstream::out);
            //                    outA1<<A1.toDense()<<std::endl;
            //
            //                }
            
            if(use_directSolver)
            {
                directSolver.compute(A1);
                assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
            }
            else
            {
                iterativeSolver.compute(A1);
                assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
            }
            
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        // Solve
        const auto t2= std::chrono::system_clock::now();
        Eigen::VectorXd b1(T.transpose()*(b-A*g));
        //            Eigen::VectorXd x(Eigen::VectorXd::Zero(gSize));
        Eigen::VectorXd x(Eigen::VectorXd::Zero(b1.rows()));
        
        
        if(use_directSolver)
        {
            //#ifdef _MODEL_PARDISO_SOLVER_
            //            std::cout<<"PardisoLLT: solving..."<<std::flush;
            //#else
            //            std::cout<<"SimplicialLLT: solving..."<<std::flush;
            //#endif
            std::cout<<"DirectSolver: solving..."<<std::flush;
            
            //                x=T*directSolver.solve(b1)+g;
            x=directSolver.solve(b1);
            assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
            const double b1Norm(b1.norm());
            if(b1Norm>0.0)
            {
                const double axbNorm((A1*x-b1).norm());
                if(axbNorm/b1Norm>tolerance)
                {
                    std::cout<<"norm(A*x-b)/norm(b)="<<axbNorm/b1Norm<<std::endl;
                    std::cout<<"tolerance="<<tolerance<<std::endl;
                    assert(0 && "SOLVER FAILED");
                }
            }
        }
        else
        {
            std::cout<<"IterativeSolver: solving..."<<std::flush;
            iterativeSolver.setTolerance(tolerance);
            //                x=T*iterativeSolver.solveWithGuess(b1,guess)+g;
            if(useGuess)
            {
                x=iterativeSolver.solveWithGuess(b1,guess);
            }
            else
            {
                x=iterativeSolver.solve(b1);
            }
            
            std::cout<<" (relative error ="<<iterativeSolver.error()<<", tolerance="<<iterativeSolver.tolerance()<<")";
            assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
        }
        
        
        
        
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<defaultColor<<std::endl;
        return T*x+g;
    }
    
};
    
}
#endif

