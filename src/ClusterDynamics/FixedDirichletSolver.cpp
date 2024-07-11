#ifndef model_FixedDirichletSolver_cpp_
#define model_FixedDirichletSolver_cpp_

#include<FixedDirichletSolver.h>

namespace model
{

    const typename FixedDirichletSolver::SparseMatrixType& FixedDirichletSolver::getA() const
    {
        return A;
    }

    const typename FixedDirichletSolver::SparseMatrixType& FixedDirichletSolver::getT() const
    {
        return T;
    }

    
FixedDirichletSolver::FixedDirichletSolver(const bool& use_directSolver_in,const double& tol) :
    /* init */ dirichletConditions(nullptr)
    /* init */,dofVector(nullptr)
    /* init */,gSize(0)
    /* init */,cSize(0)
    /* init */,tSize(0)
//    /* init */,bWF(bWF_in)
    /* init */,use_directSolver(use_directSolver_in)
//    /* init */,tolerance(tol)
    {
        iterativeSolver.setTolerance(tol);
    }
        
    Eigen::VectorXd FixedDirichletSolver::solve(const Eigen::VectorXd& b) const
    {/*!@param[in] b the rhs of A*x=b
      * @param[in] y the guess for x
      */
        
        Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
        Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
        
        size_t startRow=0;
        size_t col=0;
        
        for (const auto& cIter : *dirichletConditions)
        {
            const size_t& endRow = cIter.first;
            g(endRow)= cIter.second;
//            std::cout<< g(endRow)<<std::endl;
            for (size_t row=startRow;row!=endRow;++row)
            {
                guess(col)=(*dofVector)(row);
                col++;
            }
            startRow=endRow+1;
        }
        for (size_t row=startRow;row!=gSize;++row)
        {
            guess(col)=(*dofVector)(row);
            col++;
        }
                
        Eigen::VectorXd b1(T.transpose()*(b-A*g));
        Eigen::VectorXd x(Eigen::VectorXd::Zero(b1.rows()));
        if(use_directSolver)
        {
            x=directSolver.solve(b1);
            if(directSolver.info()!=Eigen::Success)
            {
                throw std::runtime_error("Direct FixedDirichletSolver failed.");
            }
//            std::cout<<(A1*x-b1).norm()<<std::endl;
//            const double b1Norm(b1.norm());
//            if(b1Norm>0.0)
//            {
//                const double axbNorm((A1*x-b1).norm());
//                if(axbNorm/b1Norm>tolerance)
//                {
//                    std::cout<<"norm(A*x-b)/norm(b)="<<axbNorm/b1Norm<<std::endl;
//                    std::cout<<"tolerance="<<tolerance<<std::endl;
//                    assert(0 && "SOLVER FAILED");
//                }
//            }
        }
        else
        {
            x=iterativeSolver.solveWithGuess(b1,guess);
            if(iterativeSolver.info()!=Eigen::Success)
            {
                throw std::runtime_error("Iterative FixedDirichletSolver failed.");
            }
        }
        return T*x+g;
    }


} // namespace model
#endif
