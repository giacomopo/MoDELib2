#ifndef model_FixedDirichletSolver_H_
#define model_FixedDirichletSolver_H_

#include<Eigen/SparseCore>
#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>

#include<TrialBase.h>

namespace model
{
    class FixedDirichletSolver
    {// solves Ax=b given x=T*x1+g. Only b, and g can change
        
    public:
        typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMatrixType;
        
    private:

        SparseMatrixType A1;
        
        Eigen::SparseLU<SparseMatrixType> directSolver;
        Eigen::BiCGSTAB<SparseMatrixType> iterativeSolver;

        const std::map<size_t,double>* dirichletConditions;
        const Eigen::Matrix<double,Eigen::Dynamic,1>* dofVector;
        size_t gSize;
        size_t cSize;
        size_t tSize;
        
        SparseMatrixType A;
        SparseMatrixType T;

        public:
        
        const bool use_directSolver;
        const SparseMatrixType& getA() const;
        const SparseMatrixType& getT() const;
        FixedDirichletSolver(const bool& use_directSolver_in,const double& tol);
        Eigen::VectorXd solve(const Eigen::VectorXd& b) const;

        template<typename BilinearWeakFormType>
        void compute(const BilinearWeakFormType& bWF)
        {
            dirichletConditions=&TrialBase<typename BilinearWeakFormType::TrialFunctionType>::dirichletConditions();
            dofVector=&TrialBase<typename BilinearWeakFormType::TrialFunctionType>::dofVector();
            gSize=TrialBase<typename BilinearWeakFormType::TrialFunctionType>::gSize();
            cSize=dirichletConditions->size();
            tSize = gSize-cSize;

            
            A.resize(gSize,gSize);
            
//            const auto bWFT(bWF.globalTriplets());
//            
//            double triMin(0.0);
//            double triMax(0.0);
//            
//            for(const auto& tri : bWFT)
//            {
//                triMin=std::min(triMin,tri.value());
//                triMax=std::max(triMin,tri.value());
//            }
//            std::cout<<"trimMin="<<triMin<<std::endl;
//            std::cout<<"triMax="<<triMax<<std::endl;
            
            const auto aTriplets(bWF.globalTriplets());
            A.setFromTriplets(aTriplets.begin(),aTriplets.end());
            
            std::vector<Eigen::Triplet<double> > tTriplets;
            tTriplets.reserve(tSize);
            size_t startRow=0;
            size_t col=0;
            for (const auto& cIter : *dirichletConditions)
            {
                const size_t& endRow = cIter.first;
                for (size_t row=startRow;row!=endRow;++row)
                {
                    tTriplets.emplace_back(row,col,1.0);
                    col++;
                }
                startRow=endRow+1;
            }
            for (size_t row=startRow;row!=gSize;++row)
            {
                tTriplets.emplace_back(row,col,1.0);
                col++;
            }
            T.resize(gSize,tSize);
            T.setFromTriplets(tTriplets.begin(),tTriplets.end());
            
            A1=T.transpose()*A*T; // store new A1
            
            if(use_directSolver)
            {
                directSolver.compute(A1);
                if(directSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("FixedDirichletSolver failed.");
                }
            }
            else
            {
                iterativeSolver.compute(A1);
                if(iterativeSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("FixedDirichletSolver failed.");
                }
            }
        }
    };

} // namespace model
#endif
