/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * Additional contributors:
 *  Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationDynamicsBase_CPP_
#define model_DislocationDynamicsBase_CPP_

#include <DislocationDynamicsBase.h>
#include <CTM.h>
namespace model
{

    template <int dim>
    double getEwaldLength(const Eigen::Matrix<double,dim,Eigen::Dynamic>& B,const double& EwaldLengthFactor)
    {
        
        if(B.cols())
        {
            const double vol(sqrt((B.transpose()*B).determinant())/CTM::factorial(B.cols()));
            return EwaldLengthFactor*std::pow(vol,1.0/B.cols());
        }
        else
        {
            return 0.0;
        }
    }

    template <int dim>
    bool checkIfFullyPeriodicDomain(const SimplicialMesh<dim>& mesh)
    {
        bool temp(true);
        for(const auto& region1 : mesh.regions())
        {
            for(const auto& face1 : region1.second->faces())
            {
                if(face1.second->isExternal())
                {
                    temp = (temp && face1.second->periodicFacePair.second);
                    if(!temp)
                    {
                        break;
                    }
                }
            }
            if(!temp)
            {
                break;
            }
        }
        return temp;
    }

    template <int dim>
    std::vector<Eigen::Matrix<double,dim,1>> getPeriodicShifts(const Eigen::Matrix<double,dim,Eigen::Dynamic>& shiftVectors,const std::vector<int>& periodicImageSize)
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        std::vector<VectorDim> temp;
        temp.push_back(VectorDim::Zero());
        if(periodicImageSize.size())
        {
            std::cout<<"Box periodicity basis ("<<shiftVectors.size()<<") in columns:"<<std::endl<<shiftVectors<<std::endl;
            if(shiftVectors.cols()!=int(periodicImageSize.size()))
            {
                std::cout<<"shiftVectors.size()="<<shiftVectors.size()<<std::endl;
                std::cout<<"periodicImageSize.size()="<<periodicImageSize.size()<<std::endl;
                throw std::runtime_error("shiftVectors.size() must equal periodicImageSize.size()");
            }
            
            for(int k=0;k<shiftVectors.cols();++k)
            {
                const auto shiftVector(shiftVectors.col(k));
                const int imageSize(std::abs(periodicImageSize[k]));
                std::vector<VectorDim> newTemp;
                for(const auto& v:temp)
                {// grab existing shift vectors
                    for(int i=-imageSize;i<=imageSize;++i)
                    {// add current shift times corresponding image size
                        newTemp.push_back(v+i*shiftVector);
                    }
                }
                temp.swap(newTemp);
            }
        }
        return temp;
    }

    template <int _dim>
    DislocationDynamicsBase<_dim>::DislocationDynamicsBase( const std::string& folderName) :
    /* init */ simulationParameters( folderName)
    /* init */,voigtTraits((typename SymmetricVoigtTraits<3>::VoigtSizeMatrixType()<<0,0,1,1,2,2,0,1,1,2,0,2).finished())
    /* init */,mesh(simulationParameters.traitsIO.meshFile,
                    TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("F",_dim,_dim,true),
                    TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("X0",1,_dim,true).transpose(),
                    simulationParameters.periodicFaceIDs)
    /* init */,isPeriodicDomain(checkIfFullyPeriodicDomain(mesh))
    /* init */,periodicImageSize(isPeriodicDomain? TextFileParser(simulationParameters.traitsIO.ddFile).readArray<int>("periodicImageSize",true) : std::vector<int>())
    /* init */,periodicLatticeBasis(mesh.periodicBasis())
    /* init */,periodicLatticeReciprocalBasis(periodicLatticeBasis*(periodicLatticeBasis.transpose()*periodicLatticeBasis).inverse())
    /* init */,periodicShifts(getPeriodicShifts(periodicLatticeBasis,periodicImageSize))
    /* init */,poly(simulationParameters.traitsIO.polyFile,mesh)
    /* init */,fe((!isPeriodicDomain && simulationParameters.useFEM) ? new FiniteElement<ElementType>(mesh) : nullptr)
    /* init */,glidePlaneFactory(poly)
    /* init */,periodicGlidePlaneFactory(poly,glidePlaneFactory)
    /* init */,EwaldLength(isPeriodicDomain? getEwaldLength(periodicLatticeBasis,TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<double>("EwaldLengthFactor",true)) : 0.0)
    {
        if(!mesh.simplices().size())
        {
            throw std::runtime_error("Mesh is empty");
        }
        
        std::cout<<"isPeriodicDomain="<<isPeriodicDomain<<std::endl;
        std::cout<<"periodicLatticeBasis="<<periodicLatticeBasis<<std::endl;
        std::cout<<"periodicLatticeReciprocalBasis="<<periodicLatticeReciprocalBasis<<std::endl;
        std::cout<<"periodicShifts("<<periodicShifts.size()<<")="<<std::endl;
        for(const auto& shift : periodicShifts)
        {
            std::cout<<shift.transpose()<<std::endl;
        }
        std::cout<<"EwaldLength="<<EwaldLength<<std::endl;
    }

    template <int _dim>
    Eigen::VectorXd DislocationDynamicsBase<_dim>::periodicCoordinates(const VectorDim& x) const
    {
        const VectorDim dx(x-mesh.xCenter());
        return periodicLatticeReciprocalBasis.transpose()*dx;
    }

    template struct DislocationDynamicsBase<3>;

} // namespace model
#endif
