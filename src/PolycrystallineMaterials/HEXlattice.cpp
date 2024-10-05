/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_cpp_
#define model_HEXlattice_cpp_

#include <HEXlattice.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXpyramidal.h>


namespace model
{
    
        
        HEXlattice<3>::HEXlattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
        /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
        /* init */,PlaneNormalContainerType(getPlaneNormals(material,polyFile))
        /* init */,SlipSystemContainerType(getSlipSystems(material,*this))
        /* init */,SecondPhaseContainerType(getSecondPhases(material,*this))
        {
            
        }
        
        Eigen::Matrix<double,3,3> HEXlattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 1.0, 0.5,           0.0,
            /*   */ 0.0, 0.5*sqrt(3.0), 0.0,
            /*   */ 0.0, 0.0,           sqrt(8.0/3.0);
            
            return temp;
        }

    const typename HEXlattice<3>::PlaneNormalContainerType& HEXlattice<3>::planeNormals() const
    {
        return *this;
    }

    const typename HEXlattice<3>::SlipSystemContainerType& HEXlattice<3>::slipSystems() const
    {
        return *this;
    }

    const typename HEXlattice<3>::SecondPhaseContainerType& HEXlattice<3>::secondPhases() const
    {
        return *this;
    }


//    const typename HEXlattice<3>::DislocationMobilityContainerType& dislocationMobilities() const
//    {
//        return *this;
//    }


        
        std::vector<std::shared_ptr<GlidePlaneBase>> HEXlattice<3>::getPlaneNormals(const PolycrystallineMaterialBase& material,
                                                                                    const std::string& ) const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the HEX lattice
          */
            
            const bool enableBasalPlanes(material.enabledSlipSystems.find("fullBasal")!=material.enabledSlipSystems.end() || material.enabledSlipSystems.find("ShockleyBasal")!=material.enabledSlipSystems.end());
            const bool enablePrismaticPlanes(material.enabledSlipSystems.find("fullPrismatic")!=material.enabledSlipSystems.end());
            const bool enablePyramidalPlanes(material.enabledSlipSystems.find("fullPyramidal")!=material.enabledSlipSystems.end());

            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),*this);

            std::vector<std::shared_ptr<GlidePlaneBase>> temp;
            if(enableBasalPlanes)
            {
                const double ISF(TextFileParser(material.materialFile).readScalar<double>("ISF_SI",true)/(material.mu_SI*material.b_SI));
                const double USF(TextFileParser(material.materialFile).readScalar<double>("USF_SI",true)/(material.mu_SI*material.b_SI));
                const double MSF(TextFileParser(material.materialFile).readScalar<double>("MSF_SI",true)/(material.mu_SI*material.b_SI));
                
                const Eigen::Matrix<double,3,2> waveVectors((Eigen::Matrix<double,3,2>()<<0.0, 0.0,
                                                             /*                        */ 0.0, 1.0,
    //                                                         /*                        */ 1.0,-1.0
                                                             /*                        */ 1.0,1.0
                                                             ).finished());
                
                const Eigen::Matrix<double,4,3> f((Eigen::Matrix<double,4,3>()<<0.00,0.0, 0.0,
                                                   /*                        */ 0.50,sqrt(3.0)/6.0, ISF,
                                                   /*                        */ 0.25,sqrt(3.0)/12.0,USF,
                                                   /*                        */ 1.00,sqrt(3.0)/3.0, MSF).finished());
                
                const int rotSymm(3);
                const std::vector<Eigen::Matrix<double,2,1>> mirSymm;
                const Eigen::Matrix<double,2,2> A((Eigen::Matrix<double,2,2>()<< 1.0,-0.5,
                                                                                 0.0,0.5*std::sqrt(3.0)).finished());

                std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(A,waveVectors,f,rotSymm,mirSymm));
                temp.emplace_back(new GlidePlaneBase(a1,a2,gammaSurface));           // basal plane
            }

            if(enablePrismaticPlanes)
            {
                temp.emplace_back(new GlidePlaneBase(a1     ,c,nullptr));           // prismatic plane
                temp.emplace_back(new GlidePlaneBase(a3     ,c,nullptr));           // prismatic plane
                temp.emplace_back(new GlidePlaneBase(a2*(-1),c,nullptr));           // prismatic plane
            }
            
            if(enablePyramidalPlanes)
            {
                temp.emplace_back(new GlidePlaneBase(a1     ,a2+c,nullptr));         // pyramidal plane
                temp.emplace_back(new GlidePlaneBase(a2     ,a3+c,nullptr));         // pyramidal plane
                temp.emplace_back(new GlidePlaneBase(a3     ,c-a1,nullptr));        // pyramidal plane
                temp.emplace_back(new GlidePlaneBase(a1*(-1),c-a2,nullptr));       // pyramidal plane
                temp.emplace_back(new GlidePlaneBase(a2*(-1),c-a3,nullptr));       // pyramidal plane
                temp.emplace_back(new GlidePlaneBase(a3*(-1),a1+c,nullptr));        // pyramidal plane
            }
            
            return temp;
        }
        
        std::vector<std::shared_ptr<SlipSystem>> HEXlattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                               const PlaneNormalContainerType& plN) const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            
            
            const std::string dislocationMobilityTypeBasal(TextFileParser(material.materialFile).readString("dislocationMobilityTypeBasal",true));
            DislocationMobilitySelector mobilitySelectorBasal("HEXbasal");
            const std::shared_ptr<DislocationMobilityBase> hexMobilityBasal(mobilitySelectorBasal.getMobility(dislocationMobilityTypeBasal,material));

            const std::string dislocationMobilityTypePrismatic(TextFileParser(material.materialFile).readString("dislocationMobilityTypePrismatic",true));
            DislocationMobilitySelector mobilitySelectorPrismatic("HEXprismatic");
            const std::shared_ptr<DislocationMobilityBase> hexMobilityPrismatic(mobilitySelectorPrismatic.getMobility(dislocationMobilityTypePrismatic,material));

            const std::string dislocationMobilityTypePyramidal(TextFileParser(material.materialFile).readString("dislocationMobilityTypePyramidal",true));
            DislocationMobilitySelector mobilitySelectorPyramidal("HEXpyramidal");
            const std::shared_ptr<DislocationMobilityBase> hexMobilityPyramidal(mobilitySelectorPyramidal.getMobility(dislocationMobilityTypePyramidal,material));

            const int solidSolutionNoiseMode(TextFileParser(material.materialFile).readScalar<int>("solidSolutionNoiseMode",true));
            const int stackingFaultNoiseMode(TextFileParser(material.materialFile).readScalar<int>("stackingFaultNoiseMode",true));
            std::shared_ptr<GlidePlaneNoise> planeNoise((solidSolutionNoiseMode||stackingFaultNoiseMode)? new GlidePlaneNoise(material) : nullptr);

            
            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            typedef LatticeVector<dim> LatticeVectorType;
            typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
//            LatticeVectorType a3(a2-a1);
            LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),*this);

            

            
            const double     dBasal(ReciprocalLatticeDirectionType(a1.cross(a2)).planeSpacing());
            const double dPrismatic(ReciprocalLatticeDirectionType(a1.cross(c)).planeSpacing());
            const double dPyramidal(ReciprocalLatticeDirectionType(a1.cross(a2+c)).planeSpacing());
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase->planeSpacing()-dBasal)<FLT_EPSILON)
                {
                    const auto& a1(planeBase->primitiveVectors.first);
                    const auto& a2(planeBase->primitiveVectors.second);

                    const auto b1(a1);
                    const auto b2(a2);
                    const auto b3(b2-b1);

                    std::vector<RationalLatticeDirection<3>> slipDirs;

                    if(material.enabledSlipSystems.find("fullBasal")!=material.enabledSlipSystems.end())
                    {
                        // Full slip systems
                        slipDirs.emplace_back(Rational( 1,1),b1);
                        slipDirs.emplace_back(Rational(-1,1),b1);
                        slipDirs.emplace_back(Rational( 1,1),b2);
                        slipDirs.emplace_back(Rational(-1,1),b2);
                        slipDirs.emplace_back(Rational( 1,1),b3);
                        slipDirs.emplace_back(Rational(-1,1),b3);
                    }

                    if(material.enabledSlipSystems.find("ShockleyBasal")!=material.enabledSlipSystems.end())
                    {
                        // Shockley partials
                        slipDirs.emplace_back(Rational(1,3),b1-b3);
                        slipDirs.emplace_back(Rational(1,3),b1-b2);
                        slipDirs.emplace_back(Rational(1,3),b2-b1);
                        slipDirs.emplace_back(Rational(1,3),b2-b3);
                        slipDirs.emplace_back(Rational(1,3),b3-b2);
                        slipDirs.emplace_back(Rational(1,3),b3-b1);
                    }
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace_back(new SlipSystem(*planeBase, slipDir,hexMobilityBasal,planeNoise));
                    }
                    
                }
                
                if(std::fabs(planeBase->planeSpacing()-dPrismatic)<FLT_EPSILON)
                {
                    const auto& b(planeBase->primitiveVectors.first);
                    std::vector<RationalLatticeDirection<3>> slipDirs;

                    if(material.enabledSlipSystems.find("fullPrismatic")!=material.enabledSlipSystems.end())
                    {
                        // Full slip systems
                        slipDirs.emplace_back(Rational( 1,1),b);
                        slipDirs.emplace_back(Rational(-1,1),b);
                    }
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace_back(new SlipSystem(*planeBase, slipDir,hexMobilityPrismatic,planeNoise));
                    }
                }
                
                if(std::fabs(planeBase->planeSpacing()-dPyramidal)<FLT_EPSILON)
                {
                    const auto& b(planeBase->primitiveVectors.first);
                    std::vector<RationalLatticeDirection<3>> slipDirs;

                    if(material.enabledSlipSystems.find("fullPyramidal")!=material.enabledSlipSystems.end())
                    {
                        // Full slip systems
                        slipDirs.emplace_back(Rational( 1,1),b);
                        slipDirs.emplace_back(Rational(-1,1),b);
                    }
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace_back(new SlipSystem(*planeBase, slipDir,hexMobilityPyramidal,planeNoise));
                    }
                }
                
            }
                
            
            return temp;
        }
        
        

typename HEXlattice<3>::SecondPhaseContainerType HEXlattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                            const PlaneNormalContainerType& ) const
    {
        SecondPhaseContainerType temp;

        for(const std::string& sp : material.enabledSecondPhases)
        {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in HEX crystals.");
        }
        return temp;
    }
        
} // namespace model
#endif

