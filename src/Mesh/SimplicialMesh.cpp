/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMesh_cpp_
#define model_SimplicialMesh_cpp_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>
#include <map>

#include <Eigen/LU>

//#include <VertexReader.h>
#include <TerminalColors.h>
#include <SimplexTraits.h>
#include <Simplex.h>
#include <SimplexReader.h>
//#include <MeshStats.h>
#include <MeshRegionObserver.h>
// defines mode::cout
#include <MeshRegionBoundary.h>
#include <SimplexObserver.h>
//#include <SimplicialMeshFace.h>
#include <GmshReader.h>
#include <SimplicialMesh.h>


namespace model
{

    template<int dim>
    void SimplicialMesh<dim>::createMesh(const std::set<int>& periodicFaceIDs)
    {/*!
      */
        
        const auto t0= std::chrono::system_clock::now();
        vol0=0.0;
        size_t eleConter=0;
        for (const auto& eIter : this->simplexReader().elements())
        {
            insertSimplex(eIter.second.first,eIter.second.second);
            eleConter++;
            std::cout<<greenBoldColor<<"\r"<<"Creating mesh "<<eleConter*100/this->simplexReader().elements().size()<<"%"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::flush;
        }
        std::cout<<defaultColor<<std::endl;
        
        this->info(); // print mesh info
        
        if(simplices().size())
        {
            _xMin=this->template observer<0>().begin()->second->P0;
            _xMax=this->template observer<0>().begin()->second->P0;
            
            for (const auto& nIter : this->template observer<0>())
            {
                for(int d=0;d<dim;++d)
                {
                    if (nIter.second->P0(d)<_xMin(d))
                    {
                        _xMin(d)=nIter.second->P0(d);
                    }
                    if (nIter.second->P0(d)>_xMax(d))
                    {
                        _xMax(d)=nIter.second->P0(d);
                    }
                }
            }
            _xC=0.5*(_xMax+_xMin);
        }
        else
        {
            std::cout<<"Mesh is empty."<<std::endl;
        }
        std::cout<<"  xMin="<< std::setprecision(15)<<std::scientific<<_xMin.transpose()<<std::endl;
        std::cout<<"  xMax="<< std::setprecision(15)<<std::scientific<<_xMax.transpose()<<std::endl;
        std::cout<<"  xCenter="<< std::setprecision(15)<<std::scientific<<_xC.transpose()<<std::endl;
        
        // Populate typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType
        regionBoundaries().clear();
        
        size_t bndSimplexCount=0;
        size_t rgnBndSimplexCount=0;
        for (const auto& simpl : this->template observer<dim-1>())
        {
            
            if(simpl.second->isBoundarySimplex())
            {// count number of bonudary simplices for later check
                bndSimplexCount++;
            }
            
            if(simpl.second->isRegionBoundarySimplex())
            {
                const auto regionIDset=simpl.second->regionIDs();
                std::pair<size_t,size_t> regionIDs(std::make_pair(*regionIDset.begin(),*regionIDset.rbegin()));
                const auto regionBndIter=regionBoundaries().find(regionIDs);
                if(regionBndIter!=regionBoundaries().end())
                {
                    regionBndIter->second.simplices().insert(simpl.second);
                }
                else
                {
                    regionBoundaries().emplace(regionIDs,regionIDs).first->second.simplices().insert(simpl.second);
                }
                rgnBndSimplexCount++;
            }
        }
        
        updateRegions();
        updateRegionBoundaries();
        identifyParallelFaces(periodicFaceIDs);
        
        size_t bndFaceSimplexSum=0;
        for(auto region : MeshRegionObserverType::regions())
        {// Sum number of external faces for final check
            std::cout<<magentaColor<<"MeshRegion "<<region.second->regionID<<defaultColor<<std::endl;
            std::cout<<"    simplices: "<<region.second->simplices().size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            for(auto& face : region.second->faces())
            {
                std::cout<<"    face "<<face.second->sID<<": size="<<face.second->size()<<",hullPts="<<face.second->convexHull().size()<<", outNormal "<<face.second->outNormal().transpose()<<std::endl;
                bndFaceSimplexSum+=face.second->size();
            }
            std::cout<<"    parallel faces:"<<std::endl;
            for(const auto& pair : region.second->parallelFaces())
            {
                std::cout<<"      "<<pair.first<<"<->"<<pair.second<<std::endl;
            }
        }
        
        size_t rgnBndFaceSimplexSum=0;
        for(const auto& rgnBnd : regionBoundaries())
        {// Sum number of internal faces for final check
            std::cout<<magentaColor<<"MeshRegionBoundary ("<<rgnBnd.second.regionBndID.first<<","<<rgnBnd.second.regionBndID.second<<")"<<defaultColor<<std::endl;
            std::cout<<"    simplices: "<<rgnBnd.second.simplices().size()<<" Simplex<"<<dim<<","<<dim-1<<">"<<std::endl;
            for(auto& face : rgnBnd.second.faces())
            {
                std::cout<<"    face "<<face.second->sID<<": hullPts="<<face.second->convexHull().size()<<", outNormal "<<face.second->outNormal().transpose()<<std::endl;
                rgnBndFaceSimplexSum+=face.second->size();
                bndFaceSimplexSum-=2*face.second->size(); // each region boundary face was added to two regions
            }
        }
        
        if(bndFaceSimplexSum!=bndSimplexCount)
        {
            std::cout<<"boundary simplices="<<bndSimplexCount<<std::endl;
            std::cout<<"simplices in external faces="<<bndFaceSimplexSum<<std::endl;
            throw std::runtime_error("WRONG NUMBER OF BOUNDAY FACE SIMPLICES");
        }
        
        if(rgnBndFaceSimplexSum!=rgnBndSimplexCount)
        {
            std::cout<<"region-boundary simplices="<<rgnBndSimplexCount<<std::endl;
            std::cout<<"simplices in internal faces="<<rgnBndFaceSimplexSum<<std::endl;
            throw std::runtime_error("WRONG NUMBER OF REGION-BOUNDARY FACE SIMPLICES");
        }
    }

    template<int dim>
    SimplicialMesh<dim>::SimplicialMesh() :
    /* init */ _xMin(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xMax(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xC(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,vol0(0.0)
    {
    }

    template<int dim>
    SimplicialMesh<dim>::SimplicialMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>& periodicFaceIDs) :
    /* init */ _xMin(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xMax(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xC(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,vol0(0.0)
    {
        this->read(meshFileName,A,x0);
        createMesh(periodicFaceIDs);
        this->simplexReader().clear();
    }

    template<int dim>
    void SimplicialMesh<dim>::readMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>& periodicFaceIDs)
    {
        simplices().clear();
        this->read(meshFileName,A,x0);
        createMesh(periodicFaceIDs);
        this->simplexReader().clear();
    }

    template<int dim>
    const typename SimplicialMesh<dim>::SimplexMapType& SimplicialMesh<dim>::simplices() const
    {
        return *this;
    }

    template<int dim>
    typename SimplicialMesh<dim>::SimplexMapType& SimplicialMesh<dim>::simplices()
    {
        return *this;
    }

    template<int dim>
    void SimplicialMesh<dim>::updateRegions()
    {
        for(auto region : MeshRegionObserverType::regions())
        {
            region.second->update();
        }
    }

    template<int dim>
    void SimplicialMesh<dim>::updateRegionBoundaries()
    {
        for(auto& rgnBnd : regionBoundaries())
        {
            rgnBnd.second.update();
            MeshRegionType* const region1(this->region(rgnBnd.second.regionBndID.first));
            MeshRegionType* const region2(this->region(rgnBnd.second.regionBndID.second));
            for(const auto& face : rgnBnd.second.faces())
            {// add each face to the region boundary to the corresponding two regions
                region1->faces().emplace(face.second->sID,face.second);
                region2->faces().emplace(face.second->sID,face.second);
                
            }
        }
    }

    template<int dim>
    void SimplicialMesh<dim>::identifyParallelFaces(const std::set<int>& periodicFaceIDs)
    {
        
        typedef std::map<VectorDim,PlanarMeshFace<dim>*,CompareVectorsByComponent<double,dim,float>> MultipleShiftMapType;
        
        std::map<PlanarMeshFace<dim>*,MultipleShiftMapType> faceShiftsMap;
        for(const auto& region1 : this->regions())
        {
            for(const auto& face1 : region1.second->faces())
            {
                if(face1.second->isExternal())
                {
                    for(const auto& region2 : this->regions())
                    {
                        for(const auto& face2 : region2.second->faces())
                        {
                            if(face2.second->isExternal() && face2.second!=face1.second && face1.second->convexHull().size()==face2.second->convexHull().size())
                            {
                                //                                std::cout<<"face pair "<<face1.second->sID<<" ("<<face1.second->convexHull().size()<<"), "<<face2.second->sID<<" ("<<face2.second->convexHull().size()<<")"<<std::endl;
                                if(abs(face1.second->outNormal().dot(face2.second->outNormal())+1.0)<FLT_EPSILON)
                                {// distinct parallel external faces
                                    const VectorDim shift(face1.second->center()-face2.second->center());
                                    std::set<const Simplex<dim,0>*> secondFaceVertices;
                                    
                                    for(const auto& face2Vertex : face2.second->convexHull())
                                    {
                                        secondFaceVertices.insert(face2Vertex);
                                    }
                                    for(const auto& face1Vertex : face1.second->convexHull())
                                    {
                                        for(const auto& face2Vertex : secondFaceVertices)
                                        {
                                            //                                            std::cout<<(face1Vertex->P0-face2Vertex->P0-shift).norm()<<std::endl;
                                            if((face1Vertex->P0-face2Vertex->P0-shift).norm()<FLT_EPSILON)
                                            {
                                                secondFaceVertices.erase(face2Vertex);
                                                break;
                                            }
                                        }
                                    }
                                    //                                    std::cout<<secondFaceVertices.size()<<std::endl;
                                    if(secondFaceVertices.size()==0)
                                    {//all vertices of face1 have been paired to a single vertex of face2 in a translation by shift
                                        const bool success(faceShiftsMap[face1.second.get()].emplace(-shift,face2.second.get()).second);
                                        if(!success)
                                        {
                                            throw std::runtime_error("Cannot insert shift for face "+std::to_string(face1.second->sID));
                                        }
                                        //                                    std::cout<<"Detected periodic face pair "<<pair.first<<"-"<<pair.second<<std::endl;
                                        //                                    faces()[pair.first] ->periodicFacePair=std::make_pair(-shift,faces()[pair.second].get());
                                        //                                    faces()[pair.second]->periodicFacePair=std::make_pair( shift,faces()[pair.first ].get());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        std::map<PlanarMeshFace<dim>*,std::pair<VectorDim,PlanarMeshFace<dim>*>> uniqueFaceShiftMap;
        for(const auto& pair : faceShiftsMap)
        {
            //            std::cout<<"Face "<<pair.first->sID<<" has "<<pair.second.size()<<" parallel faces"<<std::endl;
            //            std::cout<<pair.second.begin()->first.transpose()<<std::endl;
            if(pair.second.size()==1)
            {
                uniqueFaceShiftMap.emplace(pair.first,std::make_pair(pair.second.begin()->first,pair.second.begin()->second));
            }
            else
            {
                MultipleShiftMapType allowedShifts; // allowedShifts for current face
                for(const auto& currentShift : pair.second)
                {
                    bool shiftAllowed(true);
                    for(const auto& otherPair : faceShiftsMap)
                    {
                        if(pair.first!=otherPair.first)
                        {
                            const bool shiftIsParallel(std::fabs(otherPair.first->outNormal().dot(currentShift.first))<FLT_EPSILON);
                            const bool shiftIsSame(otherPair.second.find(currentShift.first)!=otherPair.second.end() || otherPair.second.find(-currentShift.first)!=otherPair.second.end());
                            //                        std::cout<<" face "<<otherPair.first->sID<<" "<<shiftIsParallel<<" "<<shiftIsSame<<std::endl;
                            shiftAllowed= shiftAllowed && (shiftIsParallel || shiftIsSame);
                        }
                    }
                    if(shiftAllowed)
                    {
                        allowedShifts.emplace(currentShift.first,currentShift.second);
                    }
                }
                if(allowedShifts.size())
                {
                    if(allowedShifts.size()==1)
                    {
                        uniqueFaceShiftMap.emplace(pair.first,std::make_pair(allowedShifts.begin()->first,allowedShifts.begin()->second));
                    }
                    else
                    {
                        throw std::runtime_error("Multiple periodic shift for face "+std::to_string(pair.first->sID));
                    }
                }
                else
                {
                    //            std::cout<<"Face non periodic"
                }
            }
            
        }
        
        for(auto& pair : uniqueFaceShiftMap)
        {
            const auto parallelIter(uniqueFaceShiftMap.find(pair.second.second));
            if(parallelIter!=uniqueFaceShiftMap.end())
            {// the parallel face is also in the map
                if((parallelIter->second.first+pair.second.first).norm()<FLT_EPSILON)
                {
                    if(   periodicFaceIDs.find(pair.first->sID)!=periodicFaceIDs.end()
                       || periodicFaceIDs.find(parallelIter->first->sID)!=periodicFaceIDs.end()
                       || periodicFaceIDs.find(-1)!=periodicFaceIDs.end()
                       )
                    {
                        pair.first->periodicFacePair=std::make_pair(pair.second.first,pair.second.second);
                        std::cout<<"Face "<<pair.first->sID<<" found periodic face "<<pair.first->periodicFacePair.second->sID<<", shift="<<pair.first->periodicFacePair.first.transpose()<<std::endl;
                    }
                    else
                    {
                        std::cout<<"Face not labelled periodic"<<std::endl;
                    }
                }
                else
                {
                    std::cout<<"face "<< pair.first->sID<<" shift="<<pair.second.first.transpose()<<std::endl;
                    std::cout<<"parallel face "<<parallelIter->first->sID<<" shift="<<parallelIter->second.first.transpose()<<std::endl;
                    throw std::runtime_error("Parallel Face shift not opposite");
                }
            }
            else
            {
                std::cout<<"face "<< pair.first->sID<<" parallel to "<< pair.second.second->sID<<", shift="<<pair.second.first.transpose()<<std::endl;
                throw std::runtime_error("Parallel Face not in map");
            }
        }
        
        // Check that periodic faces have been correclty identified
        for(const auto& region1 : this->regions())
        {
            for(const auto& face1 : region1.second->faces())
            {
                if(face1.second->isExternal())
                {
                    if(face1.second->periodicFacePair.second)
                    {// a parallel face found
                        if(   periodicFaceIDs.find(face1.second->sID)!=periodicFaceIDs.end()
                           || periodicFaceIDs.find(face1.second->periodicFacePair.second->sID)!=periodicFaceIDs.end()
                           || periodicFaceIDs.find(-1)!=periodicFaceIDs.end())
                        {// all good, ID found
                            
                        }
                        else
                        {
                            std::cout<<"Face "<<face1.second->sID<<std::endl;
                            std::cout<<"Parallel pace "<<face1.second->periodicFacePair.second->sID<<std::endl;
                            throw std::runtime_error("Found parallel faces but their IDs are not in periodicFaceIDs");
                        }
                    }
                    else
                    {// a parallel face not found
                        if(   periodicFaceIDs.find(face1.second->sID)!=periodicFaceIDs.end()
                           || periodicFaceIDs.find(-1)!=periodicFaceIDs.end())
                        {// all good, ID found
                            std::cout<<"Face "<<face1.second->sID<<std::endl;
                            throw std::runtime_error("Face in periodicFaceIDs but no parallel face found");
                        }
                        else
                        {// all good, not a periodic face
                            
                        }
                    }
                }
            }
        }
    }

    template<int dim>
    typename SimplicialMesh<dim>::PeriodicBasisType SimplicialMesh<dim>::periodicBasis() const
    {
        const VectorDim meshSize(xMax()-xMin());
        PeriodicBasisType basis(PeriodicBasisType::Zero(dim,0));
        size_t basisRank(0);
        
        for(const auto& region : this->regions())
        {
            for(const auto& face : region.second->faces())
            {
                if(face.second->periodicFacePair.second)
                {// face is a periodic face
                    PeriodicBasisType newBasis(PeriodicBasisType::Zero(dim,basis.cols()+1));
                    newBasis.block(0,0,dim,basis.cols())=basis;
                    if(meshSize.dot(face.second->periodicFacePair.first)>0.0)
                    {
                        newBasis.col(basis.cols())=face.second->periodicFacePair.first;
                    }
                    else
                    {
                        newBasis.col(basis.cols())=-face.second->periodicFacePair.first;
                    }
                    Eigen::FullPivLU<Eigen::Matrix<double,dim,Eigen::Dynamic>> lu(newBasis);
                    const size_t newBasisRank(lu.rank());
                    if(newBasisRank>basisRank)
                    {
                        basis=newBasis;
                        basisRank=newBasisRank;
                    }
                }
            }
        }
        return basis;
    }

    template<int dim>
    void SimplicialMesh<dim>::insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID)
    {/*!@param[in] xIN the (unsorted) array of mesh node IDs defining a simplex
      *\brief Inserts a Simplex<dim> into the mesh, with ID equal to the
      * sorted array xIN.
      */
        const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
        const auto pair=simplices().emplace(std::piecewise_construct,
                                            std::make_tuple(xID),
                                            std::make_tuple(this,xID, regionID)
                                            );
        assert(pair.second);
        vol0+=pair.first->second.vol0;
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::search(const Eigen::Matrix<double,dim,1>& P) const
    {/*!@param[in] P position to search for
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      *
      * By default the search starts at this->begin()->second
      */
        return searchWithGuess(P,&(simplices().begin()->second));
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                                 const Simplex<dim,dim>* const guess) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        
        std::set<const Simplex<dim,dim>*> searchSet;
        if(guess)
        {
            return searchWithGuess(true,P,guess,searchSet);
        }
        else
        {
            return searchWithGuess(true,P,&(simplices().begin()->second),searchSet);
        }
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchRegion(const int& regionID,
                                                                              const Eigen::Matrix<double,dim,1>& P) const
    {/*!@param[in] P position to search for
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      *
      * By default the search starts at this->begin()->second
      */
        return searchRegionWithGuess(P,*this->region(regionID)->simplices().begin());
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchRegionWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                                       const Simplex<dim,dim>* const guess) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        std::set<const Simplex<dim,dim>*> searchSet;
        return searchWithGuess(false,P,guess,searchSet);
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchWithGuess(const bool& searchAllRegions,
                                                                                 const Eigen::Matrix<double,dim,1>& P,
                                                                                 const Simplex<dim,dim>* const guess,
                                                                                 std::set<const Simplex<dim,dim>*>& searchSet) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts. If searchAllRegions=false, only the region of guess is searched
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        
        searchSet.clear();
        std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
        guess->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
        checkSearch(searchAllRegions,P,guess,lastSearched);
        
        if(!lastSearched.first)
        {// if search was not successful, try to search neighbors of last searched Simplex
            const std::pair<bool,const Simplex<dim,dim>*> temp(lastSearched);
            const Eigen::Matrix<double,dim+1,1> bary=temp.second->pos2bary(P);
            for(int k=0;k<dim+1;++k)
            {
                if(bary(k)<=FLT_EPSILON)
                {
                    //                        for(typename Simplex<dim,dim-1>::ParentContainerType::const_iterator pIter=temp.second->child(k).parentBegin();
                    //                            /*                                                            */ pIter!=temp.second->child(k).parentEnd();++pIter)
                    //                        {
                    //                            if((*pIter)->region->regionID==temp.second->region->regionID || searchAllRegions)
                    //                            {
                    //                                (*pIter)->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                    //                                if (lastSearched.first)
                    //                                {
                    //                                    break;
                    //                                }
                    //                            }
                    //                        }
                    for(const auto& pIter : temp.second->child(k).parents())
                    {
                        if(pIter.second->region->regionID==temp.second->region->regionID || searchAllRegions)
                        {
                            pIter.second->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                            if (lastSearched.first)
                            {
                                break;
                            }
                        }
                    }
                }
            }
            
            if(searchAllRegions)
            {// if search is still unsuccessful, reset lastSearched to temp
                if(!lastSearched.first && !lastSearched.second->isBoundarySimplex())
                {
                    lastSearched=temp;
                }
            }
            else
            {// if search is still unsuccessful, reset lastSearched to temp
                if(!lastSearched.first && !lastSearched.second->isBoundarySimplex() && !lastSearched.second->isRegionBoundarySimplex())
                {
                    lastSearched=temp;
                }
            }
        }
        
        checkSearch(searchAllRegions,P,guess,lastSearched);
        
        return lastSearched;
    }

    template<int dim>
    void SimplicialMesh<dim>::checkSearch(const bool& searchAllRegions,
                                          const Eigen::Matrix<double,dim,1>& P,
                                          const Simplex<dim,dim>* const guess,
                                          const std::pair<bool,const Simplex<dim,dim>*>& lastSearched) const
    {
        if(searchAllRegions)
        {// Check that search was successful, or that it ended on a boundary (point outside)
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY SIMPLEX");
            }
        }
        else
        {// Check that search was successful, or that it ended on a boundary or region-boundary
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex() || lastSearched.second->isRegionBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY/REGION-BOUNDARY SIMPLEX");
            }
        }
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::isStrictlyInsideMesh(const Eigen::Matrix<double,dim,1>& P,
                                                                                      const Simplex<dim,dim>* const guess,
                                                                                      const double& tol) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
        if(temp.first)
        {
            const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
            int kMin;
            const double baryMin(bary.minCoeff(&kMin));
            if (std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex()) // on a boundary face
            {
                temp.first=false;
            }
        }
        return temp;
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::isOnMeshBoundary(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess, const double& tol) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
        if(temp.first)
        {
            const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
            int kMin;
            const double baryMin(bary.minCoeff(&kMin));
            if (!(std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex())) // not on a boundary face
            {
                temp.first=false;
            }
        }
        return temp;
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xMin() const
    {
        return _xMin;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::xMin(const int& k) const
    {
        return _xMin(k);
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xMax() const
    {
        return _xMax;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::xMax(const int& k) const
    {
        return _xMax(k);
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xCenter() const
    {
        return _xC;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::volume() const
    {
        return vol0;
    }

    template<int dim>
    const typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType& SimplicialMesh<dim>::regionBoundaries() const
    {
        return *this;
    }

    template<int dim>
    typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType& SimplicialMesh<dim>::regionBoundaries()
    {
        return *this;
    }

    template<int dim>
    const typename SimplicialMesh<dim>::MeshRegionBoundaryType& SimplicialMesh<dim>::regionBoundary(const int& i,const int& j) const
    {
        return regionBoundaries().at(std::make_pair(std::min(i,j),std::max(i,j)));
    }

    template class SimplicialMesh<3>;
}
#endif
