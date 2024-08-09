/* This file is part of GreatWhite, a library for integrating
 * MOOSE and MoDeLib.
 *
 * (c) 2018 Multiscale Materials Solutions, LLC
 * ALL RIGTHS REVERVED
 *   Prepared by 2018 Multiscale Materials Solutions
 *     under contract N. 126180
 *   with the NAVAL NUCLEAR LABORATORY
 *
 *  See COPYRIGHT for full restriction
*/


#ifndef model_ClusterDynamicsParameters_H_
#define model_ClusterDynamicsParameters_H_

#include <TerminalColors.h>
//#include <DislocatedMaterialBase.h>
#include <Polycrystal.h>
#include <DislocationDynamicsBase.h>
//#include <TrialBase.h>
//#include <EvalFunction.h>
//#include <EvalExpression.h>
//#include <DislocationStress.h>

namespace model
{

template<int dim>
struct ClusterDynamicsParameters
{
    static constexpr int mSize=4;     // e.g. Cv, Ci, C2i, C3i
    static constexpr int iSize=0;  // e.g. Nc, Na1, Na2, Na3, cv, ca1, ca2, ca3

    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;

    // Materials parameters
    const double kB;
    const double T; // Studied temperature in K
    const double omega; // Atomic volume
    const double b; // Burgers vector
    const double G0; // dpa/s, Dose rate

    // Mobile Species
    const Eigen::Array<double,1,mSize> msVector;
    const Eigen::Array<double,1,mSize> msRelRelaxVol;
    const Eigen::Array<double,1,mSize> msEf;
    const Eigen::Matrix<double,mSize,dim*(dim+1)/2> msEm; // for each species (row) em11 em12 em13 em22 em23 em33
    const Eigen::Matrix<double,mSize,dim*(dim+1)/2> msD0; // diffusion pre-exponential coeff. for each species (row) D011 D012 D013 D022 D023 D033
    const std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> D;   // map<grain_ID, vector of diffusion coeff for each species>
    const std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> invD;
    const std::map<size_t,Eigen::Array<double,mSize,1>> detD;
    const Eigen::Matrix<double,1,mSize> msCascadeFractions;
    const double msSurvivingEfficiency; // Surviving effciency
//    const Eigen::Matrix<double,1,mSize> G; // dpa/s, Effective dose rate
    const Eigen::Matrix<double,mSize,1> G; // dpa/s, Effective dose rate

    // First-order reaction
    const Eigen::Array<double,1,mSize> otherSinks;
//    const Eigen::Array<double,1,iSize/2> dislocationSinks;
//    const Eigen::Array<double,1,iSize> initloopSinks;
    const Eigen::Matrix<double,mSize,mSize> R1;
    const Eigen::Matrix<double,mSize,mSize> R1cd;
    // Second-order reaction
    const std::map<std::pair<int,int>,double> reactionMap;
    const std::vector<Eigen::Matrix<double,mSize,mSize>> R2;


    // Immobile Species
    const Eigen::Array<double,1,iSize/2> immobileSpeciesVector;
    const Eigen::Array<double,1,iSize/2> immobileSpeciesRelRelaxVol;
    const Eigen::Matrix<double,dim,iSize/2> immobileSpeciesBurgers;
    const Eigen::Array<double,1,iSize/2> immobileSpeciesBurgersMagnitude;
    const double a_bp; // bi-pyramid sink strength coefficient
    const double delVPyramid;
    const double w0;
    const double n_s;
    const Eigen::Array<double,1,mSize> Eb;

    // Irradiation Production
    const double evc; // Vacancy cluster generation efficiency
    const double Nvmax; // Satuatrion number density for vacancy loops in m^-3
    const Eigen::Array<double,1,iSize/2> nmin; // Critical size for <c> pyramid -> loop
    const Eigen::Array<double,1,iSize/2> nmax;
    const Eigen::Array<double,1,iSize/2> r_min; // minimal loop sizes

    // Reaction map (types: parameters)
    const bool computeReactions;
    const int use0DsinkStrength;

    // Bias factors
    const Eigen::Array<double,1,dim> Zv;
    const Eigen::Array<double,1,dim> Zi;
    const Eigen::Array<double,1,mSize> ZVec;
    const Eigen::Array<double,1,dim> rc_il;


    // Discrete Loop Generation
    const double discreteDistanceFactor;

    
    ClusterDynamicsParameters(const DislocationDynamicsBase<dim>& ddBase /*const Polycrystal<dim>& poly*/);
    std::map<std::pair<int,int>,double> getMap(const Eigen::Array<double,mSize*(mSize+1)/2,3> matrix_in) const;
    Eigen::Matrix<double,mSize,mSize> getR1() const;
    std::vector<Eigen::Matrix<double,mSize,mSize>> getR2() const;
    Eigen::Array<double,1,iSize/2> getImmobileSpeciesBurgersMagnitude(const std::map<size_t,Grain<dim>>& grains) const;
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> getD(const std::map<size_t,Grain<dim>>& grains) const;
    std::vector<Eigen::Matrix<double,dim,dim>> getDlocal() const;
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> getInvD() const;
    std::map<size_t,Eigen::Array<double,mSize,1>> getDetD() const;
    Eigen::Array<double,1,iSize> getInitLoopSinks(const Eigen::Array<double,1,iSize> initloopSinks_SI, const double b_SI) const;
    Eigen::Array<double,1,mSize> equilibriumMobileConcentration(const double& stressTrace) const;
    Eigen::Array<double,1,mSize> dislocationMobileConcentration(const VectorDim& b,const VectorDim& t,const VectorDim& fPK,const MatrixDim& stress) const;
    Eigen::Array<double,1,mSize> boundaryMobileConcentration(const double& stressTrace,const double& normalTraction) const;
    Eigen::Array<double,1,iSize/2> sigmoid(const Eigen::Array<double,1,iSize/2>& n) const;
    Eigen::Array<double,1,iSize/2> rpyr(const Eigen::Array<double,1,iSize/2>& n) const;
    Eigen::Array<double,1,iSize/2> rloop(const Eigen::Array<double,1,iSize/2>& n) const;
    Eigen::Array<double,1,iSize/2> sigmoidalVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const;
    Eigen::Array<double,1,iSize/2> clusterRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const;
    Eigen::Array<double,1,iSize/2> clusterDensity(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const;
    Eigen::Array<double,dim,dim> sigmoidalMatrixInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& index) const;
    Eigen::Array<double,1,iSize/2> sigmoidalPlotVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const;
    Eigen::Array<double,1,iSize/2> clusterPlotRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const;
};

}
#endif
