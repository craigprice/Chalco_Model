
#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <string>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "ClassicalSpin3D.h"
class ClassicalSpin3D;

const static unsigned int PR = 14;//string precision
const static unsigned int NUMVECTORS = 5;
const static unsigned int SPINDIMEN = 3;//Components to be manipulated at each site
const static unsigned int POWERSOP = 4;
const static unsigned int TYPESOP = 4;//Operations done on the spin over the entire lattice
const static unsigned int NUMCORRSPINSBIG = 210*2;
const static unsigned int FTTOUPDATES = 4000;
const static unsigned int CORRTOUPDATES = 20;
const static unsigned int NUMSUBLATTICES = 4;

struct FourierTransformOfSpin{
    double ReXComponent;
    double ReYComponent;
    double ReZComponent;
    double ImXComponent;
    double ImYComponent;
    double ImZComponent;
    double magnitude;
};

struct MCParameters{
    double  j1;
    double  j2;
    double  j3;
    double  k1;
    double  k2;
    double  k3;
    
    double  cubicD;
    
    bool    isBField;
    double  bField_x;
    double  bField_y;
    double  bField_z;
    double  bFieldMag;
    
    int     cellsA;
    int     cellsB;
    int     cellsC;
    
    double  KbT;
    double  estimatedTc;
    int     numSweepsToPerformTotal;
    
    double  bField_xNorm;
    double  bField_yNorm;
    double  bField_zNorm;
    std::string path;
};

struct MetropolisDetails{
    double  flipAttempts;
    double  successfulFlips;
    double  range;
    double  rangeOld;
    int     numConfigsDone;
    bool    isThermal;
    int     numSweepsPerformedToCheckIfThermalized;
    int     numSweepsPerformed;
    int     numSweepsUsedToThermalize;
    clock_t timeOfInitialization;
};


const static double B_VECTOR_1[2] = {
    2.0 * PI * 1,
    2.0 * PI * 0};

const static double B_VECTOR_2[2] = {
    2.0 * PI * 0,
    2.0 * PI * 1};

const static double A_VECTOR_1[2] = {
    1,
    0};

const static double A_VECTOR_2[2] = {
    0,
    1};



//Running sums of the Order Parameters in powers of 1, 2, 3, 4.
struct LatticeCharacteristics{
    
    double sumE[POWERSOP];//components are pow ^ 1/2/3/4
    double sumPX[POWERSOP];
    double sumPY[POWERSOP];
    double sumQ[POWERSOP];
    
    //Sum over the Characteristics of the Order Parameters.
    //These are Scalar order parameters.
    double sumOP[NUMVECTORS][TYPESOP][POWERSOP];
    //M[4][4];
    //N[4][4];
    //S1[4][4];
    //S2[4][4];
    //S3[4][4];
    //S[4][4];
    //N3D[4][4];
    //ZZ1[4][4];
    //ZZ2[4][4];
    //ZZ3[4][4];
    //ZZ[4][4];
    //components of 2nd array: M, Mb1, Mb2, Mb, M_x, M_y, M_z, M_xy, M_yz, M_zx
    //Mb = sqrt(Mb1*Mb1 + Mb2*Mb2)
    //components of 3rd array: pow(OP, (1/2/3/4))
    
    double sumOPCorrFunc[NUMVECTORS][TYPESOP][NUMCORRSPINSBIG];
    double sumECorrFunc[NUMCORRSPINSBIG];
    double sumPXCorrFunc[NUMCORRSPINSBIG];
    double sumPYCorrFunc[NUMCORRSPINSBIG];
    double sumQCorrFunc[NUMCORRSPINSBIG];
    
};


struct OrderParameters{
    
    //Total E, E not per site
    double E;
    
    double PX;
    double PY;
    double Q;
    
    //All Order parameters are per site
    
    //Characteristics of the Order Parameters. These are Scalar order parameters.
    double charOP[NUMVECTORS][TYPESOP];
    //M[4];
    //N[4];
    //S1[4];
    //S2[4];
    //S3[4];
    //S[4];
    //N3D[4];
    //ZZ1[4];
    //ZZ2[4];
    //ZZ3[4];
    //ZZ[4];
    //components of 2nd array: M, Mb1, Mb2, Mb.
    //Mb = sqrt(Mb1*Mb1 + Mb2*Mb2)
    
    double corrFunc[NUMVECTORS][TYPESOP][NUMCORRSPINSBIG];
    double ECorrFunc[NUMCORRSPINSBIG];
    double PXCorrFunc[NUMCORRSPINSBIG];
    double PYCorrFunc[NUMCORRSPINSBIG];
    double QCorrFunc[NUMCORRSPINSBIG];
    
};

class MonteCarlo {
public:
    MonteCarlo();
    MonteCarlo(bool reThermalize, std::string filename);
    MonteCarlo(MCParameters input_parameters);
    
    /**
     *This updates the MC's information so that the data about the MC is saved.
     */
    void addToMagStats();
    
    void sweepSnapshot(int numSweeps, int numIndep, double hours);
    
    void sweep(int numIndep_, double durationOfSingleRun);
    
    void printSnapshot();
    
    void updateOrderParameters();
    
    void updateCorrelationFunction();
    void updateFourierTransformOnRecipLattice();
    
    /**
     *This prints the average of the previous "configurations" number of
     *statistics from the metropolis sweep. Then it resets the values back
     *to zero.
     */
    void finalPrint();
    
    /**
     *Main working part of the MC. This performs the spin flipping (and
     *checking to see if it accepts) for each one of atoms in the lattice once.
     */
    void metropolisSweep();
    void init();
    void setLCToZero();
    void setOPToZero();
    void setMCPToZero();
    void setMDToZero();
    void setParamNames();
    void setRecipToZero();
    bool thermalize(double KbT_, int numSweepsToDo, double durationOfSingleRun);
    bool isThermalized (double durationOfSingleRun);
    
    
    std::string outputMag;
    MCParameters MCP;
    MetropolisDetails MD;
private:
    //MonteCarlo(const MonteCarlo& MC);
    //MonteCarlo & operator = (const MonteCarlo& MC);
    
    
    LatticeCharacteristics LC;
    OrderParameters OP;
    std::vector< std::vector< std::vector< std::vector<ClassicalSpin3D> > > > lattice;
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > recipLattice;
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > sumOverRecipLattice;
    
    void adjustRange();
    
    std::string toString();
    
    double dotProd(const double arr1[], const double arr2[]);
    
    int findParameterInt(int currentValue, std::string str, std::string str2) const;
    
    double findParameterDbl(double currentValue,
                            std::string search,
                            std::string str) const;
    
    void setSpinComponent(std::string component_);
    void setRecipLatComponent(std::string component_);
    void setSumOverRecipLatComponent(std::string component_);
    
    
    
    void initializeLattice();
    
    /**
     *This returns the energy of the argument's atom interacting with all of
     *it's neighboring atoms. This includes next nearest neighbors as well if
     *the lattice is keeping track of those.
     */
    double getLocalEnergy(const ClassicalSpin3D& a, bool fromSpinFlip) const;
    double getLocalEnergyRotated(const ClassicalSpin3D& a) const;
    
    
    void rotateAllSpins(double phi, double theta);
    
    const ClassicalSpin3D& getSpinInTheXDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheYDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheMinusXDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheMinusYDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheUpDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheDownDir(const ClassicalSpin3D& a) const;
    
};



inline double MonteCarlo::findParameterDbl(double currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atof(l.substr(str.size()).c_str())<< std::endl;
        return atof(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
}

inline int MonteCarlo::findParameterInt(int currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atoi(l.substr(str.size()).c_str())<< std::endl;
        return atoi(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
}

inline double MonteCarlo::dotProd(const double arr1[], const double arr2[]){
    return (arr1[0] * arr2[0] +
            arr1[1] * arr2[1] +
            arr1[2] * arr2[2]);
}

#endif
