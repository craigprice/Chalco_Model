
/**
 MonteCarlo.h
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of Chalcogenides. Specifically, exploring the
 magnetic phases at finite temperature.
 
 @author Craig Price
 @version 1.0 2015/03/04
 */

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



//The object that contains all the information about the simulation and contains
//the spin objects.
class MonteCarlo {
public:
    
    /**************************************************************************/
    //Global Variables
    /**************************************************************************/
    
    
    const double PI = 3.14159265358979323846;
    
    //Numerical precision for printing strings
    const static unsigned int PR = 14;
    
    //Number of Vector order parameters
    const static unsigned int NUMVECTORS = 5;
    
    //dimension of spin. Components to be manipulated at each site.
    const static unsigned int SPINDIMEN = 3;
    
    //Number of powers of the order parameter kept e.g. to the fourth power
    const static unsigned int POWERSOP = 4;
    
    //Operations done on the vector OP over the entire lattice.
    //e.g. Order Parameter, X component of the OP, Y component of the OP...
    const static unsigned int TYPESOP = 4;
    
    //A large number that is used to initialize arrays that hold information about
    //correlations.
    const static unsigned int NUMCORRSPINSBIG = 210*2;
    
    //Number of configurations performed before calculating the Fourie Transform.
    const static unsigned int FTTOUPDATES = 4000*10;
    
    //Number of configurations performed before calculating the correlation function
    const static unsigned int CORRTOUPDATES = 20;
    
    //Number of sublattices
    const static unsigned int NUMSUBLATTICES = 4;
    
    //Real space lattice vector 1st
    const double A_VECTOR_1[2] = {
        1,
        0};
    
    //Real space lattice vecotor 2nd
    const double A_VECTOR_2[2] = {
        0,
        1};
    
    //Reciprocal space lattice vector 1st
    const double B_VECTOR_1[2] = {
        2.0 * PI * 1,
        2.0 * PI * 0};
    
    //Reciprocal space lattice vector 2nd
    const double B_VECTOR_2[2] = {
        2.0 * PI * 0,
        2.0 * PI * 1};
    
    //String label of each of the vector order parameters
    std::string orderParameter[NUMVECTORS] = {
        "Magnetization",
        "Neel",
        "Stripy-x",
        "Stripy-y",
        "Stripy"
    };
    
    //Small string label for each component of the order parameters
    std::string orderParameterSmall[NUMVECTORS] = {
        "M",
        "N",
        "SX",
        "SY",
        "S"
    };
    
    //String label for each component of the order parameters
    std::string spinComponent[SPINDIMEN] = {
        "X direction",
        "Y direction",
        "Z direction",
    };
    
    //A small string label for each component of the order parameters - used when
    //printing many repititions of the order parameters - like histogramming the
    //data
    std::string spinComponentSmall[SPINDIMEN] = {
        "X",
        "Y",
        "Z",
    };
    
    //String label for each power of the order parameters saved.
    std::string powerOP[POWERSOP] = {
        " to the 1st power",
        " to the 2nd power",
        " to the 3rd power",
        " to the 4th power"
    };
    
    //String label for each of the operations saved that are performed on the vector
    //order parameters.
    std::string typeOP[TYPESOP] = {
        "Order Parameter",
        "X Component of the Spin",
        "Y Component of the Spin",
        "Z Component of the Spin",
    };
    
    
    
    /**************************************************************************/
    //Groups of constants for Monte Carlo
    /**************************************************************************/
    
    //Arrays to contain the string names for the sum over order paramets
    std::string sumOPNames       [NUMVECTORS][TYPESOP][POWERSOP] = {""};
    std::string charOPNames      [NUMVECTORS][TYPESOP] = {""};
    
    //The information about the FT at a site in reciprocal space. Analoge to the
    //"ClassicalSpin3D" object.
    struct FourierTransformOfSpin{
        double ReXComponent;
        double ReYComponent;
        double ReZComponent;
        double ImXComponent;
        double ImYComponent;
        double ImZComponent;
        double magnitude;
    };
    
    //Input parameters for Monte Carlo simulation. Used in the program as "MCP"
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
        
        unsigned long monte_carlo_seed;
    };
    
    //Characteristics of the current simulation. Used in the program as "MD"
    struct MetropolisDetails{
        
        //Number of all spin flip attempts
        double  flipAttempts;
        
        //Number of successful spin flips. Want this to be about 65% of total.
        double  successfulFlips;
        
        //the polar angle theta that defines the polar cap within which the spin
        //flips within. At low temperature, the theta becomes small, and at
        //high temperature, it goes to PI.
        double  range;
        double  rangeOld;
        
        //The number of configurations for which the magnetization was summed over.
        int     numConfigsDone;
        
        //whether the system is thermal or not.
        bool    isThermal;
        
        //A few sweeps are done to see whether our system has converged so that the
        //change of the order parameters is ~< 1%
        int     numSweepsPerformedToCheckIfThermalized;
        
        //Recording how many sweeps we've done on the lattice, but these will not
        //be considered for our sum over magnetic order parameters.
        int     numSweepsUsedToThermalize;
        
        //Total number of monte carlo sweeps over the whole lattice.
        int     numSweepsPerformed;
        
        clock_t timeOfInitialization;
        
        clock_t totalTimeRunning;
    };
    
    //Running sums of the Order Parameters in powers of 1, 2, 3, 4.
    struct LatticeCharacteristics{
        
        double sumE[POWERSOP];//components are pow ^ 1/2/3/4
        double sumPX[POWERSOP];
        double sumPY[POWERSOP];
        double sumQ[POWERSOP];
        
        //Sum over the Characteristics of the Order Parameters.
        double sumOP[NUMVECTORS][TYPESOP][POWERSOP];
        
        //Sum over the correlation function of the the different OP's
        double sumOPCorrFunc[NUMVECTORS][TYPESOP][NUMCORRSPINSBIG];
        double sumECorrFunc[NUMCORRSPINSBIG];
        double sumPXCorrFunc[NUMCORRSPINSBIG];
        double sumPYCorrFunc[NUMCORRSPINSBIG];
        double sumQCorrFunc[NUMCORRSPINSBIG];
        
    };
    
    //Snapshot of the current order parameters
    //All Order parameters are per site
    struct OrderParameters{
        
        //Total E, E not per site
        double E;
        
        //dot product of spin in the x direction
        double PX;
        
        //dot product of spin in the y direction
        double PY;
        
        //Nematic order parameter. <PX-PY>
        double Q;
        
        //Characteristics of the Order Parameters. These are vector order parameters.
        double charOP[NUMVECTORS][TYPESOP];
        
        //correlation functions of each of the order parameters
        double corrFunc[NUMVECTORS][TYPESOP][NUMCORRSPINSBIG];
        double ECorrFunc[NUMCORRSPINSBIG];
        double PXCorrFunc[NUMCORRSPINSBIG];
        double PYCorrFunc[NUMCORRSPINSBIG];
        double QCorrFunc[NUMCORRSPINSBIG];
        
    };
    
    
    /**************************************************************************/
    //All information about the simulation
    /**************************************************************************/
    
    MCParameters MCP;
    MetropolisDetails MD;
    OrderParameters OP;
    LatticeCharacteristics LC;
    
    //The path to the output file.
    std::string outputMag;
    
    /**************************************************************************/
    //Helper functions
    /**************************************************************************/
    
    /**
     Constructing the string that suscintly describes the data that was being
     saved by the program. Used in initialization of the program.
     */
    void setParamNames();
    
    /**
     Sets the spin's components at a particular lattice site.
     
     @param component_ the string that contains the information for the spin.
     */
    void setSpinComponent(std::string component_);
    
    /**
     Sets the Reciprocal lattice site's components at a particular lattice site.
     
     @param component_ the string that contains the information for the site.
     */
    void setRecipLatComponent(std::string component_);
    
    /**
     Sets the sum over the Reciprocal lattice site's components at a particular
     lattice site.
     
     @param component_ the string that contains the information for the site.
     */
    void setSumOverRecipLatComponent(std::string component_);
    
    /**
     Sets the Reciprocal lattice site's components at a particular lattice site.
     
     @return a string that contain's the simulation's order parameter
     information.
     */
    std::string toString();
    
    /**
     Prints out all data of the simulation to the output file.
     */
    void finalPrint();
    
    /**
     Prints out minimal information about the current configuration's order
     parameters.
     */
    void printSnapshot();
    
    /**
     Rotates all spins of the lattice by an amount phi, theta with respect
     to the lattice coordinates.
     
     @param phi angle to the lattice - phi ~ atan2(y,x)
     @param theta angle to the "z" direction of the spin/lattice
     */
    void rotateAllSpins(double phi, double theta);
    
    /**
     Returns the dot product between two 3D arrays
     
     @param arr1 a 3-component array that will be treated as a 3-vector
     @param arr2 another 3-component array to be treated as a 3-vector
     @return the dot product of the given arguments.
     */
    double dotProd(const double arr1[], const double arr2[]);
    
    /**
     If the argument string str contains str2, then it returns the int cast
     of the end of the str2. For example, if str = "M: 2.3" then this would
     return (int) 2
     
     @param str the string to search within
     @param str2 the string to find.
     @param currentValue return this if str doesn't contain str2
     @return currentValue if str doesn't contain str2 and return the integer
     value that comes after str2 in str if it exists.
     */
    int findParameterInt(int currentValue, std::string str, std::string str2) const;
    
    /**
     If the argument string search contains str, then it returns the int cast
     of the end of the str. For example, if search = "M: 2.3" then this would
     return (int) 2
     
     @param search the string to search within
     @param str the string to find.
     @param currentValue return this if search doesn't contain str
     @return currentValue if search doesn't contain str and return the integer
     value that comes after str in search if it exists.
     */
    double findParameterDbl(double currentValue,
                            std::string search,
                            std::string str) const;
    
    
    /**************************************************************************/
    //Functions for initialization
    /**************************************************************************/
    
    /**
     Initializes all the arrays that store the running sum over the order
     parameters to zero. This is setting the entire global variable, "LC", which
     is the "LatticeCharacteristics" instantiation to zero/default.
     */
    void setLCToZero();
    
    /**
     Initializes all the arrays that store the order parameters to zero. This is
     setting the entire global variable, "OP", which
     is the "OrderParameters" instantiation to zero/default.
     */
    void setOPToZero();
    
    /**
     Initializes all the Monte Carlo input parameters to zero/default. This is
     setting the entire global variable, "MCP", which
     is the "MCParameters" instantiation to zero/default.
     */
    void setMCPToZero();
    
    /**
     Initializes all the Monte Carlo information about the simulation to zero
     or defaults. This is
     setting the entire global variable, "MD", which
     is the "MetropolisDetails" instantiation to zero/default.
     */
    void setMDToZero();
    
    /**
     Initializes all the Monte Carlo information about the reciprocal lattice
     and the running-sum-over-the-reciprocal lattice to zero/defaults.
     */
    void setRecipToZero();
    
    /**
     Constructs the array (implemented as a c++ vector) that contains all the
     spins and the FT of the spins. Then initiallizes all the spin information
     to zero.
     */
    void initializeLattice();
    
    //sublattice structure is:
    //2,3
    //0,1
    
    /**
     the lattice-geometry path from a current spin to the next spin to the x.
     
     @param a the spin to find the nearest neighbor from.
     @return the spin in the direction of the nearest neighbor.
     */
    const ClassicalSpin3D& getSpinInTheXDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheYDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheMinusXDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheMinusYDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheUpDir(const ClassicalSpin3D& a) const;
    
    const ClassicalSpin3D& getSpinInTheDownDir(const ClassicalSpin3D& a) const;
    
    void setNearestNeighbors();
    
    /**
     Create and initialize all variables used in the simulation.
     */
    void init();
    
    /**
     Default constructor for the simulation object.
     */
    MonteCarlo();
    
    /**
     Monte Carlo constructor that begins a new simulation.
     
     @param input_parameters is the instantiation of the global variable,
     "MCParameters" that are all the input parameters for the simulation to run.
     */
    MonteCarlo(MCParameters input_parameters);
    
    /**
     Monte Carlo constructor that continues a simulation from a previously saved
     output file.
     
     @param reThermalize boolean for whether you want to raise the temperature
     of the simulation (that you read in from a previously saved output file)
     and clear all simulation details in order to thermalize the simulation,
     and re-anneal the spins of the lattice. This effectively starts a brand new
     simulation.
     @param filename the path to the output file to read in the characteristics
     of the old simulation and continue the simulation from it.
     */
    MonteCarlo(bool reThermalize, std::string filename);
    
    /**************************************************************************/
    //Thermalization
    /**************************************************************************/
    
    /**
     From the total number of spin flips, and the number of successful spin
     flips, this function reduces the range of the allowed polar cap space
     to flip the spin into in order to on average have 65% successful spin flips
     */
    void adjustRange();
    
    /**
     Run sweeps over the lattice at a specific temperature until the order
     parameters converge to some equilibrium value.
     
     @param KbT_ the temperature at which we want to thermalize.
     @param numSweepsToDo The minimum number of sweeps to perform before
     declaring the system is thermalized.
     @param durationOfSingleRun how long this function should continue in
     real time before exiting and saving the characteristics.
     */
    bool thermalize(double KbT_, int numSweepsToDo, double durationOfSingleRun);
    
    /**
     query the system to see if the order parameters have converged.
     
     @param durationOfSingleRun how long this function should continue in real time
     before exiting and saving the characteristics.
     */
    bool isThermalized (double durationOfSingleRun);
    
    
    /**************************************************************************/
    //Updating order parameters
    /**************************************************************************/
    
    /**
     update the global variable, "OP", which holds a snapshot of all the order
     parameters of the lattice
     */
    void updateOrderParameters();
    
    /**
     update the global variables of the correlation of the order parameters.
     */
    void updateCorrelationFunction();
    
    /**
     calculates the fourier transform of the lattice and updates the current
     values of the FT in the latticeFT.
     */
    void updateFourierTransformOnRecipLattice();
    
    /**
     *This updates the MC's information so that the data about the MC is saved.
     */
    void addToMagStats();
    
    /**************************************************************************/
    //Core Functions
    /**************************************************************************/
    
    
    /**
     This returns the energy of the argument's atom interacting with all of
     it's neighboring atoms. This includes next nearest neighbors as well if
     the lattice is keeping track of those.
     
     @param a the spin to calculate the local energy around.
     @param fromSpinFlip true if calculating the total energy of the lattice.
     In that case, one needs to divide by two on all bonds to avoid double
     counting.
     @return energy from all interactions with the spin a
     */
    double getLocalEnergy(const ClassicalSpin3D& a, bool fromSpinFlip) const;
    
    /**
     Performs the spin flipping over the each lattice site once.
     It for each flip it checks to see if it accepts the new configuration or
     not.
     */
    void metropolisSweep();
    
    /**
     Calls metropolisSweep continuously until it reaches the desired number, or
     runs out of time.
     
     @param numIndep_ the number of sweeps to perform before recording the
     order parameters from the current configuration and updating the running
     average.
     @param durationOfSingleRun how long this job will run before automatically
     ending itself and saving the output to a output file.
     */
    void sweep(int numIndep_, double durationOfSingleRun);
    
    
    /**
     Loops over Monte Carlo sweeps, printing out information for each sweep.
     
     @param numSweeps The number of sweeps that you want to run.
     @param numIndep The number of Monte Carlo steps to take between calculating
     and saving the order parameter, magnetization data.
     @param hours The number of hours to run before exiting the program.
     */
    void sweepSnapshot(int numSweeps, int numIndep, double hours);
    
    //The array that holds all the spins in the simulation
    std::vector< std::vector< std::vector< std::vector<ClassicalSpin3D> > > > lattice;
    
    //Nearest Neighbor arrays
    //The array that holds the pointer to the first nearesst neighbor to the x.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN_x;
    
    //The array that holds the pointer to the first nearesst neighbor to the -x.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN_mx;
    
    //The array that holds the pointer to the first nearesst neighbor to the y.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN_y;
    
    //The array that holds the pointer to the first nearesst neighbor to the -y.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN_my;
    
    //The array that holds the pointer to the second nearesst neighbor along the x then y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN2_xy;
    
    //The array that holds the pointer to the second nearesst neighbor along the -x then y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN2_mxy;
    
    //The array that holds the pointer to the second nearesst neighbor along the -x then -y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN2_mxmy;
    
    //The array that holds the pointer to the second nearesst neighbor along the x then -y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN2_xmy;
    
    //The array that holds the pointer to the thrid nearesst neighbor along the x then x bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN3_xx;
    
    //The array that holds the pointer to the thrid nearesst neighbor along the -x then -x bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN3_mxmx;
    
    //The array that holds the pointer to the thrid nearesst neighbor along the y then y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN3_yy;
    
    //The array that holds the pointer to the thrid nearesst neighbor along the -y then -y bonds.
    std::vector< std::vector< std::vector< std::vector<const ClassicalSpin3D*> > > > NN3_mymy;
    
    //The array that holds the current fourier transform of the lattice.
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > recipLattice;
    
    //The array that holds the running sum over the fourier transform of the
    //lattice.
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > sumOverRecipLattice;
    
    
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



#endif
