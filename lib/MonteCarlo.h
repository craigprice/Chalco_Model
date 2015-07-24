
/**
 MonteCarlo.h
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of magnetic models in a regular geometry at 
 finite temperature.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/07/10
 */

#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <string>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>
const static double PI = 3.14159265358979323846;
#include "Chalco_Square.h"
//#include "KH_Honeycomb.h"

//uint _fast8_t does not work with stringstream - or at least cout.
//Must cast to int
struct SpinCoordinates{
    uint_fast8_t c; //position of the spin along the "c" axis
    uint_fast8_t a; //position of the spin along the "a" axis
    uint_fast8_t b; //position of the spin along the "b" axis
    uint_fast8_t s; //sublatice that the spin belongs within
};

/**************************************************************************/
//Global Variables
/**************************************************************************/

//Numerical precision for printing strings
const static int PR = 14;

//A large number that is used to initialize arrays that hold information
//about correlations.
const static int NUMCORRSPINSBIG = 204*2;

//minimum number of sweeps done for thermalization
const static int MINTHERMALSWEEPS = 100;// * 1000;

//factor that is multiplied to the temperature on the thermalization ramp
const static double THERMALIZATIONRATE = 0.001;

//time that a job is allowed to run before quitting
const static double LIMITONHOURSOFJOBRUN = 0.1*1.0 / 60.0;//20.0/60.0;

//Number of sweeps to perform between adding to the running sum.
int NUMSWEEPSBETWEENCONFIGS = 6;

/**************************************************************************/
//Groups of constants for Monte Carlo
/**************************************************************************/

//Arrays to contain the string names for the sum over order paramets
std::string sumOverVectorOPNames    [NUMVECTOROP][TYPESOP][POWERSOP] = {""};
std::string sumOverScalarOPNames    [NUMSCALAROP][POWERSOP] = {""};
std::string vectorOPAndType         [NUMVECTOROP][TYPESOP] = {""};

struct ClassicalSpin3D{
    double x; //x component of the 3D unit vector that is the spin
    double y; //y component of the 3D unit vector that is the spin
    double z; //z component of the 3D unit vector that is the spin
};

//The information about the FT at a site in reciprocal space. Analog to the
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

//Characteristics of the current simulation. Used in the program as "MD"
struct MetropolisDetails{
    
    //Number of all spin flip attempts
    long int  flipAttempts;
    
    //Number of successful spin flips. Want this to be about 65% of
    //flipAttempts.
    long int  successfulFlips;
    
    //the polar angle theta that defines the polar cap within which the spin
    //flips. At low temperature, the theta becomes small, and at
    //high temperature, it goes to PI. It is varied to maintain
    //successfulFlips/flipAttempts ~ 65%
    double  range;
    double  rangeOld;
    
    //The number of configurations for which the order parameters were
    //summed over.
    int     numConfigsDone;
    
    //whether the system is thermalized.
    bool    isThermal;
    
    //A few sweeps are done to see whether our system has converged so
    //that the change of the order parameters is < ~1%
    int     numSweepsPerformedToCheckIfThermalized;
    
    /*
     Recording how many sweeps we've done on the lattice while thermalizing.
     These sweeps will not be considered for determining the average order
     parameter.
     */
    int     numSweepsUsedToThermalize;
    
    //Total number of all monte carlo sweeps performed over the lattice.
    int     numSweepsPerformed;
    
    //time of initialization for this particular job.
    time_t timeAtThisJobInitialization;
    
    //Total time running, in seconds. Updated at end of simmulation.cpp
    int totalTime;
};

//Running sums of the Order Parameters in powers of: OP^(1/2/3/4).
struct RunningSumOverOP{
    
    //Sum over the scalar Order Parameters.
    double sumOverScalarOP[NUMSCALAROP][POWERSOP];
    
    //Sum over the vector Order Parameters.
    double sumOverVectorOP[NUMVECTOROP][TYPESOP][POWERSOP];
    
    //Sum over the correlation function of the scalar OP's
    double sumOverScalarOPCorrFunc[NUMSCALAROP][POWERSOP][NUMCORRSPINSBIG];
    
    //Sum over the correlation function of the vector OP's
    double sumOverVectorOPCorrFunc[NUMVECTOROP][TYPESOP][POWERSOP][NUMCORRSPINSBIG];
};

//Snapshot of the current order parameters
//All Order parameters are per site
struct OrderParameters{
    
    //Characteristics of the Order Parameters.
    //These are scalar order parameters.
    double scalarOP[NUMSCALAROP];
    
    //correlation functions of each of the scalar order parameters
    double scalarOPCorrFunc[NUMSCALAROP][NUMCORRSPINSBIG];
    
    //Characteristics of the Order Parameters.
    //These are vector order parameters.
    double vectorOP[NUMVECTOROP][TYPESOP];
    
    //correlation functions of each of the vector order parameters
    double vectorOPCorrFunc[NUMVECTOROP][TYPESOP][NUMCORRSPINSBIG];
};



//The object that contains all the information about the simulation and contains
//the spin objects.
class MonteCarlo {
public:
    
    
    /**************************************************************************/
    //All information about the simulation
    /**************************************************************************/
    
    MCParameters MCP;
    MetropolisDetails MD;
    OrderParameters OP;
    RunningSumOverOP RS;
    int numSites;//cellsA * cellsB * cellsC * NUMSUBLATTICES
    
    //The path to the output file.
    std::string outputFileName;
    
    
    /**************************************************************************/
    // Functions in FileIO.cpp
    /**************************************************************************/
    
    /**
     Constructing the string that succinctly describes the data that was being
     saved by the program. Used in initialization of the program.
     */
    void setParamNames();
    
    /**
     If the argument string str contains str2, then it returns the int cast
     of the end of the str2. For example, if str = "M: 2.3" then this would
     return (int) 2.3
     
     @param str the string to search within
     @param str2 the string to find.
     @param currentValue return this if str doesn't contain str2
     @return currentValue if str doesn't contain str2 and return the integer
     value that comes after str2 in str if it exists.
     */
    int findParameterInt(int currentValue,
                         std::string str,
                         std::string str2) const;
    unsigned long int findParameterULongInt(unsigned long int currentValue,
                                            std::string str,
                                            std::string str2) const;
    long int findParameterLongInt(long int currentValue,
                                  std::string str,
                                  std::string str2) const;
    double findParameterDbl(double currentValue,
                            std::string search,
                            std::string str) const;
    
    /**
     Functions to print the data about the simulation to a string.
     
     @return a string that contain's the simulation's information.
     */
    std::string toString();
    std::string toStringMCP();
    std::string toStringMD();
    std::string toStringRS();
    std::string toStringRSCorr();
    std::string toStringOP();
    std::string toStringSpinSnapshot();
    std::string toStringFourierTransform();
    std::string toStringSumOverFourierTransform();
    
    /**
     Functions to read in the data about the simulation and sets the input to
     the appropriate local variables.
     
     @param readline the input file stream object that the input data should
     be being read from.
     */
    void readStringMCP(std::ifstream &readline);
    void readStringMD(std::ifstream &readline);
    void readStringRS(std::ifstream &readline);
    void readStringRSCorr(std::ifstream &readline);
    void readStringSpinSnapshot(std::ifstream &readline);
    void readStringFourierTransform(std::ifstream &readline);
    void readStringSumOverFourierTransform(std::ifstream &readline);
    
    /**
     Prints out all data of the simulation to the output file. Mostly calls
     the above "toString" functions
     */
    void finalPrint();
    
    /**
     Prints out minimal information about the current configuration's order
     parameters.
     */
    void printSnapshot();
    
    /**************************************************************************/
    //Helper functions
    /**************************************************************************/
    
    /**
     @return true if the job is out of time false if not out of time.
     */
    bool isOutOfTime();
    
    /**
     Not tested
     */
    void setDomainWall();
    
    /**
     Not tested
     */
    void setDomainWallSpins(double spin1X, double spin1Y, double spin1Z,
                            int radius1, double spin2X, double spin2Y,
                            double spin2Z, int radius2);
    
    /**
     Not tested
     */
    bool isDomain(SpinCoordinates sc, int radius);
    
    /**
     Not tested
     */
    double getPolarAngleTB(const double spinX, const double spinY,
                        const double spinZ);
    
    /**
     Not tested
     */
    double getAziAngleAB(const double spinX, const double spinY,
                         const double spinZ);
    
    /**************************************************************************/
    //Functions for initialization
    /**************************************************************************/
    
    /**
     Initializes all the arrays that store the running sum over the order
     parameters to zero. This is setting the entire global variable, "RS", which
     is the "RunningSumOverOP" instantiation to zero/default.
     */
    void setRSToZero();
    
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
    
    /**
     the lattice-geometry path from a current spin to the next spin along the a
     lattice axis, or the X direction.
     
     @param sc the lattice-coordinates of the present spin
     to find the nearest neighbor from.
     @return the spin coordinates of the spin in the direction of
     the nearest neighbor.
     */
    SpinCoordinates getSpinCoordsInTheXDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheYDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheZDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheMinusXDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheMinusYDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheUpDir(SpinCoordinates sc) const;
    
    SpinCoordinates getSpinCoordsInTheDownDir(SpinCoordinates sc) const;
    
    /**
     Sets the lookup table / array that holds the spincoordinates of all the
     saved nearest neighbors to the spin at location, (c,a,b,s).
     */
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
     */
    void thermalize();
    
    /**
     query the system to see if the order parameters have converged. If not,
     continue thermalizing in sets of 2000 until the order paramters converge.
     */
    void checkThermalized();
    
    
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
     *This updates the running sum over the order parameters.
     */
    void updateRunningSumOP();
    
    
    /**************************************************************************/
    //Spin Functions
    /**************************************************************************/
    
    /**
     This makes sure that the value stored in the x and y and z values are
     appropriate for a 3D classical spin with magnitude of 1.
     
     @param c the z-axis coordinate of the spin
     @param a the x-axis coordinate of the spin
     @param b the y-axis coordinate of the spin
     @param s the sublattice of the spin
     */
    void checkSpin(uint_fast8_t c, uint_fast8_t a,
                   uint_fast8_t b, uint_fast8_t s);
    
    /**
     Sets this spin to be oriented randomly in the unit spherical surface.
     */
    void setRandomOrientation(uint_fast8_t c, uint_fast8_t a,
                              uint_fast8_t b, uint_fast8_t s);
    
    /**
     Rotate the spin.
     */
    void specifyRotation(uint_fast8_t c, uint_fast8_t a,
                         uint_fast8_t b, uint_fast8_t s,
                         double rotPhi, double rotTheta);
    
    /**
     flips this spin from it's current orientation to a random new orientation
     Where the value of the polar angular rotation is at most, range.
     
     @param range the polar angle theta that defines the polar spherical cap
     within which lies the possible orientations for the new spin to point in.
     */
    void flip(uint_fast8_t c, uint_fast8_t a,
              uint_fast8_t b, uint_fast8_t s, double range);
    
    /**
     This function rotates the spin about an axis (specified by 'theta' and
     'phi') by an angle of 'angle'.
     
     @param theta
     @param phi
     @param angle
     */
    void rotate(uint_fast8_t c, uint_fast8_t a,
                uint_fast8_t b, uint_fast8_t s,
                double theta_, double phi_, double Beta);
    
    
    /**
     Takes the dot product between two spins.
     
     @param sc the coordinates of the other spin-vector to take the
     dot product with.
     @return the value of the dot product.
     */
    double dotProd(uint_fast8_t c, uint_fast8_t a,
                   uint_fast8_t b, uint_fast8_t s,
                   SpinCoordinates sc1) const;
    
    /**************************************************************************/
    //Core Functions
    /**************************************************************************/
    
    
    /**
     This returns the energy of the argument's atom interacting with all of
     it's neighboring atoms. This includes next nearest neighbors as well if
     the lattice is keeping track of those. This function does *not* use the
     lookup table to determine what its nearest neighbors are.
     
     @param c the z-axis coordinate of the spin
     @param a the x-axis coordinate of the spin
     @param b the y-axis coordinate of the spin
     @param s the sublattice of the spin
     @param fromSpinFlip true if calculating the total energy of the lattice.
     In that case, one needs to divide by two on all bonds to avoid double
     counting.
     @return energy from all interactions with the spin at (c,a,b,s)
     */
    double getLocalEnergyFromCalc(uint_fast8_t c, uint_fast8_t a,
                                  uint_fast8_t b, uint_fast8_t s,
                                  bool fromSpinFlip) const;
    
    /**
     This returns the energy of the argument's atom interacting with all of
     it's neighboring atoms. This includes next nearest neighbors as well if
     the lattice is keeping track of those. This function *does* use the
     lookup table to determine what its nearest neighbors are.
     
     @param c the z-axis coordinate of the spin
     @param a the x-axis coordinate of the spin
     @param b the y-axis coordinate of the spin
     @param s the sublattice of the spin
     @param fromSpinFlip true if calculating the total energy of the lattice.
     In that case, one needs to divide by two on all bonds to avoid double
     counting.
     @return energy from all interactions with the spin at (c,a,b,s)
     */
    double getLocalEnergy(uint_fast8_t c, uint_fast8_t a,
                          uint_fast8_t b, uint_fast8_t s,
                          bool fromSpinFlip) const;
    
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
    void sweep();
    
    
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
    
    //Nearest Neighbor arrays. neighbors[c][a][b][s] is a vector of the
    //coordinates of the neighbors of the spin at location: lattice[c][a][b][s]
    std::vector< std::vector< std::vector< std::vector< std::vector<SpinCoordinates> > > > > neighbors;
    
    //The array that holds the current fourier transform of the lattice.
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > recipLattice;
    
    //The array that holds the running sum over the fourier transform of the
    //lattice.
    std::vector< std::vector< std::vector< std::vector<FourierTransformOfSpin> > > > sumOverRecipLattice;
    
    
};

inline double MonteCarlo::dotProd(uint_fast8_t c, uint_fast8_t a,
                                  uint_fast8_t b, uint_fast8_t s,
                                  SpinCoordinates s1) const{
    return (lattice[s1.c][s1.a][s1.b][s1.s].x * lattice[c][a][b][s].x +
            lattice[s1.c][s1.a][s1.b][s1.s].y * lattice[c][a][b][s].y +
            lattice[s1.c][s1.a][s1.b][s1.s].z * lattice[c][a][b][s].z);
}



#endif
