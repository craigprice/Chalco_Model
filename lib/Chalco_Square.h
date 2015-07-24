
/**
 MonteCarlo.h
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of Chalcogenides. Specifically, exploring the
 magnetic phases at finite temperature.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/07/10
 */

#ifndef CHALCO_SQUARE_H
#define CHALCO_SQUARE_H
#include <string>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>

//Number of Vector order parameters
const static int NUMVECTOROP = 5;

//Number of Scalar order parameters
const static int NUMSCALAROP = 5;

//dimension of spin. Components to be manipulated at each site.
const static int SPINDIMEN = 3;

//Operations done on the vector OP over the entire lattice.
//e.g. Order Parameter, X component of the OP, Y component of the OP...
const static int TYPESOP = 4;

//Number of powers of the order parameter kept
const static int POWERSOP = 4;

//Number of configurations performed before calculating
//the Fourier Transform.
const static int UPDATESTOFT = 20*1000; //1% of the metropolis sweep is 2000

//Number of configurations performed before calculating the
//correlation function
const static int UPDATESTOCORR = 4;//This doesn't actually take that long

//Number of sublattices
const static uint_fast8_t NUMSUBLATTICES = 4;

//Number of Nearest neighbors whose position is kept in lookup table
const static uint_fast8_t NUMNEIGHBORS = 12;

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
std::string vectorOrderParameterName[NUMVECTOROP] = {
    "Magnetization",
    "Neel",
    "Stripy-x",
    "Stripy-y",
    "Stripy"
};

//String label of each of the scalar order parameters
std::string scalarOrderParameterName[NUMSCALAROP] = {
    "Energy",
    "Product-x",
    "Product-y",
    "Nematic-Q",
    "Dimer"
};

//Small string label for each component of the order parameters
std::string vectorOrderParameterNameSmall[NUMVECTOROP] = {
    "M",
    "N",
    "SX",
    "SY",
    "S"
};

//Small string label for each component of the order parameters
std::string scalarOrderParameterNameSmall[NUMSCALAROP] = {
    "E",
    "PX",
    "PY",
    "Q",
    "Dim"
};

//String label for each component of the order parameters
std::string spinComponent[SPINDIMEN] = {
    "X direction",
    "Y direction",
    "Z direction",
};

//A small string label for each component of the order parameters - used
//when printing many repititions of the order parameters - like
//histogramming the data
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

//String label for each of the operations saved that are performed on the
//vector order parameters.
std::string typeOP[TYPESOP] = {
    "Order Parameter",
    "X Component of the Spin",
    "Y Component of the Spin",
    "Z Component of the Spin",
};

//sublattice structure is:
//0,1,0,1
//2,3,2,3
//0,1,0,1

//Input parameters for Monte Carlo simulation. Used in the program as "MCP"
struct MCParameters{
    double  j1;
    double  j2;
    double  j3;
    double  k1;
    double  k2;
    double  k3;
    
    double param1;
    
    double  cubicD;
    
    bool    isBField;
    double  bField_x;
    double  bField_y;
    double  bField_z;
    double  bFieldMag;
    
    uint_fast8_t     cellsA;
    uint_fast8_t     cellsB;
    uint_fast8_t     cellsC;
    
    double  KbT;
    double  estimatedTc;
    int     numSweepsToPerformTotal;
    
    double  bField_xNorm;
    double  bField_yNorm;
    double  bField_zNorm;
    std::string path;
    
    unsigned long monte_carlo_seed;
};

//Nearest Neighbor arrays. neighbors[c][a][b][s] is a vector of the
//coordinates of the neighbors of the spin at location: lattice[c][a][b][s]
//The order of neighbors are:
//[0]: first NN in +x direction
//[1]: first NN in -x direction
//[2]: first NN in +y direction
//[3]: first NN in -y direction
//[4]: second NN, going first in the +x direction, then the +y direction
//[5]: second NN, going first in the -x direction, then the +y direction
//[6]: second NN, going first in the -x direction, then the -y direction
//[7]: second NN, going first in the +x direction, then the -y direction
//[8]:  third NN, going first in the +x direction, then the +x direction
//[9]:  third NN, going first in the -x direction, then the -x direction
//[10]: third NN, going first in the +y direction, then the +y direction
//[11]: third NN, going first in the -y direction, then the -y direction


#endif
