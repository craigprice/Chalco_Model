
/**
 KH_Honeycomb.h
 Purpose: Contains information to implement the Kitaev Heisenberg hamiltonian
 on the honeycomb lattice.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/07/10
 */

#ifndef KH_HONEYCOMB_H
#define KH_HONEYCOMB_H
#include <string>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <fstream>

//#define CUBIC_HAMILTONIAN
//#define ROTATED_HAMILTONIAN
//#define ALPHA_HAMILTONIAN
//#define EXTENDED_HAMILTONIAN
//#define BFIELD_HAMILTONIAN
#define J23K23_HAMILTONIAN

//Number of Vector order parameters
const static int NUMVECTOROP = 11;

//Number of Scalar order parameters
const static int NUMSCALAROP = 1;

//dimension of spin. Components to be manipulated at each site.
const static int SPINDIMEN = 3;

//Operations done on the vector OP over the entire lattice.
//e.g. Order Parameter, X component of the OP, Y component of the OP...
const static int TYPESOP = 10;

//Number of powers of the order parameter kept
const static int POWERSOP = 4;

//Number of configurations performed before calculating
//the Fourier Transform.
const static int UPDATESTOFT = 20*1000; //1% of the metropolis sweep is 2000

//Number of configurations performed before calculating the
//correlation function
const static int UPDATESTOCORR = 4;//This doesn't actually take that long

//Number of sublattices
const static uint_fast8_t NUMSUBLATTICES = 2;

//Number of Nearest neighbors whose position is kept in lookup table
const static uint_fast8_t NUMNEIGHBORS = 12;

//Real space lattice vector 1st
const double A_VECTOR_1[2] = {
    1,
    0};

//Real space lattice vecotor 2nd
const double A_VECTOR_2[2] = {
    0.5,
    sqrt(3.0) / 2.0};


//Reciprocal space lattice vector 1st
const double B_VECTOR_1[2] = {
    (4 * PI / sqrt(3.0)) * sqrt(3.0) / 2.0,
    (4 * PI / sqrt(3.0)) * (-1) / 2.0};

//Reciprocal space lattice vector 2nd
const double B_VECTOR_2[2] = {
    0,
    (4 * PI / sqrt(3.0)) * 1};



//2BZ lattice vector - for a triangular lattice with 1 = the width of the hexagon
//or the altitude of two triangles
const double A2_VECTOR_1[2] = {
    1.0 / 2.0,
    1.0 / (2.0 * sqrt(3.0) )};

const double A2_VECTOR_2[2] = {
    0,
    1.0 / sqrt(3.0)};

const double B2_VECTOR_1[2] = {
    (4 * PI * sqrt(3.0)) * (1 / sqrt(3.0) ),
    0};

const double B2_VECTOR_2[2] = {
    (4 * PI * sqrt(3.0)) * (-1.0 / (2.0 * sqrt(3.0) ) ),
    (4 * PI * sqrt(3.0)) * (1.0 / 2.0) };


//String label of each of the vector order parameters
std::string vectorOrderParameterName[NUMVECTOROP] = {
    "Magnetization",
    "Neel",
    "Stripy-1",
    "Stripy-2",
    "Stripy-3",
    "Stripy",
    "Neel-3D",
    "ZigZag-1",
    "ZigZag-2",
    "ZigZag-3",
    "ZigZag"
};

//String label of each of the scalar order parameters
std::string scalarOrderParameterName[NUMSCALAROP] = {
    "Energy"
};

//Small string label for each component of the order parameters
std::string vectorOrderParameterNameSmall[NUMVECTOROP] = {
    "M",
    "N",
    "S1",
    "S2",
    "S3",
    "S",
    "N3D",
    "ZZ1",
    "ZZ2",
    "ZZ3",
    "ZZ"
};

//Small string label for each component of the order parameters
std::string scalarOrderParameterNameSmall[NUMSCALAROP] = {
    "E"
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
    "Cubic Order Parameter, b1",
    "Cubic Order Parameter, b2",
    "Magnitude of Cubic Order Parameter, |b|",
    "Order Parameter, m1",
    "Order Parameter, m2",
    "Magnitude of Order Parameter, |m|",
    "X Component of the Spin",
    "Y Component of the Spin",
    "Z Component of the Spin"
};

//sublattice structure is:

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
//[0]: first NN in x direction
//[1]: first NN in y direction
//[2]: first NN in z direction
//[3]: second NN, going first x direction, then y
//[4]: second NN, going first y direction, then x
//[5]: second NN, going first x direction, then z
//[6]: second NN, going first z direction, then x
//[7]: second NN, going first y direction, then z
//[8]: second NN, going first z direction, then y
//[9]:  third NN, going first x direction, then y, then z
//[10]: third NN, going first y direction, then z, then x
//[11]: third NN, going first z direction, then x, then y

#endif
