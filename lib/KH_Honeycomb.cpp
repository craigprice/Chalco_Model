/**
 Honeycomb.cpp
 Purpose: Contains information to implement the Kitaev Heisenberg hamiltonian
 on the honeycomb lattice.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/07/10
 */
#include <ctime>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include "MonteCarlo.h"
#include "KH_Honeycomb.h"

using namespace std;

/******************************************************************************/
//Initialization
/******************************************************************************/

void MonteCarlo::setMCPToZero(){
    cout<<"Begin setMCPToZero()"<<endl;
    MCP.j1          = 0;
    MCP.j2          = 0;
    MCP.j3          = 0;
    MCP.k1          = 0;
    MCP.k2          = 0;
    MCP.k3          = 0;
    MCP.cubicD      = 0;
    MCP.param1      = 0;
    
    MCP.isBField    = false;
    MCP.bField_x    = 0;
    MCP.bField_y    = 0;
    MCP.bField_z    = 0;
    MCP.bFieldMag   = 0;
    
    MCP.cellsA      = 2;
    MCP.cellsB      = 2;
    MCP.cellsC      = 1;
    MCP.KbT         = 0;
    
    MCP.estimatedTc = 0.2;
    MCP.numSweepsToPerformTotal = 1;
    MCP.bField_xNorm = 0;
    MCP.bField_yNorm = 0;
    MCP.bField_zNorm = 0;
}

string MonteCarlo::toStringMCP(){
    string function = "MonteCarlo::toStringMCP()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "KbT: " << MCP.KbT << endl;
    stream << "cellsA: " << (int)MCP.cellsA << endl;
    stream << "cellsB: " << (int)MCP.cellsB << endl;
    stream << "cellsC: " << (int)MCP.cellsC << endl;
    stream << "cubicD: " << MCP.cubicD << endl;
    stream << "j1: " << MCP.j1 << endl;
    stream << "j2: " << MCP.j2 << endl;
    stream << "j3: " << MCP.j3 << endl;
    stream << "k1: " << MCP.k1 << endl;
    stream << "k2: " << MCP.k2 << endl;
    stream << "k3: " << MCP.k3 << endl;
    stream << "param1: " << MCP.param1 << endl;
    stream << "isBField: " << MCP.isBField << endl;
    stream << "bField_x: "<< MCP.bField_x << endl;
    stream << "bField_y: "<< MCP.bField_y << endl;
    stream << "bField_z: "<< MCP.bField_z << endl;
    stream << "bFieldMag: "<< MCP.bFieldMag << endl;
    stream << "estimatedTc: " << MCP.estimatedTc << endl;
    
    stream << "numSweepsToPerformTotal: " <<
    MCP.numSweepsToPerformTotal << endl;
    
    stream << "monte_carlo_seed: " << MCP.monte_carlo_seed << endl;
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringMCP(ifstream &readline){
    string function = "MonteCarlo::readStringMCP(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    string line = "";
    while(readline.good()){
        getline(readline, line);
        //cout << line << endl;
        
        MCP.j1 = findParameterDbl(MCP.j1, line, "j1: ");
        MCP.j2 = findParameterDbl(MCP.j2, line, "j2: ");
        MCP.j3 = findParameterDbl(MCP.j3, line, "j3: ");
        MCP.k1 = findParameterDbl(MCP.k1, line, "k1: ");
        MCP.k2 = findParameterDbl(MCP.k2, line, "k2: ");
        MCP.k3 = findParameterDbl(MCP.k3, line, "k3: ");
        MCP.cubicD = findParameterDbl(MCP.cubicD, line, "cubicD: ");
        MCP.param1 = findParameterDbl(MCP.param1, line, "param1: ");
        
        MCP.isBField = findParameterInt(MCP.isBField, line, "isBField: ");
        MCP.bField_x = findParameterDbl(MCP.bField_x, line, "bField_x: ");
        MCP.bField_y = findParameterDbl(MCP.bField_y, line, "bField_y: ");
        MCP.bField_z = findParameterDbl(MCP.bField_z, line, "bField_z: ");
        MCP.bFieldMag = findParameterDbl(MCP.bFieldMag, line, "bFieldMag: ");
        
        MCP.cellsA = (uint_fast8_t) findParameterInt(MCP.cellsA, line, "cellsA: ");
        MCP.cellsB = (uint_fast8_t) findParameterInt(MCP.cellsB, line, "cellsB: ");
        MCP.cellsC = (uint_fast8_t) findParameterInt(MCP.cellsC, line, "cellsC: ");
        
        MCP.KbT = findParameterDbl(MCP.KbT, line, "KbT: ");
        MCP.estimatedTc = findParameterDbl(MCP.estimatedTc,
                                           line,
                                           "estimatedTc: ");
        MCP.numSweepsToPerformTotal =
        findParameterInt(MCP.numSweepsToPerformTotal,
                         line,
                         "numSweepsToPerformTotal: ");
        MCP.monte_carlo_seed =
        findParameterULongInt(MCP.monte_carlo_seed,
                              line,
                              "monte_carlo_seed: ");
        
        
        if(line.compare("Finish MonteCarlo::toStringMCP()") == 0){
            break;
        }
        
    }
    
    //cout << toStringMCP() << endl;
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheXDir(SpinCoordinates sc) const{
    
    switch (sc.s) {
        case 0:
            sc.s = 1;
            if((sc.a == MCP.cellsA - 1) && (sc.b == MCP.cellsB - 1)){
                sc.b = sc.b - 1;
            }else if((sc.a == MCP.cellsA - 1) && (sc.b == 0)){
                sc.b = MCP.cellsB - 1;
            }else if((sc.a == 0) && (sc.b == MCP.cellsB - 1)){
                sc.b = sc.b - 1;
            }else if((sc.a == 0) && (sc.b == 0)){
                sc.b = MCP.cellsB - 1;
            }else if(sc.a == 0){
                sc.b = sc.b - 1;
            }else if(sc.b == 0){
                sc.b = MCP.cellsB - 1;
            }else if(sc.a == MCP.cellsA - 1){
                sc.b = sc.b - 1;
            }else if(sc.b == MCP.cellsB - 1){
                sc.b = sc.b - 1;
            }else{
                sc.b = sc.b - 1;
            }
            return sc;
            break;
            
        case 1:
            sc.s = 0;
            if((sc.a == MCP.cellsA - 1) && (sc.b == MCP.cellsB - 1)){
                sc.b = 0;
            }else if((sc.a == MCP.cellsA - 1) && (sc.b == 0)){
                sc.b = sc.b + 1;
            }else if((sc.a == 0) && (sc.b == MCP.cellsB - 1)){
                sc.b = 0;
            }else if((sc.a == 0) && (sc.b == 0)){
                sc.b = sc.b + 1;
            }else if(sc.a == 0){
                sc.b = sc.b + 1;
            }else if(sc.b == 0){
                sc.b = sc.b + 1;
            }else if(sc.a == MCP.cellsA - 1){
                sc.b = sc.b + 1;
            }else if(sc.b == MCP.cellsB - 1){
                sc.b = 0;
            }else{
                sc.b = sc.b + 1;
            }
            return sc;
            break;
            
            
        default:
            exit(1);
            break;
    }
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheYDir(SpinCoordinates sc) const{
    
    switch (sc.s) {
        case 0:
            sc.s = 1;
            if((sc.a == MCP.cellsA - 1) && (sc.b == MCP.cellsB - 1)){
                sc.a = 0;
                sc.b = sc.b - 1;
            }else if((sc.a == MCP.cellsA - 1) && (sc.b == 0)){
                sc.a = 0;
                sc.b = MCP.cellsB - 1;
            }else if((sc.a == 0) && (sc.b == MCP.cellsB - 1)){
                sc.a = sc.a + 1;
                sc.b = sc.b - 1;
            }else if((sc.a == 0) && (sc.b == 0)){
                sc.a = sc.a + 1;
                sc.b = MCP.cellsB - 1;
            }else if(sc.a == 0){
                sc.a = sc.a + 1;
                sc.b = sc.b - 1;
            }else if(sc.b == 0){
                sc.a = sc.a + 1;
                sc.b = MCP.cellsB - 1;
            }else if(sc.a == MCP.cellsA - 1){
                sc.a = 0;
                sc.b = sc.b - 1;
            }else if(sc.b == MCP.cellsB - 1){
                sc.a = sc.a + 1;
                sc.b = sc.b - 1;
            }else{
                sc.a = sc.a + 1;
                sc.b = sc.b - 1;
            }
            return sc;
            break;
            
        case 1:
            sc.s = 0;
            if((sc.a == MCP.cellsA - 1) && (sc.b == MCP.cellsB - 1)){
                sc.a = sc.a - 1;
                sc.b = 0;
            }else if((sc.a == MCP.cellsA - 1) && (sc.b == 0)){
                sc.a = sc.a - 1;
                sc.b = sc.b + 1;
            }else if((sc.a == 0) && (sc.b == MCP.cellsB - 1)){
                sc.a = MCP.cellsA - 1;
                sc.b = 0;
            }else if((sc.a == 0) && (sc.b == 0)){
                sc.a = MCP.cellsA - 1;
                sc.b = sc.b + 1;
            }else if(sc.a == 0){
                sc.a = MCP.cellsA - 1;
                sc.b = sc.b + 1;
            }else if(sc.b == 0){
                sc.a = sc.a - 1;
                sc.b = sc.b + 1;
            }else if(sc.a == MCP.cellsA - 1){
                sc.a = sc.a - 1;
                sc.b = sc.b + 1;
            }else if(sc.b == MCP.cellsB - 1){
                sc.a = sc.a - 1;
                sc.b = 0;
            }else{
                sc.a = sc.a - 1;
                sc.b = sc.b + 1;
            }
            return sc;
            break;
            
        default:
            exit(1);
            break;
    }
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheUpDir(SpinCoordinates sc) const{
    
    
    if(MCP.cellsC == 1){
        return sc;
    }else{
        if(sc.c == MCP.cellsC - 1){
            sc.c = 0;
        }else{
            sc.c = sc.c + 1;
        }
    }
    return sc;
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheDownDir(SpinCoordinates sc) const{
    
    if(MCP.cellsC == 1){
        return sc;
    }else{
        if(sc.c == 0){
            sc.c = MCP.cellsC - 1;
        }else{
            sc.c = sc.c - 1;
        }
    }
    return sc;
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheZDir(SpinCoordinates sc) const{
    
    switch (sc.s) {
        case 0:
            sc.s = 1;
            return sc;
            break;
            
        case 1:
            sc.s = 0;
            return sc;
            break;
            
        default:
            exit(1);
            break;
    }
    
}


void MonteCarlo::setNearestNeighbors(){
    cout<<"Begin setNearestNeighbors"<<endl;
    SpinCoordinates sc;
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    
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
                    
                    //
                    SpinCoordinates sc_x = getSpinCoordsInTheXDir(sc);
                    SpinCoordinates sc_y = getSpinCoordsInTheYDir(sc);
                    SpinCoordinates sc_z = getSpinCoordsInTheZDir(sc);
                    
                    SpinCoordinates sc2_xy = getSpinCoordsInTheYDir(sc_x);
                    SpinCoordinates sc2_yx = getSpinCoordsInTheXDir(sc_y);
                    
                    SpinCoordinates sc2_xz = getSpinCoordsInTheZDir(sc_x);
                    SpinCoordinates sc2_zx = getSpinCoordsInTheXDir(sc_z);
                    
                    SpinCoordinates sc2_yz = getSpinCoordsInTheZDir(sc_y);
                    SpinCoordinates sc2_zy = getSpinCoordsInTheYDir(sc_z);
                    
                    SpinCoordinates sc3_1 = getSpinCoordsInTheZDir(sc2_xy);
                    SpinCoordinates sc3_2 = getSpinCoordsInTheXDir(sc2_yz);
                    SpinCoordinates sc3_3 = getSpinCoordsInTheYDir(sc2_zx);
                    
                    sc.a = a;
                    sc.b = b;
                    sc.c = c;
                    sc.s = s;
                    //first NN
                    neighbors[c][a][b][s][0] = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][1] = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][2] = getSpinCoordsInTheZDir(sc);
                    
                    //second NN
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][3] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][4] = getSpinCoordsInTheXDir(sc);
                    
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][5] = getSpinCoordsInTheZDir(sc);
                    
                    sc = getSpinCoordsInTheZDir(sc);
                    neighbors[c][a][b][s][6] = getSpinCoordsInTheXDir(sc);
                    
                    sc = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][7] = getSpinCoordsInTheZDir(sc);
                    
                    sc = getSpinCoordsInTheZDir(sc);
                    neighbors[c][a][b][s][8] = getSpinCoordsInTheYDir(sc);
                    
                    //third NN
                    sc = getSpinCoordsInTheXDir(sc);
                    sc = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][9] = getSpinCoordsInTheZDir(sc);
                    
                    sc = getSpinCoordsInTheYDir(sc);
                    sc = getSpinCoordsInTheZDir(sc);
                    neighbors[c][a][b][s][10] = getSpinCoordsInTheXDir(sc);
                    
                    sc = getSpinCoordsInTheZDir(sc);
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][11] = getSpinCoordsInTheYDir(sc);
                    
                }
            }
        }
    }
    
}

MonteCarlo::MonteCarlo(MCParameters input_parameters){
    cout<<"Begin MonteCarlo"<<endl;
    
    init();
    
    //initialize private variables
    MCP = input_parameters;
    
    int bigNum  = 10 * 1000 * 1000;
    MCP.j1          = ((int) (bigNum * MCP.j1)) / (1.0 * bigNum);
    MCP.j2          = ((int) (bigNum * MCP.j2)) / (1.0 * bigNum);
    MCP.j3          = ((int) (bigNum * MCP.j3)) / (1.0 * bigNum);
    MCP.k1          = ((int) (bigNum * MCP.k1)) / (1.0 * bigNum);
    MCP.k2          = ((int) (bigNum * MCP.k2)) / (1.0 * bigNum);
    MCP.k3          = ((int) (bigNum * MCP.k3)) / (1.0 * bigNum);
    MCP.cubicD      = ((int) (bigNum * MCP.cubicD)) / (1.0 * bigNum);
    MCP.param1      = ((int) (bigNum * MCP.param1)) / (1.0 * bigNum);
    MCP.bFieldMag   = ((int) (bigNum * MCP.bFieldMag)) / (1.0 * bigNum);
    MCP.KbT         = ((int) (bigNum * MCP.KbT)) / (1.0 * bigNum);
    
    MD.flipAttempts    = 0;
    MD.successfulFlips = 0;
    MD.range           = PI-0.1;
    MD.rangeOld        = 0.1;
    MD.numConfigsDone  = 0;
    MD.isThermal        = false;
    MD.numSweepsPerformedToCheckIfThermalized = 0;
    MD.numSweepsPerformed           = 0;
    MD.numSweepsUsedToThermalize    = 0;
    
    outputFileName           = MCP.path;
    
    if(MCP.isBField){
        double size = sqrt(MCP.bField_x * MCP.bField_x +
                           MCP.bField_y * MCP.bField_y +
                           MCP.bField_z * MCP.bField_z);
        
        MCP.bField_xNorm = MCP.bFieldMag * MCP.bField_x / (1.0*size);
        MCP.bField_yNorm = MCP.bFieldMag * MCP.bField_y / (1.0*size);
        MCP.bField_zNorm = MCP.bFieldMag * MCP.bField_z / (1.0*size);
    }
    
    initializeLattice();
    numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES;
    
    MD.timeAtThisJobInitialization = time(NULL);
    
    cout<<"Finished with (new) MC constructor"<<endl;
    
}

/******************************************************************************/
//Updating order parameter
/******************************************************************************/


void MonteCarlo::updateOrderParameters(){
    
    double E = 0;
    double X = 0;
    double Y = 0;
    double Z = 0;
    double M[SPINDIMEN] = {0};
    double N[SPINDIMEN] = {0};
    double S1[SPINDIMEN] = {0};
    double S2[SPINDIMEN] = {0};
    double S3[SPINDIMEN] = {0};
    double N3D[SPINDIMEN] = {0};
    double ZZ1[SPINDIMEN] = {0};
    double ZZ2[SPINDIMEN] = {0};
    double ZZ3[SPINDIMEN] = {0};
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    
                    X = lattice[c][a][b][s].x;
                    Y = lattice[c][a][b][s].y;
                    Z = lattice[c][a][b][s].z;
                    
                    //OP: E
                    E += getLocalEnergy(c, a, b, s, false);
                    //End OP: E
                    
                    //OP: M
                    M[0] += X;
                    M[1] += Y;
                    M[2] += Z;
                    //End OP: M
                    
                    //OP: N
                    //if s == 0 - count negative M, s == 1 - positive M
                    N[0] += (2 * s - 1) * X;
                    N[1] += (2 * s - 1) * Y;
                    N[2] += (2 * s - 1) * Z;
                    //End OP: N
                    
                    //OP: S1
                    //along perp-y axis (ferro exchange along y direction)
                    if((a % 2 == 0) && (b % 2 == 0) && (s == 1)){
                        S1[0] += X;
                        S1[1] += Y;
                        S1[2] += Z;
                    }else if((a % 2 == 1) && (b % 2 == 1) && (s == 1)){
                        S1[0] += X;
                        S1[1] += Y;
                        S1[2] += Z;
                    }else if((a % 2 == 0) && (b % 2 == 1) && (s == 1)){
                        S1[0] += (-1) * X;
                        S1[1] += (-1) * Y;
                        S1[2] += (-1) * Z;
                    }else if((a % 2 == 1) && (b % 2 == 0) && (s == 1)){
                        S1[0] += (-1) * X;
                        S1[1] += (-1) * Y;
                        S1[2] += (-1) * Z;
                    }else if((a % 2 == 1) && (b % 2 == 0) && (s == 0)){
                        S1[0] += X;
                        S1[1] += Y;
                        S1[2] += Z;
                    }else if((a % 2 == 0) && (b % 2 == 1) && (s == 0)){
                        S1[0] += X;
                        S1[1] += Y;
                        S1[2] += Z;
                    }else if((a % 2 == 0) && (b % 2 == 0) && (s == 0)){
                        S1[0] += (-1) * X;
                        S1[1] += (-1) * Y;
                        S1[2] += (-1) * Z;
                    }else if((a % 2 == 1) && (b % 2 == 1) && (s == 0)){
                        S1[0] += (-1) * X;
                        S1[1] += (-1) * Y;
                        S1[2] += (-1) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: S1
                    
                    //OP: S2
                    //along y axis (ferro exchange along x direction)
                    if((a % 2 == 0) && (s == 1)){
                        S2[0] += X;
                        S2[1] += Y;
                        S2[2] += Z;
                    }else if((a % 2 == 1) && (s == 1)){
                        S2[0] += (-1) * X;
                        S2[1] += (-1) * Y;
                        S2[2] += (-1) * Z;
                    }else if((a % 2 == 0) && (s == 0)){
                        S2[0] += (-1) * X;
                        S2[1] += (-1) * Y;
                        S2[2] += (-1) * Z;
                    }else if((a % 2 == 1) && (s == 0)){
                        S2[0] += X;
                        S2[1] += Y;
                        S2[2] += Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: S2
                    
                    //OP: S3
                    //along x axis (ferro exchange along z direction)
                    if(b % 2 == 0){
                        S3[0] += X;
                        S3[1] += Y;
                        S3[2] += Z;
                    }else if(b % 2 == 1){
                        S3[0] += (-1) * X;
                        S3[1] += (-1) * Y;
                        S3[2] += (-1) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: S3
                    
                    //OP: N3D
                    //if k == 0 - count negative M, k == 1 - positive M
                    if(c % 2 == 0){
                        N3D[0] += (2 * s - 1) * X;
                        N3D[1] += (2 * s - 1) * Y;
                        N3D[2] += (2 * s - 1) * Z;
                    }else if(c % 2 == 1){
                        //Alternatively vertically, count negative M.
                        N3D[0] += (-1) * (2 * s - 1) * X;
                        N3D[1] += (-1) * (2 * s - 1) * Y;
                        N3D[2] += (-1) * (2 * s - 1) * Z;
                    }else{
                        cerr << "shouldn't happen" << endl;
                        exit(3161990);
                    }
                    //End OP: N3D
                    
                    //OP: ZZ1 (Zig Zag) along x axis
                    if((b % 2 == 0) && (s == 1)){
                        ZZ1[0] += X;
                        ZZ1[1] += Y;
                        ZZ1[2] += Z;
                    }else if((b % 2 == 1) && (s == 1)){
                        ZZ1[0] += (-1) * X;
                        ZZ1[1] += (-1) * Y;
                        ZZ1[2] += (-1) * Z;
                    }else if((b % 2 == 0) && (s == 0)){
                        ZZ1[0] += (-1) * X;
                        ZZ1[1] += (-1) * Y;
                        ZZ1[2] += (-1) * Z;
                    }else if((b % 2 == 1) && (s == 0)){
                        ZZ1[0] += X;
                        ZZ1[1] += Y;
                        ZZ1[2] += Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: ZZ1
                    
                    //OP: ZZ2 (Zig Zag) - along Y axis
                    if((a % 2 == 0))
                    {
                        ZZ2[0] += (-1) * X;
                        ZZ2[1] += (-1) * Y;
                        ZZ2[2] += (-1) * Z;
                    }
                    else if((a % 2 == 1))
                    {
                        ZZ2[0] += X;
                        ZZ2[1] += Y;
                        ZZ2[2] += Z;
                    }
                    else
                    {
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: ZZ2
                    
                    //OP: ZZ3 (Zig Zag) - perp to y axis (along
                    if((a % 2 == 0) && (b % 2 == 0))
                    {
                        ZZ3[0] += X;
                        ZZ3[1] += Y;
                        ZZ3[2] += Z;
                    }
                    else if((a % 2 == 0) && (b % 2 == 1))
                    {
                        ZZ3[0] += (-1) * X;
                        ZZ3[1] += (-1) * Y;
                        ZZ3[2] += (-1) * Z;
                    }
                    else if((a % 2 == 1) && (b % 2 == 0))
                    {
                        ZZ3[0] += (-1) * X;
                        ZZ3[1] += (-1) * Y;
                        ZZ3[2] += (-1) * Z;
                    }
                    else if((a % 2 == 1) && (b % 2 == 1))
                    {
                        ZZ3[0] += X;
                        ZZ3[1] += Y;
                        ZZ3[2] += Z;
                    }
                    else
                    {
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: ZZ3
                    
                    
                    
                }
            }
        }
    }
    //End OP
    
    
    OP.scalarOP[0] = E; //Total E, E not per site
    
    //All Order parameters are per site
    
    for(int i = 0; i < SPINDIMEN; i++){
        M[i] = M[i] / numSites;
        N[i] = N[i] / numSites;
        S1[i] = S1[i] / numSites;
        S2[i] = S2[i] / numSites;
        S3[i] = S3[i] / numSites;
        N3D[i] = N3D[i] / numSites;
        ZZ1[i] = ZZ1[i] / numSites;
        ZZ2[i] = ZZ2[i] / numSites;
        ZZ3[i] = ZZ3[i] / numSites;
    }
    
    //Now to calculate the characteristics of the order parameters.
    double OPArr[NUMVECTOROP][SPINDIMEN] = {0};
    for(int i = 0; i < SPINDIMEN; i++){
        OPArr[0][i] = M[i];
        OPArr[1][i] = N[i];
        OPArr[2][i] = S1[i];
        OPArr[3][i] = S2[i];
        OPArr[4][i] = S3[i];
        OPArr[5][i] = sqrt(S1[i]*S1[i] + S2[i]*S2[i] + S3[i]*S3[i]);
        OPArr[6][i] = N3D[i];
        OPArr[7][i] = ZZ1[i];
        OPArr[8][i] = ZZ2[i];
        OPArr[9][i] = ZZ3[i];
        OPArr[10][i] = sqrt(ZZ1[i]*ZZ1[i] + ZZ2[i]*ZZ2[i] + ZZ3[i]*ZZ3[i]);
    }
    
    const double a1 = 0;
    const double a2 = PI/3.0;
    const double a3 = 2*PI/3.0;
    const double a4 = PI;
    const double a5 = 4*PI/3.0;
    const double a6 = 5*PI/3.0;
    
    for(int i = 0; i < NUMVECTOROP; i++){
        
        
        //Scalar quantity eg. M
        OP.vectorOP[i][0] = sqrt(OPArr[i][0]*OPArr[i][0] +
                                 OPArr[i][1]*OPArr[i][1] +
                                 OPArr[i][2]*OPArr[i][2]);
        
        /*
         double tempB1 = (OP.snapOP[i][0]*OP.snapOP[i][0] +
         OP.snapOP[i][1]*OP.snapOP[i][1] -
         2*OP.snapOP[i][2]*OP.snapOP[i][2])/sqrt(6.0);
         
         double tempB2 = (OP.snapOP[i][0]*OP.snapOP[i][0] -
         OP.snapOP[i][1]*OP.snapOP[i][1])/sqrt(2.0);
         */
        
        OP.vectorOP[i][1] = fabs(OPArr[i][0])*cos(a1);
        OP.vectorOP[i][2] = fabs(OPArr[i][0])*sin(a1);
        
        OP.vectorOP[i][1] += fabs(OPArr[i][1])*cos(a3);
        OP.vectorOP[i][2] += fabs(OPArr[i][1])*sin(a3);
        
        OP.vectorOP[i][1] += fabs(OPArr[i][2])*cos(a5);
        OP.vectorOP[i][2] += fabs(OPArr[i][2])*sin(a5);
        
        //sqrt(b1*b1 + b2*b2)
        OP.vectorOP[i][3] = sqrt(OP.vectorOP[i][1]*OP.vectorOP[i][1] +
                                 OP.vectorOP[i][2]*OP.vectorOP[i][2]);
        
        //m
        //m1
        OP.vectorOP[i][4] = 0;
        
        //m2
        OP.vectorOP[i][5] = 0;
        
        if(OPArr[i][0] > 0){
            OP.vectorOP[i][4] += fabs(OPArr[i][0])*cos(a2);
            OP.vectorOP[i][5] += fabs(OPArr[i][0])*sin(a2);
        }else{
            OP.vectorOP[i][4] += fabs(OPArr[i][0])*cos(a5);
            OP.vectorOP[i][5] += fabs(OPArr[i][0])*sin(a5);
        }
        
        if(OPArr[i][1] > 0){
            OP.vectorOP[i][4] += fabs(OPArr[i][1])*cos(a1);
            OP.vectorOP[i][5] += fabs(OPArr[i][1])*sin(a1);
        }else{
            OP.vectorOP[i][4] += fabs(OPArr[i][1])*cos(a4);
            OP.vectorOP[i][5] += fabs(OPArr[i][1])*sin(a4);
        }
        
        if(OPArr[i][2] > 0){
            OP.vectorOP[i][4] += fabs(OPArr[i][2])*cos(a6);
            OP.vectorOP[i][5] += fabs(OPArr[i][2])*sin(a6);
        }else{
            OP.vectorOP[i][4] += fabs(OPArr[i][2])*cos(a3);
            OP.vectorOP[i][5] += fabs(OPArr[i][2])*sin(a3);
        }
        
        
        
        //sqrt(m1*m1 + m2*m2)
        OP.vectorOP[i][6] = sqrt(OP.vectorOP[i][4]*OP.vectorOP[i][4] +
                                 OP.vectorOP[i][5]*OP.vectorOP[i][5]);
        
        //fabs x, y, z, directions
        OP.vectorOP[i][7] = fabs(OPArr[i][0]);
        
        OP.vectorOP[i][8] = fabs(OPArr[i][1]);
        
        OP.vectorOP[i][9] = fabs(OPArr[i][2]);
        
        
    }
    
    //if((MD.numConfigsDone % UPDATESTOCORR == 0)&&(MD.numConfigsDone > 0)){
    //    updateCorrelationFunction();
    //}
    
    //if((MD.numConfigsDone % UPDATESTOFT == 0)&&(MD.numConfigsDone > 0)){
    //    updateFourierTransformOnRecipLattice();
    //}
    
}

void MonteCarlo::updateCorrelationFunction(){
    
    double E[NUMCORRSPINSBIG] = {0};
    double X = 0;
    double Y = 0;
    double Z = 0;
    double M[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double N[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double S1[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double S2[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double S3[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double N3D[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double ZZ1[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double ZZ2[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double ZZ3[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    
    int b = -1;//To start on zero for count = 0.
    int a = 0;
    int s = 0;
    int c = 0;
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        
        //This control structure picks the diagonal spins out across the lattice.
        if(count % 2 == 0 ){
            s = 0;
        }else{
            s = 1;
        }
        
        if(count % 4 == 0 ){
            a = a + 0;
        }else if(count % 4 == 1 ){
            a = a + 0;
        }else if(count % 4 == 2 ){
            a = a + 0;
        }else{
            a = a + 1;
        }
        
        if(count % 4 == 0 ){
            b = b + 1;
        }else if(count % 4 == 1 ){
            b = b + 0;
        }else if(count % 4 == 2 ){
            b = b + 1;
        }else{
            b = b - 1;
        }
        
        c = 0;
        
        //OP: E
        E[count] = getLocalEnergy(c, a, b, s, false);
        //End OP: E
        
        
        X = lattice[c][a][b][s].x;
        Y = lattice[c][a][b][s].y;
        Z = lattice[c][a][b][s].z;
        
        //OP: M
        M[0][count] += X;
        M[1][count] += Y;
        M[2][count] += Z;
        //End OP: M
        
        //OP: N
        //if k == 0 - count negative M, k == 1 - positive M
        N[0][count] += (2 * s - 1) * X;
        N[1][count] += (2 * s - 1) * Y;
        N[2][count] += (2 * s - 1) * Z;
        //End OP: N
        
        //OP: S1
        //along perp-y axis (ferro exchange along y direction)
        if((a % 2 == 0) && (b % 2 == 0) && (s == 1)){
            S1[0][count] += X;
            S1[1][count] += Y;
            S1[2][count] += Z;
        }else if((a % 2 == 1) && (b % 2 == 1) && (s == 1)){
            S1[0][count] += X;
            S1[1][count] += Y;
            S1[2][count] += Z;
        }else if((a % 2 == 0) && (b % 2 == 1) && (s == 1)){
            S1[0][count] += (-1) * X;
            S1[1][count] += (-1) * Y;
            S1[2][count] += (-1) * Z;
        }else if((a % 2 == 1) && (b % 2 == 0) && (s == 1)){
            S1[0][count] += (-1) * X;
            S1[1][count] += (-1) * Y;
            S1[2][count] += (-1) * Z;
        }else if((a % 2 == 1) && (b % 2 == 0) && (s == 0)){
            S1[0][count] += X;
            S1[1][count] += Y;
            S1[2][count] += Z;
        }else if((a % 2 == 0) && (b % 2 == 1) && (s == 0)){
            S1[0][count] += X;
            S1[1][count] += Y;
            S1[2][count] += Z;
        }else if((a % 2 == 0) && (b % 2 == 0) && (s == 0)){
            S1[0][count] += (-1) * X;
            S1[1][count] += (-1) * Y;
            S1[2][count] += (-1) * Z;
        }else if((a % 2 == 1) && (b % 2 == 1) && (s == 0)){
            S1[0][count] += (-1) * X;
            S1[1][count] += (-1) * Y;
            S1[2][count] += (-1) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: S1
        
        //OP: S2
        //along y axis (ferro exchange along x direction)
        if((a % 2 == 0) && (s == 1)){
            S2[0][count] += X;
            S2[1][count] += Y;
            S2[2][count] += Z;
        }else if((a % 2 == 1) && (s == 1)){
            S2[0][count] += (-1) * X;
            S2[1][count] += (-1) * Y;
            S2[2][count] += (-1) * Z;
        }else if((a % 2 == 0) && (s == 0)){
            S2[0][count] += (-1) * X;
            S2[1][count] += (-1) * Y;
            S2[2][count] += (-1) * Z;
        }else if((a % 2 == 1) && (s == 0)){
            S2[0][count] += X;
            S2[1][count] += Y;
            S2[2][count] += Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: S2
        
        //OP: S3
        //along x axis (ferro exchange along z direction)
        if(b % 2 == 0){
            S3[0][count] += X;
            S3[1][count] += Y;
            S3[2][count] += Z;
        }else if(b % 2 == 1){
            S3[0][count] += (-1) * X;
            S3[1][count] += (-1) * Y;
            S3[2][count] += (-1) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: S3
        
        //OP: N3D
        //if k == 0 - count negative M, k == 1 - positive M
        if(c % 2 == 0){
            N3D[0][count] += (2 * s - 1) * X;
            N3D[1][count] += (2 * s - 1) * Y;
            N3D[2][count] += (2 * s - 1) * Z;
        }else if(c % 2 == 1){
            //Alternatively vertically, count negative M.
            N3D[0][count] += (-1) * (2 * s - 1) * X;
            N3D[1][count] += (-1) * (2 * s - 1) * Y;
            N3D[2][count] += (-1) * (2 * s - 1) * Z;
        }else{
            cerr << "shouldn't happen" << endl;
            exit(3161990);
        }
        //End OP: N3D
        
        //OP: ZZ1 (Zig Zag) along x axis
        if((b % 2 == 0) && (s == 1)){
            ZZ1[0][count] += X;
            ZZ1[1][count] += Y;
            ZZ1[2][count] += Z;
        }else if((b % 2 == 1) && (s == 1)){
            ZZ1[0][count] += (-1) * X;
            ZZ1[1][count] += (-1) * Y;
            ZZ1[2][count] += (-1) * Z;
        }else if((b % 2 == 0) && (s == 0)){
            ZZ1[0][count] += (-1) * X;
            ZZ1[1][count] += (-1) * Y;
            ZZ1[2][count] += (-1) * Z;
        }else if((b % 2 == 1) && (s == 0)){
            ZZ1[0][count] += X;
            ZZ1[1][count] += Y;
            ZZ1[2][count] += Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: ZZ1
        
        //OP: ZZ2 (Zig Zag) - along Y axis
        if(a % 2 == 0)
        {
            ZZ2[0][count] += (-1) * X;
            ZZ2[1][count] += (-1) * Y;
            ZZ2[2][count] += (-1) * Z;
        }
        else if(a % 2 == 1)
        {
            ZZ2[0][count] += X;
            ZZ2[1][count] += Y;
            ZZ2[2][count] += Z;
        }
        else
        {
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: ZZ2
        
        //OP: ZZ3 (Zig Zag) - perp to y axis (along
        if((a % 2 == 0) && (b % 2 == 0))
        {
            ZZ3[0][count] += X;
            ZZ3[1][count] += Y;
            ZZ3[2][count] += Z;
        }
        else if((a % 2 == 0) && (b % 2 == 1))
        {
            ZZ3[0][count] += (-1) * X;
            ZZ3[1][count] += (-1) * Y;
            ZZ3[2][count] += (-1) * Z;
        }
        else if((a % 2 == 1) && (b % 2 == 0))
        {
            ZZ3[0][count] += (-1) * X;
            ZZ3[1][count] += (-1) * Y;
            ZZ3[2][count] += (-1) * Z;
        }
        else if((a % 2 == 1) && (b % 2 == 1))
        {
            ZZ3[0][count] += X;
            ZZ3[1][count] += Y;
            ZZ3[2][count] += Z;
        }
        else
        {
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: ZZ3
        
        //OP: TB
        //TB[count] = getPolarAngleTB(X,Y,Z);
        //End OP: TB
        
        //OP: AB
        //AB[count] = getAziAngleAB(X,Y,Z);
        //End OP: AB
        
    }
    
    
    double OPArr[NUMVECTOROP][SPINDIMEN][NUMCORRSPINSBIG] = {0};
    const int DIMOFB = 2;
    double bOP[NUMVECTOROP][DIMOFB][NUMCORRSPINSBIG] = {0};
    double mOP[NUMVECTOROP][DIMOFB][NUMCORRSPINSBIG] = {0};
    
    for(int i = 0; i < SPINDIMEN; i++){
        for(int j = 0; j < 2*MCP.cellsA; j++){
            OPArr[0][i][j] = M[i][j];
            OPArr[1][i][j] = N[i][j];
            OPArr[2][i][j] = S1[i][j];
            OPArr[3][i][j] = S2[i][j];
            OPArr[4][i][j] = S3[i][j];
            OPArr[5][i][j] = sqrt(S1[i][j]*S1[i][j] +
                                  S2[i][j]*S2[i][j] +
                                  S3[i][j]*S3[i][j]);
            
            OPArr[6][i][j] = N3D[i][j];
            OPArr[7][i][j] = ZZ1[i][j];
            OPArr[8][i][j] = ZZ2[i][j];
            OPArr[9][i][j] = ZZ3[i][j];
            OPArr[10][i][j] = sqrt(ZZ1[i][j]*ZZ1[i][j] +
                                   ZZ2[i][j]*ZZ2[i][j] +
                                   ZZ3[i][j]*ZZ3[i][j]);
        }
    }
    
    const double a1 = 0;
    const double a2 = PI/3.0;
    const double a3 = 2*PI/3.0;
    const double a4 = PI;
    const double a5 = 4*PI/3.0;
    const double a6 = 5*PI/3.0;
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int k = 0; k < 2*MCP.cellsA; k++){
            
            
            double tempB1 = 0;
            double tempB2 = 0;
            
            tempB1 += fabs(OPArr[i][0][k])*cos(a1);
            tempB2 += fabs(OPArr[i][0][k])*sin(a1);
            
            tempB1 += fabs(OPArr[i][1][k])*cos(a3);
            tempB2 += fabs(OPArr[i][1][k])*sin(a3);
            
            tempB1 += fabs(OPArr[i][2][k])*cos(a5);
            tempB2 += fabs(OPArr[i][2][k])*sin(a5);
            
            //b1, b2
            bOP[i][0][k] = tempB1;
            bOP[i][1][k] = tempB2;
            
            //m
            double tempM1 = 0;
            double tempM2 = 0;
            
            if(OPArr[i][0] > 0){
                tempM1 += fabs(OPArr[i][0][k])*cos(a2);
                tempM2 += fabs(OPArr[i][0][k])*sin(a2);
            }else{
                tempM1 += fabs(OPArr[i][0][k])*cos(a5);
                tempM2 += fabs(OPArr[i][0][k])*sin(a5);
            }
            
            if(OPArr[i][1] > 0){
                tempM1 += fabs(OPArr[i][1][k])*cos(a1);
                tempM2 += fabs(OPArr[i][1][k])*sin(a1);
            }else{
                tempM1 += fabs(OPArr[i][1][k])*cos(a4);
                tempM2 += fabs(OPArr[i][1][k])*sin(a4);
            }
            
            if(OPArr[i][2] > 0){
                tempM1 += fabs(OPArr[i][2][k])*cos(a6);
                tempM2 += fabs(OPArr[i][2][k])*sin(a6);
            }else{
                tempM1 += fabs(OPArr[i][2][k])*cos(a3);
                tempM2 += fabs(OPArr[i][2][k])*sin(a3);
            }
            
            //m1, m2
            mOP[i][0][k] = tempM1;
            mOP[i][1][k] = tempM2;
            
        }
    }
    
    //Calculate correlation function, dot products of spins
    for(int i = 0; i < NUMVECTOROP; i++){
        
        //OP:
        for(int count = 0; count < 2*MCP.cellsA; count++){
            OP.vectorOPCorrFunc[i][0][count] = (OPArr[i][0][0] * OPArr[i][0][count] +
                                                OPArr[i][1][0] * OPArr[i][1][count] +
                                                OPArr[i][2][0] * OPArr[i][2][count]);
            
            //b OP:
            OP.vectorOPCorrFunc[i][1][count] = (bOP[i][0][0] * bOP[i][0][count]);
            
            OP.vectorOPCorrFunc[i][2][count] = (bOP[i][0][0] * bOP[i][0][count]);
            
            OP.vectorOPCorrFunc[i][3][count] = (bOP[i][0][0] * bOP[i][0][count] +
                                                bOP[i][1][0] * bOP[i][1][count]);
            
            //m OP:
            OP.vectorOPCorrFunc[i][4][count] = (mOP[i][0][0] * mOP[i][0][count]);
            
            OP.vectorOPCorrFunc[i][5][count] = (mOP[i][1][0] * mOP[i][1][count]);
            
            OP.vectorOPCorrFunc[i][6][count] = (mOP[i][0][0] * mOP[i][0][count] +
                                                mOP[i][1][0] * mOP[i][1][count]);
            
            //spin components
            OP.vectorOPCorrFunc[i][7][count] = OPArr[i][0][0] * OPArr[i][0][count];
            
            OP.vectorOPCorrFunc[i][8][count] = OPArr[i][1][0] * OPArr[i][1][count];
            
            OP.vectorOPCorrFunc[i][9][count] = OPArr[i][2][0] * OPArr[i][2][count];
            
            
        }
    }
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        OP.scalarOPCorrFunc[0][count] = E[0] * E[count];
        //OP.TBCorrFunc[count] = TB[0] * TB[count];
        //OP.ABCorrFunc[count] = AB[0] * AB[count];
    }
    
    if(MD.numConfigsDone > 0){
        double temp = 0;
        
        
        for(int i = 0; i < NUMSCALAROP; i++){
            for(int count = 0; count < 2*MCP.cellsA; count++){
                temp = OP.scalarOPCorrFunc[i][count] * OP.scalarOPCorrFunc[i][count];
                RS.sumOverScalarOPCorrFunc[i][0][count] += OP.scalarOPCorrFunc[i][count];
                RS.sumOverScalarOPCorrFunc[i][1][count] += temp;
                RS.sumOverScalarOPCorrFunc[i][2][count] += OP.scalarOPCorrFunc[i][count] * temp;
                RS.sumOverScalarOPCorrFunc[i][3][count] += temp * temp;
            }
        }
        
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                for(int count = 0; count < 2*MCP.cellsA; count++){
                    temp = OP.vectorOPCorrFunc[i][j][count] * OP.vectorOPCorrFunc[i][j][count];
                    RS.sumOverVectorOPCorrFunc[i][j][0][count] += OP.vectorOPCorrFunc[i][j][count];
                    RS.sumOverVectorOPCorrFunc[i][j][1][count] += temp;
                    RS.sumOverVectorOPCorrFunc[i][j][2][count] += temp * OP.vectorOPCorrFunc[i][j][count];
                    RS.sumOverVectorOPCorrFunc[i][j][3][count] += temp * temp;
                }
            }
        }
        
        
    }
    
    
};

void MonteCarlo::updateFourierTransformOnRecipLattice(){
    
    double kvector_x = 0;
    double kvector_y = 0;
    double rvector_x = 0;
    double rvector_y = 0;
    double kr = 0;
    double cos_kr = 0;
    double sin_kr = 0;
    
    //Loop over all k points
    for(uint_fast8_t k_c = 0; k_c < MCP.cellsC; k_c++){
        for(uint_fast8_t k_a = 0; k_a < MCP.cellsA; k_a++){
            for(uint_fast8_t k_b = 0; k_b < MCP.cellsB; k_b++){
                for(uint_fast8_t k_s = 0; k_s < 1; k_s++){
                    
                    //Loop over all sites in real space
                    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
                        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
                            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                                for(uint_fast8_t s = 0; s < 1; s++){
                                    
                                    //sublattice A
                                    kvector_x = k_a * B2_VECTOR_1[0] / (MCP.cellsA*1.0) +
                                    k_b * B2_VECTOR_2[0] / (MCP.cellsB*1.0);
                                    kvector_y = k_a * B2_VECTOR_1[1] / (MCP.cellsA*1.0) +
                                    k_b * B2_VECTOR_2[1] / (MCP.cellsB*1.0);
                                    rvector_x = a * A_VECTOR_1[0] + b * A_VECTOR_2[0];
                                    rvector_y = a * A_VECTOR_1[1] + b * A_VECTOR_2[1];
                                    kr = kvector_x * rvector_x + kvector_y * rvector_y;
                                    cos_kr = cos(kr);
                                    sin_kr = sin(kr);
                                    recipLattice[k_c][k_a][k_b][k_s].ReXComponent +=
                                    cos_kr * lattice[c][a][b][s].x / (numSites*1.0);
                                    recipLattice[k_c][k_a][k_b][k_s].ImXComponent +=
                                    sin_kr * lattice[c][a][b][s].x / (numSites*1.0);
                                    recipLattice[k_c][k_a][k_b][k_s].ReYComponent +=
                                    cos_kr * lattice[c][a][b][s].y / (numSites*1.0);
                                    recipLattice[k_c][k_a][k_b][k_s].ImYComponent +=
                                    sin_kr * lattice[c][a][b][s].y / (numSites*1.0);
                                    recipLattice[k_c][k_a][k_b][k_s].ReZComponent +=
                                    cos_kr * lattice[c][a][b][s].z / (numSites*1.0);
                                    recipLattice[k_c][k_a][k_b][k_s].ImZComponent +=
                                    sin_kr * lattice[c][a][b][s].z / (numSites*1.0);
                                    
                                    /*
                                     //sublattice B
                                     k_s = 1;
                                     s = 1;
                                     rvector_y = a * A_VECTOR_1[1] + b * A_VECTOR_2[1] + 1/sqrt(3.0);
                                     kr = kvector_x * rvector_x + kvector_y * rvector_y;
                                     cos_kr = cos(kr);
                                     sin_kr = sin(kr);
                                     recipLattice[k_c][k_a][k_b][k_s].ReXComponent +=
                                     cos_kr * lattice[c][a][b][s].x / (numSites*1.0);
                                     recipLattice[k_c][k_a][k_b][k_s].ImXComponent +=
                                     sin_kr * lattice[c][a][b][s].x / (numSites*1.0);
                                     recipLattice[k_c][k_a][k_b][k_s].ReYComponent +=
                                     cos_kr * lattice[c][a][b][s].y / (numSites*1.0);
                                     recipLattice[k_c][k_a][k_b][k_s].ImYComponent +=
                                     sin_kr * lattice[c][a][b][s].y / (numSites*1.0);
                                     recipLattice[k_c][k_a][k_b][k_s].ReZComponent +=
                                     cos_kr * lattice[c][a][b][s].z / (numSites*1.0);
                                     recipLattice[k_c][k_a][k_b][k_s].ImZComponent +=
                                     sin_kr * lattice[c][a][b][s].z / (numSites*1.0);
                                     */
                                    
                                }
                            }
                        }
                    }
                    
                    recipLattice[k_c][k_a][k_b][k_s].magnitude =
                    (
                     recipLattice[k_c][k_a][k_b][k_s].ReXComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReXComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ReYComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReYComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ReZComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReZComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImXComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImXComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImYComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImYComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImZComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImZComponent
                     );
                    
                    
                    k_s = 1;
                    recipLattice[k_c][k_a][k_b][k_s].magnitude = (0);
                    /*
                     recipLattice[k_c][k_a][k_b][k_s].ReXComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReXComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ReYComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReYComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ReZComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ReZComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImXComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImXComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImYComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImYComponent +
                     recipLattice[k_c][k_a][k_b][k_s].ImZComponent *
                     recipLattice[k_c][k_a][k_b][k_s].ImZComponent
                     );
                     */
                }
            }
        }
    }
    //Loop over all k points
    for(uint_fast8_t k_c = 0; k_c < MCP.cellsC; k_c++){
        for(uint_fast8_t k_a = 0; k_a < MCP.cellsA; k_a++){
            for(uint_fast8_t k_b = 0; k_b < MCP.cellsB; k_b++){
                for(uint_fast8_t k_s = 0; k_s < NUMSUBLATTICES; k_s++){
                    
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReXComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ReXComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImXComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ImXComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReYComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ReYComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImYComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ImYComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReZComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ReZComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImZComponent +=
                    recipLattice[k_c][k_a][k_b][k_s].ImZComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].magnitude    +=
                    recipLattice[k_c][k_a][k_b][k_s].magnitude;
                    
                }
            }
        }
    }
    
}


double MonteCarlo::getLocalEnergy(uint_fast8_t c, uint_fast8_t a,
                                  uint_fast8_t b, uint_fast8_t s,
                                  const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
#ifdef CUBIC_HAMILTONIAN
    double interaction1 = 0;
    if(fromSpinFlip){
        interaction1 = (-1) * MCP.cubicD * (lattice[c][a][b][s].x*lattice[c][a][b][s].x*
                                            lattice[c][a][b][s].x*lattice[c][a][b][s].x +
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y*
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y +
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z*
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z);
    }else if(!fromSpinFlip){//From Finding the total E of the lattice
        interaction1 = (-1) * MCP.cubicD * (lattice[c][a][b][s].x*lattice[c][a][b][s].x*
                                            lattice[c][a][b][s].x*lattice[c][a][b][s].x +
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y*
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y +
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z*
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z)/2.0;
    }
    totalInteraction += interaction1;
#endif
    
#ifdef ROTATED_HAMILTONIAN
    //Determine the sublattice of a
    if((a % 2 == 0)&&(b % 2 == 1)&&(s == 0)||
       (a % 2 == 1)&&(b % 2 == 1)&&(s == 1))
    {
        //Black circle sublattice = A
        //Spin is unchanged
    }
    if((a % 2 == 0)&&(b % 2 == 0)&&(s == 1)||
       (a % 2 == 1)&&(b % 2 == 0)&&(s == 0))
    {
        //Black square sublattice = B
        lattice[c][a][b][s].y = (-1.0) * lattice[c][a][b][s].y;
        lattice[c][a][b][s].z = (-1.0) * lattice[c][a][b][s].z;
    }
    if((a % 2 == 0)&&(b % 2 == 0)&&(s == 0)||
       (a % 2 == 1)&&(b % 2 == 0)&&(s == 1))
    {
        //white square sublattice = C
        lattice[c][a][b][s].x = (-1.0) * lattice[c][a][b][s].x;
        lattice[c][a][b][s].z = (-1.0) * lattice[c][a][b][s].z;
    }
    if((a % 2 == 0)&&(b % 2 == 1)&&(s == 1)||
       (a % 2 == 1)&&(b % 2 == 1)&&(s == 0))
    {
        //White circle sublattice = D
        lattice[c][a][b][s].x = (-1.0) * lattice[c][a][b][s].x;
        lattice[c][a][b][s].y = (-1.0) * lattice[c][a][b][s].y;
    }
    
    //Heisenberg interaction
    
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    double interaction1 = (MCP.alpha - 1.0) * (dotProd(c,a,b,s, neighbors[c][a][b][s][0]) +
                                               dotProd(c,a,b,s, neighbors[c][a][b][s][1]) +
                                               dotProd(c,a,b,s, neighbors[c][a][b][s][2]));
    
    //Kitaev interaction
    SpinCoordinates sc_x = neighbors[c][a][b][s][0];
    SpinCoordinates sc_y = neighbors[c][a][b][s][1];
    SpinCoordinates sc_z = neighbors[c][a][b][s][2];
    
    double interaction2 = (2.0 - 4.0*MCP.alpha) * (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                                                   lattice[c][a][b][s].x +
                                                   lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                                                   lattice[c][a][b][s].y +
                                                   lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                                                   lattice[c][a][b][s].z);
    
    
    if(fromSpinFlip)
    {
        totalInteraction += interaction1 + interaction2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += interaction1/2.0 + interaction2/2.0;
    }
#endif
    
    
#ifdef ALPHA_HAMILTONIAN
    
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    SpinCoordinates sc_x = neighbors[c][a][b][s][0];
    SpinCoordinates sc_y = neighbors[c][a][b][s][1];
    SpinCoordinates sc_z = neighbors[c][a][b][s][2];
    
    //Heisenberg interaction
    double interactionH = (dotProd(c,a,b,s, sc_x) +
                           dotProd(c,a,b,s, sc_y) +
                           dotProd(c,a,b,s, sc_z));
    
    double interactionK = (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                           lattice[c][a][b][s].x +
                           lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                           lattice[c][a][b][s].y +
                           lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                           lattice[c][a][b][s].z);
    
    if(fromSpinFlip){
        totalInteraction += (1-MCP.alpha)*(interactionH) + (-2.0*MCP.alpha) * (interactionK);
    }else{//From Finding the total E of the lattice
        totalInteraction += (1-MCP.alpha)*(interactionH/2.0) + (-2.0*MCP.alpha)*(interactionK/2.0);
    }
#endif
    
    
#ifdef EXTENDED_HAMILTONIAN
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    SpinCoordinates sc_x = neighbors[c][a][b][s][0];
    SpinCoordinates sc_y = neighbors[c][a][b][s][1];
    SpinCoordinates sc_z = neighbors[c][a][b][s][2];
    //Heisenberg interaction
    //Normal rotated coords
    double interaction1 = cos(MCP.phi) * (dotProd(c,a,b,s, sc_x) +
                                          dotProd(c,a,b,s, sc_y) +
                                          dotProd(c,a,b,s, sc_z));
    //Kitaev interaction
    double interaction2 = 2 * sin(MCP.phi) * (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                                              lattice[c][a][b][s].x +
                                              lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                                              lattice[c][a][b][s].y +
                                              lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                                              lattice[c][a][b][s].z);
    
    
    
    if(fromSpinFlip)
    {
        totalInteraction += interaction1 + interaction2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += interaction1/2.0 + interaction2/2.0;
    }
#endif
    
    
#ifdef J23K23_HAMILTONIAN
    SpinCoordinates sc_x = neighbors[c][a][b][s][0];
    SpinCoordinates sc_y = neighbors[c][a][b][s][1];
    SpinCoordinates sc_z = neighbors[c][a][b][s][2];
    
    SpinCoordinates sc2_xy = neighbors[c][a][b][s][3];
    SpinCoordinates sc2_yx = neighbors[c][a][b][s][4];
    
    SpinCoordinates sc2_xz = neighbors[c][a][b][s][5];
    SpinCoordinates sc2_zx = neighbors[c][a][b][s][6];
    
    SpinCoordinates sc2_yz = neighbors[c][a][b][s][7];
    SpinCoordinates sc2_zy = neighbors[c][a][b][s][8];
    
    SpinCoordinates sc3_1 = neighbors[c][a][b][s][9];
    SpinCoordinates sc3_2 = neighbors[c][a][b][s][10];
    SpinCoordinates sc3_3 = neighbors[c][a][b][s][11];
    
    double interactionH = MCP.j1 * (dotProd(c,a,b,s, sc_x) +
                                    dotProd(c,a,b,s, sc_y) +
                                    dotProd(c,a,b,s, sc_z));
    
    double interactionK = MCP.k1 * (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                                    lattice[c][a][b][s].x +
                                    lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                                    lattice[c][a][b][s].y +
                                    lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                                    lattice[c][a][b][s].z);
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (dotProd(c,a,b,s, sc2_xy) +
                                     dotProd(c,a,b,s, sc2_xz) +
                                     dotProd(c,a,b,s, sc2_yz) +
                                     dotProd(c,a,b,s, sc2_yx) +
                                     dotProd(c,a,b,s, sc2_zx) +
                                     dotProd(c,a,b,s, sc2_zy));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(c,a,b,s, sc3_1) +
                                     dotProd(c,a,b,s, sc3_2) +
                                     dotProd(c,a,b,s, sc3_3));
    
    //Kitaev Interaction
    
    double interactionK2 = MCP.k2 * (lattice[sc2_yz.c][sc2_yz.a][sc2_yz.b][sc2_yz.s].x *
                                     lattice[c][a][b][s].x +
                                     lattice[sc2_zy.c][sc2_zy.a][sc2_zy.b][sc2_zy.s].x *
                                     lattice[c][a][b][s].x +
                                     lattice[sc2_zx.c][sc2_zx.a][sc2_zx.b][sc2_zx.s].y *
                                     lattice[c][a][b][s].y +
                                     lattice[sc2_xz.c][sc2_xz.a][sc2_xz.b][sc2_xz.s].y *
                                     lattice[c][a][b][s].y +
                                     lattice[sc2_yx.c][sc2_yx.a][sc2_yx.b][sc2_yx.s].z *
                                     lattice[c][a][b][s].z +
                                     lattice[sc2_xy.c][sc2_xy.a][sc2_xy.b][sc2_xy.s].z *
                                     lattice[c][a][b][s].z);
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 +
        interactionK + interactionK2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 +
                             interactionK + interactionK2) / 2.0;
    }
#endif
    
    
#ifdef BFIELD_HAMILTONIAN
    double interactionB = (-1) * (MCP.bField_xNorm * lattice[c][a][b][s].x +
                                  MCP.bField_yNorm * lattice[c][a][b][s].y +
                                  MCP.bField_zNorm * lattice[c][a][b][s].z);
    totalInteraction += interactionB;
#endif
    
    return totalInteraction;
    
}


double MonteCarlo::getLocalEnergyFromCalc(uint_fast8_t c, uint_fast8_t a,
                                          uint_fast8_t b, uint_fast8_t s,
                                          const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
#ifdef CUBIC_HAMILTONIAN
    double interaction1 = 0;
    if(fromSpinFlip){
        interaction1 = (-1) * MCP.cubicD * (lattice[c][a][b][s].x*lattice[c][a][b][s].x*
                                            lattice[c][a][b][s].x*lattice[c][a][b][s].x +
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y*
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y +
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z*
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z);
    }else if(!fromSpinFlip){//From Finding the total E of the lattice
        interaction1 = (-1) * MCP.cubicD * (lattice[c][a][b][s].x*lattice[c][a][b][s].x*
                                            lattice[c][a][b][s].x*lattice[c][a][b][s].x +
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y*
                                            lattice[c][a][b][s].y*lattice[c][a][b][s].y +
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z*
                                            lattice[c][a][b][s].z*lattice[c][a][b][s].z)/2.0;
    }
    totalInteraction += interaction1;
#endif
    
#ifdef ROTATED_HAMILTONIAN
    //Determine the sublattice of a
    if((a % 2 == 0)&&(b % 2 == 1)&&(s == 0)||
       (a % 2 == 1)&&(b % 2 == 1)&&(s == 1))
    {
        //Black circle sublattice = A
        //Spin is unchanged
    }
    if((a % 2 == 0)&&(b % 2 == 0)&&(s == 1)||
       (a % 2 == 1)&&(b % 2 == 0)&&(s == 0))
    {
        //Black square sublattice = B
        lattice[c][a][b][s].y = (-1.0) * lattice[c][a][b][s].y;
        lattice[c][a][b][s].z = (-1.0) * lattice[c][a][b][s].z;
    }
    if((a % 2 == 0)&&(b % 2 == 0)&&(s == 0)||
       (a % 2 == 1)&&(b % 2 == 0)&&(s == 1))
    {
        //white square sublattice = C
        lattice[c][a][b][s].x = (-1.0) * lattice[c][a][b][s].x;
        lattice[c][a][b][s].z = (-1.0) * lattice[c][a][b][s].z;
    }
    if((a % 2 == 0)&&(b % 2 == 1)&&(s == 1)||
       (a % 2 == 1)&&(b % 2 == 1)&&(s == 0))
    {
        //White circle sublattice = D
        lattice[c][a][b][s].x = (-1.0) * lattice[c][a][b][s].x;
        lattice[c][a][b][s].y = (-1.0) * lattice[c][a][b][s].y;
    }
    
    //Heisenberg interaction
    
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    double interaction1 = (MCP.alpha - 1.0) * (dotProd(c,a,b,s,
                                                       getSpinCoordsInTheXDir(sc)) +
                                               dotProd(c,a,b,s,
                                                       getSpinCoordsInTheYDir(sc)) +
                                               dotProd(c,a,b,s,
                                                       getSpinCoordsInTheZDir(sc)));
    
    //Kitaev interaction
    SpinCoordinates sc_x = getSpinCoordsInTheXDir(sc);
    SpinCoordinates sc_y = getSpinCoordsInTheYDir(sc);
    SpinCoordinates sc_z = getSpinCoordsInTheZDir(sc);
    
    double interaction2 = (2.0 - 4.0*MCP.alpha) * (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                                                   lattice[c][a][b][s].x +
                                                   lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                                                   lattice[c][a][b][s].y +
                                                   lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                                                   lattice[c][a][b][s].z);
    
    
    if(fromSpinFlip)
    {
        totalInteraction += interaction1 + interaction2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += interaction1/2.0 + interaction2/2.0;
    }
#endif
    
    
#ifdef ALPHA_HAMILTONIAN
    
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    SpinCoordinates sc_x = getSpinCoordsInTheXDir(sc);
    SpinCoordinates sc_y = getSpinCoordsInTheYDir(sc);
    SpinCoordinates sc_z = getSpinCoordsInTheZDir(sc);
    
    //Heisenberg interaction
    double interactionH = (dotProd(c,a,b,s, sc_x) +
                           dotProd(c,a,b,s, sc_y) +
                           dotProd(c,a,b,s, sc_z));
    
    double interactionK = (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                           lattice[c][a][b][s].x +
                           lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                           lattice[c][a][b][s].y +
                           lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                           lattice[c][a][b][s].z);
    
    if(fromSpinFlip){
        totalInteraction += (1-MCP.alpha)*(interactionH) + (-2.0*MCP.alpha) * (interactionK);
    }else{//From Finding the total E of the lattice
        totalInteraction += (1-MCP.alpha)*(interactionH/2.0) + (-2.0*MCP.alpha)*(interactionK/2.0);
    }
#endif
    
    
#ifdef EXTENDED_HAMILTONIAN
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    SpinCoordinates sc_x = getSpinCoordsInTheXDir(sc);
    SpinCoordinates sc_y = getSpinCoordsInTheYDir(sc);
    SpinCoordinates sc_z = getSpinCoordsInTheZDir(sc);
    //Heisenberg interaction
    //Normal rotated coords
    double interaction1 = cos(MCP.phi) * (dotProd(c,a,b,s, sc_x) +
                                          dotProd(c,a,b,s, sc_y) +
                                          dotProd(c,a,b,s, sc_z));
    //Kitaev interaction
    double interaction2 = 2 * sin(MCP.phi) *
    (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
     lattice[c][a][b][s].x +
     lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
     lattice[c][a][b][s].y +
     lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
     lattice[c][a][b][s].z);
    
    
    
    if(fromSpinFlip)
    {
        totalInteraction += interaction1 + interaction2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += interaction1/2.0 + interaction2/2.0;
    }
#endif
    
    
#ifdef J23K23_HAMILTONIAN
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    SpinCoordinates sc_x = getSpinCoordsInTheXDir(sc);
    SpinCoordinates sc_y = getSpinCoordsInTheYDir(sc);
    SpinCoordinates sc_z = getSpinCoordsInTheZDir(sc);
    
    SpinCoordinates sc2_xy = getSpinCoordsInTheYDir(sc_x);
    SpinCoordinates sc2_yx = getSpinCoordsInTheXDir(sc_y);
    
    SpinCoordinates sc2_xz = getSpinCoordsInTheZDir(sc_x);
    SpinCoordinates sc2_zx = getSpinCoordsInTheXDir(sc_z);
    
    SpinCoordinates sc2_yz = getSpinCoordsInTheZDir(sc_y);
    SpinCoordinates sc2_zy = getSpinCoordsInTheYDir(sc_z);
    
    SpinCoordinates sc3_1 = getSpinCoordsInTheZDir(sc2_xy);
    SpinCoordinates sc3_2 = getSpinCoordsInTheXDir(sc2_yz);
    SpinCoordinates sc3_3 = getSpinCoordsInTheYDir(sc2_zx);
    
    double interactionH = MCP.j1 * (dotProd(c,a,b,s, sc_x) +
                                    dotProd(c,a,b,s, sc_y) +
                                    dotProd(c,a,b,s, sc_z));
    
    double interactionK = MCP.k1 * (lattice[sc_x.c][sc_x.a][sc_x.b][sc_x.s].x *
                                    lattice[c][a][b][s].x +
                                    lattice[sc_y.c][sc_y.a][sc_y.b][sc_y.s].y *
                                    lattice[c][a][b][s].y +
                                    lattice[sc_z.c][sc_z.a][sc_z.b][sc_z.s].z *
                                    lattice[c][a][b][s].z);
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (dotProd(c,a,b,s, sc2_xy) +
                                     dotProd(c,a,b,s, sc2_xz) +
                                     dotProd(c,a,b,s, sc2_yz) +
                                     dotProd(c,a,b,s, sc2_yx) +
                                     dotProd(c,a,b,s, sc2_zx) +
                                     dotProd(c,a,b,s, sc2_zy));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(c,a,b,s, sc3_1) +
                                     dotProd(c,a,b,s, sc3_2) +
                                     dotProd(c,a,b,s, sc3_3));
    
    //Kitaev Interaction
    
    double interactionK2 = MCP.k2 * (lattice[sc2_yz.c][sc2_yz.a][sc2_yz.b][sc2_yz.s].x *
                                     lattice[c][a][b][s].x +
                                     lattice[sc2_zy.c][sc2_zy.a][sc2_zy.b][sc2_zy.s].x *
                                     lattice[c][a][b][s].x +
                                     lattice[sc2_zx.c][sc2_zx.a][sc2_zx.b][sc2_zx.s].y *
                                     lattice[c][a][b][s].y +
                                     lattice[sc2_xz.c][sc2_xz.a][sc2_xz.b][sc2_xz.s].y *
                                     lattice[c][a][b][s].y +
                                     lattice[sc2_yx.c][sc2_yx.a][sc2_yx.b][sc2_yx.s].z *
                                     lattice[c][a][b][s].z +
                                     lattice[sc2_xy.c][sc2_xy.a][sc2_xy.b][sc2_xy.s].z *
                                     lattice[c][a][b][s].z);
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 + interactionK + interactionK2;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 + interactionK + interactionK2) / 2.0;
    }
#endif
    
    
#ifdef BFIELD_HAMILTONIAN
    double interactionB = (-1) * (MCP.bField_xNorm * lattice[c][a][b][s].x +
                                  MCP.bField_yNorm * lattice[c][a][b][s].y +
                                  MCP.bField_zNorm * lattice[c][a][b][s].z);
    totalInteraction += interactionB;
#endif
    
    return totalInteraction;
    
}



/******************************************************************************/
//Not Currently Used
/******************************************************************************/

void MonteCarlo::setDomainWall(){//Test before using
    
    
    //Neel Z, -Y
    //radius = cellsA/2 & cellsA/2 -1
    //Each set of three spins must be a normalized vector.
    
    int radius1 = MCP.cellsA/2 -1;
    double spin1X = 0;
    double spin1Y = 0;
    double spin1Z = 1;
    
    int radius2 = MCP.cellsA/2;
    double spin2X = 0;
    double spin2Y = -1;
    double spin2Z = 0;
    
    
    /*
     //Neel Z, -X
     //radius = cellsA/2 & cellsA/2 -1
     //Each set of three spins must be a normalized vector.
     
     int radius1 = cellsA/2 -1;
     double spin1X = 0;
     double spin1Y = 0;
     double spin1Z = 1;
     
     int radius2 = cellsA/2;
     double spin2X = -1;
     double spin2Y = 0;
     double spin2Z = 0;
     */
    
    /*
     //Neel Z, -Z
     //radius = cellsA/2 & cellsA/2 -1
     //Each set of three spins must be a normalized vector.
     
     int radius1 = cellsA/2 -1;
     double spin1X = 0;
     double spin1Y = 0;
     double spin1Z = 1;
     
     int radius2 = cellsA/2;
     double spin2X = 0;
     double spin2Y = 0;
     double spin2Z = -1;
     */
    
    /*
     //Neel Z, X
     //radius = cellsA/2 & cellsA/2 -1
     //Each set of three spins must be a normalized vector.
     
     int radius1 = cellsA/2 -1;
     double spin1X = 0;
     double spin1Y = 0;
     double spin1Z = 1;
     
     int radius2 = cellsA/2;
     double spin2X = 1;
     double spin2Y = 0;
     double spin2Z = 0;
     */
    
    /*
     //Neel Z, Y
     //radius = cellsA/2 & cellsA/2 -1
     //Each set of three spins must be a normalized vector.
     
     int radius1 = cellsA/2 -1;
     double spin1X = 0;
     double spin1Y = 0;
     double spin1Z = 1;
     
     int radius2 = cellsA/2;
     double spin2X = 0;
     double spin2Y = 1;
     double spin2Z = 0;
     */
    
    /*
     //Neel Z, Z
     //radius = cellsA/2 & cellsA/2 -1
     //Each set of three spins must be a normalized vector.
     
     int radius1 = cellsA/2 -1;
     double spin1X = 0;
     double spin1Y = 0;
     double spin1Z = 1;
     
     int radius2 = cellsA/2;
     double spin2X = 0;
     double spin2Y = 0;
     double spin2Z = 1;
     */
    
    
    setDomainWallSpins(spin1X,spin1Y,spin1Z,radius1,//inner wall
                       spin2X,spin2Y,spin2Z,radius2//outer wall
                       );
}

void MonteCarlo::setDomainWallSpins(double spin1X, double spin1Y,
                                    double spin1Z, int radius1,//inner wall
                                    double spin2X, double spin2Y,
                                    double spin2Z, int radius2//outer wall
){//Test before using
    //cout<<"Start Making Domain Wall"<<endl;
    //Only Implemented for Neel Order
    SpinCoordinates sc;
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    //cout<<"i = " << i << " j = " << j << " k = " << k << " z = " << z;
                    sc.a = a;
                    sc.b = b;
                    sc.c = c;
                    sc.s = s;
                    if(isDomain(sc, radius1)){
                        //if k == 0 - count negative M, k == 1 - positive M
                        lattice[c][a][b][s].x = ((2 * s - 1) * spin1X);
                        lattice[c][a][b][s].y = ((2 * s - 1) * spin1Y);
                        lattice[c][a][b][s].z = ((2 * s - 1) * spin1Z);
                        //cout<< " <- Dom1";
                    }
                    
                    if(isDomain(sc, radius2)){
                        //if k == 0 - count negative M, k == 1 - positive M
                        lattice[c][a][b][s].x = ((2 * s - 1) * spin2X);
                        lattice[c][a][b][s].y = ((2 * s - 1) * spin2Y);
                        lattice[c][a][b][s].z = ((2 * s - 1) * spin2Z);
                        //cout<< " <- Dom2";
                    }
                    //cout<<endl;
                }
            }
        }
    }
    //cout<<"End Making Domain Wall"<<endl;
    //exit(0);
}

bool MonteCarlo::isDomain(SpinCoordinates sc, int radius){//Test before using
    //indepedant of z
    
    int startPosX = MCP.cellsA/2;
    
    //left side
    if((sc.s == 1)&&(sc.a == (startPosX - radius))&&(sc.b == 0)){
        return true;
    }
    if((sc.a == (startPosX - radius))&&(sc.b > 0)&&(sc.b <= radius -1)){
        return true;
    }
    
    //top side
    if((sc.s == 0)&&
       (sc.a <= (startPosX -1))&&
       (sc.a >= (startPosX - radius))&&
       (sc.b == radius)){
        return true;
    }
    if((sc.s == 1)&&
       (sc.a <= (startPosX -1))&&
       (sc.a >= (startPosX - radius +1))&&
       (sc.b == radius-1)){
        return true;
    }
    
    //right side
    if((sc.a >= startPosX)&&
       ((sc.a - startPosX) + sc.b == radius - 1)&&
       ((sc.b != 0)||((sc.s == 1)&&(sc.b == 0)))
       ){
        return true;
    }
    
    return false;
}



double MonteCarlo::getPolarAngleTB(const double spinX,
                                   const double spinY, const double spinZ){
    //OP: TB, Polar Angle to B Field
    double Bfield[3] = {0};
    Bfield[0] = MCP.bField_x;
    Bfield[1] = MCP.bField_y;
    Bfield[2] = MCP.bField_z;
    double spin[3] = {0};
    spin[0] = spinX;
    spin[1] = spinY;
    spin[2] = spinZ;
    
    //Norming is important because in general the vectors (of spin at least) are
    //not normalized
    double normB = sqrt(Bfield[0]*Bfield[0] + Bfield[1]*Bfield[1] + Bfield[2]*Bfield[2]);
    double normS = sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
    double cosTB = (spin[0]*Bfield[0] + spin[1]*Bfield[1] + spin[2]*Bfield[2])/(normB * normS);
    
    if((cosTB > 1.0000001)||(cosTB < -1.00000001)){
        cout<<"bad TB "<<cosTB<<endl;
        cerr<<"bad TB "<<cosTB<<endl;
        exit(1);
    }
    if((cosTB > 1.0)){
        cosTB = 1;
    }
    if((cosTB < -1.0)){
        cosTB = -1;
    }
    
    return acos(cosTB);
    //End OP TB
    
}

double MonteCarlo::getAziAngleAB(const double spinX,
                                 const double spinY, const double spinZ){
    //OP: AB, Azimuthal Angle to the B Field
    
    /*
     Unlike finding the polar angle to the B field, the azimuthal angle is
     slightly more subtle. To completely specify the orientation of the
     azimuthal angle, we will construct a right handed coordinate system which
     will fully define the azimuthal angle.
     If the direction of the B field is analogous to a "z" direction,
     then define n to be analogous to the "x" direction via: n = B x z
     where B is the direction of the external field and z is the "c"crystal axis.
     if B is along z, then let n = x. where x is the "a" crystal axis.
     Finally, we define m = B x n - analogous to the "b" crystal axis.
     */
    
    //Norming is important because in general the vectors (of spin at least) are
    //not normalized
    double Bfield[3] = {0};
    Bfield[0] = MCP.bField_x;
    Bfield[1] = MCP.bField_y;
    Bfield[2] = MCP.bField_z;
    double normB = sqrt(Bfield[0]*Bfield[0] +
                        Bfield[1]*Bfield[1] +
                        Bfield[2]*Bfield[2]);
    
    Bfield[0] = Bfield[0] / normB;
    Bfield[1] = Bfield[1] / normB;
    Bfield[2] = Bfield[2] / normB;
    
    double spin[3] = {0};
    spin[0] = spinX;
    spin[1] = spinY;
    spin[2] = spinZ;
    double normS = sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
    
    spin[0] = spin[0] / normS;
    spin[1] = spin[1] / normS;
    spin[2] = spin[2] / normS;
    
    //Bfield is analogous to the z direction
    //nhat = B cross zhat - This is analogous to the x direction
    double nhat[3] = {0};
    nhat[0] = MCP.bField_y;//y->x in cross products
    nhat[1] = (-1)*MCP.bField_x;
    nhat[2] = 0;
    
    if(fabs(MCP.bField_z - 1) < 1e-5){
        nhat[0] = 1;
        nhat[1] = 0;
        nhat[2] = 0;
    }
    if(fabs(MCP.bField_z - (-1)) < 1e-5){
        nhat[0] = -1;
        nhat[1] = 0;
        nhat[2] = 0;
    }
    
    double normN = sqrt(nhat[0]*nhat[0] + nhat[1]*nhat[1] + nhat[2]*nhat[2]);
    nhat[0] = nhat[0] / normN;
    nhat[1] = nhat[1] / normN;
    nhat[2] = nhat[2] / normN;
    
    //Third axis definition, the resultant from B x n
    double vert[3] = {0};
    vert[0] =           Bfield[1] * nhat[2] - Bfield[2] * nhat[1];
    vert[1] = (-1) *   (Bfield[0] * nhat[2] - Bfield[2] * nhat[0]);
    vert[2] =           Bfield[0] * nhat[1] - Bfield[1] * nhat[0];
    double normV = sqrt(vert[0]*vert[0] + vert[1]*vert[1] + vert[2]*vert[2]);
    
    vert[0] = vert[0] / normV;
    vert[1] = vert[1] / normV;
    vert[2] = vert[2] / normV;
    
    //projecting the spin into the plane perpendicular to the B field.
    double spinProjIntoPerpPlane[3] = {0};
    
    double dotSB = (spin[0]*Bfield[0] +
                    spin[1]*Bfield[1] + spin[2]*Bfield[2]);
    normB = sqrt(Bfield[0]*Bfield[0] +
                 Bfield[1]*Bfield[1] + Bfield[2]*Bfield[2]);
    
    spinProjIntoPerpPlane[0] = spinX - (dotSB * Bfield[0]/normB);
    spinProjIntoPerpPlane[1] = spinY - (dotSB * Bfield[1]/normB);
    spinProjIntoPerpPlane[2] = spinZ - (dotSB * Bfield[2]/normB);
    
    
    double normSP = sqrt(spinProjIntoPerpPlane[0]*spinProjIntoPerpPlane[0] +
                         spinProjIntoPerpPlane[1]*spinProjIntoPerpPlane[1] +
                         spinProjIntoPerpPlane[2]*spinProjIntoPerpPlane[2]);
    spinProjIntoPerpPlane[0] = spinProjIntoPerpPlane[0] / normSP;
    spinProjIntoPerpPlane[1] = spinProjIntoPerpPlane[1] / normSP;
    spinProjIntoPerpPlane[2] = spinProjIntoPerpPlane[2] / normSP;
    
    //At this point, the dot product of spinProjIntoPerpPlane with nhat will
    //yield cos(AB). But this is sign insensitive - choose nhat to be the "x"
    //and ccw, [0,2pi) is the AB angle. This means that when the dot product of
    //vert and spinProjIntoPerpPlane are greater than 0, then 0 < AB < pi and
    //vert and spinProjIntoPerpPlane are less than 0, then pi < AB < 2pi.
    
    double cosAB = (spinProjIntoPerpPlane[0]*nhat[0] +
                    spinProjIntoPerpPlane[1]*nhat[1] +
                    spinProjIntoPerpPlane[2]*nhat[2]);
    double dotPV = (spinProjIntoPerpPlane[0]*vert[0] +
                    spinProjIntoPerpPlane[1]*vert[1] +
                    spinProjIntoPerpPlane[2]*vert[2]);
    if(dotPV == 0){
        if(cosAB > 0){
            return 0;
        }else{
            return PI;
        }
    }
    
    cosAB = (spinProjIntoPerpPlane[0]*nhat[0] +
             spinProjIntoPerpPlane[1]*nhat[1] +
             spinProjIntoPerpPlane[2]*nhat[2]);
    if((cosAB > 1.0000001)||(cosAB < -1.00000001)){
        cout<<"bad AB "<<cosAB<<endl;
        cerr<<"bad AB "<<cosAB<<endl;
        exit(1);
    }
    if((cosAB > 1.0)){
        cosAB = 1;
    }
    if((cosAB < -1.0)){
        cosAB = -1;
    }
    double OPAB = acos(cosAB);
    
    dotPV = (spinProjIntoPerpPlane[0]*vert[0] +
             spinProjIntoPerpPlane[1]*vert[1] +
             spinProjIntoPerpPlane[2]*vert[2]);
    if(dotPV < 0){
        OPAB = 2 * PI - OPAB;
    }
    
    return OPAB;
    
    //End OP: AB
}

