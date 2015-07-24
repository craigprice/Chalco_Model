/**
 chalco.cpp
 
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
#include "MonteCarlo.h"
#include "Chalco_Square.h"

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
            return sc;
            break;
            
        case 2:
            sc.s = 3;
            return sc;
            break;
            
        case 1:
            sc.s = 0;
            if(sc.a ==  MCP.cellsA - 1){
                sc.a = 0;
                return sc;
            }else{
                sc.a = sc.a + 1;
                return sc;
            }
            break;
            
        case 3:
            sc.s = 2;
            if(sc.a ==  MCP.cellsA - 1){
                sc.a = 0;
                return sc;
            }else{
                sc.a = sc.a + 1;
                return sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheMinusXDir(SpinCoordinates sc) const{
    switch (sc.s) {
        case 1:
            sc.s = 0;
            return sc;
            break;
            
        case 3:
            sc.s = 2;
            return sc;
            break;
            
        case 0:
            sc.s = 1;
            if(sc.a == 0){
                sc.a = MCP.cellsA - 1;
                return sc;
            }else{
                sc.a = sc.a - 1;
                return sc;
            }
            break;
            
        case 2:
            sc.s = 3;
            if(sc.a == 0){
                sc.a = MCP.cellsA - 1;
                return sc;
            }else{
                sc.a = sc.a - 1;
                return sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
}



SpinCoordinates MonteCarlo::getSpinCoordsInTheMinusYDir(SpinCoordinates sc) const{
    switch (sc.s) {
        case 2:
            sc.s = 0;
            return sc;
            break;
            
        case 3:
            sc.s = 1;
            return sc;
            break;
            
        case 0:
            sc.s = 2;
            if(sc.b == 0){
                sc.b = MCP.cellsB - 1;
                return sc;
            }else{
                sc.b = sc.b - 1;
                return sc;
            }
            break;
            
        case 1:
            sc.s = 3;
            if(sc.b == 0){
                sc.b = MCP.cellsB - 1;
                return sc;
            }else{
                sc.b = sc.b - 1;
                return sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheYDir(SpinCoordinates sc) const{
    switch (sc.s) {
        case 0:
            sc.s = 2;
            return sc;
            break;
            
        case 1:
            sc.s = 3;
            return sc;
            break;
            
        case 2:
            sc.s = 0;
            if(sc.b == MCP.cellsB - 1){
                sc.b = 0;
                return sc;
            }else{
                sc.b = sc.b + 1;
                return sc;
            }
            break;
            
        case 3:
            sc.s = 1;
            if(sc.b == MCP.cellsB - 1){
                sc.b = 0;
                return sc;
            }else{
                sc.b = sc.b + 1;
                return sc;
            }
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
            return sc;
        }else{
            sc.c = sc.c + 1;
            return sc;
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
            return sc;
        }else{
            sc.c = sc.c - 1;
            return sc;
        }
    }
    return sc;
    
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
                    //[1]: first NN in -x direction
                    //[2]: first NN in y direction
                    //[3]: first NN in -y direction
                    //[4]: second NN, going first +x direction, then +y direction
                    //[5]: second NN, going first -x direction, then +y direction
                    //[6]: second NN, going first -x direction, then -y direction
                    //[7]: second NN, going first +x direction, then -y direction
                    //[8]:  third NN, going first +x direction, then +x direction
                    //[9]:  third NN, going first -x direction, then -x direction
                    //[10]: third NN, going first +y direction, then +y direction
                    //[11]: third NN, going first -y direction, then -y direction
                    
                    sc.a = a;
                    sc.b = b;
                    sc.c = c;
                    sc.s = s;
                    neighbors[c][a][b][s][0] = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][1] = getSpinCoordsInTheMinusXDir(sc);
                    neighbors[c][a][b][s][2] = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][3] = getSpinCoordsInTheMinusYDir(sc);
                    
                    
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][4] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(sc);
                    neighbors[c][a][b][s][5] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(sc);
                    neighbors[c][a][b][s][6] = getSpinCoordsInTheMinusYDir(sc);
                    
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][7] = getSpinCoordsInTheMinusYDir(sc);
                    
                    
                    sc = getSpinCoordsInTheXDir(sc);
                    neighbors[c][a][b][s][8] = getSpinCoordsInTheXDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(sc);
                    neighbors[c][a][b][s][9] = getSpinCoordsInTheMinusXDir(sc);
                    
                    sc = getSpinCoordsInTheYDir(sc);
                    neighbors[c][a][b][s][10] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusYDir(sc);
                    neighbors[c][a][b][s][11] = getSpinCoordsInTheMinusYDir(sc);
                    
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


//sublattice structure is:
//2,3,2,3
//0,1,0,1
//2,3,2,3
void MonteCarlo::updateOrderParameters(){
    
    double E = 0;
    double X = 0;
    double Y = 0;
    double Z = 0;
    double PX = 0;
    double PY = 0;
    double PXmX = 0;
    double PYmY = 0;
    double Q = 0;
    double M[SPINDIMEN] = {0};
    double N[SPINDIMEN] = {0};
    double SX[SPINDIMEN] = {0};
    double SY[SPINDIMEN] = {0};
    double Dimer = 0;
    SpinCoordinates sc;
    ClassicalSpin3D a1_x;
    ClassicalSpin3D a1_mx;
    ClassicalSpin3D a1_y;
    ClassicalSpin3D a1_my;
    
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
                    //if k == 0,3 - count negative M, k == 1,2 - positive M
                    N[0] += ((s == 0 || s == 3) ? ((-1.0)*X) : (X));
                    N[1] += ((s == 0 || s == 3) ? ((-1.0)*Y) : (Y));
                    N[2] += ((s == 0 || s == 3) ? ((-1.0)*Z) : (Z));
                    //End OP: N
                    
                    //OP: SX
                    if((b % 2 == 0)){
                        SX[0] += X;
                        SX[1] += Y;
                        SX[2] += Z;
                    }else if(b % 2 == 1){
                        SX[0] += (-1.0) * X;
                        SX[1] += (-1.0) * Y;
                        SX[2] += (-1.0) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: SX
                    
                    //OP: SY
                    if((a % 2 == 0)){
                        SY[0] += X;
                        SY[1] += Y;
                        SY[2] += Z;
                    }else if((a % 2 == 1)){
                        SY[0] += (-1.0) * X;
                        SY[1] += (-1.0) * Y;
                        SY[2] += (-1.0) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: S2
                    
                    /*
                     a1_x = getSpinInTheXDir(lattice[i][j][k][c]);
                     a1_y = getSpinInTheYDir(lattice[i][j][k][c]);
                     */
                    
                    sc = neighbors[c][a][b][s][0];//NN_x
                    a1_x = lattice[sc.c][sc.a][sc.b][sc.s];
                    sc = neighbors[c][a][b][s][1];//NN_mx
                    a1_mx = lattice[sc.c][sc.a][sc.b][sc.s];
                    sc = neighbors[c][a][b][s][2];//NN_y
                    a1_y = lattice[sc.c][sc.a][sc.b][sc.s];
                    sc = neighbors[c][a][b][s][3];//NN_my
                    a1_my = lattice[sc.c][sc.a][sc.b][sc.s];
                    
                    PX += (X * a1_x.x + Y * a1_x.y + Z * a1_x.z);
                    PY += (X * a1_y.x + Y * a1_y.y + Z * a1_y.z);
                    PXmX += (a1_x.x * a1_mx.x + a1_x.y * a1_mx.y + a1_x.z * a1_mx.z);
                    PYmY += (a1_y.x * a1_my.x + a1_y.y * a1_my.y + a1_y.z * a1_my.z);
                    
                }
            }
        }
    }
    Q = PX - PY;
    
    //S[A-1,B].S[A+1,B]-S[A,B+1].S[A,B-1]
    Dimer = PXmX - PYmY;
    //End OP
    
    OP.scalarOP[0] = E; //Total E, E not per site
    OP.scalarOP[1] = PX / (1.0*numSites);
    OP.scalarOP[2] = PY / (1.0*numSites);
    OP.scalarOP[3] = Q / (1.0*numSites);
    OP.scalarOP[4] = Dimer / (1.0*numSites);
    
    //All Order parameters are per site
    for(int i = 0; i < SPINDIMEN; i++){
        M[i] = M[i] / (1.0*numSites);
        N[i] = N[i] / (1.0*numSites);
        SX[i] = SX[i] / (1.0*numSites);
        SY[i] = SY[i] / (1.0*numSites);
    }
    
    //Now to calculate the characteristics of the order parameters.
    double OPArr[NUMVECTOROP][SPINDIMEN] = {0};
    for(int i = 0; i < SPINDIMEN; i++){
        OPArr[0][i] = M[i];
        OPArr[1][i] = N[i];
        OPArr[2][i] = SX[i];
        OPArr[3][i] = SY[i];
        OPArr[4][i] = sqrt(SX[i]*SX[i] + SY[i]*SY[i]);
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        
        //Scalar quantity eg. M
        OP.vectorOP[i][0] = sqrt(OPArr[i][0]*OPArr[i][0] +
                                 OPArr[i][1]*OPArr[i][1] +
                                 OPArr[i][2]*OPArr[i][2]);
        
        
        //fabs x, y, z, directions
        OP.vectorOP[i][1] = fabs(OPArr[i][0]);
        
        OP.vectorOP[i][2] = fabs(OPArr[i][1]);
        
        OP.vectorOP[i][3] = fabs(OPArr[i][2]);
        
        
        
    }
    
}

void MonteCarlo::updateCorrelationFunction(){
    
    //Total E, E not per site
    double E[NUMCORRSPINSBIG] = {0};
    //dot product of spin in the x direction
    double PX[NUMCORRSPINSBIG] = {0};
    //dot product of spin in the y direction
    double PY[NUMCORRSPINSBIG] = {0};
    //dot product of spin in the x direction with spin in -x direction
    double PXmX[NUMCORRSPINSBIG] = {0};
    //dot product of spin in the y direction with spin in -y direction
    double PYmY[NUMCORRSPINSBIG] = {0};
    //Nematic order parameter. <PX-PY>
    double Q[NUMCORRSPINSBIG] = {0};
    //Dimer Order parameter. S[A-1,B].S[A+1,B]-S[A,B+1].S[A,B-1]
    double Dimer[NUMCORRSPINSBIG] = {0};
    double X = 0;
    double Y = 0;
    double Z = 0;
    double M[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double N[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double SX[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double SY[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    SpinCoordinates sc;
    ClassicalSpin3D a1_x;
    ClassicalSpin3D a1_mx;
    ClassicalSpin3D a1_y;
    ClassicalSpin3D a1_my;
    
    
    uint_fast8_t b = 0;
    uint_fast8_t a = 0;
    uint_fast8_t s = 0;
    uint_fast8_t c = 0;
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        
        //This control structure picks the diagonal spins out across the lattice.
        if(count % 2 == 0 ){
            s = 0;
            a = (uint_fast8_t) (count/2.0);
            b = (uint_fast8_t) (count/2.0);
        }else{
            s = 3;
            a = (uint_fast8_t) (count/2.0);
            b = (uint_fast8_t) (count/2.0);
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
        N[0][count] += ((s == 0 || s == 3) ? ((-1.0)*X) : (X));
        N[1][count] += ((s == 0 || s == 3) ? ((-1.0)*Y) : (Y));
        N[2][count] += ((s == 0 || s == 3) ? ((-1.0)*Z) : (Z));
        //End OP: N
        
        //OP: SX
        if((b % 2 == 0)){
            SX[0][count] += X;
            SX[1][count] += Y;
            SX[2][count] += Z;
        }else if((b % 2 == 1)){
            SX[0][count] += (-1.0) * X;
            SX[1][count] += (-1.0) * Y;
            SX[2][count] += (-1.0) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: SX
        
        //OP: SY
        if((a % 2 == 0)){
            SY[0][count] += X;
            SY[1][count] += Y;
            SY[2][count] += Z;
        }else if((a % 2 == 1)){
            SY[0][count] += (-1) * X;
            SY[1][count] += (-1) * Y;
            SY[2][count] += (-1) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: SY
        
        /*
         a1_x = getSpinInTheXDir(lattice[i][j][k][c]);
         a1_y = getSpinInTheYDir(lattice[i][j][k][c]);
         */
        
        
        sc = neighbors[c][a][b][s][0];//NN_x
        a1_x = lattice[sc.c][sc.a][sc.b][sc.s];
        sc = neighbors[c][a][b][s][1];//NN_mx
        a1_mx = lattice[sc.c][sc.a][sc.b][sc.s];
        sc = neighbors[c][a][b][s][2];//NN_y
        a1_y = lattice[sc.c][sc.a][sc.b][sc.s];
        sc = neighbors[c][a][b][s][3];//NN_my
        a1_my = lattice[sc.c][sc.a][sc.b][sc.s];
        
        PX[count] += (X * a1_x.x + Y * a1_x.y + Z * a1_x.z);
        PY[count] += (X * a1_y.x + Y * a1_y.y + Z * a1_y.z);
        PXmX[count] += (a1_x.x * a1_mx.x + a1_x.y * a1_mx.y + a1_x.z * a1_mx.z);
        PYmY[count] += (a1_y.x * a1_my.x + a1_y.y * a1_my.y + a1_y.z * a1_my.z);
        Q[count] += PX[count] - PY[count];
        
        //S[A-1,B].S[A+1,B]-S[A,B+1].S[A,B-1]
        Dimer[count] += PXmX[count] - PYmY[count];
        
    }
    
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        OP.scalarOPCorrFunc[0][count] = E[0] * E[count];
        OP.scalarOPCorrFunc[1][count] = PX[0] * PX[count];
        OP.scalarOPCorrFunc[2][count] = PY[0] * PY[count];
        OP.scalarOPCorrFunc[3][count] = Q[0] * Q[count];
        OP.scalarOPCorrFunc[4][count] = Dimer[0] * Dimer[count];
    }
    
    double OPArr[NUMVECTOROP][SPINDIMEN][NUMCORRSPINSBIG] = {0};
    for(int i = 0; i < SPINDIMEN; i++){
        for(int j = 0; j < 2*MCP.cellsA; j++){
            OPArr[0][i][j] = M[i][j];
            OPArr[1][i][j] = N[i][j];
            OPArr[2][i][j] = SX[i][j];
            OPArr[3][i][j] = SY[i][j];
            OPArr[4][i][j] = sqrt(SX[i][j]*SX[i][j] +
                                  SY[i][j]*SY[i][j]);
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int count = 0; count < 2*MCP.cellsA; count++){
            OP.vectorOPCorrFunc[i][0][count] = (OPArr[i][0][0] * OPArr[i][0][count] +
                                                OPArr[i][1][0] * OPArr[i][1][count] +
                                                OPArr[i][2][0] * OPArr[i][2][count]);
            
            
            
            //spin components
            OP.vectorOPCorrFunc[i][1][count] = OPArr[i][0][0] * OPArr[i][0][count];
            
            OP.vectorOPCorrFunc[i][2][count] = OPArr[i][1][0] * OPArr[i][1][count];
            
            OP.vectorOPCorrFunc[i][3][count] = OPArr[i][2][0] * OPArr[i][2][count];
            
            
            
        }
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
    
}

//sublattice structure is:
//0,1,0,1
//2,3,2,3
//0,1,0,1
void MonteCarlo::updateFourierTransformOnRecipLattice(){
    
    //cout<<MD.numConfigsDone<<endl;
    //cout<<FTTOUPDATES<<endl;
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
                for(uint_fast8_t k_s = 0; k_s < NUMSUBLATTICES; k_s++){
                    
                    //Loop over all sites in real space
                    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
                        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
                            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                                    
                                    kvector_x = (2*k_a + (k_s%2)) * B_VECTOR_1[0] / (MCP.cellsA*1.0) +
                                    (2*k_b + ((k_s-1)>0)) * B_VECTOR_2[0] / (MCP.cellsB*1.0);
                                    kvector_y = (2*k_a + (k_s%2)) * B_VECTOR_1[1] / (MCP.cellsA*1.0) +
                                    (2*k_b + ((k_s-1)>0)) * B_VECTOR_2[1] / (MCP.cellsB*1.0);
                                    rvector_x = (2*a + (s%2)) * A_VECTOR_1[0] +
                                    (2*b + ((s-1)>0)) * A_VECTOR_2[0];
                                    rvector_y = (2*a + (s%2)) * A_VECTOR_1[1] +
                                    (2*b + ((s-1)>0)) * A_VECTOR_2[1];
                                    kr = kvector_x * rvector_x + kvector_y * rvector_y;
                                    cos_kr = cos(kr);
                                    sin_kr = sin(kr);
                                    recipLattice[k_c][k_a][k_b][k_s].ReXComponent +=
                                    cos_kr * lattice[c][a][b][s].x;
                                    recipLattice[k_c][k_a][k_b][k_s].ImXComponent +=
                                    sin_kr * lattice[c][a][b][s].x;
                                    recipLattice[k_c][k_a][k_b][k_s].ReYComponent +=
                                    cos_kr * lattice[c][a][b][s].y;
                                    recipLattice[k_c][k_a][k_b][k_s].ImYComponent +=
                                    sin_kr * lattice[c][a][b][s].y;
                                    recipLattice[k_c][k_a][k_b][k_s].ReZComponent +=
                                    cos_kr * lattice[c][a][b][s].z;
                                    recipLattice[k_c][k_a][k_b][k_s].ImZComponent +=
                                    sin_kr * lattice[c][a][b][s].z;
                                    
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
    
    //The order of neighbors are:
    //[0]: first NN in x direction
    //[1]: first NN in -x direction
    //[2]: first NN in y direction
    //[3]: first NN in -y direction
    //[4]: second NN, going first in the +x direction, then the +y direction
    //[5]: second NN, going first in the -x direction, then the +y direction
    //[6]: second NN, going first in the -x direction, then the -y direction
    //[7]: second NN, going first in the +x direction, then the -y direction
    //[8]:  third NN, going first in the +x direction, then the +x direction
    //[9]:  third NN, going first in the -x direction, then the -x direction
    //[10]: third NN, going first in the +y direction, then the +y direction
    //[11]: third NN, going first in the -y direction, then the -y direction
    
    
    double temp = (dotProd(c, a, b, s, neighbors[c][a][b][s][0]) +
                   dotProd(c, a, b, s, neighbors[c][a][b][s][1]) +
                   dotProd(c, a, b, s, neighbors[c][a][b][s][2]) +
                   dotProd(c, a, b, s, neighbors[c][a][b][s][3]));
    
    double interactionH = MCP.j1 * (temp);
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (dotProd(c, a, b, s, neighbors[c][a][b][s][4]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][5]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][6]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][7]));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(c, a, b, s, neighbors[c][a][b][s][8]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][9]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][10]) +
                                     dotProd(c, a, b, s, neighbors[c][a][b][s][11]));
    
    double interactionK = MCP.k1 * temp * temp;
    
    double interactionCubic = MCP.cubicD * (-1.0) *
    lattice[c][a][b][s].z * lattice[c][a][b][s].z;
    
    if(fromSpinFlip){
        totalInteraction += interactionH + interactionH2 + interactionH3 +
        interactionK + interactionCubic;
    }else{
        //From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 +
                             interactionK + interactionCubic) / 2.0;
    }
    
    
    //B field
    if (MCP.isBField){
        double interactionB = (-1) * (MCP.bField_xNorm * lattice[c][a][b][s].x +
                                      MCP.bField_yNorm * lattice[c][a][b][s].y +
                                      MCP.bField_zNorm * lattice[c][a][b][s].z);
        totalInteraction += interactionB;
    }
    
    return totalInteraction;
    
}

double MonteCarlo::getLocalEnergyFromCalc(uint_fast8_t c, uint_fast8_t a,
                                          uint_fast8_t b, uint_fast8_t s,
                                          const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
    SpinCoordinates sc;
    sc.a = a;
    sc.b = b;
    sc.c = c;
    sc.s = s;
    double temp = (dotProd(c, a, b, s, getSpinCoordsInTheXDir(sc)) +
                   dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(sc)) +
                   dotProd(c, a, b, s, getSpinCoordsInTheYDir(sc)) +
                   dotProd(c, a, b, s, getSpinCoordsInTheMinusYDir(sc)));
    
    double interactionH = MCP.j1 * (temp);
    
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 *
    (dotProd(c, a, b, s, getSpinCoordsInTheYDir(getSpinCoordsInTheXDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheYDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusYDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheXDir(getSpinCoordsInTheMinusYDir(sc))));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 *
    (dotProd(c, a, b, s, getSpinCoordsInTheXDir(getSpinCoordsInTheXDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusXDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheYDir(getSpinCoordsInTheYDir(sc))) +
     dotProd(c, a, b, s, getSpinCoordsInTheMinusYDir(getSpinCoordsInTheMinusYDir(sc))));
    
    double interactionK = MCP.k1 * temp * temp;
    
    double interactionCubic = MCP.cubicD * (-1.0) *
    lattice[c][a][b][s].z * lattice[c][a][b][s].z;
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 +
        interactionK + interactionCubic;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 +
                             interactionK + interactionCubic) / 2.0;
    }
    
    
    //B field
    if (MCP.isBField) {
        double interactionB = (-1) * (MCP.bField_xNorm * lattice[c][a][b][s].x +
                                      MCP.bField_yNorm * lattice[c][a][b][s].y +
                                      MCP.bField_zNorm * lattice[c][a][b][s].z);
        totalInteraction += interactionB;
    }
    
    return totalInteraction;
    
}

