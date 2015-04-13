
/**
 MonteCarlo.cpp
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of Chalcogenides. Specifically, exploring the
 magnetic phases at finite temperature.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/04/12
 */

#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include "MonteCarlo.h"
#include "FileIO.cpp"
using namespace std;


/******************************************************************************/
//Functions for initialization
/******************************************************************************/



void MonteCarlo::setRSToZero(){
    cout<<"Begin setRSToZero"<<endl;
    
    //Setting Lat Chars to zero.
    for(int i = 0; i < NUMSCALAROP; i++){
        for(int k = 0; k < POWERSOP; k++){
            RS.sumOverScalarOP[i][k] = 0;
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                RS.sumOverVectorOP[i][j][k] = 0;
            }
        }
    }
    
    //Setting correlation function to zero
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        for(int i = 0; i < NUMSCALAROP; i++){
            for(int k = 0; k < POWERSOP; k++){
                RS.sumOverScalarOPCorrFunc[i][k][count] = 0;
            }
        }
    }
    
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                for(int k = 0; k < POWERSOP; k++){
                    RS.sumOverVectorOPCorrFunc[i][j][k][count] = 0;
                }
            }
        }
    }
    
}

void MonteCarlo::setOPToZero(){
    cout<<"Begin setOPToZero"<<endl;
    for(int i = 0; i < NUMVECTOROP; i++){
        OP.scalarOP[i] = 0;
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            OP.vectorOP[i][j] = 0;
        }
    }
    
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        for(int i = 0; i < NUMSCALAROP; i++){
            OP.scalarOPCorrFunc[i][count] = 0;
        }
    }
    
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                OP.vectorOPCorrFunc[i][j][count] = 0;
            }
        }
    }
    
}

void MonteCarlo::setMCPToZero(){
    cout<<"Begin setMCPToZero"<<endl;
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

void MonteCarlo::setMDToZero(){
    cout<<"Begin setMDToZero"<<endl;
    
    MD.flipAttempts    = 0;
    MD.successfulFlips = 0;
    MD.range           = PI-0.1;
    MD.rangeOld        = 0.1;
    MD.numConfigsDone  = 0;
    MD.isThermal        = false;
    MD.numSweepsPerformedToCheckIfThermalized = 0;
    MD.numSweepsPerformed           = 0;
    MD.numSweepsUsedToThermalize    = 0;
    MD.timeOfInitialization         = 0;
    
}

void MonteCarlo::setRecipToZero(){
    cout<<"Begin setRecipToZero"<<endl;
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    sumOverRecipLattice[c][a][b][s].ReXComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ReYComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ReZComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImXComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImYComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImZComponent = 0;
                    sumOverRecipLattice[c][a][b][s].magnitude = 0;
                    
                    recipLattice[c][a][b][s].ReXComponent = 0;
                    recipLattice[c][a][b][s].ReYComponent = 0;
                    recipLattice[c][a][b][s].ReZComponent = 0;
                    recipLattice[c][a][b][s].ImXComponent = 0;
                    recipLattice[c][a][b][s].ImYComponent = 0;
                    recipLattice[c][a][b][s].ImZComponent = 0;
                    recipLattice[c][a][b][s].magnitude = 0;
                }
            }
        }
    }
}

void MonteCarlo::initializeLattice(){
    cout<<"Begin initializeLattice"<<endl;
    
    ClassicalSpin3D spin;
    for(int c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<ClassicalSpin3D> > > row1;
        for(int a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<ClassicalSpin3D> > row2;
            for(int b = 0; b < MCP.cellsB; b++){
                std::vector<ClassicalSpin3D> row3;
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    row3.push_back(spin);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        lattice.push_back(row1);
    }
    
    cout << "neighbors" << endl;
    SpinCoordinates sc;
    for(int c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector< std::vector<SpinCoordinates> > > > row1;
        for(int a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector< std::vector<SpinCoordinates> > > row2;
            for(int b = 0; b < MCP.cellsB; b++){
                std::vector< std::vector<SpinCoordinates> > row3;
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    std::vector<SpinCoordinates> row4;
                    for(int n = 0; n < NUMNEIGHBORS; n++){
                        row4.push_back(sc);
                    }
                    row3.push_back(row4);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        neighbors.push_back(row1);
    }
    
    
    cout << "recipLattice" << endl;
    FourierTransformOfSpin spinFT;
    for(int c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row1;
        for(int a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<FourierTransformOfSpin> > row2;
            for(int b = 0; b < MCP.cellsB; b++){
                std::vector<FourierTransformOfSpin> row3;
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    row3.push_back(spinFT);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        recipLattice.push_back(row1);
    }
    
    cout << "sumOverRecipLattice" << endl;
    for(int c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row1;
        for(int a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<FourierTransformOfSpin> > row2;
            for(int b = 0; b < MCP.cellsB; b++){
                std::vector<FourierTransformOfSpin> row3;
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    row3.push_back(spinFT);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        sumOverRecipLattice.push_back(row1);
    }
    
    cout << "init" << endl;
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    lattice[c][a][b][s].x = 0;
                    lattice[c][a][b][s].y = 0;
                    lattice[c][a][b][s].z = 1;
                    setRandomOrientation(c, a, b, s);
                    
                    recipLattice[c][a][b][s].ReXComponent = 0;
                    recipLattice[c][a][b][s].ReYComponent = 0;
                    recipLattice[c][a][b][s].ReZComponent = 0;
                    recipLattice[c][a][b][s].ImXComponent = 0;
                    recipLattice[c][a][b][s].ImYComponent = 0;
                    recipLattice[c][a][b][s].ImZComponent = 0;
                    recipLattice[c][a][b][s].magnitude = 0;
                    
                    sumOverRecipLattice[c][a][b][s].ReXComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ReYComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ReZComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImXComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImYComponent = 0;
                    sumOverRecipLattice[c][a][b][s].ImZComponent = 0;
                    sumOverRecipLattice[c][a][b][s].magnitude = 0;
                    
                    for(uint_fast8_t n = 0; n < NUMNEIGHBORS; n++){
                        neighbors[c][a][b][s][n].a = a;
                        neighbors[c][a][b][s][n].b = b;
                        neighbors[c][a][b][s][n].s = s;
                        neighbors[c][a][b][s][n].c = c;
                    }
                    
                }
            }
        }
    }
    
    setNearestNeighbors();
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
                    //[4]: second NN, going first in the x direction, then the y direction
                    //[5]: second NN, going first in the -x direction, then the y direction
                    //[6]: second NN, going first in the -x direction, then the -y direction
                    //[7]: second NN, going first in the x direction, then the -y direction
                    //[8]: third NN, going first in the x direction, then the x direction
                    //[9]: third NN, going first in the -x direction, then the -x direction
                    //[10]: third NN, going first in the y direction, then the y direction
                    //[11]: third NN, going first in the -y direction, then the -y direction
                    
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


void MonteCarlo::init(){
    cout<<"Begin init"<<endl;
    
    //initialize private variables
    
    setRSToZero();
    
    setOPToZero();
    
    setMCPToZero();
    
    setMDToZero();
    
    setParamNames();
    
    outputFileName = "";
    
    
    cout<<"Finished with init()"<<endl;
    
}

MonteCarlo::MonteCarlo(){
    cout<<"Begin MonteCarlo"<<endl;
    init();
    
    initializeLattice();
    numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES;
    
    
    cout<<"Finished with default MC constructor"<<endl;
    MD.timeOfInitialization = clock();
    
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
    
    MD.timeOfInitialization = clock();
    
    cout<<"Finished with (new) MC constructor"<<endl;
    
}

//This function initializes the lattice through copying data from a previous run.
MonteCarlo::MonteCarlo(bool reThermalize, string filename){
    string function = "MonteCarlo::MonteCarlo(bool reThermalize, string filename)";
    cout << "Begin " << function.c_str() << endl;
    
    init();
    
    //Reading in the rest of the lattice information from the saved file
    outputFileName = filename;
    
    ifstream readline;
    readline.open(outputFileName.c_str(), ifstream::in);
    
    if (!readline){
        cerr << "Can't open file: " << endl;
        cerr << outputFileName <<endl;
        exit(1);
    }
    
    string line = "";
    stringstream ss;
    while(readline.good()){
        getline (readline,line);
        
        
        if (line.compare("Begin MonteCarlo::toStringMCP()") == 0){
            readStringMCP(readline);
        }
        
        if (line.compare("Begin MonteCarlo::toStringMD()") == 0){
            readStringMD(readline);
        }
        
        if ((!reThermalize)&&
            (MD.numSweepsPerformed >= MCP.numSweepsToPerformTotal)) {
            cout<<"!!!Finished Simulation!!!"<<endl;
            exit(0);
        }
        
        
        if (line.compare("Begin MonteCarlo::toStringRS()") == 0){
            initializeLattice();
            readStringRS(readline);
        }
        
        if (line.compare("Begin MonteCarlo::toStringSpinSnapshot()") == 0){
            readStringSpinSnapshot(readline);
        }
        
        if (line.compare("Begin MonteCarlo::toStringRSCorr()") == 0){
            readStringRSCorr(readline);
        }
        
        if (line.compare("Begin MonteCarlo::toStringSumOverFourierTransform()") == 0){
            readStringSumOverFourierTransform(readline);
            break;
        }
        
    }
    
    
    cout<<"Finished with reading in all lattice chars"<<endl;
    readline.close();
    
    
    
    if(MCP.isBField){
        double size = sqrt(MCP.bField_x * MCP.bField_x +
                           MCP.bField_y * MCP.bField_y +
                           MCP.bField_z * MCP.bField_z);
        
        MCP.bField_xNorm = MCP.bFieldMag * MCP.bField_x / (1.0*size);
        MCP.bField_yNorm = MCP.bFieldMag * MCP.bField_y / (1.0*size);
        MCP.bField_zNorm = MCP.bFieldMag * MCP.bField_z / (1.0*size);
    }
    
    
    cout<<"Finished with reading in all components"<<endl;
    
    //cout<< toString()<<endl;
    
    cout<<"Finished with MC constructor"<<endl;
    
    numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES;
    updateOrderParameters();
    
    if(reThermalize){
        cout<< "ReThermalizing. Setting All Char to Zero and setting ";
        cout<< "Thermalization data to Zero." << endl;
        setRSToZero();
        setOPToZero();
        setMDToZero();
        setRecipToZero();
        
        for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
            for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
                for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                    for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                        
                        setRandomOrientation(c, a, b, s);
                        
                        recipLattice[c][a][b][s].ReXComponent = 0;
                        recipLattice[c][a][b][s].ReYComponent = 0;
                        recipLattice[c][a][b][s].ReZComponent = 0;
                        recipLattice[c][a][b][s].ImXComponent = 0;
                        recipLattice[c][a][b][s].ImYComponent = 0;
                        recipLattice[c][a][b][s].ImZComponent = 0;
                        recipLattice[c][a][b][s].magnitude = 0;
                        sumOverRecipLattice[c][a][b][s].ReXComponent = 0;
                        sumOverRecipLattice[c][a][b][s].ReYComponent = 0;
                        sumOverRecipLattice[c][a][b][s].ReZComponent = 0;
                        sumOverRecipLattice[c][a][b][s].ImXComponent = 0;
                        sumOverRecipLattice[c][a][b][s].ImYComponent = 0;
                        sumOverRecipLattice[c][a][b][s].ImZComponent = 0;
                        sumOverRecipLattice[c][a][b][s].magnitude = 0;
                    }
                }
            }
        }
        
        MD.timeOfInitialization = clock();
        
    }
    
    
}


/******************************************************************************/
//Thermalization
/******************************************************************************/

void MonteCarlo::adjustRange(){
    
    double goal = 0.65;
    double successRate = MD.successfulFlips / (MD.flipAttempts * 1.0);
    double diff = successRate - goal;
    if(fabs(diff) < 0.01){
        return;
    }
    
    /**
     *This uses a binary search to change the range of the spin flips until
     *the rate of successful flips is ~goal.
     */
    double temp = fabs(MD.range - MD.rangeOld) / 2.0;
    if(temp == 0){temp = 0.1;}
    
    MD.rangeOld = MD.range;
    if(diff > 0){
        MD.range += temp;
        if(MD.range > PI){
            MD.range = PI;
            return;
        }
    }else{
        MD.range -= temp;
        if(MD.range < 0.01){
            MD.range = 0.01;
            return;
        }
    }
    MD.successfulFlips = 0;
    MD.flipAttempts = 0;
}

bool MonteCarlo::thermalize(double KbT_, int numSweepsToDo, double durationOfSingleRun){
    
    /*
     double totalInteraction = 0;
     
     
     clock_t t0 = clock();
     
     totalInteraction = getLocalEnergy(lattice[0][0][0][0],true);
     
     cout<< clock() - t0 << endl;
     t0 = clock();
     
     totalInteraction = getLocalEnergy0(lattice[0][0][0][0],true);
     
     cout << clock() - t0 << endl;
     totalInteraction = 10;
     cout << totalInteraction << endl;
     exit(0);
     */
    
    double oldKbT = MCP.KbT;
    MCP.KbT = KbT_;
    
    
    for(int sweeps = MD.numSweepsUsedToThermalize; sweeps <= numSweepsToDo; sweeps++){
        unsigned long clocks = clock() - MD.timeOfInitialization;
        double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60.0 * 60.0);
        if(hoursElapsed > durationOfSingleRun){
            break;
        }
        
        metropolisSweep();
        MD.numSweepsPerformed += 1;
        MD.numSweepsUsedToThermalize += 1;
        adjustRange();
        
        
        
    }
    
    MCP.KbT = oldKbT;
    
    if(MD.numSweepsUsedToThermalize >= numSweepsToDo){
        cout<<"Returning from Thermalize. Finished."<<endl;
        return false;
    }
    
    cout<<"Returning from Thermalize. Out of Time."<<endl;
    return true;//true if out of time.
}

bool MonteCarlo::isThermalized(double durationOfSingleRun){//return true if out of time.
    
    
    for(int i = 0; i < 10; i++){
        metropolisSweep();
        MD.numSweepsPerformed += 1;
        MD.numSweepsUsedToThermalize += 1;
        adjustRange();
        updateRunningSumOP();//This function is the only place this is called except sweep()
        MD.numSweepsPerformedToCheckIfThermalized += 1;
    }
    
    bool outOfTime = false;
    int thermalCycles = -1;
    while(!outOfTime){
        thermalCycles ++;
        cout<<"Thermal Cycle # " << thermalCycles<<endl;
        
        double clocks = clock() - MD.timeOfInitialization;
        double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60.0 * 60.0);
        if(hoursElapsed > durationOfSingleRun){
            outOfTime = true;
            return outOfTime;
        }
        
        RunningSumOverOP firstRS = RS;
        
        int numSweepsAveOver = (int) (2000);
        for(int i = 0; i < numSweepsAveOver; i++){
            metropolisSweep();
            MD.numSweepsPerformed += 1;
            MD.numSweepsUsedToThermalize += 1;
            adjustRange();
            updateRunningSumOP();
            MD.numSweepsPerformedToCheckIfThermalized += 1;
        }
        
        RunningSumOverOP secRS = RS;
        
        double firstBinder = 0;
        double secBinder = 0;
        bool smallDiffBinder = true;
        for(int i = 0; i < NUMVECTOROP; i++){
            double m4_1 = firstRS.sumOverVectorOP[i][0][3] / (1.0*MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m2_1 = firstRS.sumOverVectorOP[i][0][1] / (1.0*MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m4_2 = secRS.sumOverVectorOP[i][0][3] / (1.0*MD.numSweepsPerformedToCheckIfThermalized);
            double m2_2 = secRS.sumOverVectorOP[i][0][1] / (1.0*MD.numSweepsPerformedToCheckIfThermalized);
            
            firstBinder = 1.0 - m4_1/ (3.0 * m2_1 * m2_1);
            
            secBinder = 1.0 - m4_2/ (3.0 * m2_2 * m2_2);
            
            cout<<"Is Thermal? Diff is: "<<fabs(firstBinder - secBinder)<<endl;
            if(fabs(firstBinder - secBinder) > 1e2){
                smallDiffBinder = false;
                break;
            }
            
        }
        
        if(smallDiffBinder){// thermalized.
            setRSToZero();
            setOPToZero();
            setRecipToZero();
            
            double range_t = MD.range;
            double rangeOld_t =  MD.rangeOld;
            int numSweepsPerformed_t =  MD.numSweepsPerformed;
            int numSweepsUsedToThermalize_t =  MD.numSweepsUsedToThermalize;
            clock_t timeOfInitialization_t =  MD.timeOfInitialization;
            
            setMDToZero();
            
            MD.range = range_t;
            MD.rangeOld = rangeOld_t;
            MD.numSweepsPerformed = numSweepsPerformed_t;
            MD.numSweepsUsedToThermalize = numSweepsUsedToThermalize_t;
            MD.timeOfInitialization = timeOfInitialization_t;
            MD.isThermal = true;
            
            cout<<"Returning from isThermalized. Not Out of Time."<<endl;
            return outOfTime;
        }
    }
    
    
    return outOfTime;
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
                    Q += PX - PY;
                    
                    //S[A-1,B].S[A+1,B]-S[A,B+1].S[A,B-1]
                    Dimer += PXmX - PYmY;
                    
                }
            }
        }
    }
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
                                    
                                    kvector_x = (2*k_a + (k_s%2)) * B_VECTOR_1[0] / (MCP.cellsA*1.0) + (2*k_b + ((k_s-1)>0)) * B_VECTOR_2[0] / (MCP.cellsB*1.0);
                                    kvector_y = (2*k_a + (k_s%2)) * B_VECTOR_1[1] / (MCP.cellsA*1.0) + (2*k_b + ((k_s-1)>0)) * B_VECTOR_2[1] / (MCP.cellsB*1.0);
                                    rvector_x = (2*a + (s%2)) * A_VECTOR_1[0] + (2*b + ((s-1)>0)) * A_VECTOR_2[0];
                                    rvector_y = (2*a + (s%2)) * A_VECTOR_1[1] + (2*b + ((s-1)>0)) * A_VECTOR_2[1];
                                    kr = kvector_x * rvector_x + kvector_y * rvector_y;
                                    cos_kr = cos(kr);
                                    sin_kr = sin(kr);
                                    recipLattice[k_c][k_a][k_b][k_s].ReXComponent += cos_kr * lattice[c][a][b][s].x;
                                    recipLattice[k_c][k_a][k_b][k_s].ImXComponent += sin_kr * lattice[c][a][b][s].x;
                                    recipLattice[k_c][k_a][k_b][k_s].ReYComponent += cos_kr * lattice[c][a][b][s].y;
                                    recipLattice[k_c][k_a][k_b][k_s].ImYComponent += sin_kr * lattice[c][a][b][s].y;
                                    recipLattice[k_c][k_a][k_b][k_s].ReZComponent += cos_kr * lattice[c][a][b][s].z;
                                    recipLattice[k_c][k_a][k_b][k_s].ImZComponent += sin_kr * lattice[c][a][b][s].z;
                                    
                                }
                            }
                        }
                    }
                    
                    
                    recipLattice[k_c][k_a][k_b][k_s].magnitude = (
                                                                  recipLattice[k_c][k_a][k_b][k_s].ReXComponent * recipLattice[k_c][k_a][k_b][k_s].ReXComponent +
                                                                  recipLattice[k_c][k_a][k_b][k_s].ReYComponent * recipLattice[k_c][k_a][k_b][k_s].ReYComponent +
                                                                  recipLattice[k_c][k_a][k_b][k_s].ReZComponent * recipLattice[k_c][k_a][k_b][k_s].ReZComponent +
                                                                  recipLattice[k_c][k_a][k_b][k_s].ImXComponent * recipLattice[k_c][k_a][k_b][k_s].ImXComponent +
                                                                  recipLattice[k_c][k_a][k_b][k_s].ImYComponent * recipLattice[k_c][k_a][k_b][k_s].ImYComponent +
                                                                  recipLattice[k_c][k_a][k_b][k_s].ImZComponent * recipLattice[k_c][k_a][k_b][k_s].ImZComponent
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
                    
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReXComponent += recipLattice[k_c][k_a][k_b][k_s].ReXComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImXComponent += recipLattice[k_c][k_a][k_b][k_s].ImXComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReYComponent += recipLattice[k_c][k_a][k_b][k_s].ReYComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImYComponent += recipLattice[k_c][k_a][k_b][k_s].ImYComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ReZComponent += recipLattice[k_c][k_a][k_b][k_s].ReZComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].ImZComponent += recipLattice[k_c][k_a][k_b][k_s].ImZComponent;
                    sumOverRecipLattice[k_c][k_a][k_b][k_s].magnitude    += recipLattice[k_c][k_a][k_b][k_s].magnitude;
                    
                }
            }
        }
    }
}

void MonteCarlo::updateRunningSumOP(){
    
    updateOrderParameters();
    double temp;
    
    for(int i = 0; i < NUMSCALAROP; i++){
        
        temp = OP.scalarOP[i]*OP.scalarOP[i];
        
        RS.sumOverScalarOP[i][0] += OP.scalarOP[i];
        RS.sumOverScalarOP[i][1] += temp;
        RS.sumOverScalarOP[i][2] += OP.scalarOP[i]*temp;
        RS.sumOverScalarOP[i][3] += temp*temp;
    }
    
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            
            temp = OP.vectorOP[i][j]*OP.vectorOP[i][j];
            
            RS.sumOverVectorOP[i][j][0] += OP.vectorOP[i][j];
            RS.sumOverVectorOP[i][j][1] += temp;
            RS.sumOverVectorOP[i][j][2] += OP.vectorOP[i][j]*temp;
            RS.sumOverVectorOP[i][j][3] += temp*temp;
        }
    }
    
    
    MD.numConfigsDone++;
    if(MD.numConfigsDone % (10 * 1000) == 0){
        adjustRange();
    }
}


/******************************************************************************/
//Spin Functions
/******************************************************************************/

void MonteCarlo::checkSpin(uint_fast8_t c, uint_fast8_t a, uint_fast8_t b, uint_fast8_t s){
    double size = sqrt(lattice[c][a][b][s].x*
                       lattice[c][a][b][s].x +
                       lattice[c][a][b][s].y*
                       lattice[c][a][b][s].y +
                       lattice[c][a][b][s].z*
                       lattice[c][a][b][s].z);
    
    lattice[c][a][b][s].x = lattice[c][a][b][s].x / (1.0*size);
    lattice[c][a][b][s].y = lattice[c][a][b][s].y / (1.0*size);
    lattice[c][a][b][s].z = lattice[c][a][b][s].z / (1.0*size);
    
    if((lattice[c][a][b][s].x > 1 || lattice[c][a][b][s].x < -1)||
       (lattice[c][a][b][s].y > 1 || lattice[c][a][b][s].y < -1)||
       (lattice[c][a][b][s].z > 1 || lattice[c][a][b][s].z < -1)||
       (lattice[c][a][b][s].x != lattice[c][a][b][s].x)||
       (lattice[c][a][b][s].y != lattice[c][a][b][s].y)||
       (lattice[c][a][b][s].z != lattice[c][a][b][s].z)||
       ((fabs(sqrt(lattice[c][a][b][s].x *
                   lattice[c][a][b][s].x +
                   lattice[c][a][b][s].y *
                   lattice[c][a][b][s].y +
                   lattice[c][a][b][s].z *
                   lattice[c][a][b][s].z) - 1)) >= 0.000001) ){
        std::cerr << "Error: Bad Spin values: " << "x: " << lattice[c][a][b][s].x << std::endl;
        std::cerr << "Error: Bad Spin values: " << "y: " << lattice[c][a][b][s].y << std::endl;
        std::cerr << "Error: Bad Spin values: " << "z: " << lattice[c][a][b][s].z << std::endl;
        std::cerr << "Error: Bad Spin values: " << "Magnitude: " <<
        std::setprecision(15) << std::setw(15) <<
        sqrt(lattice[c][a][b][s].x *
             lattice[c][a][b][s].x +
             lattice[c][a][b][s].y *
             lattice[c][a][b][s].y +
             lattice[c][a][b][s].z *
             lattice[c][a][b][s].z) << std::endl;
        exit(1);
    }
}

void MonteCarlo::setRandomOrientation(uint_fast8_t c, uint_fast8_t a,
                                      uint_fast8_t b, uint_fast8_t s){
    lattice[c][a][b][s].x = 2 * drand48() - 1;
    lattice[c][a][b][s].y = 2 * drand48() - 1;
    lattice[c][a][b][s].z = 2 * drand48() - 1;
    double size = sqrt(lattice[c][a][b][s].x *
                       lattice[c][a][b][s].x +
                       lattice[c][a][b][s].y *
                       lattice[c][a][b][s].y +
                       lattice[c][a][b][s].z *
                       lattice[c][a][b][s].z);
    lattice[c][a][b][s].x = lattice[c][a][b][s].x / (1.0*size);
    lattice[c][a][b][s].y = lattice[c][a][b][s].y / (1.0*size);
    lattice[c][a][b][s].z = lattice[c][a][b][s].z / (1.0*size);
    
    flip(c, a, b, s, PI-0.01);
}


void MonteCarlo::specifyRotation(uint_fast8_t c, uint_fast8_t a,
                                 uint_fast8_t b, uint_fast8_t s,
                                 double rotPhi, double rotTheta){
    
    //Converts to spherical (in order to construct a perpendicular vector).
    //measured from the vertical.
    double theta = acos(lattice[c][a][b][s].z);
    double phi =  atan2(lattice[c][a][b][s].y,
                        lattice[c][a][b][s].x);
    
    theta = theta + rotTheta;
    phi = phi + rotPhi;
    
    if((theta > PI) && (theta < 2 * PI)){
        theta = PI - (theta - PI);
        phi = phi + PI;//new
    }else if (theta >= 2 * PI){
        theta = theta - 2 * PI;
    }
    
    if(phi > 2 * PI){
        phi = phi - 2 * PI;
    }
    if(phi > 2 * PI){
        phi = phi - 2 * PI;
    }
    
    lattice[c][a][b][s].x = sin(theta) * cos(phi);
    lattice[c][a][b][s].y = sin(theta) * sin(phi);
    lattice[c][a][b][s].z = cos(theta);
    
    checkSpin(c, a, b, s);
    
}

/**Does rotation using quaternions. Requires that the spin has cartesian
 *components but that the axis that it rotates around is in spherical.
 *
 *Note: this function should be replaced with a function that only uses
 *spherical coordinates or only cartesian coordinates.
 */
void MonteCarlo::rotate(uint_fast8_t c, uint_fast8_t a,
                        uint_fast8_t b, uint_fast8_t s,
                        double theta_, double phi_, double Beta) {
    double M[3][3];
    double q0, q1, q2, q3;
    double newS[3] = {0};
    q0 = cos(Beta / 2.0);
    q1 = sin(Beta / 2.0) * sin(theta_) * cos(phi_);
    q2 = sin(Beta / 2.0) * sin(theta_) * sin(phi_);
    q3 = sin(Beta / 2.0) * cos(theta_);
    M[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    M[0][1] = 2.0 * (q1 * q2 - q0 * q3);
    M[0][2] = 2.0 * (q1 * q3 + q0 * q2);
    M[1][0] = 2.0 * (q2 * q1 + q0 * q3);
    M[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    M[1][2] = 2.0 * (q2 * q3 - q0 * q1);
    M[2][0] = 2.0 * (q3 * q1 - q0 * q2);
    M[2][1] = 2.0 * (q3 * q2 + q0 * q1);
    M[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    
    /*
     for(int i = 0; i < 3; i ++){
     for(int j = 0; j < 3; j++){
     newS[i] += M[i][j] * components[j];
     }
     }
     */
    
    for(int i = 0; i < 3; i ++){
        for(int j = 0; j < 3; j++){
            if(i == 0 && j == 0){
                newS[i] += M[i][j] * lattice[c][a][b][s].x;
            }else if(i == 0 && j == 1){
                newS[i] += M[i][j] * lattice[c][a][b][s].y;
            }else if(i == 0 && j == 2){
                newS[i] += M[i][j] * lattice[c][a][b][s].z;
            }else if(i == 1 && j == 0){
                newS[i] += M[i][j] * lattice[c][a][b][s].x;
            }else if(i == 1 && j == 1){
                newS[i] += M[i][j] * lattice[c][a][b][s].y;
            }else if(i == 1 && j == 2){
                newS[i] += M[i][j] * lattice[c][a][b][s].z;
            }else if(i == 2 && j == 0){
                newS[i] += M[i][j] * lattice[c][a][b][s].x;
            }else if(i == 2 && j == 1){
                newS[i] += M[i][j] * lattice[c][a][b][s].y;
            }else if(i == 2 && j == 2){
                newS[i] += M[i][j] * lattice[c][a][b][s].z;
            }else{
                std::cerr << "Bad Rotation" << std::endl;
                exit(1);
            }
        }
    }
    
    //In Cartesian
    lattice[c][a][b][s].x = newS[0];
    lattice[c][a][b][s].y = newS[1];
    lattice[c][a][b][s].z = newS[2];
}

void MonteCarlo::flip(uint_fast8_t c, uint_fast8_t a,
                      uint_fast8_t b, uint_fast8_t s, double range){
    
    
    double u[3];
    //Cross product with a sufficiently far away "g" to find a perp vector.
    //u[0] = g2s3 - g3s2;
    //u[1] = g3s1 - g1s3;
    //u[2] = g1s2 - g2s1;
    if(lattice[c][a][b][s].x > 0.57) {//g = (0,1,0)
        u[0] = lattice[c][a][b][s].z;
        u[1] = 0;
        u[2] = (-1) * lattice[c][a][b][s].x;
    }else if(lattice[c][a][b][s].y > 0.57){//g = (0,0,1)
        u[0] = (-1) * lattice[c][a][b][s].y;
        u[1] = lattice[c][a][b][s].x;
        u[2] = 0;
    }else{//g = (1,0,0)
        u[0] = 0;
        u[1] = (-1) * lattice[c][a][b][s].z;
        u[2] = lattice[c][a][b][s].y;
    }
    
    //Cross product to find third Orthogonal vector
    double v[3] = { u[1]*lattice[c][a][b][s].z - u[2]*lattice[c][a][b][s].y,
        u[2]*lattice[c][a][b][s].x - u[0]*lattice[c][a][b][s].z,
        u[0]*lattice[c][a][b][s].y - u[1]*lattice[c][a][b][s].x};
    
    //Normalize vectors
    double size = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    u[0] = u[0] / (1.0*size);
    u[1] = u[1] / (1.0*size);
    u[2] = u[2] / (1.0*size);
    
    size = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] = v[0] / (1.0*size);
    v[1] = v[1] / (1.0*size);
    v[2] = v[2] / (1.0*size);
    
    //spin, u, v are orthonormal.
    //See these two sites:
    //http://www.math.niu.edu/~rusin/known-math/96/sph.rand
    //http://objectmix.com/graphics/314285-random-points-spherical-cap.html
    //(copied below)
    
    //Size of circle:
    if(range >= PI){range = PI;}
    double distFromOrig = cos(range);//works both for pos dist and neg dist.
    double circDist = drand48() * (1 - distFromOrig) + distFromOrig;
    double circleRadius = sqrt(1 - circDist*circDist);
    
    double r = drand48() * 2.0 * PI;
    double cosr = cos(r);
    double sinr = sin(r);
    
    lattice[c][a][b][s].x = lattice[c][a][b][s].x * circDist;
    lattice[c][a][b][s].y = lattice[c][a][b][s].y * circDist;
    lattice[c][a][b][s].z = lattice[c][a][b][s].z * circDist;
    
    lattice[c][a][b][s].x = lattice[c][a][b][s].x + circleRadius * (cosr * u[0] + sinr * v[0]);
    lattice[c][a][b][s].y = lattice[c][a][b][s].y + circleRadius * (cosr * u[1] + sinr * v[1]);
    lattice[c][a][b][s].z = lattice[c][a][b][s].z + circleRadius * (cosr * u[2] + sinr * v[2]);
    
}

/******************************************************************************/
//Core Functions
/******************************************************************************/



double MonteCarlo::getLocalEnergy(uint_fast8_t c, uint_fast8_t a,
                                  uint_fast8_t b, uint_fast8_t s,
                                  const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
    //The order of neighbors are:
    //[0]: first NN in x direction
    //[1]: first NN in -x direction
    //[2]: first NN in y direction
    //[3]: first NN in -y direction
    //[4]: second NN, going first in the x direction, then the y direction
    //[5]: second NN, going first in the -x direction, then the y direction
    //[6]: second NN, going first in the -x direction, then the -y direction
    //[7]: second NN, going first in the x direction, then the -y direction
    //[8]: third NN, going first in the x direction, then the x direction
    //[9]: third NN, going first in the -x direction, then the -x direction
    //[10]: third NN, going first in the y direction, then the y direction
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
    
    double interactionCubic = MCP.cubicD * (-1.0) * lattice[c][a][b][s].z * lattice[c][a][b][s].z;
    
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
    double interactionH2 = MCP.j2 * (dotProd(c, a, b, s, getSpinCoordsInTheYDir(getSpinCoordsInTheXDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheYDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusYDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheXDir(getSpinCoordsInTheMinusYDir(sc))));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(c, a, b, s, getSpinCoordsInTheXDir(getSpinCoordsInTheXDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusXDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheYDir(getSpinCoordsInTheYDir(sc))) +
                                     dotProd(c, a, b, s, getSpinCoordsInTheMinusYDir(getSpinCoordsInTheMinusYDir(sc))));
    
    double interactionK = MCP.k1 * temp * temp;
    
    double interactionCubic = MCP.cubicD * (-1.0) * lattice[c][a][b][s].z * lattice[c][a][b][s].z;
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic) / 2.0;
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


void MonteCarlo::metropolisSweep(){
    double   changeInEnergy  = 0;
    double   initialEnergy   = 0;
    double   finalEnergy     = 0;
    double   xOld = 0;
    double   yOld = 0;
    double   zOld = 0;
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    //getSpin(i, j, k).checkSpin();
                    initialEnergy = getLocalEnergy(c, a, b, s, true);
                    //initialEnergy = getLocalEnergyFromCalc(c, a, b, s, true);
                    xOld = lattice[c][a][b][s].x;
                    yOld = lattice[c][a][b][s].y;
                    zOld = lattice[c][a][b][s].z;
                    flip(c, a, b, s, MD.range);
                    finalEnergy = getLocalEnergy(c, a, b, s, true);
                    //finalEnergy = getLocalEnergyFromCalc(c, a, b, s, true);
                    changeInEnergy = finalEnergy - initialEnergy;
                    MD.flipAttempts += 1;
                    MD.successfulFlips += 1;
                    //The lattice automatically keeps the spin flip unless the
                    //energy is > 0
                    
                    if(changeInEnergy > 0){
                        /**
                         *Then it keeps the flip if a function f(KbT, delta_E) is
                         *less than a certain probability (0,1).
                         */
                        if(drand48() >= exp((-1) * (1.0/(1.0*MCP.KbT)) * changeInEnergy)){
                            
                            //resets spin
                            lattice[c][a][b][s].x = xOld;
                            lattice[c][a][b][s].y = yOld;
                            lattice[c][a][b][s].z = zOld;
                            MD.successfulFlips -= 1;
                        }else{
                            //Keeps Spin flip if rand < exp()
                        }
                    }else if(changeInEnergy == 0 && initialEnergy == 0){
                        //resets spin
                        lattice[c][a][b][s].x = xOld;
                        lattice[c][a][b][s].y = yOld;
                        lattice[c][a][b][s].z = zOld;
                        MD.successfulFlips -= 1;
                    }
                }
            }
        }
    }
}




void MonteCarlo::sweep(int numIndep_, double durationOfSingleRun){
    
    for(int sweeps = MD.numSweepsPerformed; sweeps <= MCP.numSweepsToPerformTotal; sweeps++){
        
        
        metropolisSweep();
        MD.numSweepsPerformed += 1;
        if(sweeps % numIndep_ == 0){
            updateRunningSumOP();
            
            
            double clocks = clock() - MD.timeOfInitialization;
            double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60.0 * 60.0);
            if(hoursElapsed > durationOfSingleRun){
                cout<<"Ran out of time, in sweep()"<<endl;
                break;
            }
        }
        
        if((MD.numConfigsDone % UPDATESTOCORR == 0)&&(MD.numConfigsDone > 10)){
            updateCorrelationFunction();
        }
        
        if((MD.numConfigsDone % UPDATESTOFT == 0)&&(MD.numConfigsDone > 10)){
            updateFourierTransformOnRecipLattice();
        }
        
        if (MD.totalTimeRunning >= OVERFLOWLIMIT){
            MD.totalTimeRunning -= OVERFLOWLIMIT;
            MD.numClockOverflows++;
        }
        
    }
    
    cout<<"Finished all sweeps. Printing. In function: sweep()"<<endl;
    
}

void MonteCarlo::sweepSnapshot(int numSweeps_, int numIndep_, double hours_){
    /*
     clock_t start = clock();
     while((clock() - start) * (1/(CLOCKS_PER_SEC*1.0)) < (hours_*60*60) )
     {
     double KbT_ = KbT;
     thermalize(KbT_ * 64);
     thermalize(KbT_ * 16);
     thermalize(KbT_ * 4);
     thermalize(KbT_ * 2);
     thermalize(KbT_);
     int swept = 0;
     for(int sweeps = 1; sweeps <= numSweeps_; sweeps++)
     {
     metropolisSweep();
     if(sweeps % numIndep_ == 0)
     {
     printSnapshot();
     }
     
     if((clock() - start) * (1/(CLOCKS_PER_SEC*1.0)) > (hours_*60*60) ){
     swept = sweeps;
     break;
     }
     
     }
     
     //
     std::stringstream stream;
     std::ofstream fileOut;
     
     stringstream out;
     out.str("");
     string temp = outputFileName.substr(0,outputFileName.size()-4);
     out << temp.c_str() << "OPSnap.txt";
     
     fileOut.open(out.str().c_str(), std::ios::out | std::ios::app);
     stream.str("");
     
     stream << "Done with one run of sweepSnapshot(), consisting of " << swept;
     stream << " runs, skipping every " << numIndep_ << " sweep." <<endl;
     
     fileOut << stream.str() << std::endl;
     fileOut.close();
     //
     
     }
     
     cout<<"finished finalPrint(), in sweepSnapshot()"<<endl;
     */
}


/*
 From: ags@seaman.cc.purdue.edu (Dave Seaman)
 Newsgroups: sci.math.num-analysis
 Subject: Re: N-dim spherical random number drawing
 Date: 20 Sep 1996 13:04:58 -0500
 
 In article <wmEgVDG00VUtMEhUwg@andrew.cmu.edu>,
 Shing-Te Li  <sl2x+@andrew.cmu.edu> wrote:
 >Does anyone know a good algorithm that can draw a vector of random
 >numbers uniformly from a N dimensional sphere?
 
 Here are four methods for generating uniformly-distributed points on a unit
 sphere:
 
 (1) The normal-deviate method.  Choose x, y, and z, each
 normally-distributed with mean 0 and variance 1.  (Actually, the
 variance does not matter, provided it remains fixed.)  Normalize
 the result to obtain a uniform density on the unit sphere.
 
 This method also generalizes nicely to n dimensions.  It works
 because the vector chosen (before normalization) has a density
 that depends only on the distance from the origin.  The method is
 mentioned in Knuth, _Seminumerical Algorithms_, 2nd. ed., pp.
 130-131.
 
 (2) The hypercube rejection method.  Choose x, y, and z, uniformly
 distributed on [-1,1].  Compute s = x^2 + y^2 + z^2.  If s > 1,
 reject the triplet and start over.  Once you have an acceptable
 vector, normalize it to obtain a point on the unit sphere.
 
 This method also generalizes to n dimensions, but as the dimension
 increases the probability of rejection becomes high and many random
 numbers are wasted.
 
 (3) The trig method.  This method works only in 3-space, but it is
 very fast.  It depends on the slightly counterintuitive fact (see
 proof below) that each of the three coordinates of a uniformly
 distributed point on S^2 is uniformly distributed on [-1,1] (but
 the three are not independent, obviously).  Therefore, it
 suffices to choose one axis (Z, say) and generate a uniformly
 distributed value on that axis.  This constrains the chosen point
 to lie on a circle parallel to the X-Y plane, and the obvious
 trig method may be used to obtain the remaining coordinates.
 
 (a) Choose z uniformly distributed in [-1,1].
 (b) Choose t uniformly distributed on [0, 2*pi).
 (c) Let r = sqrt(1-z^2).
 (d) Let x = r * cos(t).
 (e) Let y = r * sin(t).
 
 (4) A variation of method (3) that doesn't use trig may be found in
 Knuth, _loc. cit._.  The method is equivalent to the following:
 
 (a) Choose u and v uniformly distributed on [-1,1].
 (b) Let s = u^2 + v^2.
 (c) If s > 1, return to step (a).
 (d) Let a = 2 * sqrt(1-s).
 (e) The desired point is (a*u, a*v, 2*s-1).
 
 This method uses two-dimensional rejection, which gives a higher
 probability of acceptance compared to algorithm (2) in three
 dimensions.  I group this with the trig method because both
 depend on the fact that the projection on any axis is uniformly
 distributed on [-1,1], and therefore both methods are limited to
 use on S^2.
 
 I have found the trig method to be fastest in tests on an IBM RS/6000
 model 590, using the XLF 3.2 Fortran compiler, with the IMSL
 subroutines DRNUN and DRNNOA generating the random numbers.  The trig
 routines were accelerated by using the MASS library, and sqrt
 computations were performed by the hardware sqrt instruction available
 on the 590.  Timings for 1 million random points on S^2:
 
 (1) normal-deviate method                   9.60 sec.
 (2) 3-D rejection                           8.84 sec.
 (3) trig method                             2.71 sec.
 (4) polar with 2-D rejection                5.08 sec.
 
 The fact that algorithm (3) does not involve rejection turns out to be
 a huge advantage because it becomes possible to generate all the points
 in a single loop with no conditional statements, which optimizes very
 effectively and more than offsets the cost of using trig functions.
 The alternative with algorithms (2) and (4) is to generate the points
 one at a time and test each one for acceptance before proceeding to the
 next, which requires nested loops and many more subroutine calls to
 generate the random numbers.  Algorithm (1) does not explicitly use
 rejection, but the IMSL subroutine DRNNOA uses rejection in generating
 the normally-distributed numbers that are required.
 
 In higher dimensions, the normal-deviate method is preferred
 because the probability of rejection in the hypercube method
 grows very rapidly with the dimension.
 
 --------------------------------------------------------------------
 
 Here is a proof of correctness for the fact that the z-coordinate is
 uniformly distributed.  The proof also works for the x- and y-
 coordinates, but note that this works only in 3-space.
 
 The area of a surface of revolution in 3-space is given by
 
 A = 2 * pi * int_a^b f(x) * sqrt(1 + [f'(x}]^2) dx
 
 Consider a zone of a sphere of radius R.  Since we are integrating in
 the z direction, we have
 
 f(z) = sqrt(R^2 - z^2)
 f'(z) = -z / sqrt(R^2-z^2)
 
 1 + [f'(z)]^2 = r^2 / (R^2-z^2)
 
 A = 2 * pi * int_a^b sqrt(R^2-z^2) * R/sqrt(R^2-z^2) dz
 
 = 2 * pi * R int_a^b dz
 
 = 2 * pi * R * (b-a)
 
 = 2 * pi * R * h,
 
 where h = b-a is the vertical height of the zone.  Notice how the integrand
 reduces to a constant.  The density is therefore uniform.
 ==============================================================================
 */

/*
 Hi,
 
 if we have two sphere, S1 and S2, their intersection will be a circle
 C (suppose that the intersection of the two spheres exist and it is
 not a point). On S1 and S2 we will have a two spherical caps. C will
 be the border of the spherical cap.
 
 There is an efficient way to generate random points on one of the
 spherical cap ?
 
 One way that should work well, if the radius of the circle is small,
 is to generate a point p on the corresponding disk and to compute the
 intersection of the ray o-p with the sphere S1 (for example), where o
 is the center of S1.
 
 But if the radius of the circle is not small, the disk will not be a
 good approximation of the spherical cap, so the points will not be
 uniformly distributed on the cap...
 
 Thank you
 Luca
 Reply With Quote Reply With Quote
 12-11-2007 09:39 PM #2
 Default Re: Random points on a spherical cap
 <kingpin@freemail.it> wrote in message
 news:96bb8ce2-0d44-4b94-b5ec-9ee829e485e7@p69g2000hsa.googlegroups.com...
 
 To generate random points on the spherical cap for S1, let K1 be
 the sphere center and let R1 be the sphere radius. Let A be the
 circle center and let S be the circle radius. Compute
 W = (A-K1)/Length(A-K1). Choose U and V to be unit-length
 vectors so that {U,V,W} is an orthonormal set (mutually
 perpendicular, all unit-length).
 
 Select a random angle t in [0,2*pi) and a random length L in [0,S].
 The point P = A + L*(cos(t)*U + sin(t)*V) lies inside the circle.
 Define D = P - K1. The corresponding point on the spherical cap is
 Q = R1*D/Length(D).
 
 --
 Dave Eberly
 http://www.geometrictools.com
 */




