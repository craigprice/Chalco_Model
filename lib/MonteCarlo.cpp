
/**
 MonteCarlo.cpp
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of magnetic models in a regular geometry at
 finite temperature.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/07/10
 */

#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include "MonteCarlo.h"
#include "FileIO.cpp"
#include "Chalco_Square.cpp"
//#include "KH_Honeycomb.cpp"
using namespace std;


/******************************************************************************/
//Functions for initialization
/******************************************************************************/



void MonteCarlo::setRSToZero(){
    cout<<"Begin setRSToZero()"<<endl;
    
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
    cout<<"Begin setOPToZero()"<<endl;
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



void MonteCarlo::setMDToZero(){
    cout<<"Begin setMDToZero()"<<endl;
    
    MD.flipAttempts    = 0;
    MD.successfulFlips = 0;
    MD.range           = PI-0.1;
    MD.rangeOld        = 0.1;
    MD.numConfigsDone  = 0;
    MD.isThermal        = false;
    MD.numSweepsPerformedToCheckIfThermalized = 0;
    MD.numSweepsPerformed           = 0;
    MD.numSweepsUsedToThermalize    = 0;
    
}

void MonteCarlo::setRecipToZero(){
    cout<<"Begin setRecipToZero()"<<endl;
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
    cout<<"Begin initializeLattice()"<<endl;
    
    ClassicalSpin3D spin;
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<ClassicalSpin3D> > > row1;
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<ClassicalSpin3D> > row2;
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                std::vector<ClassicalSpin3D> row3;
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
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
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector< std::vector<SpinCoordinates> > > > row1;
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector< std::vector<SpinCoordinates> > > row2;
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                std::vector< std::vector<SpinCoordinates> > row3;
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    std::vector<SpinCoordinates> row4;
                    for(uint_fast8_t n = 0; n < NUMNEIGHBORS; n++){
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
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row1;
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<FourierTransformOfSpin> > row2;
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                std::vector<FourierTransformOfSpin> row3;
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    row3.push_back(spinFT);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        recipLattice.push_back(row1);
    }
    
    cout << "sumOverRecipLattice" << endl;
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row1;
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            std::vector< std::vector<FourierTransformOfSpin> > row2;
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                std::vector<FourierTransformOfSpin> row3;
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    row3.push_back(spinFT);
                }
                row2.push_back(row3);
            }
            row1.push_back(row2);
        }
        sumOverRecipLattice.push_back(row1);
    }
    
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



void MonteCarlo::init(){
    cout<<"Begin init()"<<endl;
    
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
    cout<<"Begin MonteCarlo()"<<endl;
    init();
    
    initializeLattice();
    numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES;
    
    
    cout<<"Finished with default MC constructor"<<endl;
    MD.timeAtThisJobInitialization = time(NULL);
    
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
    }
    
    MD.timeAtThisJobInitialization = time(NULL);
}


/******************************************************************************/
//Thermalization
/******************************************************************************/

void MonteCarlo::adjustRange(){
    //cout<<"Begin adjustRange()"<<endl;
    
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


void MonteCarlo::thermalize(){
    cout<<"Begin thermalize()"<<endl;
    
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
    //Ramp temperature down.
    if(MD.numSweepsUsedToThermalize < (MINTHERMALSWEEPS / 2.0) ){
        for(int sweeps = MD.numSweepsUsedToThermalize;
            sweeps <= (MINTHERMALSWEEPS / 2.0);
            sweeps++){
            
            MCP.KbT = 1000.0 *
            (MINTHERMALSWEEPS / 2.0 - MD.numSweepsUsedToThermalize) /
            (MINTHERMALSWEEPS / 2.0) + oldKbT;
            
            if(isOutOfTime()){
                break;
            }
            
            metropolisSweep();
            MD.numSweepsUsedToThermalize += 1;
            adjustRange();
            
        }
    }
    
    MCP.KbT = oldKbT;
    
    //spend time thermalizing at the given kbt
    if(MD.numSweepsUsedToThermalize >= (MINTHERMALSWEEPS / 2.0) ){
        for(int sweeps = MD.numSweepsUsedToThermalize;
            sweeps <= MINTHERMALSWEEPS;
            sweeps++){
            
            if(isOutOfTime()){
                break;
            }
            
            metropolisSweep();
            MD.numSweepsUsedToThermalize += 1;
            adjustRange();
            
        }
    }
    
    while((!isOutOfTime()) && !MD.isThermal){
        checkThermalized();
    }
    
    if(MD.numSweepsUsedToThermalize >= MINTHERMALSWEEPS){
        cout<<"Returning from Thermalize. Finished."<<endl;
    }else{
        cout<<"Returning from Thermalize. Out of Time."<<endl;
    }
    
}

void MonteCarlo::checkThermalized(){
    cout<<"Begin checkThermalized()"<<endl;
    
    for(int i = 0; i < 10; i++){
        metropolisSweep();
        MD.numSweepsUsedToThermalize += 1;
        adjustRange();
        updateRunningSumOP();//This function is the only place this is called except sweep()
        MD.numSweepsPerformedToCheckIfThermalized += 1;
    }
    
    int thermalCycles = -1;
    while(! isOutOfTime() ){
        thermalCycles ++;
        cout<<"Thermal Cycle # " << thermalCycles<<endl;
        
        RunningSumOverOP firstRS = RS;
        
        int numSweepsAveOver = (int) (2000);
        for(int i = 0; i < numSweepsAveOver; i++){
            metropolisSweep();
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
            double m4_1 = firstRS.sumOverVectorOP[i][0][3] /
            (1.0*MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m2_1 = firstRS.sumOverVectorOP[i][0][1] /
            (1.0*MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m4_2 = secRS.sumOverVectorOP[i][0][3] /
            (1.0*MD.numSweepsPerformedToCheckIfThermalized);
            double m2_2 = secRS.sumOverVectorOP[i][0][1] /
            (1.0*MD.numSweepsPerformedToCheckIfThermalized);
            
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
            MD.isThermal = true;
            cout<<"Returning from isThermalized. Not Out of Time."<<endl;
            break;
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

void MonteCarlo::checkSpin(uint_fast8_t c, uint_fast8_t a,
                           uint_fast8_t b, uint_fast8_t s){
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
        std::cerr << "Error: Bad Spin values: " << "x: " <<
        lattice[c][a][b][s].x << std::endl;
        std::cerr << "Error: Bad Spin values: " << "y: " <<
        lattice[c][a][b][s].y << std::endl;
        std::cerr << "Error: Bad Spin values: " << "z: " <<
        lattice[c][a][b][s].z << std::endl;
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
    double size = sqrt(lattice[c][a][b][s].x * lattice[c][a][b][s].x +
                       lattice[c][a][b][s].y * lattice[c][a][b][s].y +
                       lattice[c][a][b][s].z * lattice[c][a][b][s].z);
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
 *SPEEDUP: this function should be replaced with a function that only uses
 *spherical coordinates or only cartesian coordinates.
 *
 *This function is not being used. Kept for reference for quaterion rotations.
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
    double v[3] =
    {u[1]*lattice[c][a][b][s].z - u[2]*lattice[c][a][b][s].y,
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
    
    double tempx = lattice[c][a][b][s].x * circDist;
    double tempy = lattice[c][a][b][s].y * circDist;
    double tempz = lattice[c][a][b][s].z * circDist;
    
    lattice[c][a][b][s].x = tempx + circleRadius * (cosr * u[0] + sinr * v[0]);
    lattice[c][a][b][s].y = tempy + circleRadius * (cosr * u[1] + sinr * v[1]);
    lattice[c][a][b][s].z = tempz + circleRadius * (cosr * u[2] + sinr * v[2]);
    
}

/******************************************************************************/
//Core Functions
/******************************************************************************/


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
    MD.numSweepsPerformed += 1;
}


bool MonteCarlo::isOutOfTime(){
    time_t currentElapsedTime = difftime(time(NULL),
                                         MD.timeAtThisJobInitialization);
    double hoursElapsed = currentElapsedTime / (60.0 * 60.0);
    if(hoursElapsed > LIMITONHOURSOFJOBRUN){
        cout << "Is out of Time " << endl;
        return true;
    }
    return false;
}


void MonteCarlo::sweep(){
    cout<<"Begin sweep()"<<endl;
    
    for(int sweeps = MD.numSweepsPerformed;
        sweeps <= MCP.numSweepsToPerformTotal;
        sweeps++){
        
        
        metropolisSweep();
        if(sweeps % NUMSWEEPSBETWEENCONFIGS == 0){
            updateRunningSumOP();
            
            if(isOutOfTime()){
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
        
    }
    
    cout<<"Finished sweep()"<<endl;
    
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




