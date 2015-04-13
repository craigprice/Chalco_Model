/**
 FileIO.cpp
 Purpose: A collection of functions for MonteCarlo.cpp that manage writing
 the simulation out to a text file and reading in a previously written text file
 to continue where it left off.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/04/12
 */
#include <ctime>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "MonteCarlo.h"

using namespace std;

void MonteCarlo::setParamNames(){
    string function = "MonteCarlo::setParamNames()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream name;
    for(int i = 0; i < NUMVECTOROP; i ++){
        for(int j = 0; j < TYPESOP; j++){
            name.str("");
            name << vectorOrderParameterName[i] << "'s "<< typeOP[j];
            vectorOPAndType[i][j] = name.str();
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i ++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                name.str("");
                name << "Sum over Configurations, " << vectorOrderParameterName[i] << "'s "<< typeOP[j] << powerOP[k] << ": ";
                sumOverVectorOPNames[i][j][k] = name.str();
            }
        }
    }
    
    for(int i = 0; i < NUMSCALAROP; i ++){
        for(int j = 0; j < POWERSOP; j++){
            name.str("");
            name << "Sum over Configurations, " << scalarOrderParameterName[i] << powerOP[j] << ": ";
            sumOverScalarOPNames[i][j] = name.str();
        }
    }
    
    cout<<"End setParamNames"<<endl;
}



double MonteCarlo::findParameterDbl(double currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atof(l.substr(str.size()).c_str())<< std::endl;
        return atof(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
}

int MonteCarlo::findParameterInt(int currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atoi(l.substr(str.size()).c_str())<< std::endl;
        return atoi(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
}

long int MonteCarlo::findParameterLongInt(long int currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atoi(l.substr(str.size()).c_str())<< std::endl;
        return atol(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
}

unsigned long int MonteCarlo::findParameterULongInt(unsigned long int currentValue, std::string l, std::string str) const{
    if(l.find(str.c_str()) != std::string::npos){
        //std::cout<<l << " "<< atoi(l.substr(str.size()).c_str())<< std::endl;
        return stoul(l.substr(str.size()).c_str());
    }else{
        return currentValue;
    }
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
        
        MCP.cellsA = findParameterInt(MCP.cellsA, line, "cellsA: ");
        MCP.cellsB = findParameterInt(MCP.cellsB, line, "cellsB: ");
        MCP.cellsC = findParameterInt(MCP.cellsC, line, "cellsC: ");
        
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

string MonteCarlo::toStringMD(){
    string function = "MonteCarlo::toStringMD()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "flipAttempts: " << MD.flipAttempts << endl;
    stream << "successfulFlips: " << MD.successfulFlips << endl;
    stream << "range: " << MD.range << endl;
    stream << "numConfigsDone: " << MD.numConfigsDone << endl;
    stream << "isThermal: " << MD.isThermal << endl;
    
    stream << "numSweepsPerformedToCheckIfThermalized: " <<
    MD.numSweepsPerformedToCheckIfThermalized << endl;
    
    stream << "numSweepsUsedToThermalize: " <<
    MD.numSweepsUsedToThermalize << endl;
    stream << "numSweepsPerformed: " << MD.numSweepsPerformed << endl;
    stream << "timeOfInitialization: " << MD.timeOfInitialization << endl;
    stream << "totalTimeRunning: " << MD.totalTimeRunning << endl;
    stream << "numClockOverflows: " << MD.numClockOverflows << endl;
    
    stream << "Flip Success percentage: " <<
    MD.successfulFlips/(1.0*MD.flipAttempts)<<endl;
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringMD(ifstream &readline){
    string function = "MonteCarlo::readStringMD(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    
    string line = "";
    while(readline.good()){
        getline(readline, line);
        //cout << line << endl;
        
        MD.flipAttempts = findParameterInt(MD.flipAttempts,
                                           line,
                                           "flipAttempts: ");
        MD.successfulFlips = findParameterInt(MD.successfulFlips,
                                              line,
                                              "successfulFlips: ");
        MD.range = findParameterDbl(MD.range,
                                    line,
                                    "range: ");
        MD.numConfigsDone = findParameterInt(MD.numConfigsDone,
                                             line,
                                             "numConfigsDone: ");
        MD.isThermal = findParameterInt(MD.isThermal,
                                        line,
                                        "isThermal: ");
        
        MD.numSweepsPerformedToCheckIfThermalized =
        findParameterInt(MD.numSweepsPerformedToCheckIfThermalized,
                         line,
                         "numSweepsPerformedToCheckIfThermalized: ");
        
        MD.numSweepsUsedToThermalize =
        findParameterInt(MD.numSweepsUsedToThermalize,
                         line,
                         "numSweepsUsedToThermalize: ");
        
        
        MD.totalTimeRunning = findParameterInt(MD.totalTimeRunning,
                                               line,
                                               "totalTimeRunning: ");
        
        MD.numClockOverflows = findParameterInt(MD.numClockOverflows,
                                                line,
                                                "numClockOverflows: ");
        
        MD.numSweepsPerformed = findParameterInt(MD.numSweepsPerformed,
                                                 line,
                                                 "numSweepsPerformed: ");
        
        if(line.compare("Finish MonteCarlo::toStringMD()") == 0){
            break;
        }
        
    }
    
    //cout << toStringMD() << endl;
    
    
}

string MonteCarlo::toStringOP(){
    string function = "MonteCarlo::toStringOP()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Current snapshot, " << vectorOPAndType[i][j] << ": ";
            stream << std::scientific << std::setprecision(PR) <<
            OP.vectorOP[i][j] << endl;
        }
    }
    
    for(int i = 0; i < NUMSCALAROP; i++){
        stream << "Current snapshot, " << scalarOrderParameterName[i] << ": ";
        stream << std::scientific << std::setprecision(PR) <<
        OP.scalarOP[i] << endl;
    }
    
    
    //Correlation Function - Snapshot
    for(int count = 0; count < MCP.cellsA*2; count++){
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                stream << "Current Snapshot, Correlation Function - ";
                stream << "between " << count << " Nearest Neighbors - ";
                stream << vectorOPAndType[i][j] << ": " <<
                OP.vectorOPCorrFunc[i][j][count] << endl;
            }
        }
    }
    
    for(int count = 0; count < MCP.cellsA*2; count++){
        for(int i = 0; i < NUMSCALAROP; i++){
            stream << "Current Snapshot, Correlation Function - ";
            stream << "between " << count << " Nearest Neighbors - ";
            stream << scalarOrderParameterName[i] << ": " <<
            OP.scalarOPCorrFunc[i][count] << endl;
        }
    }
    
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}



string MonteCarlo::toStringRS(){
    string function = "MonteCarlo::toStringRS()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                stream << sumOverVectorOPNames[i][j][k];
                stream << std::scientific << std::setprecision(PR) <<
                RS.sumOverVectorOP[i][j][k] << endl;
            }
        }
    }
    
    for(int i = 0; i < NUMSCALAROP; i++){
        for(int j = 0; j < POWERSOP; j++){
            stream << sumOverScalarOPNames[i][j];
            stream << std::scientific << std::setprecision(PR) <<
            RS.sumOverScalarOP[i][j] << endl;
        }
    }
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringRS(ifstream &readline){
    string function = "MonteCarlo::toStringRS(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    string line = "";
    double value = 0;
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                getline(readline, line);
                //cout << line << endl;
                line = line.substr(sumOverVectorOPNames[i][j][k].size());
                value = atof(line.c_str());
                RS.sumOverVectorOP[i][j][k] = value;
                //cout << std::scientific << std::setprecision(PR) << value << endl;
            }
        }
    }
    
    for(int i = 0; i < NUMSCALAROP; i++){
        for(int j = 0; j < POWERSOP; j++){
            getline(readline, line);
            //cout << line << endl;
            line = line.substr(sumOverScalarOPNames[i][j].size());
            value = atof(line.c_str());
            RS.sumOverScalarOP[i][j] = value;
            //cout << std::scientific << std::setprecision(PR) << value << endl;
        }
    }
    
    
    //cout << toStringRS() << endl;
    
}

string MonteCarlo::toStringRSCorr(){
    string function = "MonteCarlo::toStringRSCorr()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "Sum over Configurations. Correlation function starting at " <<
    "(0,0,0,0) and incrementing across the lattice. First number is which " <<
    "number nearest neighbor, second is the correlation ([dot] product)" << endl;
    
    for(int count = 0; count < MCP.cellsA*2; count++){
        
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                for(int k = 0; k < POWERSOP; k++){
                    stream << "Corr: " << count << " " <<
                    sumOverVectorOPNames[i][j][k] << " " <<
                    std::scientific << std::setprecision(PR) <<
                    RS.sumOverVectorOPCorrFunc[i][j][k][count] << endl;
                }
            }
        }
        
        
        for(int i = 0; i < NUMSCALAROP; i++){
            for(int k = 0; k < POWERSOP; k++){
                stream << "Corr: " << count << " " <<
                sumOverScalarOPNames[i][k] << " " <<
                std::scientific << std::setprecision(PR) <<
                RS.sumOverScalarOPCorrFunc[i][k][count] << endl;
            }
        }
    }
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringRSCorr(ifstream &readline){
    string function = "MonteCarlo::readStringRSCorr(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream ss;
    ss.str("");
    
    string line = "";
    double value = 0;
    getline(readline, line);
    
    for(int count = 0; count < MCP.cellsA*2; count++){
        
        for(int i = 0; i < NUMVECTOROP; i++){
            for(int j = 0; j < TYPESOP; j++){
                for(int k = 0; k < POWERSOP; k++){
                    getline(readline, line);
                    //cout << line << endl;
                    ss.str("");
                    ss << "Corr: " << count << " " <<
                    sumOverVectorOPNames[i][j][k] << " ";
                    line = line.substr(ss.str().size());
                    value = atof(line.c_str());
                    RS.sumOverVectorOPCorrFunc[i][j][k][count] = value;
                    //cout << std::scientific << std::setprecision(PR) << value << endl;
                }
            }
        }
        
        
        for(int i = 0; i < NUMSCALAROP; i++){
            for(int k = 0; k < POWERSOP; k++){
                getline(readline, line);
                //cout << line << endl;
                ss.str("");
                ss << "Corr: " << count << " " <<
                sumOverVectorOPNames[i][k] << " ";
                line = line.substr(ss.str().size());
                value = atof(line.c_str());
                RS.sumOverScalarOPCorrFunc[i][k][count] = value;
                //cout << std::scientific << std::setprecision(PR) << value << endl;
            }
        }
    }
    
    //cout << toStringRSCorr() << endl;
    
}

string MonteCarlo::toStringSpinSnapshot(){
    string function = "MonteCarlo::toStringSpinSnapshot()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "Begin Components of Spin. Crystal Axis Coordinates: " <<
    "z (c axis), i (a axis), j (b axis), k (sublattice), " <<
    "Spin Components: x, y, z" << endl;
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    stream << "Spin: " << (int)c << " " << (int)a << " " <<
                    (int)b << " " <<
                    (int)s << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[c][a][b][s].x << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[c][a][b][s].y << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[c][a][b][s].z << endl;
                }
            }
        }
    }
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringSpinSnapshot(ifstream &readline){
    string function = "MonteCarlo::readStringSpinSnapshot(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream ss;
    ss.str("");
    
    string line = "";
    string cstr = "";
    string astr = "";
    string bstr = "";
    string sstr = "";
    string xstr = "";
    string ystr = "";
    string zstr = "";
    uint_fast8_t cint = 0;
    uint_fast8_t aint = 0;
    uint_fast8_t bint = 0;
    uint_fast8_t sint = 0;
    double xdoub = 0;
    double ydoub = 0;
    double zdoub = 0;
    getline(readline, line);
    
    for(int c = 0; c < MCP.cellsC; c++){
        for(int a = 0; a < MCP.cellsA; a++){
            for(int b = 0; b < MCP.cellsB; b++){
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    getline(readline, line);
                    //cout << line << endl;
                    ss.clear();
                    ss.str("");
                    ss << "Spin: ";
                    line = line.substr(ss.str().size());
                    ss.str(line);
                    ss >> cstr >> astr >> bstr >> sstr >> xstr >> ystr >> zstr;
                    cint = (uint_fast8_t) atoi(cstr.c_str());
                    aint = (uint_fast8_t) atoi(astr.c_str());
                    bint = (uint_fast8_t) atoi(bstr.c_str());
                    sint = (uint_fast8_t) atoi(sstr.c_str());
                    xdoub = (double) atof(xstr.c_str());
                    ydoub = (double) atof(ystr.c_str());
                    zdoub = (double) atof(zstr.c_str());
                    lattice[cint][aint][bint][sint].x = xdoub;
                    lattice[cint][aint][bint][sint].y = ydoub;
                    lattice[cint][aint][bint][sint].z = zdoub;
                    //cout << (int)cint << " " << (int)aint << " " <<
                    //(int)bint << " " << (int)sint <<
                    //" " << std::scientific << std::setprecision(PR) << xdoub <<
                    //" " << std::scientific << std::setprecision(PR) << ydoub <<
                    //" " << std::scientific << std::setprecision(PR) <<
                    //zdoub << endl;
                }
            }
        }
    }
    
    //cout << toStringSpinSnapshot() << endl;
    
}

string MonteCarlo::toStringFourierTransform(){
    string function = "MonteCarlo::toStringFourierTransform()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "Begin Reciprocal Lattice, S(q). Crystal Axis Coordinates: " <<
    "z (c axis), i (a axis), j (b axis), k (sublattice), " <<
    "Fourier Transform Components: ReX, ReY, ReZ, ImX, ImY, ImZ, " <<
    "magnitude" << endl;
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    stream << "FT: " << (int)c << " " << (int)a << " " <<
                    (int)b << " " <<
                    (int)s << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ReXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ReYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ReZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ImXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ImYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].ImZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[c][a][b][s].magnitude << endl;
                }
            }
        }
    }
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringFourierTransform(ifstream &readline){
    string function = "MonteCarlo::readStringFourierTransform(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream ss;
    ss.str("");
    
    string line = "";
    string cstr = "";
    string astr = "";
    string bstr = "";
    string sstr = "";
    string ReXstr = "";
    string ReYstr = "";
    string ReZstr = "";
    string ImXstr = "";
    string ImYstr = "";
    string ImZstr = "";
    string magstr = "";
    
    uint_fast8_t cint = 0;
    uint_fast8_t aint = 0;
    uint_fast8_t bint = 0;
    uint_fast8_t sint = 0;
    double ReXdoub = 0;
    double ReYdoub = 0;
    double ReZdoub = 0;
    double ImXdoub = 0;
    double ImYdoub = 0;
    double ImZdoub = 0;
    double magdoub = 0;
    getline(readline, line);
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    getline(readline, line);
                    //cout << line << endl;
                    ss.clear();
                    ss.str("");
                    ss << "FT: ";
                    line = line.substr(ss.str().size());
                    ss.str(line);
                    ss >> cstr >> astr >> bstr >> sstr >> ReXstr >> ReYstr >>
                    ReZstr >> ImXstr >> ImYstr >> ImZstr >> magstr;
                    cint = (uint_fast8_t) atoi(cstr.c_str());
                    aint = (uint_fast8_t) atoi(astr.c_str());
                    bint = (uint_fast8_t) atoi(bstr.c_str());
                    sint = (uint_fast8_t) atoi(sstr.c_str());
                    ReXdoub = (double) atof(ReXstr.c_str());
                    ReYdoub = (double) atof(ReYstr.c_str());
                    ReZdoub = (double) atof(ReZstr.c_str());
                    ImXdoub = (double) atof(ImXstr.c_str());
                    ImYdoub = (double) atof(ImYstr.c_str());
                    ImZdoub = (double) atof(ImZstr.c_str());
                    magdoub = (double) atof(magstr.c_str());
                    recipLattice[cint][aint][bint][sint].ReXComponent = ReXdoub;
                    recipLattice[cint][aint][bint][sint].ReYComponent = ReYdoub;
                    recipLattice[cint][aint][bint][sint].ReZComponent = ReZdoub;
                    recipLattice[cint][aint][bint][sint].ImXComponent = ImXdoub;
                    recipLattice[cint][aint][bint][sint].ImYComponent = ImYdoub;
                    recipLattice[cint][aint][bint][sint].ImZComponent = ImZdoub;
                    recipLattice[cint][aint][bint][sint].magnitude = magdoub;
                    //cout << (int)cint << " " << (int)aint << " " << (int)bint << " " <<
                    //(int)sint << " " <<
                    //std::scientific << std::setprecision(PR) << ReXdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ReYdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ReZdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImXdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImYdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImZdoub << " " <<
                    //std::scientific << std::setprecision(PR) << magdoub << endl;
                }
            }
        }
    }
    
    //cout << toStringFourierTransform() << endl;
    
}


string MonteCarlo::toStringSumOverFourierTransform(){
    string function = "MonteCarlo::toStringSumOverFourierTransform()";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream stream;
    stream.str("");
    stream << "Begin " << function.c_str() << endl;
    
    stream << "Begin Reciprocal Lattice, Sum Over S(q). " <<
    "Crystal Axis Coordinates: c (c axis), a (a axis), b (b axis), " <<
    "s (sublattice), Fourier Transform Components: ReX, ReY, ReZ, ImX, ImY, " <<
    "ImZ, magnitude" << endl;
    
    for(uint_fast8_t c = 0; c < MCP.cellsC; c++){
        for(uint_fast8_t a = 0; a < MCP.cellsA; a++){
            for(uint_fast8_t b = 0; b < MCP.cellsB; b++){
                for(uint_fast8_t s = 0; s < NUMSUBLATTICES; s++){
                    stream << "Sum over FT: " << (int)c << " " << (int)a << " " << (int)b <<
                    " " << (int)s << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ReXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ReYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ReZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ImXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ImYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].ImZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[c][a][b][s].magnitude << endl;
                }
            }
        }
    }
    
    stream << "Finish " << function.c_str() << endl;
    return stream.str();
}

void MonteCarlo::readStringSumOverFourierTransform(ifstream &readline){
    string function = "MonteCarlo::readSumOverStringFourierTransform(ifstream &readline)";
    cout << "Begin " << function.c_str() << endl;
    
    stringstream ss;
    ss.str("");
    
    string line = "";
    string cstr = "";
    string astr = "";
    string bstr = "";
    string sstr = "";
    string ReXstr = "";
    string ReYstr = "";
    string ReZstr = "";
    string ImXstr = "";
    string ImYstr = "";
    string ImZstr = "";
    string magstr = "";
    
    uint_fast8_t cint = 0;
    uint_fast8_t aint = 0;
    uint_fast8_t bint = 0;
    uint_fast8_t sint = 0;
    double ReXdoub = 0;
    double ReYdoub = 0;
    double ReZdoub = 0;
    double ImXdoub = 0;
    double ImYdoub = 0;
    double ImZdoub = 0;
    double magdoub = 0;
    getline(readline, line);
    
    for(int c = 0; c < MCP.cellsC; c++){
        for(int a = 0; a < MCP.cellsA; a++){
            for(int b = 0; b < MCP.cellsB; b++){
                for(int s = 0; s < NUMSUBLATTICES; s++){
                    getline(readline, line);
                    //cout << line << endl;
                    ss.clear();
                    ss.str("");
                    ss << "Sum over FT: ";
                    line = line.substr(ss.str().size());
                    ss.str(line);
                    ss >> cstr >> astr >> bstr >> sstr >> ReXstr >>
                    ReYstr >> ReZstr >> ImXstr >> ImYstr >> ImZstr >> magstr;
                    cint = (uint_fast8_t) atoi(cstr.c_str());
                    aint = (uint_fast8_t) atoi(astr.c_str());
                    bint = (uint_fast8_t) atoi(bstr.c_str());
                    sint = (uint_fast8_t) atoi(sstr.c_str());
                    ReXdoub = (double) atof(ReXstr.c_str());
                    ReYdoub = (double) atof(ReYstr.c_str());
                    ReZdoub = (double) atof(ReZstr.c_str());
                    ImXdoub = (double) atof(ImXstr.c_str());
                    ImYdoub = (double) atof(ImYstr.c_str());
                    ImZdoub = (double) atof(ImZstr.c_str());
                    magdoub = (double) atof(magstr.c_str());
                    sumOverRecipLattice[cint][aint][bint][sint].ReXComponent = ReXdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].ReYComponent = ReYdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].ReZComponent = ReZdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].ImXComponent = ImXdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].ImYComponent = ImYdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].ImZComponent = ImZdoub;
                    sumOverRecipLattice[cint][aint][bint][sint].magnitude = magdoub;
                    //cout << (int)cint << " " << (int)aint << " " << (int)bint << " " <<
                    //(int)sint << " " <<
                    //std::scientific << std::setprecision(PR) << ReXdoub <<
                    //" " << std::scientific << std::setprecision(PR) <<
                    //ReYdoub << " " << std::scientific << std::setprecision(PR) <<
                    //ReZdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImXdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImYdoub << " " <<
                    //std::scientific << std::setprecision(PR) << ImZdoub << " " <<
                    //std::scientific << std::setprecision(PR) << magdoub << endl;
                }
            }
        }
    }
    
    //cout << toStringSumOverFourierTransform() << endl;
    
}

string MonteCarlo::toString(){
    string function = "MonteCarlo::toString()";
    cout << "Begin " << function.c_str() << endl;
    
    //Print out the essential information about the Monte Carlo simulation.
    
    stringstream stream;
    stream.str("");
    stream << endl;
    stream << toStringMCP() << endl;
    stream << toStringMD() << endl;
    stream << toStringRS() << endl;
    stream << toStringSpinSnapshot() << endl;
    stream << toStringRSCorr() << endl;
    stream << toStringSumOverFourierTransform() << endl;
    
    return stream.str();
}


void MonteCarlo::finalPrint(){
    string function = "MonteCarlo::finalPrint()";
    cout << "Begin " << function.c_str() << endl;
    
    updateOrderParameters();
    
    //For convenience, I calculate some of the more commonly desired quantities
    //so that the most important raw data is human readable.
    
    //Averaging.
    double aveScalarOP[NUMVECTOROP][POWERSOP] = {0};
    double aveVectorOP[NUMVECTOROP][TYPESOP][POWERSOP] = {0};
    
    
    for(int i = 0; i < NUMSCALAROP; i++){
        for(int k = 0; k < POWERSOP; k++){
            aveScalarOP[i][k] = RS.sumOverScalarOP[i][k] / (MD.numConfigsDone*1.0);
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                aveVectorOP[i][j][k] = RS.sumOverVectorOP[i][j][k] / (MD.numConfigsDone*1.0);
            }
        }
    }
    
    //Calculate Lattice Quantity
    double suscepVectorOP[NUMVECTOROP][TYPESOP] = {0};
    double binderVectorOP[NUMVECTOROP][TYPESOP] = {0};
    double suscepScalarOP[NUMSCALAROP] = {0};
    double binderScalarOP[NUMSCALAROP] = {0};
    
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            suscepVectorOP[i][j] = numSites*(aveVectorOP[i][j][1] - aveVectorOP[i][j][0]*aveVectorOP[i][j][0])/(1.0*MCP.KbT);
            binderVectorOP[i][j] = 1.0 - aveVectorOP[i][j][3]/(3.0*aveVectorOP[i][j][1]*aveVectorOP[i][j][1]);
        }
    }
    
    for(int i = 0; i < NUMSCALAROP; i++){
        suscepScalarOP[i] = numSites*(aveScalarOP[i][1] - aveScalarOP[i][0]*aveScalarOP[i][0])/(1.0*MCP.KbT);
        binderScalarOP[i] = 1.0 - aveScalarOP[i][3]/(3.0*aveScalarOP[i][1]*aveScalarOP[i][1]);
    }
    
    
    double spHeat = 0;
    double energy = 0;
    for(int i = 0; i < NUMSCALAROP; i++){
        if (scalarOrderParameterName[i].compare("Energy") == 0){
            spHeat = (aveScalarOP[i][1] - (aveScalarOP[i][0]*aveScalarOP[i][0])) / (MCP.KbT*MCP.KbT*numSites*1.0);
            energy = aveScalarOP[i][0];
        }
    }
    
    //Printing out info
    std::stringstream stream;
    std::ofstream fileOut;
    time_t t = time(0);
    tm time = *localtime(&t);
    
    fileOut.open(outputFileName.c_str(), std::ios::out | std::ios::trunc);
    if (!fileOut.is_open() || !fileOut.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputFileName <<endl;
        exit(3161990);
    }
    stream.str("");
    stream << "Start characteristics for this configuration" << endl;
    stream << "Time: " << std::setw(2) << std::setfill('0') <<
    time.tm_mon + 1 << "_" << std::setw(2) << std::setfill('0') <<
    time.tm_mday << "_" << time.tm_year + 1900 << "_" << time.tm_hour + 1 << endl;
    if((MD.successfulFlips/(1.0*MD.flipAttempts)) * 100 < 50){
        cerr << "Alert! acceptance rate is too low!" << endl;
        cerr << "Flip Success percentage: " << MD.successfulFlips/(1.0*MD.flipAttempts)<<endl;
        cerr << "Range: " << MD.range <<endl;
        stream << "Alert! acceptance rate is too low!" << endl;
        stream << "Flip Success percentage: " << MD.successfulFlips/(1.0*MD.flipAttempts)<<endl;
        stream << "Range: " << MD.range <<endl;
    }
    
    stream << toString();
    
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, " << vectorOPAndType[i][j] << ": " <<
            std::scientific << std::setprecision(PR) << aveVectorOP[i][j][0] << endl;
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, Susceptibility - " << vectorOPAndType[i][j] << ": "<<
            std::scientific << std::setprecision(PR) << suscepVectorOP[i][j] << endl;
        }
    }
    
    for(int i = 0; i < NUMVECTOROP; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, Binder - " << vectorOPAndType[i][j] << ": "<<
            std::scientific << std::setprecision(PR) << binderVectorOP[i][j] << endl;
        }
    }
    
    //Scalar
    for(int i = 0; i < NUMSCALAROP; i++){
        stream << "Lattice Averaged Quantity, " << scalarOrderParameterName[i] << ": " <<
        std::scientific << std::setprecision(PR) << aveScalarOP[i][0] << endl;
    }
    
    
    for(int i = 0; i < NUMSCALAROP; i++){
        stream << "Lattice Averaged Quantity, Susceptibility - " << scalarOrderParameterName[i] << ": "<<
        std::scientific << std::setprecision(PR) << suscepScalarOP[i] << endl;
    }
    
    
    for(int i = 0; i < NUMSCALAROP; i++){
        stream << "Lattice Averaged Quantity, Binder - " << scalarOrderParameterName[i] << ": "<<
        std::scientific << std::setprecision(PR) << binderScalarOP[i] << endl;
    }
    
    
    stream << "Specific Heat: ";
    stream << std::scientific << std::setprecision(PR) << spHeat << endl;
    stream << "Energy Per Site: " << std::scientific <<
    std::setprecision(PR) << energy / (1.0*numSites) << endl;
    stream << "End information for this configuration" << endl;
    stream << endl;
    fileOut << stream.str() << std::endl;
    //std::cout << stream.str() << std::endl;
    fileOut.close();
    
}

void MonteCarlo::printSnapshot(){
    string function = "MonteCarlo::printSnapshot()";
    cout << "Begin " << function.c_str() << endl;
    /*
     updateOrderParameters();
     
     std::stringstream stream;
     std::ofstream fileOut;
     
     stringstream out;
     out.str("");
     string temp = outputFileName.substr(0,outputFileName.size()-4);
     out << temp.c_str() << "OPSnap.txt";
     
     fileOut.open(out.str().c_str(), std::ios::out | std::ios::app);
     stream.str("");
     */
    /*
     for(int i = 0; i < NUMVECTOROP; i++){
     for(int j = 0; j < SPINDIMEN; j++){
     stream << snapOPNamesSmall[i][j] << ": ";
     stream << std::scientific << std::setprecision(PR) << OP.snapOP[i][j] << endl;
     }
     }
     stream << "E: ";
     stream << std::scientific << std::setprecision(PR) << OP.E << endl;
     */
    /*
     double a1 = 0;
     double a2 = PI/3.0;
     double a3 = 2*PI/3.0;
     double a4 = PI;
     double a5 = 4*PI/3.0;
     double a6 = 5*PI/3.0;
     
     double m1 = 0;
     double m2 = 0;
     
     double N_X = OP.charOP[1][0];
     double N_Y = OP.charOP[1][1];
     double N_Z = OP.charOP[1][2];
     
     if(N_X > 0){m1 += fabs(N_X)*cos(a2); m2 += fabs(N_X)*sin(a2);}
     else       {m1 += fabs(N_X)*cos(a5); m2 += fabs(N_X)*sin(a5);}
     
     if(N_Y > 0){m1 += fabs(N_Y)*cos(a1); m2 += fabs(N_Y)*sin(a1);}
     else       {m1 += fabs(N_Y)*cos(a4); m2 += fabs(N_Y)*sin(a4);}
     
     if(N_Z > 0){m1 += fabs(N_Z)*cos(a6); m2 += fabs(N_Z)*sin(a6);}
     else       {m1 += fabs(N_Z)*cos(a3); m2 += fabs(N_Z)*sin(a3);}
     
     stream << "MN: " << std::scientific << std::setprecision(PR) << m1;
     stream << ", " << std::scientific << std::setprecision(PR) << m2 << endl;
     
     fileOut << stream.str() << std::endl;
     fileOut.close();
     */
    
}
