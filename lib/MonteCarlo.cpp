
/**
 MonteCarlo.cpp
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of Chalcogenides. Specifically, exploring the
 magnetic phases at finite temperature.
 
 @author Craig Price
 @version 2.0 2015/03/07
 */

#include <cstddef>
#include <cmath>
#include <ctime>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
#include "MonteCarlo.h"
using namespace std;


/******************************************************************************/
//Helper functions
/******************************************************************************/

void MonteCarlo::setParamNames(){
    cout<<"Set Names"<<endl;
    
    stringstream name;
    for(int i = 0; i < NUMVECTORS; i ++){
        for(int j = 0; j < TYPESOP; j++){
            name.str("");
            name << orderParameter[i] << "'s "<< typeOP[j];
            charOPNames[i][j] = name.str();
        }
    }
    
    for(int i = 0; i < NUMVECTORS; i ++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                name.str("");
                name << orderParameter[i] << "'s "<< typeOP[j] << powerOP[k];
                sumOPNames[i][j][k] = name.str();
            }
        }
    }
    
}

void MonteCarlo::setSpinComponent(string component){
    int xPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int yPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int zPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int sPos = atoi(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double xComp = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double yComp = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double zComp = atof(component.substr(0,component.find(" ")).c_str());
    lattice[xPos][yPos][sPos][zPos].x = xComp;
    lattice[xPos][yPos][sPos][zPos].y = yComp;
    lattice[xPos][yPos][sPos][zPos].z = zComp;
    //cout<< "Component: " << xPos << " " << yPos << " " << sPos << " " <<
    //zPos << " " << xComp<< " " << yComp << " " << zComp<<endl;
}

void MonteCarlo::setRecipLatComponent(string component){
    int xPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int yPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int zPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int sPos = atoi(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_x = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_y = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_z = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_x = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_y = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_z = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double magnitude = atof(component.substr(0,component.find(" ")).c_str());
    recipLattice[xPos][yPos][sPos][zPos].ReXComponent = ReSq_x;
    recipLattice[xPos][yPos][sPos][zPos].ReYComponent = ReSq_y;
    recipLattice[xPos][yPos][sPos][zPos].ReZComponent = ReSq_z;
    recipLattice[xPos][yPos][sPos][zPos].ImXComponent = ImSq_x;
    recipLattice[xPos][yPos][sPos][zPos].ImYComponent = ImSq_y;
    recipLattice[xPos][yPos][sPos][zPos].ImZComponent = ImSq_z;
    recipLattice[xPos][yPos][sPos][zPos].magnitude = magnitude;
    //cout<< "Reciprocal Lattice, S(q): " << xPos << " " << yPos << " " << sPos << " " <<
    //zPos << " " << ReSq_x << " " << ReSq_y << " " << ReSq_z << " " <<
    //ImSq_x << " " << ImSq_y << " " << ImSq_z << " " << magnitude << endl;
}

void MonteCarlo::setSumOverRecipLatComponent(string component){
    int xPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int yPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int zPos = atoi(component.substr(0,component.find(",")).c_str());
    component = component.substr(1 + component.find(","));
    int sPos = atoi(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_x = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_y = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ReSq_z = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_x = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_y = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double ImSq_z = atof(component.substr(0,component.find(" ")).c_str());
    component = component.substr(1 + component.find(" "));
    double magnitude = atof(component.substr(0,component.find(" ")).c_str());
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ReXComponent = ReSq_x;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ReYComponent = ReSq_y;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ReZComponent = ReSq_z;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ImXComponent = ImSq_x;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ImYComponent = ImSq_y;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].ImZComponent = ImSq_z;
    sumOverRecipLattice[xPos][yPos][sPos][zPos].magnitude = magnitude;
    //cout<< "Reciprocal Lattice, Sum Over S(q): " << xPos << " " << yPos << " " << sPos << " " <<
    //zPos << " " << ReSq_x << " " << ReSq_y << " " << ReSq_z << " " <<
    //ImSq_x << " " << ImSq_y << " " << ImSq_z << " " << magnitude << endl;
}

string MonteCarlo::toString(){
    
    //To print out the Lattice characteristics
    stringstream stream;
    stream.str("");
    stream << "KbT: " << MCP.KbT << endl;
    stream << "Cells-A: " << MCP.cellsA << endl;
    stream << "Cells-B: " << MCP.cellsB << endl;
    stream << "Cells-C: " << MCP.cellsC << endl;
    stream << "cubicD: " << MCP.cubicD << endl;
    stream << "j1: " << MCP.j1 << endl;
    stream << "j2: " << MCP.j2 << endl;
    stream << "j3: " << MCP.j3 << endl;
    stream << "k1: " << MCP.k1 << endl;
    stream << "k2: " << MCP.k2 << endl;
    stream << "k3: " << MCP.k3 << endl;
    stream << "Is External Magnetic Field?: " << MCP.isBField << endl;
    stream << "bField_x: "<< MCP.bField_x << endl;
    stream << "bField_y: "<< MCP.bField_y << endl;
    stream << "bField_z: "<< MCP.bField_z << endl;
    stream << "bFieldMag: "<< MCP.bFieldMag << endl;
    stream << "estimatedTc: " << MCP.estimatedTc << endl;
    stream << "numSweepsToPerformTotal: " << MCP.numSweepsToPerformTotal << endl;
    stream << "Monte Carlo Seed: " << MCP.monte_carlo_seed << endl;
    
    stream << "Attempted Flips: " << MD.flipAttempts << endl;
    stream << "Successful Flips: " << MD.successfulFlips << endl;
    stream << "Flip Range: " << MD.range << endl;
    stream << "numConfigsDone: " << MD.numConfigsDone << endl;
    stream << "isThermal: " << MD.isThermal << endl;
    stream << "numSweepsPerformedToCheckIfThermalized: " << MD.numSweepsPerformedToCheckIfThermalized << endl;
    stream << "numSweepsUsedToThermalize: " << MD.numSweepsUsedToThermalize << endl;
    stream << "numSweepsPerformed: " << MD.numSweepsPerformed << endl;
    stream << "timeOfInitialization: " << MD.timeOfInitialization << endl;
    stream << "totalTimeRunning: " << MD.totalTimeRunning << endl;
    
    stream << "Flip Success percentage: " << MD.successfulFlips/MD.flipAttempts<<endl;
    
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Current Snapshot, " << charOPNames[i][j] << ": ";
            stream << std::scientific << std::setprecision(PR) << OP.charOP[i][j] << endl;
        }
    }
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                stream << "Sum over Configurations, " << sumOPNames[i][j][k] << ": ";
                stream << std::scientific << std::setprecision(PR) << LC.sumOP[i][j][k] << endl;
            }
        }
    }
    
    
    //E
    stream << "Current Snapshot, Energy: ";
    stream << std::scientific << std::setprecision(PR) << OP.E << endl;
    stream << "Sum over Configurations, Energy to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumE[0] << endl;
    stream << "Sum over Configurations, Energy to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumE[1] << endl;
    stream << "Sum over Configurations, Energy to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumE[2] << endl;
    stream << "Sum over Configurations, Energy to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumE[3] << endl;
    //PX
    stream << "Current Snapshot, Dot Product with Spin in X direction: ";
    stream << std::scientific << std::setprecision(PR) << OP.PX << endl;
    stream << "Sum over Configurations, Dot Product with Spin in X direction to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPX[0] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in X direction to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPX[1] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in X direction to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPX[2] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in X direction to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPX[3] << endl;
    //PY
    stream << "Current Snapshot, Dot Product with Spin in Y direction: ";
    stream << std::scientific << std::setprecision(PR) << OP.PY << endl;
    stream << "Sum over Configurations, Dot Product with Spin in Y direction to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPY[0] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in Y direction to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPY[1] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in Y direction to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPY[2] << endl;
    stream << "Sum over Configurations, Dot Product with Spin in Y direction to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumPY[3] << endl;
    //Q
    stream << "Current Snapshot, Nematic Bi-Quadratic: ";
    stream << std::scientific << std::setprecision(PR) << OP.Q << endl;
    stream << "Sum over Configurations, Nematic Bi-Quadratic to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumQ[0] << endl;
    stream << "Sum over Configurations, Nematic Bi-Quadratic to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumQ[1] << endl;
    stream << "Sum over Configurations, Nematic Bi-Quadratic to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumQ[2] << endl;
    stream << "Sum over Configurations, Nematic Bi-Quadratic to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << LC.sumQ[3] << endl;
    
    //Spin snapshot
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    stream << "Components of Spin " << i << "," << j << "," <<
                    z << "," << k << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[i][j][k][z].x << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[i][j][k][z].y << " " <<
                    std::scientific << std::setprecision(PR) <<
                    lattice[i][j][k][z].z << endl;
                }
            }
        }
    }
    
    //Snapshot of Recip Lattice
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    stream << "Reciprocal Lattice, S(q): " << i << "," << j << "," <<
                    z << "," << k << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ReXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ReYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ReZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ImXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ImYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].ImZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    recipLattice[i][j][k][z].magnitude << endl;
                }
            }
        }
    }
    
    //Sum over recip lattice
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    
                    stream << "Reciprocal Lattice, Sum Over S(q): " << i << "," << j << "," <<
                    z << "," << k << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ReXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ReYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ReZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ImXComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ImYComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].ImZComponent << " " <<
                    std::scientific << std::setprecision(PR) <<
                    sumOverRecipLattice[i][j][k][z].magnitude << " " << endl;
                }
            }
        }
    }
    
    //Correlation Function - Snapshot
    for(int count = 0; count < MCP.cellsA*2; count++){
        for(int i = 0; i < NUMVECTORS; i++){
            for(int j = 0; j < TYPESOP; j++){
                stream << "Current Snapshot, Correlation Function - ";
                stream << "between " << count << " Nearest Neighbors - ";
                stream << charOPNames[i][j] << ": " << OP.corrFunc[i][j][count] << endl;
            }
        }
    }
    //E
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Current Snapshot, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Energy to the 1st power: " << OP.ECorrFunc[count] << endl;
    }
    //PX
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Current Snapshot, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Dot Product with Spin in X direction to the 1st power: " << OP.PXCorrFunc[count] << endl;
    }
    //PY
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Current Snapshot, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Dot Product with Spin in Y direction to the 1st power: " << OP.PYCorrFunc[count] << endl;
    }
    //Q
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Current Snapshot, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Nematic Bi-Quadratic to the 1st power: " << OP.QCorrFunc[count] << endl;
    }
    
    
    //Correlation function - SumChar
    for(int count = 0; count < MCP.cellsA*2; count++){
        for(int i = 0; i < NUMVECTORS; i++){
            for(int j = 0; j < TYPESOP; j++){
                stream << "Sum over Configurations, Correlation Function - ";
                stream << "between " << count << " Nearest Neighbors - ";
                stream << charOPNames[i][j] << ": " << LC.sumOPCorrFunc[i][j][count] << endl;
            }
        }
    }
    //E
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Sum over Configurations, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Energy to the 1st power: " << LC.sumECorrFunc[count] << endl;
    }
    //PX
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Sum over Configurations, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Dot Product with Spin in X direction to the 1st power: " << LC.sumPXCorrFunc[count] << endl;
    }
    //PY
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Sum over Configurations, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Dot Product with Spin in Y direction to the 1st power: " << LC.sumPYCorrFunc[count] << endl;
    }
    //Q
    for(int count = 0; count < MCP.cellsA*2; count++){
        stream << "Sum over Configurations, Correlation Function - ";
        stream << "between " << count << " Nearest Neighbors - ";
        stream << "Nematic Bi-Quadratic to the 1st power: " << LC.sumQCorrFunc[count] << endl;
    }
    
    
    return stream.str();
}



void MonteCarlo::finalPrint(){
    
    updateOrderParameters();
    
    //Averaging.
    double aveE[POWERSOP] = {0};
    double avePX[POWERSOP] = {0};
    double avePY[POWERSOP] = {0};
    double aveQ[POWERSOP] = {0};
    double aveOP[NUMVECTORS][TYPESOP][POWERSOP] = {0};
    
    //E
    aveE[0] = LC.sumE[0] / (MD.numConfigsDone*1.0);
    aveE[1] = LC.sumE[1] / (MD.numConfigsDone*1.0);
    aveE[2] = LC.sumE[2] / (MD.numConfigsDone*1.0);
    aveE[3] = LC.sumE[3] / (MD.numConfigsDone*1.0);//not divided by numSites yet.
    //PX
    avePX[0] = LC.sumPX[0] / (MD.numConfigsDone*1.0);
    avePX[1] = LC.sumPX[1] / (MD.numConfigsDone*1.0);
    avePX[2] = LC.sumPX[2] / (MD.numConfigsDone*1.0);
    avePX[3] = LC.sumPX[3] / (MD.numConfigsDone*1.0);
    //PY
    avePY[0] = LC.sumPY[0] / (MD.numConfigsDone*1.0);
    avePY[1] = LC.sumPY[1] / (MD.numConfigsDone*1.0);
    avePY[2] = LC.sumPY[2] / (MD.numConfigsDone*1.0);
    avePY[3] = LC.sumPY[3] / (MD.numConfigsDone*1.0);
    //Q
    aveQ[0] = LC.sumQ[0] / (MD.numConfigsDone*1.0);
    aveQ[1] = LC.sumQ[1] / (MD.numConfigsDone*1.0);
    aveQ[2] = LC.sumQ[2] / (MD.numConfigsDone*1.0);
    aveQ[3] = LC.sumQ[3] / (MD.numConfigsDone*1.0);
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                aveOP[i][j][k] = LC.sumOP[i][j][k] / (MD.numConfigsDone*1.0);
            }
        }
    }
    
    //Calculate Lattice Quantity
    double suscep[NUMVECTORS][TYPESOP] = {0};
    double binder[NUMVECTORS][TYPESOP] = {0};
    
    double suscepPX = 0;
    double suscepPY = 0;
    double suscepQ = 0;
    
    double binderPX = 0;
    double binderPY = 0;
    double binderQ = 0;
    
    double spHeat = 0;
    
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            suscep[i][j] = numSites*(aveOP[i][j][1] - aveOP[i][j][0]*aveOP[i][j][0])/MCP.KbT;
            binder[i][j] = 1.0 - aveOP[i][j][3]/(3.0*aveOP[i][j][1]*aveOP[i][j][1]);
        }
    }
    
    suscepPX = numSites*(avePX[1] - avePX[0]*avePX[0])/MCP.KbT;
    suscepPY = numSites*(avePY[1] - avePY[0]*avePY[0])/MCP.KbT;
    suscepQ = numSites*(aveQ[1] - aveQ[0]*aveQ[0])/MCP.KbT;
    
    binderPX = 1.0 - avePX[3]/(3.0*avePX[1]*avePX[1]);
    binderPY = 1.0 - avePY[3]/(3.0*avePY[1]*avePY[1]);
    binderQ = 1.0 - aveQ[3]/(3.0*aveQ[1]*aveQ[1]);
    
    spHeat = (aveE[1] - (aveE[0]*aveE[0])) / (MCP.KbT*MCP.KbT*numSites);
    
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
    if((MD.successfulFlips/MD.flipAttempts) * 100 < 50){
        cerr << "Alert! acceptance rate is too low!" << endl;
        cerr << "Flip Success percentage: " << MD.successfulFlips/MD.flipAttempts<<endl;
        cerr << "Range: " << MD.range <<endl;
        stream << "Alert! acceptance rate is too low!" << endl;
        stream << "Flip Success percentage: " << MD.successfulFlips/MD.flipAttempts<<endl;
        stream << "Range: " << MD.range <<endl;
    }
    
    stream << toString();
    
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, " << charOPNames[i][j] << ": " <<
            std::scientific << std::setprecision(PR) << aveOP[i][j][0] << endl;
        }
    }
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, Susceptibility - " << charOPNames[i][j] << ": "<<
            std::scientific << std::setprecision(PR) << suscep[i][j] << endl;
        }
    }
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            stream << "Lattice Averaged Quantity, Binder - " << charOPNames[i][j] << ": "<<
            std::scientific << std::setprecision(PR) << binder[i][j] << endl;
        }
    }
    
    //E
    stream << "Lattice Averaged Quantity, Energy to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << aveE[0] << endl;
    stream << "Lattice Averaged Quantity, Energy to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << aveE[1] << endl;
    stream << "Lattice Averaged Quantity, Energy to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << aveE[2] << endl;
    stream << "Lattice Averaged Quantity, Energy to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << aveE[3] << endl;
    //PX
    stream << "Lattice Averaged Quantity, Dot Product with Spin in X direction to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << avePX[0] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in X direction to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << avePX[1] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in X direction to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << avePX[2] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in X direction to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << avePX[3] << endl;
    //PY
    stream << "Lattice Averaged Quantity, Dot Product with Spin in Y direction to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << avePY[0] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in Y direction to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << avePY[1] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in Y direction to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << avePY[2] << endl;
    stream << "Lattice Averaged Quantity, Dot Product with Spin in Y direction to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << avePY[3] << endl;
    //Q
    stream << "Lattice Averaged Quantity, Nematic Bi-Quadratic to the 1st power: ";
    stream << std::scientific << std::setprecision(PR) << aveQ[0] << endl;
    stream << "Lattice Averaged Quantity, Nematic Bi-Quadratic to the 2nd power: ";
    stream << std::scientific << std::setprecision(PR) << aveQ[1] << endl;
    stream << "Lattice Averaged Quantity, Nematic Bi-Quadratic to the 3rd power: ";
    stream << std::scientific << std::setprecision(PR) << aveQ[2] << endl;
    stream << "Lattice Averaged Quantity, Nematic Bi-Quadratic to the 4th power: ";
    stream << std::scientific << std::setprecision(PR) << aveQ[3] << endl;
    //Sus
    //PX
    stream << "Lattice Averaged Quantity, Susceptibility - Dot Product with Spin in X direction: ";
    stream << std::scientific << std::setprecision(PR) << suscepPX << endl;
    //PY
    stream << "Lattice Averaged Quantity, Susceptibility - Dot Product with Spin in Y direction: ";
    stream << std::scientific << std::setprecision(PR) << suscepPY << endl;
    //Q
    stream << "Lattice Averaged Quantity, Susceptibility - Nematic Bi-Quadratic: ";
    stream << std::scientific << std::setprecision(PR) << suscepQ << endl;
    //
    //Bin
    //PX
    stream << "Lattice Averaged Quantity, Binder - Dot Product with Spin in X direction: ";
    stream << std::scientific << std::setprecision(PR) << binderPX << endl;
    //PY
    stream << "Lattice Averaged Quantity, Binder - Dot Product with Spin in Y direction: ";
    stream << std::scientific << std::setprecision(PR) << binderPY << endl;
    //Q
    stream << "Lattice Averaged Quantity, Binder - Nematic Bi-Quadratic: ";
    stream << std::scientific << std::setprecision(PR) << binderQ << endl;
    //
    
    
    stream << "Specific Heat: ";
    stream << std::scientific << std::setprecision(PR) << spHeat << endl;
    stream << "Energy Per Site: " << std::scientific <<
    std::setprecision(PR) << aveE[0] / numSites << endl;
    stream << "End information for this configuration" << endl;
    stream << endl;
    fileOut << stream.str() << std::endl;
    //std::cout << stream.str() << std::endl;
    fileOut.close();
    //printDiagram();
    
}

void MonteCarlo::printSnapshot(){
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
     for(int i = 0; i < NUMVECTORS; i++){
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



void MonteCarlo::rotateAllSpins(double phi, double theta){
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    specifyRotation(lattice[i][j][k][z].sc, phi, theta);
                }
            }
        }
    }
}


/******************************************************************************/
//Functions for initialization
/******************************************************************************/



void MonteCarlo::setLCToZero(){
    cout<<"Set Lat Chars to Zero"<<endl;
    
    //Setting Lat Chars to zero.
    for(int i = 0; i < POWERSOP; i++){
        LC.sumE[i] = 0;
        LC.sumPX[i] = 0;
        LC.sumPY[i] = 0;
        LC.sumQ[i] = 0;
    }
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            for(int k = 0; k < POWERSOP; k++){
                LC.sumOP[i][j][k] = 0;
            }
        }
    }
    
    //Setting correlation function to zero
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        LC.sumECorrFunc[count] = 0;
        LC.sumPXCorrFunc[count] = 0;
        LC.sumPYCorrFunc[count] = 0;
        LC.sumQCorrFunc[count] = 0;
        for(int i = 0; i < NUMVECTORS; i++){
            for(int j = 0; j < TYPESOP; j++){
                LC.sumOPCorrFunc[i][j][count] = 0;
            }
        }
    }
    
}

void MonteCarlo::setOPToZero(){
    cout<<"Set OP to Zero"<<endl;
    
    //Setting OP to zero.
    OP.E = 0;
    OP.PX = 0;
    OP.PY = 0;
    OP.Q = 0;
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            OP.charOP[i][j] = 0;
        }
    }
    
    //Setting correlation function to zero
    for(int count = 0; count < NUMCORRSPINSBIG; count++){
        OP.ECorrFunc[count] = 0;
        OP.PXCorrFunc[count] = 0;
        OP.PYCorrFunc[count] = 0;
        OP.QCorrFunc[count] = 0;
        for(int i = 0; i < NUMVECTORS; i++){
            for(int j = 0; j < TYPESOP; j++){
                OP.corrFunc[i][j][count] = 0;
            }
        }
    }
}

void MonteCarlo::setMCPToZero(){
    cout<<"Set MCP to Zero"<<endl;
    MCP.j1          = 0;
    MCP.j2          = 0;
    MCP.j3          = 0;
    MCP.k1          = 0;
    MCP.k2          = 0;
    MCP.k3          = 0;
    MCP.cubicD      = 0;
    
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
    cout<<"Set MD to Zero"<<endl;
    
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
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int k = 0; k < NUMSUBLATTICES; k++){
                for(int z = 0; z < MCP.cellsC; z++){
                    sumOverRecipLattice[i][j][k][z].ReXComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ReYComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ReZComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImXComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImYComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImZComponent = 0;
                    sumOverRecipLattice[i][j][k][z].magnitude = 0;
                    
                    recipLattice[i][j][k][z].ReXComponent = 0;
                    recipLattice[i][j][k][z].ReYComponent = 0;
                    recipLattice[i][j][k][z].ReZComponent = 0;
                    recipLattice[i][j][k][z].ImXComponent = 0;
                    recipLattice[i][j][k][z].ImYComponent = 0;
                    recipLattice[i][j][k][z].ImZComponent = 0;
                    recipLattice[i][j][k][z].magnitude = 0;
                }
            }
        }
    }
    
    
}

void MonteCarlo::initializeLattice(){
    //Constructs and initializes the lattice to random spins
    ClassicalSpin3D spin;
    SpinCoordinates sc;
    
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<ClassicalSpin3D> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<ClassicalSpin3D> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<ClassicalSpin3D> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(spin);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        lattice.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN_x.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN_y.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN_mx.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN_my.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN2_xy.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN2_mxy.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN2_mxmy.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN2_xmy.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN3_xx.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN3_mxmx.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN3_yy.push_back(row);
    }
    
    //array of pointers to the nearest neighbors of the site in the lattice
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<SpinCoordinates> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<SpinCoordinates> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<SpinCoordinates> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(sc);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        NN3_mymy.push_back(row);
    }
    
    FourierTransformOfSpin spinFT;
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<FourierTransformOfSpin> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<FourierTransformOfSpin> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(spinFT);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        recipLattice.push_back(row);
    }
    
    for(int i = 0; i < MCP.cellsA; i++){
        std::vector< std::vector< std::vector<FourierTransformOfSpin> > > row;
        for(int j = 0; j < MCP.cellsB; j++){
            std::vector< std::vector<FourierTransformOfSpin> > col;
            for(int k = 0; k < NUMSUBLATTICES; k++){
                std::vector<FourierTransformOfSpin> rise;
                for(int z = 0; z < MCP.cellsC; z++){
                    rise.push_back(spinFT);
                }
                col.push_back(rise);
            }
            row.push_back(col);
        }
        sumOverRecipLattice.push_back(row);
    }
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int k = 0; k < NUMSUBLATTICES; k++){
                for(int z = 0; z < MCP.cellsC; z++){
                    lattice[i][j][k][z].x = 0;
                    lattice[i][j][k][z].y = 0;
                    lattice[i][j][k][z].z = 1;
                    //lattice[i][j][k][z].setRandomOrientation();
                    lattice[i][j][k][z].sc.x = i;
                    lattice[i][j][k][z].sc.y = j;
                    lattice[i][j][k][z].sc.s = k;
                    lattice[i][j][k][z].sc.z = z;
                    
                    recipLattice[i][j][k][z].ReXComponent = 0;
                    recipLattice[i][j][k][z].ReYComponent = 0;
                    recipLattice[i][j][k][z].ReZComponent = 0;
                    recipLattice[i][j][k][z].ImXComponent = 0;
                    recipLattice[i][j][k][z].ImYComponent = 0;
                    recipLattice[i][j][k][z].ImZComponent = 0;
                    recipLattice[i][j][k][z].magnitude = 0;
                    
                    sumOverRecipLattice[i][j][k][z].ReXComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ReYComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ReZComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImXComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImYComponent = 0;
                    sumOverRecipLattice[i][j][k][z].ImZComponent = 0;
                    sumOverRecipLattice[i][j][k][z].magnitude = 0;
                }
            }
        }
    }
    
    numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES;
    setNearestNeighbors();
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheXDir(SpinCoordinates sc) const{
    uint_fast8_t ns;//next spin's sublattice
    
    switch (sc.s) {
        case 0:
            ns = 1;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 2:
            ns = 3;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 1:
            ns = 0;
            if(sc.x ==  MCP.cellsA - 1){
                return lattice[0][sc.y][ns][sc.z].sc;
            }else{
                return lattice[sc.x + 1][sc.y][ns][sc.z].sc;
            }
            break;
            
        case 3:
            ns = 2;
            if(sc.x ==  MCP.cellsA - 1){
                return lattice[0][sc.y][ns][sc.z].sc;
            }else{
                return lattice[sc.x + 1][sc.y][ns][sc.z].sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
    
    /*
     if(sc.s == 0){
     ns = 1;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 2){
     ns = 3;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 1){
     ns = 0;
     if(sc.x == MCP.cellsA - 1){
     return lattice[0][sc.y][ns][sc.z].sc;
     }else{
     return lattice[sc.x + 1][sc.y][ns][sc.z].sc;
     }
     }else if(sc.s == 3){
     ns = 2;
     if(sc.x == MCP.cellsA - 1){
     return lattice[0][sc.y][ns][sc.z].sc;
     }else{
     return lattice[sc.x + 1][sc.y][ns][sc.z].sc;
     }
     }else{
     std::cerr << "Wrong sublattice" << std::endl;
     exit(1);
     return sc;
     }
     */
    
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheMinusXDir(SpinCoordinates sc) const{
    uint_fast8_t ns;//next spin's sublattice
    
    switch (sc.s) {
        case 1:
            ns = 0;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 3:
            ns = 2;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 0:
            ns = 1;
            if(sc.x == 0){
                return lattice[MCP.cellsA - 1][sc.y][ns][sc.z].sc;
            }else{
                return lattice[sc.x - 1][sc.y][ns][sc.z].sc;
            }
            break;
            
        case 2:
            ns = 3;
            if(sc.x == 0){
                return lattice[MCP.cellsA - 1][sc.y][ns][sc.z].sc;
            }else{
                return lattice[sc.x - 1][sc.y][ns][sc.z].sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
    
    /*
     if(sc.s == 1){
     ns = 0;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 3){
     ns = 2;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 0){
     ns = 1;
     if(sc.x == 0){
     return lattice[MCP.cellsA - 1][sc.y][ns][sc.z].sc;
     }else{
     return lattice[sc.x - 1][sc.y][ns][sc.z].sc;
     }
     }else if(sc.s == 2){
     ns = 3;
     if(sc.x == 0){
     return lattice[MCP.cellsA - 1][sc.y][ns][sc.z].sc;
     }else{
     return lattice[sc.x - 1][sc.y][ns][sc.z].sc;
     }
     }else{
     std::cerr << "Wrong sublattice" << std::endl;
     exit(1);
     return sc;
     }
     */
    
}



SpinCoordinates MonteCarlo::getSpinCoordsInTheMinusYDir(SpinCoordinates sc) const{
    uint_fast8_t ns;//next spin's sublattice
    
    switch (sc.s) {
        case 2:
            ns = 0;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 3:
            ns = 1;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 0:
            ns = 2;
            if(sc.y == 0){
                return lattice[sc.x][MCP.cellsB - 1][ns][sc.z].sc;
            }else{
                return lattice[sc.x][sc.y - 1][ns][sc.z].sc;
            }
            break;
            
        case 1:
            ns = 3;
            if(sc.y == 0){
                return lattice[sc.x][MCP.cellsB - 1][ns][sc.z].sc;
            }else{
                return lattice[sc.x][sc.y - 1][ns][sc.z].sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
    
    /*
     if(sc.s == 2){
     ns = 0;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 3){
     ns = 1;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 0){
     ns = 2;
     if(sc.y == 0){
     return lattice[sc.x][MCP.cellsB - 1][ns][sc.z].sc;
     }else{
     return lattice[sc.x][sc.y - 1][ns][sc.z].sc;
     }
     }else if(sc.s == 1){
     ns = 3;
     if(sc.y == 0){
     return lattice[sc.x][MCP.cellsB - 1][ns][sc.z].sc;
     }else{
     return lattice[sc.x][sc.y - 1][ns][sc.z].sc;
     }
     }else{
     std::cerr << "Wrong sublattice" << std::endl;
     exit(1);
     return sc;
     }
     */
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheYDir(SpinCoordinates sc) const{
    uint_fast8_t ns;//next spin's sublattice
    
    switch (sc.s) {
        case 0:
            ns = 2;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 1:
            ns = 3;
            return lattice[sc.x][sc.y][ns][sc.z].sc;
            break;
            
        case 2:
            ns = 0;
            if(sc.y == MCP.cellsB - 1){
                return lattice[sc.x][0][ns][sc.z].sc;
            }else{
                return lattice[sc.x][sc.y + 1][ns][sc.z].sc;
            }
            break;
            
        case 3:
            ns = 1;
            if(sc.y == MCP.cellsB - 1){
                return lattice[sc.x][0][ns][sc.z].sc;
            }else{
                return lattice[sc.x][sc.y + 1][ns][sc.z].sc;
            }
            break;
            
        default:
            exit(1);
            break;
    }
    
    /*
     if(sc.s == 0){
     ns = 2;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 1){
     ns = 3;
     return lattice[sc.x][sc.y][ns][sc.z].sc;
     }else if(sc.s == 2){
     ns = 0;
     if(sc.y == MCP.cellsB - 1){
     return lattice[sc.x][0][ns][sc.z].sc;
     }else{
     return lattice[sc.x][sc.y + 1][ns][sc.z].sc;
     }
     }else if(sc.s == 3){
     ns = 1;
     if(sc.y == MCP.cellsB - 1){
     return lattice[sc.x][0][ns][sc.z].sc;
     }else{
     return lattice[sc.x][sc.y + 1][ns][sc.z].sc;
     }
     }else{
     std::cerr << "Wrong sublattice" << std::endl;
     exit(1);
     return sc;
     }
     */
    
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheUpDir(SpinCoordinates sc) const{
    
    if(MCP.cellsC == 1){
        return sc;
    }else{
        if(sc.z == MCP.cellsC - 1){
            return lattice[sc.x][sc.y][sc.s][0].sc;
        }else{
            return lattice[sc.x][sc.y][sc.s][sc.z + 1].sc;
        }
    }
    return sc;
}

SpinCoordinates MonteCarlo::getSpinCoordsInTheDownDir(SpinCoordinates sc) const{
    
    if(MCP.cellsC == 1){
        return sc;
    }else{
        if(sc.z == 0){
            return lattice[sc.x][sc.y][sc.s][MCP.cellsC - 1].sc;
        }else{
            return lattice[sc.x][sc.y][sc.s][sc.z - 1].sc;
        }
    }
    return sc;
    
}

void MonteCarlo::setNearestNeighbors(){
    SpinCoordinates sc;
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int k = 0; k < NUMSUBLATTICES; k++){
                for(int z = 0; z < MCP.cellsC; z++){
                    
                    NN_x[i][j][k][z] = getSpinCoordsInTheXDir(lattice[i][j][k][z].sc);
                    NN_y[i][j][k][z] = getSpinCoordsInTheYDir(lattice[i][j][k][z].sc);
                    NN_mx[i][j][k][z] = getSpinCoordsInTheMinusXDir(lattice[i][j][k][z].sc);
                    NN_my[i][j][k][z] = getSpinCoordsInTheMinusYDir(lattice[i][j][k][z].sc);
                    
                    sc = getSpinCoordsInTheXDir(lattice[i][j][k][z].sc);
                    NN2_xy[i][j][k][z] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(lattice[i][j][k][z].sc);
                    NN2_mxy[i][j][k][z] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(lattice[i][j][k][z].sc);
                    NN2_mxmy[i][j][k][z] = getSpinCoordsInTheMinusYDir(sc);
                    
                    sc = getSpinCoordsInTheXDir(lattice[i][j][k][z].sc);
                    NN2_xmy[i][j][k][z] = getSpinCoordsInTheMinusYDir(sc);
                    
                    sc = getSpinCoordsInTheXDir(lattice[i][j][k][z].sc);
                    NN3_xx[i][j][k][z] = getSpinCoordsInTheXDir(sc);
                    
                    sc = getSpinCoordsInTheMinusXDir(lattice[i][j][k][z].sc);
                    NN3_mxmx[i][j][k][z] = getSpinCoordsInTheMinusXDir(sc);
                    
                    sc = getSpinCoordsInTheYDir(lattice[i][j][k][z].sc);
                    NN3_yy[i][j][k][z] = getSpinCoordsInTheYDir(sc);
                    
                    sc = getSpinCoordsInTheMinusYDir(lattice[i][j][k][z].sc);
                    NN3_mymy[i][j][k][z] = getSpinCoordsInTheMinusYDir(sc);
                    
                }
            }
        }
    }
    
}


void MonteCarlo::init(){
    
    //initialize private variables
    
    setLCToZero();
    
    setOPToZero();
    
    setMCPToZero();
    
    setMDToZero();
    
    setParamNames();
    
    outputFileName           = "";
    
    
    cout<<"Finished with init()"<<endl;
    
}

MonteCarlo::MonteCarlo(){
    init();
    
    initializeLattice();
    
    
    cout<<"Finished with default MC constructor"<<endl;
    MD.timeOfInitialization = clock();
    
}

MonteCarlo::MonteCarlo(MCParameters input_parameters){
    
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
        
        MCP.bField_xNorm = MCP.bFieldMag * MCP.bField_x / size;
        MCP.bField_yNorm = MCP.bFieldMag * MCP.bField_y / size;
        MCP.bField_zNorm = MCP.bFieldMag * MCP.bField_z / size;
    }
    
    initializeLattice();
    
    MD.timeOfInitialization = clock();
    
    cout<<"Finished with (new) MC constructor"<<endl;
    
}

//This function initializes the lattice through copying data from a previous run.
MonteCarlo::MonteCarlo(bool reThermalize, string filename){
    init();
    
    //Reading in the rest of the lattice information from the saved file
    outputFileName = filename;
    
    string line = "";
    ifstream file;
    file.clear();
    file.open(outputFileName.c_str());
    if (!file.is_open() || !file.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputFileName <<endl;
        exit(3161990);
    }
    line = "";
    string search = "";
    while(file.good()){
        getline (file,line);
        if(line.find("End information for this configuration") !=
           string::npos){
            break;
        }
        
        //Retrieving the lattice properties from the text fil
        
        MCP.j1 = findParameterDbl(MCP.j1, line, "j1: ");
        MCP.j2 = findParameterDbl(MCP.j2, line, "j2: ");
        MCP.j3 = findParameterDbl(MCP.j3, line, "j3: ");
        MCP.k1 = findParameterDbl(MCP.k1, line, "k1: ");
        MCP.k2 = findParameterDbl(MCP.k2, line, "k2: ");
        MCP.k2 = findParameterDbl(MCP.k3, line, "k3: ");
        MCP.cubicD = findParameterDbl(MCP.cubicD, line, "cubicD: ");
        
        MCP.isBField = findParameterInt(MCP.isBField, line,
                                        "Is External Magnetic Field?: ");
        MCP.bField_x = findParameterDbl(MCP.bField_x, line, "bField_x: ");
        MCP.bField_y = findParameterDbl(MCP.bField_y, line, "bField_y: ");
        MCP.bField_z = findParameterDbl(MCP.bField_z, line, "bField_z: ");
        MCP.bFieldMag = findParameterDbl(MCP.bFieldMag, line, "bFieldMag: ");
        
        MCP.cellsA = findParameterInt(MCP.cellsA, line, "Cells-A: ");
        MCP.cellsB = findParameterInt(MCP.cellsB, line, "Cells-B: ");
        MCP.cellsC = findParameterInt(MCP.cellsC, line, "Cells-C: ");
        
        MCP.KbT = findParameterDbl(MCP.KbT, line, "KbT: ");
        MCP.estimatedTc = findParameterDbl(MCP.estimatedTc, line, "estimatedTc: ");
        MCP.numSweepsToPerformTotal = findParameterInt(MCP.numSweepsToPerformTotal,
                                                       line, "numSweepsToPerformTotal: ");
        
        MD.flipAttempts = findParameterInt(MD.flipAttempts, line, "flipAttempts: ");
        MD.successfulFlips = findParameterInt(MD.successfulFlips, line, "successfulFlips: ");
        MD.range = findParameterDbl(MD.range, line, "Flip Range: ");
        MD.numConfigsDone = findParameterInt(MD.numConfigsDone, line, "numConfigsDone: ");
        MD.isThermal = findParameterInt(MD.isThermal, line, "isThermal: ");
        MD.numSweepsPerformedToCheckIfThermalized =
        findParameterInt(MD.numSweepsPerformedToCheckIfThermalized,
                         line, "numSweepsPerformedToCheckIfThermalized: ");
        MD.numSweepsUsedToThermalize = findParameterInt(MD.numSweepsUsedToThermalize,
                                                        line, "numSweepsUsedToThermalize: ");
        MD.numSweepsPerformed = findParameterInt(MD.numSweepsPerformed,
                                                 line, "numSweepsPerformed: ");
        
        if ((!reThermalize)&&
            (MD.numSweepsPerformed >= MCP.numSweepsToPerformTotal)) {
            cout<<"!!!Finished Simulation!!!"<<endl;
            exit(0);
        }
        
        //Reading in the Lattice Characteristics
        if (line.find("Sum over Configurations, ") != string::npos) {
            stringstream ss;
            for(int i = 0; i < NUMVECTORS; i++){
                for(int j = 0; j < TYPESOP; j++){
                    for(int k = 0; k < POWERSOP; k++){
                        ss.str("");
                        ss << "Sum over Configurations, " << sumOPNames[i][j][k] << ": ";
                        LC.sumOP[i][j][k] = findParameterDbl(LC.sumOP[i][j][k], line,
                                                             ss.str());
                    }
                }
            }
            //E
            ss.str("");
            ss << "Sum over Configurations, Energy to the 1st power: ";
            LC.sumE[0] = findParameterDbl(LC.sumE[0], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Energy to the 2nd power: ";
            LC.sumE[1] = findParameterDbl(LC.sumE[1], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Energy to the 3rd power: ";
            LC.sumE[2] = findParameterDbl(LC.sumE[2], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Energy to the 4th power: ";
            LC.sumE[3] = findParameterDbl(LC.sumE[3], line, ss.str());
            //PX
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in X direction to the 1st power: ";
            LC.sumPX[0] = findParameterDbl(LC.sumPX[0], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in X direction to the 2nd power: ";
            LC.sumPX[1] = findParameterDbl(LC.sumPX[1], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in X direction to the 3rd power: ";
            LC.sumPX[2] = findParameterDbl(LC.sumPX[2], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in X direction to the 4th power: ";
            LC.sumPX[3] = findParameterDbl(LC.sumPX[3], line, ss.str());
            //PY
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in Y direction to the 1st power: ";
            LC.sumPY[0] = findParameterDbl(LC.sumPY[0], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in Y direction to the 2nd power: ";
            LC.sumPY[1] = findParameterDbl(LC.sumPY[1], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in Y direction to the 3rd power: ";
            LC.sumPY[2] = findParameterDbl(LC.sumPY[2], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Dot Product with Spin in Y direction to the 4th power: ";
            LC.sumPY[3] = findParameterDbl(LC.sumPY[3], line, ss.str());
            //Q
            ss.str("");
            ss << "Sum over Configurations, Nematic Bi-Quadratic to the 1st power: ";
            LC.sumQ[0] = findParameterDbl(LC.sumQ[0], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Nematic Bi-Quadratic to the 2nd power: ";
            LC.sumQ[1] = findParameterDbl(LC.sumQ[1], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Nematic Bi-Quadratic to the 3rd power: ";
            LC.sumQ[2] = findParameterDbl(LC.sumQ[2], line, ss.str());
            ss.str("");
            ss << "Sum over Configurations, Nematic Bi-Quadratic to the 4th power: ";
            LC.sumQ[3] = findParameterDbl(LC.sumQ[3], line, ss.str());
        }
        //End explicit labeling for scalar quantities
        
        //Correlation function
        if (line.find("Sum over Configurations, Correlation Function - ") != string::npos) {
            stringstream ss;
            for(int count = 0; count < MCP.cellsA*2; count++){
                ss.str("");
                ss << "Sum over Configurations, Correlation Function - between " << count << " Nearest Neighbors - ";
                if (line.find(ss.str()) != string::npos) {
                    for(int i = 0; i < NUMVECTORS; i++){
                        for(int j = 0; j < TYPESOP; j++){
                            ss.str("");
                            ss << "Sum over Configurations, Correlation Function - ";
                            ss << "between " << count << " Nearest Neighbors - ";
                            ss << charOPNames[i][j] << ": ";
                            LC.sumOPCorrFunc[i][j][count] = findParameterDbl(LC.sumOPCorrFunc[i][j][count],
                                                                             line, ss.str());
                        }
                    }
                }
            }
            for(int count = 0; count < MCP.cellsA*2; count++){
                //E
                ss.str("");
                ss << "Sum over Configurations, Correlation Function - ";
                ss << "between " << count << " Nearest Neighbors - ";
                ss << "Energy to the 1st power: ";
                LC.sumECorrFunc[count] = findParameterDbl(LC.sumECorrFunc[count], line, ss.str());
                //PX
                ss.str("");
                ss << "Sum over Configurations, Correlation Function - ";
                ss << "between " << count << " Nearest Neighbors - ";
                ss << "Dot Product with Spin in X direction to the 1st power: ";
                LC.sumPXCorrFunc[count] = findParameterDbl(LC.sumPXCorrFunc[count], line, ss.str());
                //PY
                ss.str("");
                ss << "Sum over Configurations, Correlation Function - ";
                ss << "between " << count << " Nearest Neighbors - ";
                ss << "Dot Product with Spin in Y direction to the 1st power: ";
                LC.sumPYCorrFunc[count] = findParameterDbl(LC.sumPYCorrFunc[count], line, ss.str());
                //Q
                ss.str("");
                ss << "Sum over Configurations, Correlation Function - ";
                ss << "between " << count << " Nearest Neighbors - ";
                ss << "Nematic Bi-Quadratic to the 1st power: ";
                LC.sumQCorrFunc[count] = findParameterDbl(LC.sumQCorrFunc[count], line, ss.str());
            }
        }
    }
    
    
    cout<<"Finished with reading in all lattice chars"<<endl;
    file.close();
    
    initializeLattice();
    
    //Copying Spin Components
    line = "";
    file.clear();
    file.open(outputFileName.c_str());
    if (!file.is_open() || !file.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputFileName <<endl;
        exit(3161990);
    }
    line = "";
    search = "";
    while(file.good()){
        getline (file,line);
        //cout<<line<<endl;
        if(line.find("End information for this configuration") !=
           string::npos){
            break;
        }
        search = "Components of Spin ";
        if(line.find(search.c_str()) != string::npos){
            setSpinComponent(line.substr(search.size()).c_str());
        }
        search = "Reciprocal Lattice, S(q): ";
        if(line.find(search.c_str()) != string::npos){
            setRecipLatComponent(line.substr(search.size()).c_str());
        }
        search = "Reciprocal Lattice, Sum Over S(q): ";
        if(line.find(search.c_str()) != string::npos){
            setSumOverRecipLatComponent(line.substr(search.size()).c_str());
        }
    }
    file.close();
    
    if(MCP.isBField){
        double size = sqrt(MCP.bField_x * MCP.bField_x +
                           MCP.bField_y * MCP.bField_y +
                           MCP.bField_z * MCP.bField_z);
        
        MCP.bField_xNorm = MCP.bFieldMag * MCP.bField_x / size;
        MCP.bField_yNorm = MCP.bFieldMag * MCP.bField_y / size;
        MCP.bField_zNorm = MCP.bFieldMag * MCP.bField_z / size;
    }
    
    
    cout<<"Finished with reading in all components"<<endl;
    
    //cout<< toString()<<endl;
    
    cout<<"Finished with MC constructor"<<endl;
    
    if(reThermalize){
        cout<< "ReThermalizing. Setting All Char to Zero and setting ";
        cout<< "Thermalization data to Zero." << endl;
        setLCToZero();
        setOPToZero();
        setMDToZero();
        setRecipToZero();
        
        for(int i = 0; i < MCP.cellsA; i++){
            for(int j = 0; j < MCP.cellsB; j++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    for(int z = 0; z < MCP.cellsC; z++){
                        setRandomOrientation(lattice[i][j][k][z].sc);
                        lattice[i][j][k][z].sc.x = i;
                        lattice[i][j][k][z].sc.y = j;
                        lattice[i][j][k][z].sc.s = k;
                        lattice[i][j][k][z].sc.z = z;
                        
                        recipLattice[i][j][k][z].ReXComponent = 0;
                        recipLattice[i][j][k][z].ReYComponent = 0;
                        recipLattice[i][j][k][z].ReZComponent = 0;
                        recipLattice[i][j][k][z].ImXComponent = 0;
                        recipLattice[i][j][k][z].ImYComponent = 0;
                        recipLattice[i][j][k][z].ImZComponent = 0;
                        recipLattice[i][j][k][z].magnitude = 0;
                        sumOverRecipLattice[i][j][k][z].ReXComponent = 0;
                        sumOverRecipLattice[i][j][k][z].ReYComponent = 0;
                        sumOverRecipLattice[i][j][k][z].ReZComponent = 0;
                        sumOverRecipLattice[i][j][k][z].ImXComponent = 0;
                        sumOverRecipLattice[i][j][k][z].ImYComponent = 0;
                        sumOverRecipLattice[i][j][k][z].ImZComponent = 0;
                        sumOverRecipLattice[i][j][k][z].magnitude = 0;
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
        double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60 * 60);
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
        addToMagStats();//This function is the only place this is called except sweep()
        MD.numSweepsPerformedToCheckIfThermalized += 1;
    }
    
    bool outOfTime = false;
    int thermalCycles = -1;
    while(!outOfTime){
        thermalCycles ++;
        cout<<"Thermal Cycle # " << thermalCycles<<endl;
        
        double clocks = clock() - MD.timeOfInitialization;
        double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60 * 60);
        if(hoursElapsed > durationOfSingleRun){
            outOfTime = true;
            return outOfTime;
        }
        
        LatticeCharacteristics firstLC = LC;
        
        int numSweepsAveOver = (int) (2000);
        for(int i = 0; i < numSweepsAveOver; i++){
            metropolisSweep();
            MD.numSweepsPerformed += 1;
            MD.numSweepsUsedToThermalize += 1;
            adjustRange();
            addToMagStats();
            MD.numSweepsPerformedToCheckIfThermalized += 1;
        }
        
        LatticeCharacteristics secLC = LC;
        
        double firstBinder = 0;
        double secBinder = 0;
        bool smallDiffBinder = true;
        for(int i = 0; i < NUMVECTORS; i++){
            double m4_1 = firstLC.sumOP[i][0][3] / (MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m2_1 = firstLC.sumOP[i][0][1] / (MD.numSweepsPerformedToCheckIfThermalized - numSweepsAveOver);
            double m4_2 = secLC.sumOP[i][0][3] / (MD.numSweepsPerformedToCheckIfThermalized);
            double m2_2 = secLC.sumOP[i][0][1] / (MD.numSweepsPerformedToCheckIfThermalized);
            
            firstBinder = 1.0 - m4_1/ (3.0 * m2_1 * m2_1);
            
            secBinder = 1.0 - m4_2/ (3.0 * m2_2 * m2_2);
            
            cout<<"Is Thermal? Diff is: "<<fabs(firstBinder - secBinder)<<endl;
            if(fabs(firstBinder - secBinder) > 1e-2){
                smallDiffBinder = false;
                break;
            }
            
        }
        
        if(smallDiffBinder){// thermalized.
            setLCToZero();
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
    double Q = 0;
    double M[SPINDIMEN] = {0};
    double N[SPINDIMEN] = {0};
    double SX[SPINDIMEN] = {0};
    double SY[SPINDIMEN] = {0};
    double Dimer = 0;
    SpinCoordinates sc;
    ClassicalSpin3D a1_x;
    ClassicalSpin3D a1_y;
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    
                    X = lattice[i][j][k][z].x;
                    Y = lattice[i][j][k][z].y;
                    Z = lattice[i][j][k][z].z;
                    
                    //OP: E
                    E += getLocalEnergy(lattice[i][j][k][z].sc, false);
                    //End OP: E
                    
                    //OP: M
                    M[0] += X;
                    M[1] += Y;
                    M[2] += Z;
                    //End OP: M
                    
                    //OP: N
                    //if k == 0,3 - count negative M, k == 1,2 - positive M
                    N[0] += ((k == 0 || k == 3) ? ((-1.0)*X) : (X));
                    N[1] += ((k == 0 || k == 3) ? ((-1.0)*Y) : (Y));
                    N[2] += ((k == 0 || k == 3) ? ((-1.0)*Z) : (Z));
                    //End OP: N
                    
                    //OP: SX
                    if((j % 2 == 0)){
                        SX[0] += X;
                        SX[1] += Y;
                        SX[2] += Z;
                    }else if(j % 2 == 1){
                        SX[0] += (-1.0) * X;
                        SX[1] += (-1.0) * Y;
                        SX[2] += (-1.0) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: SX
                    
                    //OP: SY
                    if((i % 2 == 0)){
                        SY[0] += X;
                        SY[1] += Y;
                        SY[2] += Z;
                    }else if((i % 2 == 1)){
                        SY[0] += (-1.0) * X;
                        SY[1] += (-1.0) * Y;
                        SY[2] += (-1.0) * Z;
                    }else{
                        std::cerr << "Shouldn't occur" << std::endl;
                        exit(1);
                    }
                    //End OP: S2
                    
                    /*
                     a1_x = getSpinInTheXDir(lattice[i][j][k][z]);
                     a1_y = getSpinInTheYDir(lattice[i][j][k][z]);
                     */
                    
                    sc = NN_x[i][j][k][z];
                    a1_x = lattice[sc.x][sc.y][sc.s][sc.z];
                    sc = NN_y[i][j][k][z];
                    a1_y = lattice[sc.x][sc.y][sc.s][sc.z];
                    
                    PX += (X * a1_x.x + Y * a1_x.y + Z * a1_x.z);
                    PY += (X * a1_y.x + Y * a1_y.y + Z * a1_y.z);
                    Q += PX - PY;
                    
                    Dimer += 0;
                    
                }
            }
        }
    }
    //End OP
    
    OP.E = E; //Total E, E not per site
    
    OP.PX = PX;
    OP.PY = PY;
    OP.Q = Q;
    
    //All Order parameters are per site
    const double numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES * 1.0;
    
    for(int i = 0; i < SPINDIMEN; i++){
        M[i] = M[i] / numSites;
        N[i] = N[i] / numSites;
        SX[i] = SX[i] / numSites;
        SY[i] = SY[i] / numSites;
    }
    
    OP.PX = OP.PX / numSites;
    OP.PY = OP.PY / numSites;
    OP.Q = OP.Q / numSites;
    
    //Now to calculate the characteristics of the order parameters.
    double OPArr[NUMVECTORS][SPINDIMEN] = {0};
    for(int i = 0; i < SPINDIMEN; i++){
        OPArr[0][i] = M[i];
        OPArr[1][i] = N[i];
        OPArr[2][i] = SX[i];
        OPArr[3][i] = SY[i];
        OPArr[4][i] = sqrt(SX[i]*SX[i] + SY[i]*SY[i]);
    }
    
    for(int i = 0; i < NUMVECTORS; i++){
        
        //Scalar quantity eg. M
        OP.charOP[i][0] = sqrt(OPArr[i][0]*OPArr[i][0] +
                               OPArr[i][1]*OPArr[i][1] +
                               OPArr[i][2]*OPArr[i][2]);
        
        
        //fabs x, y, z, directions
        OP.charOP[i][1] = fabs(OPArr[i][0]);
        
        OP.charOP[i][2] = fabs(OPArr[i][1]);
        
        OP.charOP[i][3] = fabs(OPArr[i][2]);
        
        
        
    }
    
    if((MD.numConfigsDone % CORRTOUPDATES == 0)&&(MD.numConfigsDone > 10)){
        updateCorrelationFunction();
    }
    
    if((MD.numConfigsDone % FTTOUPDATES == 0)&&(MD.numConfigsDone > 10)){
        updateFourierTransformOnRecipLattice();
    }
    
}

void MonteCarlo::updateCorrelationFunction(){
    
    double E[NUMCORRSPINSBIG] = {0};
    double PX[NUMCORRSPINSBIG] = {0};
    double PY[NUMCORRSPINSBIG] = {0};
    double Q[NUMCORRSPINSBIG] = {0};
    double X = 0;
    double Y = 0;
    double Z = 0;
    double M[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double N[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double SX[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    double SY[SPINDIMEN][NUMCORRSPINSBIG] = {0};
    SpinCoordinates sc;
    ClassicalSpin3D a1_x;
    ClassicalSpin3D a1_y;
    
    int j = 0;
    int i = 0;
    int k = 0;
    int z = 0;
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        
        //This control structure picks the diagonal spins out across the lattice.
        if(count % 2 == 0 ){
            k = 0;
            i = count/2;
            j = count/2;
        }else{
            k = 3;
            i = (int) (count/2.0);
            j = (int) (count/2.0);
        }
        
        
        z = 0;
        
        //OP: E
        E[count] = getLocalEnergy(lattice[i][j][k][z].sc, false);
        //End OP: E
        
        
        X = lattice[i][j][k][z].x;
        Y = lattice[i][j][k][z].y;
        Z = lattice[i][j][k][z].z;
        
        //OP: M
        M[0][count] += X;
        M[1][count] += Y;
        M[2][count] += Z;
        //End OP: M
        
        //OP: N
        //if k == 0 - count negative M, k == 1 - positive M
        N[0][count] += ((k == 0 || k == 3) ? ((-1.0)*X) : (X));
        N[1][count] += ((k == 0 || k == 3) ? ((-1.0)*Y) : (Y));
        N[2][count] += ((k == 0 || k == 3) ? ((-1.0)*Z) : (Z));
        //End OP: N
        
        //OP: SX
        if((j % 2 == 0)){
            SX[0][count] += X;
            SX[1][count] += Y;
            SX[2][count] += Z;
        }else if((j % 2 == 1)){
            SX[0][count] += (-1.0) * X;
            SX[1][count] += (-1.0) * Y;
            SX[2][count] += (-1.0) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: SX
        
        //OP: SY
        if((i % 2 == 0)){
            SY[0][count] += X;
            SY[1][count] += Y;
            SY[2][count] += Z;
        }else if((i % 2 == 1)){
            SY[0][count] += (-1) * X;
            SY[1][count] += (-1) * Y;
            SY[2][count] += (-1) * Z;
        }else{
            std::cerr << "Shouldn't occur" << std::endl;
            exit(1);
        }
        //End OP: SY
        
        /*
         a1_x = getSpinInTheXDir(lattice[i][j][k][z]);
         a1_y = getSpinInTheYDir(lattice[i][j][k][z]);
         */
        
        sc = NN_x[i][j][k][z];
        a1_x = lattice[sc.x][sc.y][sc.s][sc.z];
        sc = NN_y[i][j][k][z];
        a1_y = lattice[sc.x][sc.y][sc.s][sc.z];
        
        PX[count] += (X * a1_x.x + Y * a1_x.y + Z * a1_x.z);
        PY[count] += (X * a1_y.x + Y * a1_y.y + Z * a1_y.z);
        Q[count] += PX[count] - PY[count];
        
        
    }
    
    
    double OPArr[NUMVECTORS][SPINDIMEN][NUMCORRSPINSBIG] = {0};
    
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
    
    
    
    //Calculate correlation function, dot products of spins
    for(int i = 0; i < NUMVECTORS; i++){
        
        //OP:
        for(int count = 0; count < 2*MCP.cellsA; count++){
            OP.corrFunc[i][0][count] = (OPArr[i][0][0] * OPArr[i][0][count] +
                                        OPArr[i][1][0] * OPArr[i][1][count] +
                                        OPArr[i][2][0] * OPArr[i][2][count]);
            
            
            
            //spin components
            OP.corrFunc[i][1][count] = OPArr[i][0][0] * OPArr[i][0][count];
            
            OP.corrFunc[i][2][count] = OPArr[i][1][0] * OPArr[i][1][count];
            
            OP.corrFunc[i][3][count] = OPArr[i][2][0] * OPArr[i][2][count];
            
            
            
        }
    }
    
    for(int count = 0; count < 2*MCP.cellsA; count++){
        OP.ECorrFunc[count] = E[0] * E[count];
        OP.PXCorrFunc[count] = PX[0] * PX[count];
        OP.PYCorrFunc[count] = PY[0] * PY[count];
        OP.QCorrFunc[count] = Q[0] * Q[count];
    }
    
    
};

//Change this for sublattices
void MonteCarlo::updateFourierTransformOnRecipLattice(){
    
    cout<<MD.numConfigsDone<<endl;
    cout<<FTTOUPDATES<<endl;
    double kvector_x = 0;
    double kvector_y = 0;
    double rvector_x = 0;
    double rvector_y = 0;
    double kr = 0;
    double cos_kr = 0;
    double sin_kr = 0;
    //Loop over all k points
    for(int k_i = 0; k_i < MCP.cellsA; k_i++){
        for(int k_j = 0; k_j < MCP.cellsB; k_j++){
            
            //Loop over all sites in real space
            for(int i = 0; i < MCP.cellsA; i++){
                for(int j = 0; j < MCP.cellsB; j++){
                    //sublattice A
                    kvector_x = k_i * B_VECTOR_1[0] / (MCP.cellsA*1.0) + k_j * B_VECTOR_2[0] / (MCP.cellsB*1.0);
                    kvector_y = k_i * B_VECTOR_1[1] / (MCP.cellsA*1.0) + k_j * B_VECTOR_2[1] / (MCP.cellsB*1.0);
                    rvector_x = i * A_VECTOR_1[0] + j * A_VECTOR_2[0];
                    rvector_y = i * A_VECTOR_1[1] + j * A_VECTOR_2[1];
                    kr = kvector_x * rvector_x + kvector_y * rvector_y;
                    cos_kr = cos(kr);
                    sin_kr = sin(kr);
                    recipLattice[k_i][k_j][0][0].ReXComponent += cos(kr) * lattice[i][j][0][0].x / (numSites*1.0);//A,B,s,C,x/y/z/ix/iy/iz
                    recipLattice[k_i][k_j][0][0].ImXComponent += sin_kr * lattice[i][j][0][0].x / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ReYComponent += cos_kr * lattice[i][j][0][0].y / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ImYComponent += sin_kr * lattice[i][j][0][0].y / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ReZComponent += cos_kr * lattice[i][j][0][0].z / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ImZComponent += sin_kr * lattice[i][j][0][0].z / (numSites*1.0);
                    
                }
            }
            
            
            recipLattice[k_i][k_j][0][0].magnitude = (
                                                      recipLattice[k_i][k_j][0][0].ReXComponent * recipLattice[k_i][k_j][0][0].ReXComponent +
                                                      recipLattice[k_i][k_j][0][0].ReYComponent * recipLattice[k_i][k_j][0][0].ReYComponent +
                                                      recipLattice[k_i][k_j][0][0].ReZComponent * recipLattice[k_i][k_j][0][0].ReZComponent +
                                                      recipLattice[k_i][k_j][0][0].ImXComponent * recipLattice[k_i][k_j][0][0].ImXComponent +
                                                      recipLattice[k_i][k_j][0][0].ImYComponent * recipLattice[k_i][k_j][0][0].ImYComponent +
                                                      recipLattice[k_i][k_j][0][0].ImZComponent * recipLattice[k_i][k_j][0][0].ImZComponent
                                                      );
            
            recipLattice[k_i][k_j][1][0].magnitude = (0);
            
        }
    }
    
}

void MonteCarlo::addToMagStats(){
    
    updateOrderParameters();
    double temp;
    
    //E
    LC.sumE[0] += OP.E;
    temp = OP.E*OP.E;
    LC.sumE[1] += temp;
    LC.sumE[2] += temp * OP.E;
    LC.sumE[3] += temp * temp;
    
    //PX
    LC.sumPX[0] += OP.PX;
    temp = OP.PX*OP.PX;
    LC.sumPX[1] += temp;
    LC.sumPX[2] += temp * OP.PX;
    LC.sumPX[3] += temp * temp;
    
    //PY
    LC.sumPY[0] += OP.PY;
    temp = OP.PY*OP.PY;
    LC.sumPY[1] += temp;
    LC.sumPY[2] += temp * OP.PY;
    LC.sumPY[3] += temp * temp;
    
    //Q
    LC.sumQ[0] += OP.Q;
    temp = OP.Q*OP.Q;
    LC.sumQ[1] += temp;
    LC.sumQ[2] += temp * OP.Q;
    LC.sumQ[3] += temp * temp;
    
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            
            temp = OP.charOP[i][j]*OP.charOP[i][j];
            
            LC.sumOP[i][j][0] += OP.charOP[i][j];
            LC.sumOP[i][j][1] += temp;
            LC.sumOP[i][j][2] += OP.charOP[i][j]*temp;
            LC.sumOP[i][j][3] += temp*temp;
        }
    }
    
    if((MD.numConfigsDone % CORRTOUPDATES == 0)&&(MD.numConfigsDone > 0)){
        //Correlation Function
        for(int count = 0; count < 2*MCP.cellsA; count++){
            LC.sumECorrFunc[count] += OP.ECorrFunc[count];
        }
        
        for(int count = 0; count < 2*MCP.cellsA; count++){
            LC.sumPXCorrFunc[count] += OP.PXCorrFunc[count];
        }
        for(int count = 0; count < 2*MCP.cellsA; count++){
            LC.sumPYCorrFunc[count] += OP.PYCorrFunc[count];
        }
        for(int count = 0; count < 2*MCP.cellsA; count++){
            LC.sumQCorrFunc[count] += OP.QCorrFunc[count];
        }
        
        for(int count = 0; count < 2*MCP.cellsA; count++){
            for(int i = 0; i < NUMVECTORS; i++){
                for(int j = 0; j < TYPESOP; j++){
                    LC.sumOPCorrFunc[i][j][count] += OP.corrFunc[i][j][count];
                }
            }
        }
    }
    
    if((MD.numConfigsDone % FTTOUPDATES == 0)&&(MD.numConfigsDone > 0)){
        //Loop over all k points
        for(int k_i = 0; k_i < MCP.cellsA; k_i++){
            for(int k_j = 0; k_j < MCP.cellsB; k_j++){
                
                //sublattice A
                sumOverRecipLattice[k_i][k_j][0][0].ReXComponent += recipLattice[k_i][k_j][0][0].ReXComponent;//A,B,s,C,x/y/z/ix/iy/iz
                sumOverRecipLattice[k_i][k_j][0][0].ImXComponent += recipLattice[k_i][k_j][0][0].ReYComponent;
                sumOverRecipLattice[k_i][k_j][0][0].ReYComponent += recipLattice[k_i][k_j][0][0].ReZComponent;
                sumOverRecipLattice[k_i][k_j][0][0].ImYComponent += recipLattice[k_i][k_j][0][0].ImXComponent;
                sumOverRecipLattice[k_i][k_j][0][0].ReZComponent += recipLattice[k_i][k_j][0][0].ImYComponent;
                sumOverRecipLattice[k_i][k_j][0][0].ImZComponent += recipLattice[k_i][k_j][0][0].ImZComponent;
                sumOverRecipLattice[k_i][k_j][0][0].magnitude += recipLattice[k_i][k_j][0][0].magnitude;
                
                //sublattice B
                sumOverRecipLattice[k_i][k_j][1][0].ReXComponent += recipLattice[k_i][k_j][1][0].ReXComponent;//A,B,s,C,x/y/z/ix/iy/iz
                sumOverRecipLattice[k_i][k_j][1][0].ImXComponent += recipLattice[k_i][k_j][1][0].ReYComponent;
                sumOverRecipLattice[k_i][k_j][1][0].ReYComponent += recipLattice[k_i][k_j][1][0].ReZComponent;
                sumOverRecipLattice[k_i][k_j][1][0].ImYComponent += recipLattice[k_i][k_j][1][0].ImXComponent;
                sumOverRecipLattice[k_i][k_j][1][0].ReZComponent += recipLattice[k_i][k_j][1][0].ImYComponent;
                sumOverRecipLattice[k_i][k_j][1][0].ImZComponent += recipLattice[k_i][k_j][1][0].ImZComponent;
                sumOverRecipLattice[k_i][k_j][1][0].magnitude += recipLattice[k_i][k_j][1][0].magnitude;
                
            }
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

void MonteCarlo::setSpinToZero(SpinCoordinates sc){
    lattice[sc.x][sc.y][sc.s][sc.z].x = 0;
    lattice[sc.x][sc.y][sc.s][sc.z].y = 0;
    lattice[sc.x][sc.y][sc.s][sc.z].z = 0;
}

void MonteCarlo::setSpin(SpinCoordinates sc, double x_, double y_, double z_){
    lattice[sc.x][sc.y][sc.s][sc.z].x = x_;
    lattice[sc.x][sc.y][sc.s][sc.z].y = y_;
    lattice[sc.x][sc.y][sc.s][sc.z].z = z_;
}

void MonteCarlo::checkSpin(SpinCoordinates sc){
    double size = sqrt(lattice[sc.x][sc.y][sc.s][sc.z].x*
                       lattice[sc.x][sc.y][sc.s][sc.z].x +
                       lattice[sc.x][sc.y][sc.s][sc.z].y*
                       lattice[sc.x][sc.y][sc.s][sc.z].y +
                       lattice[sc.x][sc.y][sc.s][sc.z].z*
                       lattice[sc.x][sc.y][sc.s][sc.z].z);
    lattice[sc.x][sc.y][sc.s][sc.z].x = lattice[sc.x][sc.y][sc.s][sc.z].x / size;
    lattice[sc.x][sc.y][sc.s][sc.z].y = lattice[sc.x][sc.y][sc.s][sc.z].y / size;
    lattice[sc.x][sc.y][sc.s][sc.z].z = lattice[sc.x][sc.y][sc.s][sc.z].z / size;
    if((lattice[sc.x][sc.y][sc.s][sc.z].x > 1 || lattice[sc.x][sc.y][sc.s][sc.z].x < -1)||
       (lattice[sc.x][sc.y][sc.s][sc.z].y > 1 || lattice[sc.x][sc.y][sc.s][sc.z].y < -1)||
       (lattice[sc.x][sc.y][sc.s][sc.z].z > 1 || lattice[sc.x][sc.y][sc.s][sc.z].z < -1)||
       (lattice[sc.x][sc.y][sc.s][sc.z].x != lattice[sc.x][sc.y][sc.s][sc.z].x)||
       (lattice[sc.x][sc.y][sc.s][sc.z].y != lattice[sc.x][sc.y][sc.s][sc.z].y)||
       (lattice[sc.x][sc.y][sc.s][sc.z].z != lattice[sc.x][sc.y][sc.s][sc.z].z)||
       ((fabs(sqrt(lattice[sc.x][sc.y][sc.s][sc.z].x *
                   lattice[sc.x][sc.y][sc.s][sc.z].x +
                   lattice[sc.x][sc.y][sc.s][sc.z].y *
                   lattice[sc.x][sc.y][sc.s][sc.z].y +
                   lattice[sc.x][sc.y][sc.s][sc.z].z *
                   lattice[sc.x][sc.y][sc.s][sc.z].z) - 1)) >= 0.000001) ){
        std::cerr << "Error: Bad Spin values: " << "x: " << lattice[sc.x][sc.y][sc.s][sc.z].x << std::endl;
        std::cerr << "Error: Bad Spin values: " << "y: " << lattice[sc.x][sc.y][sc.s][sc.z].y << std::endl;
        std::cerr << "Error: Bad Spin values: " << "z: " << lattice[sc.x][sc.y][sc.s][sc.z].z << std::endl;
        std::cerr << "Error: Bad Spin values: " << "Magnitude: " <<
        std::setprecision(15) << std::setw(15) <<
        sqrt(lattice[sc.x][sc.y][sc.s][sc.z].x *
             lattice[sc.x][sc.y][sc.s][sc.z].x +
             lattice[sc.x][sc.y][sc.s][sc.z].y *
             lattice[sc.x][sc.y][sc.s][sc.z].y +
             lattice[sc.x][sc.y][sc.s][sc.z].z *
             lattice[sc.x][sc.y][sc.s][sc.z].z) << std::endl;
        exit(1);
    }
}

void MonteCarlo::setRandomOrientation(SpinCoordinates sc){
    lattice[sc.x][sc.y][sc.s][sc.z].x = 2 * drand48() - 1;
    lattice[sc.x][sc.y][sc.s][sc.z].y = 2 * drand48() - 1;
    lattice[sc.x][sc.y][sc.s][sc.z].z = 2 * drand48() - 1;
    double size = sqrt(lattice[sc.x][sc.y][sc.s][sc.z].x *
                       lattice[sc.x][sc.y][sc.s][sc.z].x +
                       lattice[sc.x][sc.y][sc.s][sc.z].y *
                       lattice[sc.x][sc.y][sc.s][sc.z].y +
                       lattice[sc.x][sc.y][sc.s][sc.z].z *
                       lattice[sc.x][sc.y][sc.s][sc.z].z);
    lattice[sc.x][sc.y][sc.s][sc.z].x = lattice[sc.x][sc.y][sc.s][sc.z].x / size;
    lattice[sc.x][sc.y][sc.s][sc.z].y = lattice[sc.x][sc.y][sc.s][sc.z].y / size;
    lattice[sc.x][sc.y][sc.s][sc.z].z = lattice[sc.x][sc.y][sc.s][sc.z].z / size;
    flip(lattice[sc.x][sc.y][sc.s][sc.z].sc, PI-0.01);
}


void MonteCarlo::specifyRotation(SpinCoordinates sc, double rotPhi, double rotTheta){
    
    //Converts to spherical (in order to construct a perpendicular vector).
    //measured from the vertical.
    double theta = acos(lattice[sc.x][sc.y][sc.s][sc.z].z);
    double phi = atan2(lattice[sc.x][sc.y][sc.s][sc.z].y,
                       lattice[sc.x][sc.y][sc.s][sc.z].x);
    
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
    
    lattice[sc.x][sc.y][sc.s][sc.z].x = sin(theta) * cos(phi);
    lattice[sc.x][sc.y][sc.s][sc.z].y = sin(theta) * sin(phi);
    lattice[sc.x][sc.y][sc.s][sc.z].z = cos(theta);
    
    checkSpin(lattice[sc.x][sc.y][sc.s][sc.z].sc);
    
}
void MonteCarlo::reset(SpinCoordinates sc){
    //lattice[sc.x][sc.y][sc.s][sc.z].x = lattice[sc.x][sc.y][sc.s][sc.z].xOld;
    //lattice[sc.x][sc.y][sc.s][sc.z].y = lattice[sc.x][sc.y][sc.s][sc.z].yOld;
    //lattice[sc.x][sc.y][sc.s][sc.z].z = lattice[sc.x][sc.y][sc.s][sc.z].zOld;
}
void MonteCarlo::clear(SpinCoordinates sc){
    lattice[sc.x][sc.y][sc.s][sc.z].x = 0;
    lattice[sc.x][sc.y][sc.s][sc.z].y = 0;
    lattice[sc.x][sc.y][sc.s][sc.z].z = 0;
}
/**Does rotation using quaternions. Requires that the spin has cartesian
 *components but that the axis that it rotates around is in spherical.
 *
 *Note: this function should be replaced with a function that only uses
 *spherical coordinates or only cartesian coordinates.
 */
void MonteCarlo::rotate(SpinCoordinates sc, double theta_, double phi_, double Beta) {
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
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].x;
            }else if(i == 0 && j == 1){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].y;
            }else if(i == 0 && j == 2){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].z;
            }else if(i == 1 && j == 0){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].x;
            }else if(i == 1 && j == 1){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].y;
            }else if(i == 1 && j == 2){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].z;
            }else if(i == 2 && j == 0){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].x;
            }else if(i == 2 && j == 1){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].y;
            }else if(i == 2 && j == 2){
                newS[i] += M[i][j] * lattice[sc.x][sc.y][sc.s][sc.z].z;
            }else{
                std::cerr << "Bad Rotation" << std::endl;
                exit(1);
            }
        }
    }
    
    //In Cartesian
    lattice[sc.x][sc.y][sc.s][sc.z].x = newS[0];
    lattice[sc.x][sc.y][sc.s][sc.z].y = newS[1];
    lattice[sc.x][sc.y][sc.s][sc.z].z = newS[2];
}

/**
 *First we construct a new spin perpendicular to the initial configuration by
 *adding PI/2 to the zenith angle. Then, this new vector is rotated, within the
 *same plane, anywhere in 2*PI radians. Lastly, the origial vector is rotated
 *around the perp spin by the angle "range"
 */
void MonteCarlo::flip(SpinCoordinates sc, double range){
    
    
    double u[3];
    //Cross product with a sufficiently far away "g" to find a perp vector.
    //u[0] = g2s3 - g3s2;
    //u[1] = g3s1 - g1s3;
    //u[2] = g1s2 - g2s1;
    if(lattice[sc.x][sc.y][sc.s][sc.z].x > 0.57) {//g = (0,1,0)
        u[0] = lattice[sc.x][sc.y][sc.s][sc.z].z;
        u[1] = 0;
        u[2] = (-1) * lattice[sc.x][sc.y][sc.s][sc.z].x;
    }else if(lattice[sc.x][sc.y][sc.s][sc.z].y > 0.57){//g = (0,0,1)
        u[0] = (-1) * lattice[sc.x][sc.y][sc.s][sc.z].y;
        u[1] = lattice[sc.x][sc.y][sc.s][sc.z].x;
        u[2] = 0;
    }
    else{//g = (1,0,0)
        u[0] = 0;
        u[1] = (-1) * lattice[sc.x][sc.y][sc.s][sc.z].z;
        u[2] = lattice[sc.x][sc.y][sc.s][sc.z].y;
    }
    
    //Cross product to find third Orthogonal vector
    double v[3] = { u[1]*
        lattice[sc.x][sc.y][sc.s][sc.z].z -
        u[2]*lattice[sc.x][sc.y][sc.s][sc.z].y,
        u[2]*lattice[sc.x][sc.y][sc.s][sc.z].x -
        u[0]*lattice[sc.x][sc.y][sc.s][sc.z].z,
        u[0]*lattice[sc.x][sc.y][sc.s][sc.z].y -
        u[1]*lattice[sc.x][sc.y][sc.s][sc.z].x};
    
    //Normalize vectors
    double size = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    u[0] = u[0] / size;
    u[1] = u[1] / size;
    u[2] = u[2] / size;
    
    size = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] = v[0] / size;
    v[1] = v[1] / size;
    v[2] = v[2] / size;
    
    //spin, u, v are orthonormal.
    //See these two sites:
    //http://www.math.niu.edu/~rusin/known-math/96/sph.rand
    //http://objectmix.com/graphics/314285-random-points-spherical-cap.html
    //(copied below)
    //Except that the second link needs to be modified in a number of ways to
    //apply correctly to this situation. Most notably, before adding the radial
    //vector (from a vector in the plane perp to the spin) to the spin vector,
    //we need to shorten the spin vector such that the addition of the radial
    //vector lands on the surface of the sphere. We cannot simply normalize
    //at the end like they say.
    
    //Size of circle:
    if(range >= PI){range = PI;}
    double r = drand48() * 2.0 * PI;
    double distFromOrig = cos(range);//works both for pos dist and neg dist.
    double circDist = drand48() * (1 - distFromOrig) + distFromOrig;
    double circleRadius = sqrt(1 - circDist*circDist);
    
    double cosr = cos(r);
    double sinr = sin(r);
    
    lattice[sc.x][sc.y][sc.s][sc.z].x = lattice[sc.x][sc.y][sc.s][sc.z].x * circDist;
    lattice[sc.x][sc.y][sc.s][sc.z].y = lattice[sc.x][sc.y][sc.s][sc.z].y * circDist;
    lattice[sc.x][sc.y][sc.s][sc.z].z = lattice[sc.x][sc.y][sc.s][sc.z].z * circDist;
    
    lattice[sc.x][sc.y][sc.s][sc.z].x = lattice[sc.x][sc.y][sc.s][sc.z].x + circleRadius * (cosr * u[0] + sinr * v[0]);
    lattice[sc.x][sc.y][sc.s][sc.z].y = lattice[sc.x][sc.y][sc.s][sc.z].y + circleRadius * (cosr * u[1] + sinr * v[1]);
    lattice[sc.x][sc.y][sc.s][sc.z].z = lattice[sc.x][sc.y][sc.s][sc.z].z + circleRadius * (cosr * u[2] + sinr * v[2]);
    
    /*
     //Should not be needed.
     size = sqrt(x*x + y*y + z*z);
     x = x / size;
     y = y / size;
     z = z / size;
     */
    
}

/******************************************************************************/
//Core Functions
/******************************************************************************/



double MonteCarlo::getLocalEnergy(SpinCoordinates a, const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
    
    double temp = (dotProd(a, NN_x[a.x][a.y][a.s][a.z]) +
                   dotProd(a, NN_mx[a.x][a.y][a.s][a.z]) +
                   dotProd(a, NN_y[a.x][a.y][a.s][a.z]) +
                   dotProd(a, NN_my[a.x][a.y][a.s][a.z]));
    
    double interactionH = MCP.j1 * (temp);
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (dotProd(a, NN2_xy[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN2_mxy[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN2_mxmy[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN2_xmy[a.x][a.y][a.s][a.z]));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(a, NN3_xx[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN3_mxmx[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN3_yy[a.x][a.y][a.s][a.z]) +
                                     dotProd(a, NN3_mymy[a.x][a.y][a.s][a.z]));
    
    double interactionK = MCP.k1 * temp * temp;
    
    double interactionCubic = MCP.cubicD * (-1.0) * a.z * a.z;
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic) / 2.0;
    }
    
    
    //B field
    if (MCP.isBField) {
        double interactionB = (-1) * (MCP.bField_xNorm * a.x +
                                      MCP.bField_yNorm * a.y +
                                      MCP.bField_zNorm * a.z);
        totalInteraction += interactionB;
    }
    
    return totalInteraction;
    
}

double MonteCarlo::getLocalEnergyFromCalc(SpinCoordinates a, const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
    double temp = (dotProd(a, getSpinCoordsInTheXDir(a)) +
                   dotProd(a, getSpinCoordsInTheMinusXDir(a)) +
                   dotProd(a, getSpinCoordsInTheYDir(a)) +
                   dotProd(a, getSpinCoordsInTheMinusYDir(a)));
    
    double interactionH = MCP.j1 * (temp);
    
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (dotProd(a, getSpinCoordsInTheYDir(getSpinCoordsInTheXDir(a))) +
                                     dotProd(a, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheYDir(a))) +
                                     dotProd(a, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusYDir(a))) +
                                     dotProd(a, getSpinCoordsInTheXDir(getSpinCoordsInTheMinusYDir(a))));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (dotProd(a, getSpinCoordsInTheXDir(getSpinCoordsInTheXDir(a))) +
                                     dotProd(a, getSpinCoordsInTheMinusXDir(getSpinCoordsInTheMinusXDir(a))) +
                                     dotProd(a, getSpinCoordsInTheYDir(getSpinCoordsInTheYDir(a))) +
                                     dotProd(a, getSpinCoordsInTheMinusYDir(getSpinCoordsInTheMinusYDir(a))));
    
    double interactionK = MCP.k1 * temp * temp;
    
    double interactionCubic = MCP.cubicD * (-1.0) * a.z * a.z;
    
    if(fromSpinFlip)
    {
        totalInteraction += interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic;
    }else
    {//From Finding the total E of the lattice
        totalInteraction += (interactionH + interactionH2 + interactionH3 + interactionK + interactionCubic) / 2.0;
    }
    
    
    //B field
    if (MCP.isBField) {
        double interactionB = (-1) * (MCP.bField_xNorm * a.x +
                                      MCP.bField_yNorm * a.y +
                                      MCP.bField_zNorm * a.z);
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
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int k = 0; k < NUMSUBLATTICES; k++){
                for(int z = 0; z < MCP.cellsC; z++)
                {
                    //getSpin(i, j, k).checkSpin();
                    initialEnergy = getLocalEnergy(lattice[i][j][k][z].sc , true);
                    xOld = lattice[i][j][k][z].x;
                    yOld = lattice[i][j][k][z].y;
                    zOld = lattice[i][j][k][z].z;
                    flip(lattice[i][j][k][z].sc, MD.range);
                    finalEnergy = getLocalEnergy(lattice[i][j][k][z].sc, true);
                    //finalEnergy = getLocalEnergyFromCalc(lattice[i][j][k][z].sc, true);
                    changeInEnergy = finalEnergy - initialEnergy;
                    MD.flipAttempts += 1;
                    MD.successfulFlips += 1;
                    //The lattice automatically keeps the spin flip unless the
                    //energy is > 0
                    
                    if(changeInEnergy > 0)
                    {
                        /**
                         *Then it keeps the flip if a function f(KbT, delta_E) is
                         *less than a certain probability (0,1).
                         */
                        if(drand48() >= exp((-1) * (1.0/MCP.KbT) * changeInEnergy))
                        {
                            
                            //resets spin
                            //reset(lattice[i][j][k][z].sc);
                            lattice[i][j][k][z].x = xOld;
                            lattice[i][j][k][z].y = yOld;
                            lattice[i][j][k][z].z = zOld;
                            MD.successfulFlips -= 1;
                        }
                        else
                        {
                            //Keeps Spin flip if rand < exp()
                        }
                    }
                    else if(changeInEnergy == 0 && initialEnergy == 0)
                    {
                        //resets spin
                        //reset(lattice[i][j][k][z].sc);
                        lattice[i][j][k][z].x = xOld;
                        lattice[i][j][k][z].y = yOld;
                        lattice[i][j][k][z].z = zOld;
                        MD.successfulFlips -= 1;
                    }
                }
            }
        }
    }
}




void MonteCarlo::sweep(int numIndep_, double durationOfSingleRun){
    
    for(int sweeps = MD.numSweepsPerformed; sweeps <= MCP.numSweepsToPerformTotal; sweeps++)
    {
        
        
        metropolisSweep();
        MD.numSweepsPerformed += 1;
        if(sweeps % numIndep_ == 0)
        {
            addToMagStats();
            
            
            double clocks = clock() - MD.timeOfInitialization;
            double hoursElapsed = clocks * (1.0/(CLOCKS_PER_SEC*1.0)) / (60 * 60);
            if(hoursElapsed > durationOfSingleRun)
            {
                cout<<"Ran out of time, in sweep()"<<endl;
                break;
            }
            
            
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
 
 --
 Dave Seaman			dseaman@purdue.edu
 ++++ stop the execution of Mumia Abu-Jamal ++++
 ++++ if you agree copy these lines to your sig ++++
 ++++ see http://www.xs4all.nl/~tank/spg-l/sigaction.htm ++++
 ==============================================================================
 
 From: George Marsaglia <geo@stat.fsu.edu>
 Subject: Random points on a sphere.
 Date: Tue, 15 Jun 1999 18:03:07 -0400
 Newsgroups: sci.math.num-analysis,sci.math,sci.stat.math
 
 Questions on methods for random sampling from the surface
 of a sphere keep coming up, and seemingly get lost as
 new generations of computer users encounter the problem.
 
 I described the methods that I had used for years in
 "Sampling from the surface of a sphere", Ann Math Stat, v43, 1972,
 recommending the two that seemed the fastest:
 
 1) Choose standard normal variates x1,x2,x3,
 put R=1/sqrt(x1^2+x2^2+x3^2) and return the point
 (x1*R,x2*R,x3*R).
 
 2) Generate u and v, uniform in [-1,1]  until S=u^2+v^2 < 1,
 then return (1-2*S,u*r,v*r), with r=2*sqrt(1-S).
 
 The second method is fast and easiest to program---unless
 one already has a fast normal generator in his library.
 
 Method 1) seems by far the fastest (and easily applies to higher
 dimensions) if my latest normal generator is used.
 It produces normal variates at the rate of seven million
 per second in a 300MHz Pc.
 
 George Marsaglia
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




