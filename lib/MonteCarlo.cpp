
/**
 MonteCarlo.cpp
 Purpose: Manages Classical Monte Carlo Program and contains the lattice of
 spins. It is for the simulation of Chalcogenides. Specifically, exploring the
 magnetic phases at finite temperature.
 
 @author Craig Price
 @version 1.0 2015/03/04
 */

#include "MonteCarlo.h"
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
    
    double numSites = MCP.cellsA * MCP.cellsB * MCP.cellsC * NUMSUBLATTICES * 1.0;
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
    
    fileOut.open(outputMag.c_str(), std::ios::out | std::ios::trunc);
    if (!fileOut.is_open() || !fileOut.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputMag <<endl;
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
     string temp = outputMag.substr(0,outputMag.size()-4);
     out << temp.c_str() << "OPSnap.txt";
     
     fileOut.open(out.str().c_str(), std::ios::out | std::ios::app);
     stream.str("");
     
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
                    lattice[i][j][k][z].specifyRotation(phi, theta);
                }
            }
        }
    }
}

double MonteCarlo::dotProd(const double arr1[], const double arr2[]){
    return (arr1[0] * arr2[0] +
            arr1[1] * arr2[1] +
            arr1[2] * arr2[2]);
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
                    lattice[i][j][k][z].xPos = i;
                    lattice[i][j][k][z].yPos = j;
                    lattice[i][j][k][z].sPos = k;
                    lattice[i][j][k][z].zPos = z;
                    
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
}



void MonteCarlo::init(){
    
    //initialize private variables
    
    setLCToZero();
    
    setOPToZero();
    
    setMCPToZero();
    
    setMDToZero();
    
    setParamNames();
    
    outputMag           = "";
    
    
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
    
    outputMag           = MCP.path;
    
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
    outputMag = filename;
    
    string line = "";
    ifstream file;
    file.clear();
    file.open(outputMag.c_str());
    if (!file.is_open() || !file.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputMag <<endl;
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
        
        MCP.KbT = findParameterDbl(MCP.KbT, line, "KbT: ");
        MCP.cellsA = findParameterInt(MCP.cellsA, line, "Cells-A: ");
        MCP.cellsB = findParameterInt(MCP.cellsB, line, "Cells-B: ");
        MCP.cellsC = findParameterInt(MCP.cellsC, line, "Cells-C: ");
        MCP.cubicD = findParameterDbl(MCP.cubicD, line, "cubicD: ");
        MCP.j1 = findParameterDbl(MCP.j1, line, "j1: ");
        MCP.j2 = findParameterDbl(MCP.j2, line, "j2: ");
        MCP.j3 = findParameterDbl(MCP.j3, line, "j3: ");
        MCP.k1 = findParameterDbl(MCP.k1, line, "k1: ");
        MCP.k2 = findParameterDbl(MCP.k2, line, "k2: ");
        MCP.k2 = findParameterDbl(MCP.k3, line, "k3: ");
        MCP.isBField = findParameterInt(MCP.isBField, line, "Is External Magnetic Field?: ");
        MCP.bField_x = findParameterDbl(MCP.bField_x, line, "bField_x: ");
        MCP.bField_y = findParameterDbl(MCP.bField_y, line, "bField_y: ");
        MCP.bField_z = findParameterDbl(MCP.bField_z, line, "bField_z: ");
        MCP.bFieldMag = findParameterDbl(MCP.bFieldMag, line, "bFieldMag: ");
        MCP.estimatedTc = findParameterDbl(MCP.estimatedTc, line, "estimatedTc: ");
        MCP.numSweepsToPerformTotal = findParameterInt(MCP.numSweepsToPerformTotal, line, "numSweepsToPerformTotal: ");
        
        MD.range = findParameterDbl(MD.range, line, "Flip Range: ");
        MD.numConfigsDone = findParameterInt(MD.numConfigsDone, line, "numConfigsDone: ");
        MD.isThermal = findParameterInt(MD.isThermal, line, "isThermal: ");
        MD.numSweepsPerformedToCheckIfThermalized =
        findParameterInt(MD.numSweepsPerformedToCheckIfThermalized,
                         line,
                         "numSweepsPerformedToCheckIfThermalized: ");
        MD.numSweepsPerformed = findParameterInt(MD.numSweepsPerformed, line, "numSweepsPerformed: ");
        MD.numSweepsUsedToThermalize = findParameterInt(MD.numSweepsUsedToThermalize, line, "numSweepsUsedToThermalize: ");
        
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
    file.open(outputMag.c_str());
    if (!file.is_open() || !file.good()){
        cerr << "Can't open file: " << endl;
        cerr << outputMag <<endl;
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
                        lattice[i][j][k][z].setRandomOrientation();
                        lattice[i][j][k][z].xPos = i;
                        lattice[i][j][k][z].yPos = j;
                        lattice[i][j][k][z].sPos = k;
                        lattice[i][j][k][z].zPos = z;
                        
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
    
    double oldKbT = MCP.KbT;
    MCP.KbT = KbT_;
    
    
    for(int sweeps = MD.numSweepsUsedToThermalize; sweeps <= numSweepsToDo; sweeps++){
        int clocks = clock() - MD.timeOfInitialization;
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
    double M[SPINDIMEN] = {0};
    double N[SPINDIMEN] = {0};
    double SX[SPINDIMEN] = {0};
    double SY[SPINDIMEN] = {0};
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int z = 0; z < MCP.cellsC; z++){
                for(int k = 0; k < NUMSUBLATTICES; k++){
                    
                    X = lattice[i][j][k][z].x;
                    Y = lattice[i][j][k][z].y;
                    Z = lattice[i][j][k][z].z;
                    
                    //OP: E
                    E += getLocalEnergy(lattice[i][j][k][z], false);
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
                    
                    ClassicalSpin3D a1_x = getSpinInTheXDir(lattice[i][j][k][z]);
                    ClassicalSpin3D a1_y = getSpinInTheYDir(lattice[i][j][k][z]);
                    
                    PX += (X * a1_x.x + Y * a1_x.y + Z * a1_x.z);
                    PY += (X * a1_y.x + Y * a1_y.y + Z * a1_y.z);
                    
                }
            }
        }
    }
    //End OP
    
    OP.E = E; //Total E, E not per site
    
    OP.PX = PX;
    OP.PY = PY;
    OP.Q = PX - PY;
    
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
    
    if((MD.numConfigsDone % CORRTOUPDATES == 0)&&(MD.numConfigsDone > 0)){
        updateCorrelationFunction();
    }
    
    if((MD.numConfigsDone % FTTOUPDATES == 0)&&(MD.numConfigsDone > 0)){
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
        E[count] = getLocalEnergy(lattice[i][j][k][z], false);
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
        
        
        ClassicalSpin3D a1_x = getSpinInTheXDir(lattice[i][j][k][z]);
        ClassicalSpin3D a1_y = getSpinInTheYDir(lattice[i][j][k][z]);
        
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
    const double numSites = MCP.cellsA *MCP.cellsB * MCP.cellsC * NUMSUBLATTICES * 1.0;
    double kvector_x = 0;
    double kvector_y = 0;
    double rvector_x = 0;
    double rvector_y = 0;
    double kr = 0;
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
                    recipLattice[k_i][k_j][0][0].ReXComponent += cos(kr) * lattice[i][j][0][0].x / (numSites*1.0);//A,B,s,C,x/y/z/ix/iy/iz
                    recipLattice[k_i][k_j][0][0].ImXComponent += sin(kr) * lattice[i][j][0][0].x / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ReYComponent += cos(kr) * lattice[i][j][0][0].y / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ImYComponent += sin(kr) * lattice[i][j][0][0].y / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ReZComponent += cos(kr) * lattice[i][j][0][0].z / (numSites*1.0);
                    recipLattice[k_i][k_j][0][0].ImZComponent += sin(kr) * lattice[i][j][0][0].z / (numSites*1.0);
                    
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
    
    //E
    LC.sumE[0] += OP.E;
    LC.sumE[1] += OP.E*OP.E;
    LC.sumE[2] += (OP.E*OP.E) * OP.E;
    LC.sumE[3] += (OP.E*OP.E) * (OP.E*OP.E);
    
    //PX
    LC.sumPX[0] += OP.PX;
    LC.sumPX[1] += OP.PX*OP.PX;
    LC.sumPX[2] += (OP.PX*OP.PX) * OP.PX;
    LC.sumPX[3] += (OP.PX*OP.PX) * (OP.PX*OP.PX);
    
    //PY
    LC.sumPY[0] += OP.PY;
    LC.sumPY[1] += OP.PY*OP.PY;
    LC.sumPY[2] += (OP.PY*OP.PY) * OP.PY;
    LC.sumPY[3] += (OP.PY*OP.PY) * (OP.PY*OP.PY);
    
    //Q
    LC.sumQ[0] += OP.Q;
    LC.sumQ[1] += OP.Q*OP.Q;
    LC.sumQ[2] += (OP.Q*OP.Q) * OP.Q;
    LC.sumQ[3] += (OP.Q*OP.Q) * (OP.Q*OP.Q);
    
    
    for(int i = 0; i < NUMVECTORS; i++){
        for(int j = 0; j < TYPESOP; j++){
            
            double temp2 = OP.charOP[i][j]*OP.charOP[i][j];
            
            LC.sumOP[i][j][0] += OP.charOP[i][j];
            LC.sumOP[i][j][1] += temp2;
            LC.sumOP[i][j][2] += OP.charOP[i][j]*temp2;
            LC.sumOP[i][j][3] += temp2*temp2;
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
//Core Functions
/******************************************************************************/

const ClassicalSpin3D& MonteCarlo::getSpinInTheXDir(const ClassicalSpin3D& a) const{
    int aX = a.xPos;
    int aS = a.sPos;
    int aY = a.yPos;
    int aZ = a.zPos;
    int ns = 0;//next spin's sublattice
    
    if(aS == 0){
        ns = 1;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 2){
        ns = 3;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 1){
        ns = 0;
        if(aX == MCP.cellsA - 1){
            return lattice[0][aY][ns][aZ];
        }else{
            return lattice[aX + 1][aY][ns][aZ];
        }
    }else if(aS == 3){
        ns = 2;
        if(aX == MCP.cellsA - 1){
            return lattice[0][aY][ns][aZ];
        }else{
            return lattice[aX + 1][aY][ns][aZ];
        }
    }else{
        std::cerr << "Wrong sublattice" << std::endl;
        exit(1);
        return a;
    }
}

const ClassicalSpin3D& MonteCarlo::getSpinInTheMinusXDir(const ClassicalSpin3D& a) const{
    int aX = a.xPos;
    int aS = a.sPos;
    int aY = a.yPos;
    int aZ = a.zPos;
    int ns = 0;//next spin's sublattice
    
    if(aS == 1){
        ns = 0;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 3){
        ns = 2;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 0){
        ns = 1;
        if(aX == 0){
            return lattice[MCP.cellsA - 1][aY][ns][aZ];
        }else{
            return lattice[aX - 1][aY][ns][aZ];
        }
    }else if(aS == 2){
        ns = 3;
        if(aX == 0){
            return lattice[MCP.cellsA - 1][aY][ns][aZ];
        }else{
            return lattice[aX - 1][aY][ns][aZ];
        }
    }else{
        std::cerr << "Wrong sublattice" << std::endl;
        exit(1);
        return a;
    }
}



const ClassicalSpin3D& MonteCarlo::getSpinInTheMinusYDir(const ClassicalSpin3D& a) const{
    int aY = a.yPos;
    int aS = a.sPos;
    int aX = a.xPos;
    int aZ = a.zPos;
    int ns = 0;//next spin's sublattice
    
    if(aS == 2){
        ns = 0;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 3){
        ns = 1;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 0){
        ns = 2;
        if(aY == 0){
            return lattice[aX][MCP.cellsB - 1][ns][aZ];
        }else{
            return lattice[aX][aY - 1][ns][aZ];
        }
    }else if(aS == 1){
        ns = 3;
        if(aY == 0){
            return lattice[aX][MCP.cellsB - 1][ns][aZ];
        }else{
            return lattice[aX][aY - 1][ns][aZ];
        }
    }else{
        std::cerr << "Wrong sublattice" << std::endl;
        exit(1);
        return a;
    }
    
}

const ClassicalSpin3D& MonteCarlo::getSpinInTheYDir(const ClassicalSpin3D& a) const{
    int aY = a.yPos;
    int aS = a.sPos;
    int aX = a.xPos;
    int aZ = a.zPos;
    int ns = 0;//next spin's sublattice
    
    if(aS == 0){
        ns = 2;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 1){
        ns = 3;
        return lattice[aX][aY][ns][aZ];
    }else if(aS == 2){
        ns = 0;
        if(aY == MCP.cellsB - 1){
            return lattice[aX][0][ns][aZ];
        }else{
            return lattice[aX][aY + 1][ns][aZ];
        }
    }else if(aS == 3){
        ns = 1;
        if(aY == MCP.cellsB - 1){
            return lattice[aX][0][ns][aZ];
        }else{
            return lattice[aX][aY + 1][ns][aZ];
        }
    }else{
        std::cerr << "Wrong sublattice" << std::endl;
        exit(1);
        return a;
    }
    
}

const ClassicalSpin3D& MonteCarlo::getSpinInTheUpDir(const ClassicalSpin3D& a) const{
    int aY = a.yPos;
    int aS = a.sPos;
    int aX = a.xPos;
    int aZ = a.zPos;
    
    if(MCP.cellsC == 1){
        return a;
    }else{
        if(aZ == MCP.cellsC - 1){
            return lattice[aX][aY][a.sPos][0];
        }else{
            return lattice[aX][aY][a.sPos][aZ + 1];
        }
    }
    return a;
}

const ClassicalSpin3D& MonteCarlo::getSpinInTheDownDir(const ClassicalSpin3D& a) const{
    int aY = a.yPos;
    int aS = a.sPos;
    int aX = a.xPos;
    int aZ = a.zPos;
    
    if(MCP.cellsC == 1){
        return a;
    }else{
        if(aZ == 0){
            return lattice[aX][aY][a.sPos][MCP.cellsC - 1];
        }else{
            return lattice[aX][aY][a.sPos][aZ - 1];
        }
    }
    return a;
    
}

double MonteCarlo::getLocalEnergy(const ClassicalSpin3D& a, const bool fromSpinFlip) const{
    double totalInteraction = 0;
    
    ClassicalSpin3D a1_x = getSpinInTheXDir(a);
    ClassicalSpin3D a1_mx = getSpinInTheMinusXDir(a);
    ClassicalSpin3D a1_y = getSpinInTheYDir(a);
    ClassicalSpin3D a1_my = getSpinInTheMinusYDir(a);
    
    ClassicalSpin3D a2_xy =   getSpinInTheYDir(a1_x);
    ClassicalSpin3D a2_mxy =  getSpinInTheYDir(a1_mx);
    ClassicalSpin3D a2_mxmy = getSpinInTheMinusYDir(a1_mx);
    ClassicalSpin3D a2_xmy =  getSpinInTheMinusYDir(a1_x);
    
    ClassicalSpin3D a3_xx =   getSpinInTheXDir(a1_x);
    ClassicalSpin3D a3_mxmx = getSpinInTheMinusXDir(a1_mx);
    ClassicalSpin3D a3_yy =   getSpinInTheYDir(a1_y);
    ClassicalSpin3D a3_mymy = getSpinInTheMinusYDir(a1_my);
    
    double interactionH = MCP.j1 * (a.dotProd(a1_x) +
                                    a.dotProd(a1_mx) +
                                    a.dotProd(a1_y) +
                                    a.dotProd(a1_my));
    
    //Heisenberg interaction
    //2nd Neighbor interaction
    double interactionH2 = MCP.j2 * (a.dotProd(a2_xy) + a.dotProd(a2_mxy) +
                                     a.dotProd(a2_mxmy) + a.dotProd(a2_xmy));
    
    //3rd Neighbor interaction
    double interactionH3 = MCP.j3 * (a.dotProd(a3_xx) +
                                     a.dotProd(a3_mxmx) +
                                     a.dotProd(a3_yy) +
                                     a.dotProd(a3_mymy));
    
    double temp = (a.dotProd(a1_x) +
                   a.dotProd(a1_mx) +
                   a.dotProd(a1_y) +
                   a.dotProd(a1_my));
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
    
    for(int i = 0; i < MCP.cellsA; i++){
        for(int j = 0; j < MCP.cellsB; j++){
            for(int k = 0; k < NUMSUBLATTICES; k++){
                for(int z = 0; z < MCP.cellsC; z++)
                {
                    //getSpin(i, j, k).checkSpin();
                    initialEnergy = getLocalEnergy(lattice[i][j][k][z] , true);
                    lattice[i][j][k][z].flip(MD.range);
                    finalEnergy = getLocalEnergy(lattice[i][j][k][z], true);
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
                            lattice[i][j][k][z].reset();
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
                        lattice[i][j][k][z].reset();
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
     string temp = outputMag.substr(0,outputMag.size()-4);
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





