/**
 simulation.cpp
 This file is the entrance portal to running the simulation. In it we specify
 the input parameters to the simulation, create the MonteCarlo object that
 is really the lattice of spins, thermalize the system, perform enough MC
 updates to achieve our statistical uncertainty requirements, and then print
 out the result of our simulation to a text file.
 
 @author Craig Price (ccp134@psu.edu)
 @version 1.0 2015/03/04
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "ClassicalSpin3D.h"
#include "MonteCarlo.h"
#include "ClassicalSpin3D.cpp"
#include "MonteCarlo.cpp"

int main(int argc, char*argv[]){
    if(argc != 22 && argc != 4){
        std::cerr << "Need correct inputs" << std::endl;
        return 1;
    }
    
    MonteCarlo::MCParameters input_parameters;
    
    MonteCarlo* MC;// = new MonteCarlo();
    //delete MC;
    if(argc == 22){
        input_parameters.j1     =           atof(argv[2]);
        input_parameters.j2     =           atof(argv[3]);
        input_parameters.j3     =           atof(argv[4]);
        input_parameters.k1     =           atof(argv[5]);
        input_parameters.k2     =           atof(argv[6]);
        input_parameters.k3     =           atof(argv[7]);
        
        input_parameters.cubicD =           atof(argv[8]);
        
        // B field Not necessarily normalized
        input_parameters.isBField    = (bool)   (atoi(argv[9]));
        input_parameters.bField_x    =           atof(argv[10]);
        input_parameters.bField_y    =           atof(argv[11]);
        input_parameters.bField_z    =           atof(argv[12]);
        input_parameters.bFieldMag   =           atof(argv[13]);
        
        input_parameters.cellsA      =           atoi(argv[14]);
        input_parameters.cellsB      =           atoi(argv[15]);
        input_parameters.cellsC      =           atoi(argv[16]);
        
        input_parameters.KbT         =           atof(argv[17]);
        
        input_parameters.estimatedTc =           atof(argv[18]);
        input_parameters.numSweepsToPerformTotal = atoi(argv[19]);
        input_parameters.path        =           argv[20];
        
        if (atoi(argv[21]) == -1){
            //Set the random seed.
            input_parameters.monte_carlo_seed = (unsigned long) time(0);
        }else{
            //Set the seed from the command line
            input_parameters.monte_carlo_seed = atoi(argv[21]);
        }
        
        srand48( input_parameters.monte_carlo_seed );
        
        MC = new MonteCarlo(input_parameters);
        
    }else if(argc == 4){
        bool reThermalize = (bool) (atoi(argv[2]));
        std::string filename = argv[3];
        MC = new MonteCarlo(reThermalize, filename);
    }else{
        std::cerr<<"inputs"<<std::endl;
        exit(3161990);
    }
    //int numSweeps             =  500 * 1000;
    double durationOfSingleRun = 1/60.0;
    int minSweepsToThermalize = 100 * 1000;
    int numSweepsBetwnConfigs = 6;
    bool outOfTime = false;
    //MC -> outputMag = "out.txt";
    //MC -> metropolisSweep();
    //MC -> updateFourierTransformOnRecipLattice();
    
    //
    //Runs the lattice without taking data at these temperatures to get the
    //spins pointed in approximately the correct direction.
    std::cout << "Thermalizing" << std::endl;
    if(MC -> MD.numSweepsUsedToThermalize < (int) (minSweepsToThermalize / 2)){
        outOfTime = MC -> thermalize(3.0 * MC -> MCP.estimatedTc,
                                     (int) (minSweepsToThermalize / 2),
                                     durationOfSingleRun);
    }
    if((!outOfTime)&&
       (MC -> MD.numSweepsUsedToThermalize < minSweepsToThermalize)){
        outOfTime = MC -> thermalize(MC -> MCP.KbT,
                                     minSweepsToThermalize,
                                     durationOfSingleRun);
    }
    
    while (!(MC -> MD.isThermal) && !outOfTime) {
        outOfTime = MC -> isThermalized(durationOfSingleRun);
    }
    
    std::cout << "Sweeping and summing over configurations" << std::endl;
    if(!outOfTime && MC -> MD.isThermal){
        //To run a snapshot, change numSweeps, sweep().
        //MC -> sweepSnapshot(numSweeps, 6, 9);//This hasn't been edited - check function
        
        MC -> sweep(numSweepsBetwnConfigs,
                    durationOfSingleRun);
    }
    
    MC -> MD.totalTimeRunning = (clock() - MC -> MD.timeOfInitialization);
    double hours_elapsed = MC -> MD.totalTimeRunning * (1.0/(CLOCKS_PER_SEC*1.0)) / (60.0 * 60.0);
    std::cout << hours_elapsed << std::endl;
    MC->finalPrint();
    std::cout << "done!" << std::endl;
    return 0;
}
