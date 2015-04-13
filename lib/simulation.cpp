/**
 simulation.cpp
 This file is the entrance portal to running the simulation. In it we specify
 the input parameters to the simulation, create the MonteCarlo object that
 is really the lattice of spins, thermalize the system, perform enough MC
 updates to achieve our statistical uncertainty requirements, and then print
 out the result of our simulation to a text file.
 
 @author Craig Price  (ccp134@psu.edu)
 @version 3.0 2015/04/12
 */

#include <sstream>
#include <ctime>
#include "MonteCarlo.h"
#include "MonteCarlo.cpp"

int main(int argc, char*argv[]){
    if(argc != 23 && argc != 4){
        std::cerr << "Need correct inputs" << std::endl;
        return 1;
    }
    
    MonteCarlo::MCParameters input_parameters;
    
    MonteCarlo* MC;// = new MonteCarlo();
    
    //delete MC;
    if(argc == 23){
        
        if((atoi(argv[2]) > 255)||//cellsA
           (atoi(argv[3]) > 255)||//cellsB
           (atoi(argv[4]) > 255)){//cellsC
            std::cout << "Lattice Index > 255. Make larger sublattice." << std::endl;
            exit(1);
        }
        
        input_parameters.cellsA      =           atoi(argv[2]);
        input_parameters.cellsB      =           atoi(argv[3]);
        input_parameters.cellsC      =           atoi(argv[4]);
        
        input_parameters.KbT         =           atof(argv[5]);
        
        input_parameters.estimatedTc =           atof(argv[6]);
        input_parameters.numSweepsToPerformTotal = atoi(argv[7]);
        input_parameters.path        =           argv[8];
        
        
        if (atoi(argv[9]) == -1){
            //Set the random seed.
            input_parameters.monte_carlo_seed = (unsigned long) time(0);
        }else{
            //Set the seed from the command line
            input_parameters.monte_carlo_seed = atoi(argv[9]);
        }
        
        input_parameters.monte_carlo_seed = 10;
        srand48( input_parameters.monte_carlo_seed );
        
        // B field Not necessarily normalized
        input_parameters.isBField    = (bool)   (atoi(argv[10]));
        input_parameters.bField_x    =           atof(argv[11]);
        input_parameters.bField_y    =           atof(argv[12]);
        input_parameters.bField_z    =           atof(argv[13]);
        input_parameters.bFieldMag   =           atof(argv[14]);
        
        input_parameters.j1     =           atof(argv[15]);
        input_parameters.j2     =           atof(argv[16]);
        input_parameters.j3     =           atof(argv[17]);
        input_parameters.k1     =           atof(argv[18]);
        input_parameters.k2     =           atof(argv[19]);
        input_parameters.k3     =           atof(argv[20]);
        
        input_parameters.cubicD =           atof(argv[21]);
        input_parameters.param1 =           atof(argv[22]);
        
        
        MC = new MonteCarlo(input_parameters);
        
    }else if(argc == 4){
        bool reThermalize = (bool) (atoi(argv[2]));
        std::string filename = argv[3];
        MC = new MonteCarlo(reThermalize, filename);
    }else{
        std::cerr<<"inputs"<<std::endl;
        exit(3161990);
    }
    double durationOfSingleRun = 20/60.0;
    int minSweepsToThermalize = 100 * 1000;
    int numSweepsBetwnConfigs = 6;
    bool outOfTime = false;
    /*
    //int numSweeps             =  500 * 1000;
    clock_t c0, c1;
    MC -> metropolisSweep();
    MC -> updateOrderParameters();
    MC -> updateRunningSumOP();
    
    MC -> metropolisSweep();
    MC -> updateOrderParameters();
    MC -> updateRunningSumOP();
    
    c0 = clock();
    for (int i = 0; i < 20*1000; i++){
        //MC -> metropolisSweep();
        //MC -> updateRunningSumOP();
    }
    
    
    for (int i = 0; i < 20*1000; i++){
        MC -> metropolisSweep();
        MC -> updateCorrelationFunction();
        //MC -> updateFourierTransformOnRecipLattice();
    }
    */
    /*
    double temp = 100;
    int count = 0;
    while(temp > 0.01){
        temp = temp * 0.99;
        count = count + 1;
        cout << count << " " << temp << endl;
    }
     */
    /*
    c1 = clock();
    std::cout<<"Clocks: " << (c1-c0)/(1000*1000.0) << std::endl;
    MC -> finalPrint();
    exit(0);
    */
    //
    //Runs the lattice without taking data at these temperatures to get the
    //spins pointed in approximately the correct direction.
    std::cout << "Thermalizing" << std::endl;
    if(MC -> MD.numSweepsUsedToThermalize < (int) (minSweepsToThermalize / 2.0)){
        outOfTime = MC -> thermalize(3.0 * MC -> MCP.estimatedTc,
                                     (int) (minSweepsToThermalize / 2.0),
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
