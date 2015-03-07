#include "header.h"


int main(){
    cout<<"start!"<<endl;
    
    bool submitJobs = true;
    bool reThermalize = false;
    //reThermalize = true;
    
    paramSet KbTPar;    KbTPar.name = "KbT";
    KbTPar.min = 0.01;
    KbTPar.max = 1;
    KbTPar.inc = 0.5;  KbTPar.isRunParam = true;
    double estimatedTc = 0.25;
    int numSweepsToPerformTotal = 10 * 1000 * 1000;
    
    paramSet cellsPar;  cellsPar.name = "cells";
    cellsPar.min = 0;
    cellsPar.max = 3;
    cellsPar.inc = 1;   cellsPar.isRunParam = true;
    //0==36, 1==48, 2==60, 3==84,
    //4==96, 6==108, 7==120, 8==144, 9==168, 10==204
    
    paramSet cubePar;   cubePar.name = "cube";
    cubePar.min = 0.1;
    cubePar.max = 0.1;
    cubePar.inc = 0.1;  cubePar.isRunParam = true; cubePar.isHam = cubePar.isRunParam;
    
    paramSet bFieldPar;    bFieldPar.name = "bField";
    bFieldPar.min = 0;
    bFieldPar.max = 7.3;
    bFieldPar.inc = 0.1;
    bFieldPar.bFieldx = 0;
    bFieldPar.bFieldy = 0;
    bFieldPar.bFieldz = 1;  bFieldPar.isRunParam = false;
    
    
    paramSet j2Par;    j2Par.name = "j2";
    j2Par.min = 3;
    j2Par.max = 3;
    j2Par.inc = 0.1;   j2Par.isRunParam = true; j2Par.isHam = j2Par.isRunParam;
    
    paramSet j3Par;    j3Par.name = "j3";
    j3Par.min = 0;
    j3Par.max = 1;
    j3Par.inc = 0.1;   j3Par.isRunParam = false; j3Par.isHam = j3Par.isRunParam;
    
    paramSet k1Par;    k1Par.name = "k1";
    k1Par.min = 0;
    k1Par.max = 1;
    k1Par.inc = 0.1;   k1Par.isRunParam = false; k1Par.isHam = k1Par.isRunParam;
    
    paramSet k2Par;    k2Par.name = "k2";
    k2Par.min = 0;
    k2Par.max = 1;
    k2Par.inc = 0.1;   k2Par.isRunParam = false; k2Par.isHam = k2Par.isRunParam;
    
    paramSet k3Par;    k3Par.name = "k3";
    k3Par.min = 0;
    k3Par.max = 1;
    k3Par.inc = 0.1;   k3Par.isRunParam = false; k3Par.isHam = k3Par.isRunParam;
    
    //Begin Common Strings for Condor File//////////////////////////////////////
    stringstream ss;
    ss.str("");
    ss << "requirements = "<<
    "(TARGET.OSglibc_major == 2 && TARGET.OSglibc_minor >= 9)&&"<<
    //"(IS_GLIDEIN)&&"<<
    "(ARCH == \"INTEL\" || ARCH == \"X86_64\")";
    //    "(Machine != \"glow-c070.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c071.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c072.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c073.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c074.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c075.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c076.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c077.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c080.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c081.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c082.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c083.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c084.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c085.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c086.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c087.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c089.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c093.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c094.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c141.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c150.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c170.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c186.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c188.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c242.cs.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c088.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"g20n02.hep.wisc.edu\")&&"<<
    //    "(Machine != \"glow-c078.medphysics.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-01.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-02.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-03.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-04.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-05.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-06.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-07.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-08.ep.wisc.edu\")&&"<<
    //    "(Machine != \"cnerg-09.ep.wisc.edu\")";
    string requirementString = ss.str();
    
    const int numCondorCommands = 10;
    string condorCommands[numCondorCommands] = {
        "executable = simulation",
        "universe = vanilla",
        "when_to_transfer_output = ON_EXIT",
        "+WantGlidein = true",
        "+estimated_run_hours = 0.25",
        "should_transfer_files = YES",
        "getenv = True",
        // "nice_user = True",
        "notify_user = ccprice@uwalumni.com",
        "notification = Error",
        "priority = 0"
    };
    //End Common Strings for Condor File////////////////////////////////////////
    
    char pathTocwd [512];
    getcwd(pathTocwd, 512);
    string folderName = pathTocwd;
    folderName = folderName.substr(folderName.find_last_of("/")+1);
    
    while(true){
        //Check if sim file alrady exists.
        vector<simParameters> allSimFiles;
        allSimFiles.clear();
        getAllSimFiles(allSimFiles, "dataFiles/");
        double icellsIt = cellsPar.max;
        do{
            double ikbt = KbTPar.min;
            do{
                double ij2 = j2Par.min;
                do{
                    double ij3 = cubePar.min;
                    do{
                        double ik1 = k1Par.min;
                        do{
                            double ik2 = k2Par.min;
                            do{
                                double ik3 = k3Par.min;
                                do{
                                    double icube = cubePar.min;
                                    do{
                                        double ibfield = bFieldPar.min;
                                        do{
                                            
                                            int icells = 0;
                                            if(cellsPar.isRunParam == true){
                                                if(fabs(icellsIt - 0) < 1e-8){
                                                    icells = 36;
                                                }else if(fabs(icellsIt - 1) < 1e-8){
                                                    icells = 48;
                                                }else if(fabs(icellsIt - 2) < 1e-8){
                                                    icells = 60;
                                                }else if(fabs(icellsIt - 3) < 1e-8){
                                                    icells = 84;
                                                }else if(fabs(icellsIt - 4) < 1e-8){
                                                    icells = 96;
                                                }else if(fabs(icellsIt - 5) < 1e-8){
                                                    icells = 108;
                                                }else if(fabs(icellsIt - 6) < 1e-8){
                                                    icells = 120;
                                                }else if(fabs(icellsIt - 7) < 1e-8){
                                                    icells = 144;
                                                }else if(fabs(icellsIt - 8) < 1e-8){
                                                    icells = 168;
                                                }else if(fabs(icellsIt - 9) < 1e-8){
                                                    icells = 204;
                                                }else{
                                                    cerr << "Stop" << endl;
                                                    exit(1);
                                                }
                                            }
                                            
                                            
                                            int cellsA = icells;
                                            int cellsB = icells;
                                            int cellsC = 1;
                                            if(!j2Par.isRunParam)       ij2                 = 0;
                                            if(!j3Par.isRunParam)       ij3                 = 0;
                                            if(!k1Par.isRunParam)       ik1                 = 0;
                                            if(!k2Par.isRunParam)       ik2                 = 0;
                                            if(!k3Par.isRunParam)       ik3                 = 0;
                                            if(!cubePar.isRunParam)     icube               = 0;
                                            if(!bFieldPar.isRunParam)
                                            {
                                                bFieldPar.bFieldx                           = 0;
                                                bFieldPar.bFieldy                           = 0;
                                                bFieldPar.bFieldz                           = 0;
                                                ibfield                                     = 0;
                                            }
                                            if(!cellsPar.isRunParam)    icells              = 0;
                                            if(!KbTPar.isRunParam)      ikbt                = 0;
                                            
                                            if (fabs(ij2) < 1e-8)        {j2Par.isHam       = false;} else {j2Par.isHam = true;}
                                            if (fabs(ij3) < 1e-8)        {j3Par.isHam       = false;} else {j3Par.isHam = true;}
                                            if (fabs(ik1) < 1e-8)        {k1Par.isHam       = false;} else {k1Par.isHam = true;}
                                            if (fabs(ik2) < 1e-8)        {k2Par.isHam       = false;} else {k2Par.isHam = true;}
                                            if (fabs(ik3) < 1e-8)        {k3Par.isHam       = false;} else {k3Par.isHam = true;}
                                            if (fabs(icube) < 1e-8)      {cubePar.isHam     = false;} else {cubePar.isHam = true;}
                                            if (fabs(ibfield) < 1e-8)    {bFieldPar.isHam   = false;} else {bFieldPar.isHam = true;}
                                            
                                            int precision  = 1000;
                                            ij2         = floorf(ij2 * precision + 0.5) / precision + 0.0;//Get rid of -0.0
                                            ij3         = floorf(ij3 * precision + 0.5) / precision + 0.0;
                                            ik1         = floorf(ik1 * precision + 0.5) / precision + 0.0;
                                            ik2         = floorf(ik2 * precision + 0.5) / precision + 0.0;
                                            ik3         = floorf(ik3 * precision + 0.5) / precision + 0.0;
                                            icube       = floorf(icube * precision + 0.5) / precision + 0.0;
                                            ibfield     = floorf(ibfield * precision + 0.5) / precision + 0.0;
                                            ikbt        = floorf(ikbt * precision + 0.5) / precision + 0.0;
                                            
                                            
                                            double j1 = 1;
                                            
                                            //cout<< "icells: "<< icells << " ikbt: " << ikbt << " iphi: " <<
                                            //iphi << " ibfield: " << ibfield << " icube: " << icube << " ij2: " <<
                                            //ij2 << " ij3: " << ij3 << " ik2: " << k2 << " ialpha: " << ialpha << endl;
                                            
                                            MCParameters input_parameters;
                                            input_parameters.j1     =           1;
                                            input_parameters.j2     =           ij2;
                                            input_parameters.j3     =           ij3;
                                            input_parameters.k1     =           ik1;
                                            input_parameters.k2     =           ik2;
                                            input_parameters.k3     =           ik3;
                                            
                                            input_parameters.cubicD =           icube;
                                            input_parameters.isBField    = (bool)   (bFieldPar.isHam);
                                            input_parameters.bField_x    =           bFieldPar.bFieldx;
                                            input_parameters.bField_y    =           bFieldPar.bFieldy;
                                            input_parameters.bField_z    =           bFieldPar.bFieldz;
                                            input_parameters.bFieldMag   =           ibfield;
                                            
                                            input_parameters.cellsA      =           icells;
                                            input_parameters.cellsB      =           icells;
                                            input_parameters.cellsC      =           1;
                                            
                                            input_parameters.KbT         =           ikbt;
                                            
                                            input_parameters.estimatedTc =           1;
                                            input_parameters.numSweepsToPerformTotal = 1;
                                            input_parameters.path    =           "";
                                            
                                            std::stringstream ss;
                                            ss.str("");
                                            ss << "sim";
                                            ss << "_" << "J1_"<<        input_parameters.j1;
                                            ss << "_" << "J2_" <<       input_parameters.j2;
                                            ss << "_" << "J3_" <<       input_parameters.j3;
                                            ss << "_" << "K1_"<<        input_parameters.k1;
                                            ss << "_" << "K2_" <<       input_parameters.k2;
                                            ss << "_" << "K3_"<<        input_parameters.k3;
                                            ss << "_" << "D_" <<        input_parameters.cubicD;
                                            ss << "_" << "isB_" <<      input_parameters.isBField;
                                            ss << "_" << "bx_" <<       input_parameters.bField_x;
                                            ss << "_" << "by_" <<       input_parameters.bField_y;
                                            ss << "_" << "bz_" <<       input_parameters.bField_z;
                                            ss << "_" << "bMag_" <<     input_parameters.bFieldMag;
                                            ss << "_" << "celA_" <<     input_parameters.cellsA;
                                            ss << "_" << "celB_" <<     input_parameters.cellsB;
                                            ss << "_" << "celC_" <<     input_parameters.cellsC;
                                            ss << "_" << "KbT_" <<      input_parameters.KbT;
                                            string bareFileName = ss.str();
                                            
                                            //cout<< fileName.str()<<endl;
                                            simParameters sp(bareFileName, input_parameters);
                                            
                                            
                                            if (isSimDone(allSimFiles,
                                                          sp,
                                                          reThermalize,
                                                          numSweepsToPerformTotal)) continue;
                                            
                                            cout<<"Sim File Not Done"<<endl;
                                            
                                            if (isRunningJob(bareFileName))  continue;
                                            
                                            cout<<"Sim File Not Running"<<endl;
                                            
                                            string reSubFileName = reSubmitFileName(allSimFiles,
                                                                                    sp,
                                                                                    reThermalize,
                                                                                    numSweepsToPerformTotal);
                                            
                                            //Begin Writing Arguments of Simulation/////
                                            ss.str("");
                                            ss << "arguments = 21" << " " <<
                                            input_parameters.j1 << " " <<
                                            input_parameters.j2 << " " <<
                                            input_parameters.j3 << " " <<
                                            input_parameters.k1 << " " <<
                                            input_parameters.k2 << " " <<
                                            input_parameters.k3 << " " <<
                                            input_parameters.cubicD << " " <<
                                            input_parameters.isBField << " " <<
                                            input_parameters.bField_x << " " <<
                                            input_parameters.bField_y << " " <<
                                            input_parameters.bField_z << " " <<
                                            input_parameters.bFieldMag << " " <<
                                            input_parameters.cellsA << " " <<
                                            input_parameters.cellsB << " " <<
                                            input_parameters.cellsC << " " <<
                                            input_parameters.KbT << " " <<
                                            estimatedTc << " " <<
                                            numSweepsToPerformTotal << " " <<
                                            //"/scratch/cprice/" << folderName << "/dataFiles/" << fileName.str() << ".txt ";
                                            bareFileName << ".txt ";
                                            string argumentsString = ss.str();
                                            
                                            if (reSubFileName.compare("0") != 0)
                                            {
                                                ss.str("");
                                                ss << "arguments = 4" << " "<<
                                                reThermalize << " " <<
                                                reSubFileName.substr(reSubFileName.find("sim_"));
                                                argumentsString = ss.str();
                                                cout<<"Sim File Not Beginning"<<endl;
                                            }
                                            //End Writing Arguments of Simulation////////
                                            
                                            char path [512];
                                            getcwd(path, 512);
                                            string pathName = path;
                                            
                                            stringstream runJob;
                                            runJob.str("");
                                            runJob << pathName <<
                                            "/condor/runJob/runJob-" <<
                                            bareFileName << ".condor";
                                            //cout << runJob.str().c_str() << endl;
                                            
                                            //Begin Writing Condor Submit File//////////
                                            ofstream myfile;
                                            myfile.open(runJob.str().c_str());
                                            for(int k = 0; k < numCondorCommands; k++){
                                                myfile << condorCommands[k] << endl;
                                            }
                                            myfile << "executable = " <<
                                            "/afs/hep.wisc.edu/home/nperkins/" <<
                                            folderName << "/bin/simulation" << endl;
                                            myfile << requirementString << endl;
                                            myfile << argumentsString << endl;
                                            //cout<<" "<< argumentsString<<endl;
                                            myfile << "log = "<< pathName << "/condor/log/log-";
                                            myfile << bareFileName << ".condor" << endl;
                                            myfile << "output = " << pathName << "/condor/output/output-";
                                            myfile << bareFileName << ".condor" << endl;
                                            myfile << "error = " << pathName << "/condor/error/error-";
                                            myfile << bareFileName << ".condor" << endl;
                                            
                                            //Transfer Files
                                            if(reSubFileName.compare("0") != 0)
                                            {
                                                myfile << "transfer_input_files = " <<
                                                pathName << "/" <<
                                                reSubFileName << endl;//.substr(1) << endl;
                                            }
                                            //cout<<folderName<<endl;
                                            //End Transfer Files
                                            
                                            myfile << "queue";
                                            myfile.close();
                                            //End Writing Condor Submit File////////////
                                            
                                            stringstream condor;
                                            condor.str("");
                                            condor<< "condor_submit " << runJob.str() << endl;
                                            //cout<<runJob.str()<<endl;
                                            if(submitJobs){system(condor.str().c_str());}
                                            
                                            //
                                            ibfield += bFieldPar.inc;
                                        } while (((ibfield - bFieldPar.max) < 1e-8) && (bFieldPar.isRunParam));
                                        icube += cubePar.inc;
                                    } while (((icube - cubePar.max) < 1e-8) && (cubePar.isRunParam));
                                    ik3 += k3Par.inc;
                                } while (((ik3 - k3Par.max) < 1e-8) && (k3Par.isRunParam));
                                ik2 += k2Par.inc;
                            }while (((ik2 - k2Par.max) < 1e-8) && (k2Par.isRunParam));
                            ik1 += k1Par.inc;
                        } while (((ik1 - k1Par.max) < 1e-8) && (k1Par.isRunParam));
                        ij3 += j3Par.inc;
                    } while (((ij3 - j3Par.max) < 1e-8) && (j3Par.isRunParam));
                    ij2 += j2Par.inc;
                } while (((ij2 - j2Par.max) < 1e-8) && (j2Par.isRunParam));
                ikbt += KbTPar.inc;
            } while (((ikbt - KbTPar.max) < 1e-8) && (KbTPar.isRunParam));
            icellsIt -= cellsPar.inc;
        } while (((cellsPar.min - icellsIt) < 1e-8) && (cellsPar.isRunParam));
        
        
        cout<<"OK to Quit - Sleeping..."<<endl;
        sleep(5*60);
        cout<<"Don't Quit! - Processing..."<<endl;
        system("condor_release nperkins");
        sleep(2);
        system("ls sim_isDom_*.txt > _tempSubmit.txt");
        sleep(2);
        ifstream tempFile;
        tempFile.open("_tempSubmit.txt");
        if (!tempFile.is_open() || !tempFile.good()){
            cout << "Can't open file" << endl;
            cout << "_tempSubmit.txt" << endl;
            exit(1);
        }
        string line = "";
        while(tempFile.good()){
            getline (tempFile,line);
            if(line.find("sim_isDom") == string::npos) continue;
            stringstream ss;
            ss.str("");
            //ss << "cat " << line << " >> dataFiles/" << line;
            ss << "mv -f " << line << " dataFiles/" << line;
            system(ss.str().c_str());
            //cout<<ss.str()<<endl;
        }
        tempFile.close();
        system("rm -f _tempSubmit.txt");
        system("rm -f sim_isDom_*.txt");
        sleep(2);
    }
    
    cout<<"done!"<<endl;
    return 0;
    
}


