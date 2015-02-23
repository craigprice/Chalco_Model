#include "header.h"


int main(){
    cout<<"start!"<<endl;
    
    bool submitJobs = true;
    bool reThermalize = false;
    //reThermalize = true;
    
    paramSet KbTPar;    KbTPar.name = "KbT";
    KbTPar.min = 0.01;
    KbTPar.max = 0.01;
    KbTPar.inc = 0.08;  KbTPar.isRunParam = true;
    double estimatedTc = 0.25;
    int numSweepsToPerformTotal = 10 * 1000 * 1000;
    
    paramSet cellsPar;  cellsPar.name = "cells";
    cellsPar.min = 0;
    cellsPar.max = 0;
    cellsPar.inc = 1;   cellsPar.isRunParam = true;
    //0==36, 1==48, 2==60, 3==84,
    //4==96, 6==108, 7==120, 8==144, 9==168, 10==204
    
    paramSet cubePar;   cubePar.name = "cube";
    cubePar.min = 0.5;
    cubePar.max = 0.5;
    cubePar.inc = 0.1;  cubePar.isRunParam = false; cubePar.isHam = cubePar.isRunParam;
    
    paramSet bFieldPar;    bFieldPar.name = "bField";
    bFieldPar.min = 0;
    bFieldPar.max = 7.3;
    bFieldPar.inc = 0.1;
    bFieldPar.bFieldx = 0;
    bFieldPar.bFieldy = 0;
    bFieldPar.bFieldz = 1;  bFieldPar.isRunParam = false;
    
    paramSet phiPar;    phiPar.name = "phi";
    phiPar.min = 0;
    phiPar.max = 6.2;
    phiPar.inc = 0.1;   phiPar.isRunParam = false; phiPar.isHam = phiPar.isRunParam;
    /*
     phi = 0 => 0.8; N
     phi = 1 => 2.1; Z
     phi = 2 => 2.5; Z
     phi = 3 => 2.9; F
     phi = 4 => 4;   F
     phi = 5 => 5.2; S
     phi = 6 => 5.7; S
     phi = 7 => 6;   N
     */
    
    paramSet j2Par;    j2Par.name = "j2";
    j2Par.min = 0;
    j2Par.max = 1;
    j2Par.inc = 0.1;   j2Par.isRunParam = false; j2Par.isHam = j2Par.isRunParam;
    
    paramSet j3Par;    j3Par.name = "j3";
    j3Par.min = 0;
    j3Par.max = 1;
    j3Par.inc = 0.1;   j3Par.isRunParam = false; j3Par.isHam = j3Par.isRunParam;
    
    paramSet k2Par;    k2Par.name = "k2";
    k2Par.min = 0;
    k2Par.max = 1;
    k2Par.inc = 0.1;   k2Par.isRunParam = false; k2Par.isHam = k2Par.isRunParam;
    
    bool settingK2 = false;
    
    paramSet alphaPar;  alphaPar.name = "alpha";
    alphaPar.min = 0.3;
    alphaPar.max = 0.3;
    alphaPar.inc = 0.1; alphaPar.isRunParam = false; alphaPar.isHam = alphaPar.isRunParam;
    
    bool isDomainWalls = false;
    int  directionOfDomain = 0;
    bool isRotHam = false;
    
    paramSet phiK12Par; phiK12Par.name = "phiK12";
    phiK12Par.min = 0;
    phiK12Par.max = 1;
    phiK12Par.inc = 0.01; phiK12Par.isRunParam = true; phiK12Par.isHam = phiK12Par.isRunParam;
    
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
                double iphiIt = phiPar.min;
                do{
                    double ibfield = bFieldPar.min;
                    do{
                        double icube = cubePar.min;
                        do{
                            double ik2 = k2Par.min;
                            do{
                                double ij3 = j3Par.min;
                                do{
                                    double ij2 = j2Par.min;
                                    do{
                                        double ik12 = phiK12Par.min;
                                        do{
                                            double ialpha = alphaPar.min;
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
                                                
                                                double iphi = iphiIt;//0;
                                                /*
                                                 if(phiPar.isRunParam == true){
                                                 if(fabs(iphiIt - 0) < 1e-8){
                                                 iphi = 0.8;
                                                 }else if(fabs(iphiIt - 1) < 1e-8){
                                                 iphi = 2.1;
                                                 }else if(fabs(iphiIt - 2) < 1e-8){
                                                 iphi = 2.5;
                                                 }else if(fabs(iphiIt - 3) < 1e-8){
                                                 iphi = 2.9;
                                                 }else if(fabs(iphiIt - 4) < 1e-8){
                                                 iphi = 4;
                                                 }else if(fabs(iphiIt - 5) < 1e-8){
                                                 iphi = 5.2;
                                                 }else if(fabs(iphiIt - 6) < 1e-8){
                                                 iphi = 5.7;
                                                 }else if(fabs(iphiIt - 7) < 1e-8){
                                                 iphi = 6;
                                                 }else{
                                                 cerr << "Stop" << endl;
                                                 exit(1);
                                                 }
                                                 }
                                                 /**/
                                                
                                                if(alphaPar.isRunParam){
                                                    /*double alpha = 0;
                                                     if(counter == 0){
                                                     alpha = 0.8;
                                                     }else if(counter == 1){
                                                     alpha = 2.1;
                                                     }else if(counter == 2){
                                                     alpha = 2.5;
                                                     }else if(counter == 3){
                                                     alpha = 4;
                                                     }else if(counter == 4){
                                                     alpha = 4.7;
                                                     }else if(counter == 5){
                                                     alpha = 5.2;
                                                     }else if(counter == 6){
                                                     alpha = 5.7;
                                                     }else if(counter == 7){
                                                     alpha = 6;
                                                     }else{
                                                     cerr << "Stop" << endl;
                                                     exit(1);
                                                     }
                                                     //if(fabs(alpha - 4.7) <= 1e-7) continue;
                                                     //if(fabs(alpha - 4) <= 1e-7) continue;
                                                     //alpha = 0.33;
                                                     /*
                                                     if(
                                                     ((fabs(extMagField - 0.8)<0.0001)&&(cells == 60))||
                                                     ((fabs(extMagField - 6)<0.0001)&&(cells == 60))||
                                                     ((fabs(extMagField - 7)<0.0001)&&(cells == 60))
                                                     ){continue;}
                                                     */
                                                }
                                                
                                                
                                                int cellsA = icells;
                                                int cellsB = icells;
                                                int cellsC = 1;
                                                if(!isDomainWalls)          directionOfDomain   = 0;
                                                if(!alphaPar.isRunParam)    ialpha              = 0;
                                                if(!j2Par.isRunParam)       ij2                 = 0;
                                                if(!j3Par.isRunParam)       ij3                 = 0;
                                                if(!k2Par.isRunParam)       ik2                 = 0;
                                                if(!cubePar.isRunParam)     icube               = 0;
                                                if(!bFieldPar.isRunParam)
                                                {
                                                    bFieldPar.bFieldx                           = 0;
                                                    bFieldPar.bFieldy                           = 0;
                                                    bFieldPar.bFieldz                           = 0;
                                                    ibfield                                     = 0;
                                                }
                                                if(!phiPar.isRunParam)      iphi                = 0;
                                                if(!cellsPar.isRunParam)    icells              = 0;
                                                if(!KbTPar.isRunParam)      ikbt                = 0;
                                                
                                                if (fabs(ialpha) < 1e-8)     {alphaPar.isHam    = false;} else {alphaPar.isHam = true;}
                                                if (fabs(ij2) < 1e-8)        {j2Par.isHam       = false;} else {j2Par.isHam = true;}
                                                if (fabs(ij3) < 1e-8)        {j3Par.isHam       = false;} else {j3Par.isHam = true;}
                                                if (fabs(ik2) < 1e-8)        {k2Par.isHam       = false;} else {k2Par.isHam = true;}
                                                if (fabs(icube) < 1e-8)      {cubePar.isHam     = false;} else {cubePar.isHam = true;}
                                                if (fabs(ibfield) < 1e-8)    {bFieldPar.isHam   = false;} else {bFieldPar.isHam = true;}
                                                if (fabs(iphi) < 1e-8)       {phiPar.isHam      = false;} else {phiPar.isHam = true;}
                                                
                                                int precision  = 1000;
                                                ialpha      = floorf(ialpha * precision + 0.5) / precision + 0.0;//Get rid of -0.0
                                                ij2         = floorf(ij2 * precision + 0.5) / precision + 0.0;
                                                ij3         = floorf(ij3 * precision + 0.5) / precision + 0.0;
                                                ik2         = floorf(ik2 * precision + 0.5) / precision + 0.0;
                                                icube       = floorf(icube * precision + 0.5) / precision + 0.0;
                                                ibfield     = floorf(ibfield * precision + 0.5) / precision + 0.0;
                                                iphi        = floorf(iphi * precision + 0.5) / precision + 0.0;
                                                ikbt        = floorf(ikbt * precision + 0.5) / precision + 0.0;
                                                
                                                double k2 = ik2;
                                                if(settingK2) {k2 = -2*ij2 + 0.0;}
                                                
                                                double j1 = 0;
                                                double k1 = 0;
                                                
                                                if (phiK12Par.isRunParam) {
                                                    j1 = 0;
                                                    ij2 = 0;
                                                    ij3 = 0;
                                                    k1 = cos((2.0 * PI) * ik12);
                                                    k2 = sin((2.0 * PI) * ik12);
                                                }
                                                
                                                cout<<k1<< " "<<k2<<endl;
                                                
                                                k1          = floorf(k1 * precision + 0.5) / precision + 0.0;
                                                k2          = floorf(k2 * precision + 0.5) / precision + 0.0;
                                                cout<<k1<<" " <<k2<<endl;
                                                //cout<< "icells: "<< icells << " ikbt: " << ikbt << " iphi: " <<
                                                //iphi << " ibfield: " << ibfield << " icube: " << icube << " ij2: " <<
                                                //ij2 << " ij3: " << ij3 << " ik2: " << k2 << " ialpha: " << ialpha << endl;
                                                
                                                if (isContinue(iphi,ibfield,ikbt)) continue;
                                                
                                                std::stringstream ss;
                                                ss.str("");
                                                ss << "sim";
                                                ss << "_" << "isDom_" <<    isDomainWalls;
                                                ss << "_" << "DomDir_" <<   directionOfDomain;
                                                ss << "_" << "isRot_" <<    isRotHam;
                                                ss << "_" << "isA_" <<      alphaPar.isHam;
                                                ss << "_" << "a_" <<        ialpha;
                                                ss << "_" << "isJ23_" <<    (j2Par.isHam ||
                                                                             j3Par.isHam ||
                                                                             k2Par.isHam);
                                                ss << "_" << "isK12"<<      phiK12Par.isRunParam;
                                                ss << "_" << "J1_"<<        j1;
                                                ss << "_" << "K1_"<<        k1;
                                                ss << "_" << "J2_" <<       ij2;
                                                ss << "_" << "J3_" <<       ij3;
                                                ss << "_" << "K2_" <<       k2;
                                                ss << "_" << "isCub_" <<    cubePar.isHam;
                                                ss << "_" << "D_" <<        icube;
                                                ss << "_" << "isB_" <<      bFieldPar.isRunParam;
                                                ss << "_" << "bx_" <<       bFieldPar.bFieldx;
                                                ss << "_" << "by_" <<       bFieldPar.bFieldy;
                                                ss << "_" << "bz_" <<       bFieldPar.bFieldz;
                                                ss << "_" << "bMag_" <<     ibfield;
                                                ss << "_" << "isPhi_" <<    phiPar.isHam;
                                                ss << "_" << "phi_" <<      iphi;
                                                ss << "_" << "celA_" <<     icells;
                                                ss << "_" << "celB_" <<     icells;
                                                ss << "_" << "celC_" <<     1;
                                                ss << "_" << "KbT_" <<      ikbt;
                                                string bareFileName = ss.str();
                                                
                                                //cout<< fileName.str()<<endl;
                                                simParameters sp(bareFileName,
                                                                 isDomainWalls,
                                                                 directionOfDomain,
                                                                 isRotHam,
                                                                 alphaPar.isHam,
                                                                 ialpha,
                                                                 (j2Par.isHam || j3Par.isHam || k2Par.isHam),
                                                                 phiK12Par.isRunParam,
                                                                 j1,
                                                                 k1,
                                                                 ij2,
                                                                 ij3,
                                                                 k2,
                                                                 cubePar.isHam,
                                                                 icube,
                                                                 bFieldPar.isRunParam,
                                                                 bFieldPar.bFieldx,
                                                                 bFieldPar.bFieldy,
                                                                 bFieldPar.bFieldz,
                                                                 ibfield,
                                                                 phiPar.isHam,
                                                                 iphi,
                                                                 cellsA,
                                                                 cellsB,
                                                                 cellsC,
                                                                 ikbt);
                                                
                                                
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
                                                ss << "arguments = 27" << " " <<
                                                isDomainWalls << " " <<
                                                directionOfDomain << " " <<
                                                isRotHam << " " <<
                                                alphaPar.isHam << " " <<
                                                ialpha << " " <<
                                                (j2Par.isHam || j3Par.isHam || k2Par.isHam) << " " <<
                                                phiK12Par.isRunParam << " " <<
                                                j1 << " " <<
                                                k1 << " " <<
                                                ij2 << " " <<
                                                ij3 << " " <<
                                                k2 << " " <<
                                                cubePar.isHam << " " <<
                                                icube << " " <<
                                                bFieldPar.isRunParam << " " <<
                                                bFieldPar.bFieldx << " " <<
                                                bFieldPar.bFieldy << " " <<
                                                bFieldPar.bFieldz << " " <<
                                                ibfield << " " <<
                                                phiPar.isHam << " " <<
                                                iphi << " " <<
                                                cellsA << " " <<
                                                cellsB << " " <<
                                                cellsC << " " <<
                                                ikbt << " " <<
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
                                                "/afs/hep.wisc.edu/home/cprice/perkins/" <<
                                                folderName << "/simulation" << endl;
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
                                                ialpha += alphaPar.inc;
                                            } while (((ialpha - alphaPar.max) < 1e-8) && (alphaPar.isRunParam));
                                            ik12 += phiK12Par.inc;
                                        } while (((ik12 - phiK12Par.max) < 1e-8) && (phiK12Par.isRunParam));
                                        ij2 += j2Par.inc;
                                    } while (((ij2 - j2Par.max) < 1e-8) && (j2Par.isRunParam));
                                    ij3 += j3Par.inc;
                                } while (((ij3 - j3Par.max) < 1e-8) && (j3Par.isRunParam));
                                ik2 += k2Par.inc;
                            }while (((ik2 - k2Par.max) < 1e-8) && (k2Par.isRunParam));
                            icube += cubePar.inc;
                        } while (((icube - cubePar.max) < 1e-8) && (cubePar.isRunParam));
                        ibfield += bFieldPar.inc;
                    } while (((ibfield - bFieldPar.max) < 1e-8) && (bFieldPar.isRunParam));
                    iphiIt += phiPar.inc;
                } while (((iphiIt - phiPar.max) < 1e-8) && (phiPar.isRunParam));
                ikbt += KbTPar.inc;
            } while (((ikbt - KbTPar.max) < 1e-8) && (KbTPar.isRunParam));
            icellsIt -= cellsPar.inc;
        } while (((cellsPar.min - icellsIt) < 1e-8) && (cellsPar.isRunParam));
        
        
        cout<<"OK to Quit - Sleeping..."<<endl;
        sleep(5*60);
        cout<<"Don't Quit! - Processing..."<<endl;
        system("condor_release cprice");
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


