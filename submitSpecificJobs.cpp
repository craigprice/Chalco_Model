#include "header.h"

simParameters setCells(int icellsIt, simParameters sp)
{
    int icells = 0;
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
    
    sp.setCel(icells,icells,1);
    return sp;
}

void setKbT(double startRange, double endRange, double incRange,
            simParameters sp, vector<simParameters>& v)
{
    for (double i = startRange; i < endRange; i+=incRange) {
        sp.setKbT(i);
        v.push_back(sp);
    }
}

int main(){
    cout<<"start!"<<endl;
    
    
    vector<simParameters> jobsToSubmit;
    jobsToSubmit.clear();
    
    simParameters sp;
    sp.setDom(false, 0);
    sp.setRot(false);
    sp.setA(false, 0);
    sp.setJ23(false, 0,0,0);
    sp.setCub(false, 0);
    
    for(int i = 3; i>= 0; i--){
        sp.setPhi(true, 2.1);
        sp.setB(true, 1,1,1,0.05);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 2.1);
        sp.setB(true, 1,1,1,0.5);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 2.1);
        sp.setB(true, 1,1,1,1.25);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 2.1);
        sp.setB(true, 1,1,1,1.61);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 2.1);
        sp.setB(true, 1,1,1,1.75);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        
        
        sp.setPhi(true, 5.7);
        sp.setB(true, 1,1,1,0.05);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 5.7);
        sp.setB(true, 1,1,1,0.5);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 5.7);
        sp.setB(true, 1,1,1,2);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 5.7);
        sp.setB(true, 1,1,1,2.25);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
        
        sp.setPhi(true, 5.7);
        sp.setB(true, 1,1,1,2.5);
        sp = setCells(i,sp);
        setKbT(0.04,0.25,0.04,sp,jobsToSubmit);
    }
    
    
    bool submitJobs = false;
    bool reThermalize = false;
    double estimatedTc = 0.25;
    int numSweepsToPerformTotal = 10 * 1000 * 1000;
    
    //reThermalize = true;
    
    
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
        
        for (int i = 0; i < jobsToSubmit.size(); i++)
        {
            simParameters sJ;
            sJ = jobsToSubmit[i];
            
            
            
            
            
            
            std::stringstream ss;
            ss.str("");
            ss << "sim";
            ss << "_" << "isDom_" <<    sJ.isDomainWalls;
            ss << "_" << "DomDir_" <<   sJ.directionOfDomain;
            ss << "_" << "isRot_" <<    sJ.isRotHam;
            ss << "_" << "isA_" <<      sJ.isA;
            ss << "_" << "a_" <<        sJ.a;
            ss << "_" << "isJ23_" <<    sJ.isJ23;
            ss << "_" << "J2_" <<       sJ.J2;
            ss << "_" << "J3_" <<       sJ.J3;
            ss << "_" << "K2_" <<       sJ.J3;
            ss << "_" << "isCub_" <<    sJ.isCub;
            ss << "_" << "D_" <<        sJ.D;
            ss << "_" << "isB_" <<      sJ.isB;
            ss << "_" << "bx_" <<       sJ.bx;
            ss << "_" << "by_" <<       sJ.by;
            ss << "_" << "bz_" <<       sJ.bz;
            ss << "_" << "bMag_" <<     sJ.bMag;
            ss << "_" << "isPhi_" <<    sJ.isPhi;
            ss << "_" << "phi_" <<      sJ.phi;
            ss << "_" << "celA_" <<     sJ.celA;
            ss << "_" << "celB_" <<     sJ.celB;
            ss << "_" << "celC_" <<     sJ.celC;
            ss << "_" << "KbT_" <<      sJ.KbT;
            string bareFileName = ss.str();
            
            //cout<< fileName.str()<<endl;
            simParameters sp(bareFileName,
                             sJ.isDomainWalls,
                             sJ.directionOfDomain,
                             sJ.isRotHam,
                             sJ.isA,
                             sJ.a,
                             sJ.isJ23,
                             sJ.J2,
                             sJ.J3,
                             sJ.J3,
                             sJ.isCub,
                             sJ.D,
                             sJ.isB,
                             sJ.bx,
                             sJ.by,
                             sJ.bz,
                             sJ.bMag,
                             sJ.isPhi,
                             sJ.phi,
                             sJ.celA,
                             sJ.celB,
                             sJ.celC,
                             sJ.KbT);
            
            
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
            sJ.isDomainWalls << " " <<
            sJ.directionOfDomain << " " <<
            sJ.isRotHam << " " <<
            sJ.isA << " " <<
            sJ.a << " " <<
            sJ.isJ23 << " " <<
            sJ.J2 << " " <<
            sJ.J3 << " " <<
            sJ.J3 << " " <<
            sJ.isCub << " " <<
            sJ.D << " " <<
            sJ.isB << " " <<
            sJ.bx << " " <<
            sJ.by << " " <<
            sJ.bz << " " <<
            sJ.bMag << " " <<
            sJ.isPhi << " " <<
            sJ.phi << " " <<
            sJ.celA << " " <<
            sJ.celB << " " <<
            sJ.celC << " " <<
            sJ.KbT << " " <<
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
            cout<<" "<< argumentsString<<endl;
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
        }
        
        
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
            cout<<ss.str()<<endl;
        }
        tempFile.close();
        system("rm -f _tempSubmit.txt");
        system("rm -f sim_isDom_*.txt");
        sleep(2);
    }
    
    cout<<"done!"<<endl;
    return 0;
    
}


