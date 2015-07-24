
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <ostream>
#include <iomanip>
#include <map>
#include <math.h>
#include <cstdlib>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <unistd.h>
//const static double PI = 3.14159265358979323846;
#include "lib/MonteCarlo.h"
using namespace std;



class simParameters{
public:
    string  fileName;
    bool    isDomainWalls;
    int     directionOfDomain;
    bool    isRotHam;
    bool    isA;//alpha Hamiltonian
    double  a;
    bool    isJ23;//J23 Hamiltonian
    bool    isK12;
    double  J1;
    double  J2;
    double  J3;
    double  K1;
    double  K2;
    double  K3;
    double  param1;
    bool    isCub;
    double  D;
    bool    isB;
    double  bx;
    double  by;
    double  bz;
    double  bMag;
    bool    isPhi;
    double  phi;
    int     celA;
    int     celB;
    int     celC;
    double  KbT;
    simParameters(){
        fileName            = "";
        isDomainWalls       = false;
        directionOfDomain   = 0;
        isRotHam            = false;
        isA                 = false;
        a                   = 0;
        isJ23               = false;
        isK12               = false;
        J1                  = 1;
        J2                  = 0;
        J3                  = 0;
        K1                  = 0;
        K2                  = 0;
        K3                  = 0;
        param1              = 0;
        isCub               = false;
        D                   = 0;
        isB                 = false;
        bx                  = 0;
        by                  = 0;
        bz                  = 0;
        bMag                = 0;
        isPhi               = false;
        phi                 = 0;
        celA                = 0;
        celB                = 0;
        celC                = 0;
        KbT                 = 0;
    }
    simParameters(string fileName_, MCParameters params){
        fileName            = fileName_;
        isDomainWalls       = 0;
        directionOfDomain   = 0;
        isRotHam            = 0;
        isA                 = 0;
        a                   = 0;
        isJ23               = 0;
        isK12               = 0;
        J1                  = params.j1;
        J2                  = params.j2;
        J3                  = params.j3;
        K1                  = params.k1;
        K2                  = params.k2;
        K3                  = params.k3;
        isCub               = 0;
        D                   = params.cubicD;
        isB                 = params.isBField;
        bx                  = params.bField_x;
        by                  = params.bField_y;
        bz                  = params.bField_z;
        bMag                = params.bFieldMag;
        isPhi               = 0;
        phi                 = 0;
        celA                = params.cellsA;
        celB                = params.cellsB;
        celC                = params.cellsC;
        KbT                 = params.KbT;
    }
    simParameters(string fileName_,
                  bool isDomainWalls_,
                  int directionOfDomain_,
                  bool isRotHam_,
                  bool isA_,
                  double a_,
                  bool isJ23_,
                  bool isK12_,
                  double J1_,
                  double K1_,
                  double J2_,
                  double J3_,
                  double K2_,
                  bool isCub_,
                  double D_,
                  bool isB_,
                  double bx_,
                  double by_,
                  double bz_,
                  double bMag_,
                  bool isPhi_,
                  double phi_,
                  int celA_,
                  int celB_,
                  int celC_,
                  double KbT_){
        fileName            = fileName_;
        isDomainWalls       = isDomainWalls_;
        directionOfDomain   = directionOfDomain_;
        isRotHam            = isRotHam_;
        isA                 = isA_;
        a                   = a_;
        isJ23               = isJ23_;
        isK12               = isK12_;
        J1                  = J1_;
        K1                  = K1_;
        J2                  = J2_;
        J3                  = J3_;
        K2                  = K2_;
        isCub               = isCub_;
        D                   = D_;
        isB                 = isB_;
        bx                  = bx_;
        by                  = by_;
        bz                  = bz_;
        bMag                = bMag_;
        isPhi               = isPhi_;
        phi                 = phi_;
        celA                = celA_;
        celB                = celB_;
        celC                = celC_;
        KbT                 = KbT_;
    }
    void setDom(bool isDomainWalls_, int directionOfDomain_){
        isDomainWalls = isDomainWalls_;
        directionOfDomain = directionOfDomain_;
    }
    void setRot(bool isRotHam_){
        isRotHam = isRotHam_;
    }
    void setA(bool isA_, double a_){
        isA = isA_;
        a = a_;
    }
    void setJ23(bool isJ23_, double J2_, double J3_, double K2_){
        isJ23 = isJ23_;
        J2 = J2_;
        J3 = J3_;
        K2 = K2_;
        J1 = 1;
    }
    void setK12(bool isK12_, double J1_, double K1_, double J2_, double K2_, double J3_){
        isK12 = isK12_;
        J1 = J1_;
        K1 = K1_;
        K2 = K2_;
        J2 = J2_;
        J3 = J3_;
    }
    void setCub(bool isCub_, double D_){
        isCub_ = isCub_;
        D = D_;
    }
    void setB(bool isB_, double bx_, double by_, double bz_, double bMag_){
        isB = isB_;
        bx = bx_;
        by = by_;
        bz = bz_;
        bMag = bMag_;
    }
    void setPhi(bool isPhi_, double phi_){
        isPhi = isPhi_;
        phi = phi_;
    }
    void setCel(int celA_, int celB_, int celC_){
        celA = celA_;
        celB = celB_;
        celC = celC_;
    }
    void setKbT(double KbT_){
        KbT = KbT_;
    }
    void toString() const {
        cout << endl;
        cout << "simParameter Object, toString(): " << endl;
        cout << "fileName:" << endl;
        cout << fileName << endl;
        cout << "J1:..........................." << J1 << endl;
        cout << "J2:..........................." << J2 << endl;
        cout << "J3:..........................." << J3 << endl;
        cout << "K1:..........................." << K1 << endl;
        cout << "K2:..........................." << K2 << endl;
        cout << "K3:..........................." << K3 << endl;
        cout << "P1:..........................." << param1 << endl;
        cout << "D:............................" << D << endl;
        cout << "isB:.........................." << isB << endl;
        cout << "bx:..........................." << bx << endl;
        cout << "by:..........................." << by << endl;
        cout << "bz:..........................." << bz << endl;
        cout << "bMag:........................." << bMag << endl;
        cout << "celA:........................." << celA << endl;
        cout << "celB:........................." << celB << endl;
        cout << "celC:........................." << celC << endl;
        cout << "KbT:.........................." << KbT << endl;
    }
};

string extractParamValueFromTitle(string line, string flagString){
    if (line.find(flagString) == string::npos) {
        cout<< "Error in File Title Parsing"<<endl;
        exit(1);
    }
    
    int index = line.find(flagString);
    string trimFront = line.substr(index + flagString.size());
    int endIndex = 0;
    if(flagString.compare("KbT_") == 0)
    {
        endIndex = trimFront.find(".txt");
    }
    else
    {
        endIndex = trimFront.find("_");
    }
    string value = line.substr(index + flagString.size(),endIndex);
    cout<<line<<endl;
    cout<<flagString<<endl;
    cout<<value<<endl;
    return value;
}

void getAllSimFiles(vector<simParameters>& v, string path){
    stringstream ss;
    ss.str("");
    ss << "find " << path << " -name \"sim_*.txt\" > _temp.txt";
    system(ss.str().c_str());
    ifstream tempFile;
    tempFile.open("_temp.txt");
    if (!tempFile.is_open() || !tempFile.good()){
        cout << "Can't open temp file" << endl;
        exit(1);
    }
    
    string line = "";
    while(getline (tempFile,line)){
        cout<<line<<endl;
        simParameters sp;
        sp.fileName = line.erase(line.find_last_not_of(" \n\r\t")+1);
        sp.J1 =                 (double) atof(extractParamValueFromTitle(line, "J1_").c_str());
        sp.J2 =                 (double) atof(extractParamValueFromTitle(line, "J2_").c_str());
        sp.J3 =                 (double) atof(extractParamValueFromTitle(line, "J3_").c_str());
        sp.K1 =                 (double) atof(extractParamValueFromTitle(line, "K1_").c_str());
        sp.K2 =                 (double) atof(extractParamValueFromTitle(line, "K2_").c_str());
        sp.K3 =                 (double) atof(extractParamValueFromTitle(line, "K3_").c_str());
        sp.D =                  (double) atof(extractParamValueFromTitle(line, "D_").c_str());
        sp.isB =                (bool)   atoi(extractParamValueFromTitle(line, "isB_").c_str());
        sp.bx =                 (double) atof(extractParamValueFromTitle(line, "bx_").c_str());
        sp.by =                 (double) atof(extractParamValueFromTitle(line, "by_").c_str());
        sp.bz =                 (double) atof(extractParamValueFromTitle(line, "bz_").c_str());
        sp.bMag =               (double) atof(extractParamValueFromTitle(line, "bMag_").c_str());
        sp.celA =               (int)    atoi(extractParamValueFromTitle(line, "celA_").c_str());
        sp.celB =               (int)    atoi(extractParamValueFromTitle(line, "celB_").c_str());
        sp.celC =               (int)    atoi(extractParamValueFromTitle(line, "celC_").c_str());
        sp.KbT =                (double) atof(extractParamValueFromTitle(line, "KbT_").c_str());
        v.push_back(sp);
        cout<<line<<endl;
    }
    tempFile.close();
    system("rm -f _temp.txt");
}

bool isMatchParameters(const simParameters sp, const simParameters sp2) {
    //DO not check fileName identicalness because might containe folder names
    double eps = 1e-5;
    //sp2.toString();
    if(
       (fabs(sp.J1              -  sp2.J1) < eps)&&
       (fabs(sp.J2              -  sp2.J2) < eps)&&
       (fabs(sp.J3              -  sp2.J3) < eps)&&
       (fabs(sp.K1              -  sp2.K1) < eps)&&
       (fabs(sp.K2              -  sp2.K2) < eps)&&
       (fabs(sp.K3              -  sp2.K3) < eps)&&
       (fabs(sp.D               -  sp2.D) < eps)&&
       (sp.isB                  == sp2.isB)&&
       (fabs(sp.bx              -  sp2.bx) < eps)&&
       (fabs(sp.by              -  sp2.by) < eps)&&
       (fabs(sp.bz              -  sp2.bz) < eps)&&
       (fabs(sp.bMag            -  sp2.bMag) < eps)&&
       (sp.celA                 == sp2.celA)&&
       (sp.celB                 == sp2.celB)&&
       (sp.celC                 == sp2.celC)&&
       (fabs(sp.KbT             -  sp2.KbT) < eps)
       ){
        return true;
    }
    else
    {
        return false;
    }
}

bool isMatchAnyParameters(simParameters sp, const vector<simParameters> &v){
    //cout<<"DF"<<endl;
    for (int i = 0; i < v.size(); i++) {
        //if(v[i].bMag < 1)
        //{
            //v[i].toString();
            //sp.toString();
        //}
        if(isMatchParameters(sp, v[i])){
            return true;
        }
    }
    return false;
}

class parametersToSearch{
public:
    string paramName;
    string orderParameter;
    bool isFoundInFile;
    parametersToSearch(){
        paramName = "";
        isFoundInFile = false;
        orderParameter = "";
    }
    parametersToSearch(string str, string opStr, bool boo){
        paramName = str;
        orderParameter = opStr;
        isFoundInFile = boo;
    }
};

bool isAllParametersFound(const vector<parametersToSearch> &v){
    for(int i = 0; i < v.size(); i++){
        if (v[i].isFoundInFile == false)
        {
            return false;
        }
    }
    return true;
}

class extractedValueWithSimParam{
public:
    simParameters sp;
    string paramName;
    string originalLineInFile;
    string orderParameter;
    double valueDbl;
    int valueInt;
    bool valueBool;
    bool isExtracted;
    extractedValueWithSimParam(){
        paramName = "";
        originalLineInFile = "";
        valueDbl = 0;
        valueInt = 0;
        valueBool = 0;
        isExtracted = false;
    }
    void toString(){
        cout << endl;
        cout << "extractedValueWithSimParam, toString():" << endl;
        cout << "simParameters:" << endl;
        sp.toString();
        cout << "originalLineInFile:..........." << originalLineInFile << endl;
        cout << "paramName:...................." << paramName << endl;
        cout << "valueDbl:....................." << valueDbl << endl;
        cout << "valueInt:....................." << valueInt << endl;
        cout << "valueBool:...................." << valueBool << endl;
    }
    string toFileName(string xaxis){
        stringstream ss;
        ss.str("");
        if(xaxis.compare("KbT") == 0)
        {
            ss << "_" << "isDom_" <<    sp.isDomainWalls;
            ss << "_" << "DomDir_" <<   sp.directionOfDomain;
            ss << "_" << "isRot_" <<    sp.isRotHam;
            ss << "_" << "isA_" <<      sp.isA;
            ss << "_" << "a_" <<        sp.a;
            ss << "_" << "isJ23_" <<    sp.isJ23;
            ss << "_" << "isK12_" <<    sp.isK12;
            ss << "_" << "J1_" <<       sp.J1;
            ss << "_" << "K1_" <<       sp.K1;
            ss << "_" << "J2_" <<       sp.J2;
            ss << "_" << "J3_" <<       sp.J3;
            ss << "_" << "K2_" <<       sp.K2;
            ss << "_" << "isCub_" <<    sp.isCub;
            ss << "_" << "D_" <<        sp.D;
            ss << "_" << "isB_" <<      sp.isB;
            ss << "_" << "bx_" <<       sp.bx;
            ss << "_" << "by_" <<       sp.by;
            ss << "_" << "bz_" <<       sp.bz;
            ss << "_" << "bMag_" <<     sp.bMag;
            ss << "_" << "isPhi_" <<    sp.isPhi;
            ss << "_" << "phi_" <<      sp.phi;
            ss << "_" << "celA_" <<     sp.celA;
            ss << "_" << "celB_" <<     sp.celB;
            ss << "_" << "celC_" <<     sp.celC;
        }
        if(xaxis.compare("B") == 0)
        {
            ss << "_" << "isDom_" <<    sp.isDomainWalls;
            ss << "_" << "DomDir_" <<   sp.directionOfDomain;
            ss << "_" << "isRot_" <<    sp.isRotHam;
            ss << "_" << "isA_" <<      sp.isA;
            ss << "_" << "a_" <<        sp.a;
            ss << "_" << "isJ23_" <<    sp.isJ23;
            ss << "_" << "isK12_" <<    sp.isK12;
            ss << "_" << "J1_" <<       sp.J1;
            ss << "_" << "K1_" <<       sp.K1;
            ss << "_" << "J2_" <<       sp.J2;
            ss << "_" << "J3_" <<       sp.J3;
            ss << "_" << "K2_" <<       sp.K2;
            ss << "_" << "isCub_" <<    sp.isCub;
            ss << "_" << "D_" <<        sp.D;
            ss << "_" << "isB_" <<      sp.isB;
            ss << "_" << "bx_" <<       sp.bx;
            ss << "_" << "by_" <<       sp.by;
            ss << "_" << "bz_" <<       sp.bz;
            ss << "_" << "isPhi_" <<    sp.isPhi;
            ss << "_" << "phi_" <<      sp.phi;
            ss << "_" << "KbT_" <<      sp.KbT;
            ss << "_" << "celA_" <<     sp.celA;
            ss << "_" << "celB_" <<     sp.celB;
            ss << "_" << "celC_" <<     sp.celC;
        }
        
        ss << "_" << "OP_" << orderParameter;
        ss << ".txt";
        return ss.str();
    }
};

void extractValue(simParameters sp_,
                  vector<parametersToSearch> &vSearch,
                  vector<extractedValueWithSimParam> &vExtract){
    
    
    //Opening File
    ifstream dataFile;
    dataFile.clear();
    dataFile.open(sp_.fileName.c_str());
    if (!dataFile.is_open() || !dataFile.good())
    {
        cerr << "Can't open file: " <<  sp_.fileName <<endl;
        exit(1);
    }
    
    
    //Iterating through the file
    string line = "";
    while(getline (dataFile,line)){
        if (isAllParametersFound(vSearch)) {
            break;
        }
        //Checking each line to see if it's the parameter we want
        
        
        
        for (int i = 0; i < vSearch.size(); i++)
        {
            //if so, extract the value from the file
            //cout<<vSearch[i].paramName.c_str()<<endl;
            //cout<<line<<endl;
            if(line.find(vSearch[i].paramName.c_str()) == string::npos)
            {
                continue;
            }
            vSearch[i].isFoundInFile = true;
            extractedValueWithSimParam exSP;
            exSP.originalLineInFile = line;
            exSP.sp = sp_;
            exSP.paramName = vSearch[i].paramName;
            exSP.orderParameter = vSearch[i].orderParameter;
            exSP.valueDbl  = (double) atof(line.substr(vSearch[i].paramName.size()).c_str());
            exSP.valueInt  = (int) atoi(line.substr(vSearch[i].paramName.size()).c_str());
            exSP.valueBool = (bool) atoi(line.substr(vSearch[i].paramName.size()).c_str());
            vExtract.push_back(exSP);
            //exSP.toString();
        }
        
        
        
        
    }
    
    for(int j = 0; j < vSearch.size(); j++)
    {
        vSearch[j].isFoundInFile = false;
    }
    
    dataFile.close();
    
    
}

void printProgress(clock_t startClock, int finishedNumExtractions, int totalExtractions){
    cout << "Files Finished " << 100 * finishedNumExtractions / totalExtractions  << "\%" << endl;
}

struct paramSet{
    string name;
    double min;
    double max;
    double inc;
    double bFieldx;
    double bFieldy;
    double bFieldz;
    bool isRunParam;
    bool isHam;
};

bool isSimFileDone(string simFileName_, int numSweepsToPerformTotal_){
    
    //Opening File
    ifstream dataFile;
    dataFile.clear();
    dataFile.open(simFileName_.c_str());
    if (!dataFile.is_open() || !dataFile.good())
    {
        cerr << "Can't open file: " <<  simFileName_ <<endl;
        exit(1);
    }
    //Iterating through the file
    string line = "";
    int numSweepsPerformed = 1000* 1000* 1000;
    bool isThermal = false;
    
    int lineCount = 0;
    while(getline (dataFile,line)){
        //Checking each line to see if it's the parameter we want
        //if so, extract the value from the file
        string sweepsTag = "numSweepsPerformed: ";
        string isThemTag = "isThermal: ";
        if((line.find(sweepsTag.c_str()) == string::npos)&&
           (line.find(isThemTag.c_str()) == string::npos))
        {
            continue;
        }
        
        
        if(line.find(sweepsTag.c_str()) != string::npos)
        {
            numSweepsPerformed = (int) atoi(line.substr(sweepsTag.size()).c_str());
        }
        
        if(line.find(isThemTag.c_str()) != string::npos)
        {
            isThermal = (bool) atoi(line.substr(isThemTag.size()).c_str());
            if(isThermal)
            {
                cout << "Is Thermal:      " << simFileName_ << endl;
            }
            else
            {
                cout << "Not Thermal:     " << simFileName_ << endl;
            }
        }
        
        lineCount++;
        if(lineCount > 100) {break;}
    }
    dataFile.close();
    
    if((numSweepsPerformed >= numSweepsToPerformTotal_)&&
       (isThermal))
    {
        return true;
    }
    else
    {
        return false;
    }
}

const double maxBoundaryTc2Phi[63] = {
    0.335, // phi = 0
    0.335, // phi = 0.1
    0.341, // phi = 0.2
    0.345, // phi = 0.3
    0.345, // phi = 0.4
    0.345, // phi = 0.5
    0.345, // phi = 0.6
    0.345, // phi = 0.7
    0.345, // phi = 0.8
    0.29, // phi = 0.9
    0.29, // phi = 1.0
    0.29, // phi = 1.1
    0.29, // phi = 1.2
    0.157, // phi = 1.3
    0.157, // phi = 1.4
    0.157, // phi = 1.5
    0.02, // phi = 1.6
    0.11, // phi = 1.7
    0.19, // phi = 1.8
    0.19, // phi = 1.9
    0.19, // phi = 2.0
    0.23, // phi = 2.1
    0.23, // phi = 2.2
    0.23, // phi = 2.3
    0.23, // phi = 2.4
    0.21, // phi = 2.5
    0.2, // phi = 2.6
    0.13, // phi = 2.7
    0.198, // phi = 2.8
    0.248, // phi = 2.9
    0.282, // phi = 3.0
    0.31, // phi = 3.1
    0.33, // phi = 3.2
    0.33, // phi = 3.3
    0.342, // phi = 3.4
    0.35, // phi = 3.5
    0.35, // phi = 3.6
    0.35, // phi = 3.7
    0.35, // phi = 3.8
    0.35, // phi = 3.9
    0.298, // phi = 4.0
    0.298, // phi = 4.1
    0.298, // phi = 4.2
    0.214, // phi = 4.3
    0.214, // phi = 4.4
    0.214, // phi = 4.5
    0.214, // phi = 4.6
    0.02, // phi = 4.7
    0.08, // phi = 4.8
    0.12, // phi = 4.9
    0.13, // phi = 5.0
    0.14, // phi = 5.1
    0.15, // phi = 5.2
    0.156, // phi = 5.3
    0.156, // phi = 5.4
    0.152, // phi = 5.5
    0.145, // phi = 5.6
    0.131, // phi = 5.7
    0.098, // phi = 5.8
    0.162, // phi = 5.9
    0.195, // phi = 6.0
    0.225, // phi = 6.1
    0.285 // phi = 6.2
};

const double satBField[63] = {
    6.9, // phi = 0
    6.4, // phi = 0.1
    6.7, // phi = 0.2
    7.0, // phi = 0.3
    7.1, // phi = 0.4
    7.2, // phi = 0.5
    7.3, // phi = 0.6
    7.2, // phi = 0.7
    7.1, // phi = 0.8
    6.9, // phi = 0.9
    6.7, // phi = 1.0
    6.4, // phi = 1.1
    6.0, // phi = 1.2
    5.5, // phi = 1.3
    5.0, // phi = 1.4
    4.5, // phi = 1.5
    4.0, // phi = 1.6
    3.8, // phi = 1.7
    3.5, // phi = 1.8
    3.2, // phi = 1.9
    2.9, // phi = 2.0
    2.5, // phi = 2.1
    2.2, // phi = 2.2
    1.7, // phi = 2.3
    1.3, // phi = 2.4
    0.9, // phi = 2.5
    0.4, // phi = 2.6
    0.01, // phi = 2.7
    0.01, // phi = 2.8
    0.01, // phi = 2.9
    0.01, // phi = 3.0
    0.01, // phi = 3.1
    0.01, // phi = 3.2
    0.01, // phi = 3.3
    0.01, // phi = 3.4
    0.01, // phi = 3.5
    0.01, // phi = 3.6
    0.01, // phi = 3.7
    0.01, // phi = 3.8
    0.01, // phi = 3.9
    0.01, // phi = 4.0
    0.01, // phi = 4.1
    0.01, // phi = 4.2
    0.01, // phi = 4.3
    0.01, // phi = 4.4
    0.01, // phi = 4.5
    0.01, // phi = 4.6
    0.04, // phi = 4.7
    0.35, // phi = 4.8
    0.76, // phi = 4.9
    1.13, // phi = 5.0
    1.52, // phi = 5.1
    1.87, // phi = 5.2
    2.3, // phi = 5.3
    2.6, // phi = 5.4
    2.9, // phi = 5.5
    3.2, // phi = 5.6
    3.4, // phi = 5.7
    3.6, // phi = 5.8
    4.1, // phi = 5.9
    4.7, // phi = 6.0
    5.2, // phi = 6.1
    5.7 // phi = 6.2
};

bool isContinue(double iphi, double ibfield, double ikbt){
    for(int i = 0; i <= 63; i++)
    {
        if(fabs(iphi - ((double)(i/10.0)) ) < 0.09)
        {
            if( ibfield > satBField[i]  ) {return true;}
            if(ikbt > 2*maxBoundaryTc2Phi[i]  ) {return true;}
        }
    }
    return false;
}

string reSubmitFileName(const vector<simParameters> &allSimFiles_, simParameters sp_, bool reThermalize_, int numSweepsToPerformTotal_){
    string simFileName = "";
    for (int i = 0; i < allSimFiles_.size(); i++)
    {
        cout<< allSimFiles_[i].fileName<<endl;
        if(isMatchParameters(sp_, allSimFiles_[i]))
        {
            return allSimFiles_[i].fileName;
        }
    }
    cout<<"Sim File Does Not Exist " << sp_.fileName <<endl;
    sp_.toString();
    //exit(0);
    return "0";
    
    //Done checking if sim parameter exists
}

bool isSimDone(const vector<simParameters> &allSimFiles_, simParameters sp_, bool reThermalize_, int numSweepsToPerformTotal_){
    string simFileName = "";
    for (int i = 0; i < allSimFiles_.size(); i++)
    {
        //cout<< allSimFiles[i].fileName<<endl;
        if(isMatchParameters(sp_, allSimFiles_[i]))
        {
            simFileName = allSimFiles_[i].fileName;
            //cout<<simFileName<<endl;
            if (isSimFileDone(simFileName, numSweepsToPerformTotal_)&&
                (!reThermalize_))
            {
                cout<< "continue:        " << simFileName << endl;
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    
    return false;
    
    //Done checking if sim parameter exists
}

bool isRunningJob(string fileName_){
    //Begin Checking if not currently running sim Job
    stringstream ss;
    ss.str("");
    ss << "condor/log/log-" <<
    fileName_ << ".condor";
    string condorFile = ss.str();
    ifstream logFile;
    logFile.open(condorFile.c_str());
    if (!logFile.is_open() || !logFile.good())
    {
        cout << "Can't open file: " << condorFile << endl;
    }
    else
    {
        string line = "";
        string line000Sub = "";
        string line005Ter = "";
        while(getline (logFile,line)){
            if(line.find("000 (") != string::npos){
                line000Sub = line;
                line005Ter = "";
            }
            if(line.find("005 (") != string::npos){
                line005Ter = line;
            }
            
        }
        
        if(line005Ter.compare("") == 0)
        {
            return true;
        }
    }
    //End checking if not currently running sim Job
    return false;
}

