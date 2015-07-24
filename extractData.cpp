#include "header.h"

int main(){
    cout<<"start"<<endl;
    
    vector<parametersToSearch> parametersToSearchVec;
    parametersToSearchVec.clear();
    
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Magnetization's Order Parameter: ", "Magnetization_M", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Neel's Order Parameter: ", "Magnetization_N", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Stripy-x's Order Parameter: ", "Magnetization_SX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Stripy-y's Order Parameter: ", "Magnetization_SY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Stripy's Order Parameter: ", "Magnetization_S", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Product-x: ", "Magnetization_PX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Product-y: ", "Magnetization_PY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Nematic-Q: ", "Magnetization_Q", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Dimer: ", "Magnetization_Dimer", false));
    
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Magnetization's Order Parameter: ", "Binder_M", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Neel's Order Parameter: ", "Binder_N", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Stripy-x's Order Parameter: ", "Binder_SX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Stripy-y's Order Parameter: ", "Binder_SY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Stripy's Order Parameter: ", "Binder_S", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Product-x: ", "Binder_PX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Product-y: ", "Binder_PY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Nematic-Q: ", "Binder_Q", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Susceptibility - Dimer: ", "Binder_Dimer", false));
    
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Magnetization's Order Parameter: ", "Susceptibility_M", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Neel's Order Parameter: ", "Susceptibility_N", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Stripy-x's Order Parameter: ", "Susceptibility_SX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Stripy-y's Order Parameter: ", "Susceptibility_SY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Stripy's Order Parameter: ", "Susceptibility_S", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Product-x: ", "Susceptibility_PX", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Product-y: ", "Susceptibility_PY", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Nematic-Q: ", "Susceptibility_Q", false));
    parametersToSearchVec.push_back(parametersToSearch("Lattice Averaged Quantity, Binder - Dimer: ", "Susceptibility_Dimer", false));
    
    parametersToSearchVec.push_back(parametersToSearch("Energy Per Site: ", "Energy", false));
    //parametersToSearchVec.push_back(parametersToSearch("numSweepsPerformed: ", "numSweepsPerformed", false));
    parametersToSearchVec.push_back(parametersToSearch("Specific Heat: ", "SpecificHeat", false));
    
    const int bigNum = 10 * 1000 * 1000;
    
    //Remember that each of these vectors must have at least one element in it - ie. "0".
    //Interested Phis
    vector<double> phiCombinations;
    phiCombinations.clear();
    int startPhi = 0;
    int   endPhi = 0;//63;
    int   incPhi = 1;
    for(int i = startPhi; i <= endPhi; i+=incPhi){
        phiCombinations.push_back(i/10.0);
    }
    
    /*
     phiCombinations.push_back(0);
     phiCombinations.push_back(0.8);
     phiCombinations.push_back(2.1);
     phiCombinations.push_back(2.5);
     phiCombinations.push_back(2.9);
     phiCombinations.push_back(4);
     phiCombinations.push_back(5.2);
     phiCombinations.push_back(5.7);
     phiCombinations.push_back(6);
     */
    
    //Interested Cells
    vector<int> cellCombinations;
    cellCombinations.clear();
    cellCombinations.push_back(36);
    cellCombinations.push_back(48);
    cellCombinations.push_back(60);
    cellCombinations.push_back(84);
    /*cellCombinations.push_back(96);
     cellCombinations.push_back(108);
     cellCombinations.push_back(120);
     cellCombinations.push_back(144);
     cellCombinations.push_back(168);
     cellCombinations.push_back(204);*/
    
    
    //KbT
    vector<double> KbTCombinations;
    KbTCombinations.clear();
    int startKbT = 1;
    int   endKbT = 300;//80;
    int   incKbT = 1;
    for(int i = startKbT; i <= endKbT; i+=incKbT){
        //double j = ((int) (bigNum * i)) / (1.0 * bigNum);
        KbTCombinations.push_back(i/100.0);
    }
    
    //Interested bField Directions
    struct bFieldStruct{ int bFieldX; int bFieldY; int bFieldZ;};
    vector<bFieldStruct> bFieldCombinations;
    bFieldStruct bField000;
    bField000.bFieldX = 0; bField000.bFieldY = 0; bField000.bFieldZ = 0;
    bFieldStruct bField001;
    bField001.bFieldX = 0; bField001.bFieldY = 0; bField001.bFieldZ = 1;
    bFieldStruct bField110;
    bField110.bFieldX = 1; bField110.bFieldY = 1; bField110.bFieldZ = 0;
    bFieldStruct bField111;
    bField111.bFieldX = 1; bField111.bFieldY = 1; bField111.bFieldZ = 1;
    bFieldStruct bField5711;
    bField5711.bFieldX = 5; bField5711.bFieldY = 7; bField5711.bFieldZ = 11;
    bFieldCombinations.clear();
    bFieldCombinations.push_back(bField000);
    //bFieldCombinations.push_back(bField001);
    //bFieldCombinations.push_back(bField110);
    //bFieldCombinations.push_back(bField111);
    //bFieldCombinations.push_back(bField5711);
    
    //Interested bMag
    vector<double> bMagCombinations;
    bMagCombinations.clear();
    int startbMag = 0;
    int   endbMag = 0;//750;
    int   incbMag = 1;
    for(int i = startbMag; i <= endbMag; i+=incbMag){
        bMagCombinations.push_back(i/100.0);
    }
    
    //Interested J2
    vector<double> J2Combinations;
    J2Combinations.clear();
    int startJ2 = 300;//-100;
    int   endJ2 = 300;//100;
    int   incJ2 = 1;
    for(int i = startJ2; i <= endJ2; i+=incJ2){
        J2Combinations.push_back(i/100.0);
    }
    
    //Interested J3
    vector<double> J3Combinations;
    J3Combinations.clear();
    int startJ3 = 0;
    int   endJ3 = 0;//100
    int   incJ3 = 1;
    for(int i = startJ3; i <= endJ3; i+=incJ3){
        //double j = ((int) (bigNum * i)) / (1.0 * bigNum);
        J3Combinations.push_back(i/100.0);
    }
    
    //Interested K1
    vector<double> K1Combinations;
    K1Combinations.clear();
    int startK1 = -1000;//200;
    int   endK1 = 1000;//200;
    int   incK1 = 1;
    for(int i = startK1; i <= endK1; i+=incK1){
        //double j = ((int) (bigNum * i)) / (1.0 * bigNum);
        //if(i == 0){continue;}
        //K1Combinations.push_back(i/1000.0);
    }
    //K1Combinations.push_back(0.707);
    //K1Combinations.push_back(-0.707);
    //K1Combinations.push_back(1);
    //K1Combinations.push_back(-1);
    K1Combinations.push_back(0);
    
    //Interested K2
    vector<double> K2Combinations;
    K2Combinations.clear();
    int startK2 = -1000;//200;
    int   endK2 = 1000;//200;
    int   incK2 = 1;
    for(int i = startK2; i <= endK2; i+=incK2){
        //double j = ((int) (bigNum * i)) / (1.0 * bigNum);
        //if(i == 0){continue;}
        //K2Combinations.push_back(i/1000.0);
    }
    //K2Combinations.push_back(0.707);
    //K2Combinations.push_back(-0.707);
    //K2Combinations.push_back(1);
    //K2Combinations.push_back(-1);
    K2Combinations.push_back(0);
    
    //cubicD
    vector<double> cubicDCombinations;
    cubicDCombinations.clear();
    int startD = 1;
    int   endD = 1;
    int   incD = 1;
    for(int i = startD; i <= endD; i+=incD){
        //double j = ((int) (bigNum * i)) / (1.0 * bigNum);
        cubicDCombinations.push_back(i/10.0);
    }
    
    cout<<"Done Making Combination Vectors"<<endl;
    
    vector<simParameters> allSimFiles;
    allSimFiles.clear();
    //string filesToExtract = "/Volumes/Toshi/dataFromSimulations/honV9K12/hep2p/dataFiles";
    string filesToExtract = "/Users/Work/Documents/perkins_research/chalco_model/Chalco_Model/dataFiles";
    getAllSimFiles(allSimFiles, filesToExtract);
    cout<< "Done Getting All Files" << endl;
    
    vector<extractedValueWithSimParam> extractedValuesVec;
    extractedValuesVec.clear();
    
    //Making a vector of all of the parameters we want to plot versus KbT
    int size = phiCombinations.size() * cellCombinations.size() *
    KbTCombinations.size() * J2Combinations.size() * J3Combinations.size() *
    K2Combinations.size() * K1Combinations.size() * cubicDCombinations.size() * bFieldCombinations.size() *
    bMagCombinations.size();
    int count = 0;
    
    cout<<phiCombinations.size()<<endl;
    cout<<cellCombinations.size()<<endl;
    cout<<KbTCombinations.size()<<endl;
    cout<<J2Combinations.size()<<endl;
    cout<<J3Combinations.size()<<endl;
    cout<<cubicDCombinations.size()<<endl;
    cout<<K1Combinations.size()<<endl;
    cout<<K2Combinations.size()<<endl;
    cout<<bMagCombinations.size()<<endl;
    cout<<bFieldCombinations.size()<<endl;
    //exit(0);
    
    for (int i1 = 0; i1 < phiCombinations.size(); i1++) {
        //cout << "i1: " << i1 << endl;
        for (int i2 = 0; i2 < cellCombinations.size(); i2++) {
            //cout << "i2: " << i2 << endl;
            for (int i3 = 0; i3 < KbTCombinations.size(); i3++) {
                //cout << "i3: " << i3 << endl;
                for (int i4 = 0; i4 < J2Combinations.size(); i4++) {
                    //cout << "i4: " << i4 << endl;
                    for (int i5 = 0; i5 < J3Combinations.size(); i5++) {
                        //cout << "i5: " << i5 << endl;
                        for (int i6 = 0; i6 < cubicDCombinations.size(); i6++) {
                            //cout << "i6: " << i6 << endl;
                            for (int i10 = 0; i10 < K1Combinations.size(); i10++) {
                                //cout << "i10: " << i10 << endl;
                                for (int i7 = 0; i7 < K2Combinations.size(); i7++) {
                                    //cout << "i7: " << i7 << endl;
                                    for (int i8 = 0; i8 < bMagCombinations.size(); i8++) {
                                        //cout << "i8: " << i8 << endl;
                                        for (int i9 = 0; i9 < bFieldCombinations.size(); i9++) {
                                            //cout << "i9: " << i9 << endl;
                                            
                                            if(count % 10000 == 0) printProgress(0, count, size);
                                            count++;
                                            simParameters sp;
                                            sp.setCel(cellCombinations[i2],cellCombinations[i2],1);
                                            sp.setKbT(KbTCombinations[i3]);
                                            
                                            if(phiCombinations[i1] > 1e-8)
                                            {
                                                sp.setPhi(true, phiCombinations[i1]);
                                            }
                                            if((fabs(J2Combinations[i4]) > 1e-8)||(fabs(J3Combinations[i5]) > 1e-8)||(fabs(K2Combinations[i7]) > 1e-8))
                                            {
                                                //sp.setJ23(true, J2Combinations[i4], J3Combinations[i5], K2Combinations[i7]);//(-2*J2Combinations[i4]));//K2Combinations[i7]);
                                            }
                                            if((fabs(J2Combinations[i4]) > 1e-8)||(fabs(J3Combinations[i5]) > 1e-8)||(fabs(K2Combinations[i7]) > 1e-8)||(fabs(K1Combinations[i10]) > 1e-8))
                                            {
                                                sp.setK12(true, 1/*J1*/, K1Combinations[i10],  J2Combinations[i4],K2Combinations[i7], J3Combinations[i5]);//(-2*J2Combinations[i4]));//K2Combinations[i7]);
                                            }
                                            if(cubicDCombinations[i6] > 1e-8)
                                            {
                                                sp.setCub(false,cubicDCombinations[i6]);
                                            }
                                            
                                            if(bMagCombinations[i8] > 1e-8)
                                            {
                                                sp.setB(true,
                                                        bFieldCombinations[i9].bFieldX,
                                                        bFieldCombinations[i9].bFieldY,
                                                        bFieldCombinations[i9].bFieldZ,
                                                        bMagCombinations[i8]);
                                            }
                                            
                                            //cout<<"df"<<endl;
                                            //cout<<bMagCombinations[i8] << " " << phiCombinations[i1] << " " <<cellCombinations[i2] << " "<<KbTCombinations[i3]<<" "<<J2Combinations[i4] << " "<<J3Combinations[i5]<<" "<<AlphaCombinations[i6]<< K1Combinations[i10]<< " " << K2Combinations[i7] <<endl;
                                            
                                            //
                                            //Skip all files that we aren't looking for with the above parameters.
                                            //sp.toString();
                                            if (!isMatchAnyParameters(sp, allSimFiles))
                                            {
                                                //sp.toString();
                                                //exit(0);
                                                continue;
                                            }
                                            //Store values in a vector.
                                            //Extracting the values of the parameters
                                            //cout<<"df"<<endl;
                                            //cout<<allSimFiles.size()<<endl;
                                            for (int i = 0; i < allSimFiles.size(); i++) {
                                                //sp.toString();
                                                //allSimFiles[i].toString();
                                                if(isMatchParameters(sp, allSimFiles[i])){
                                                    //sp.toString();
                                                    extractValue(allSimFiles[i],
                                                                 parametersToSearchVec,
                                                                 extractedValuesVec);
                                                }
                                            }
                                            //exit(0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    cout << "Done extracting all data" << endl;
    
    //Opening clean graph files.
    string graphTitle = "readyToPlot/graphByKbT";
    //string graphTitle = "readyToPlot/graphByB";
    
    ofstream graphFile;
    stringstream ss;
    
    
    //Appending to the data files.
    for (int i = 0; i < extractedValuesVec.size(); i++) {
        if (extractedValuesVec[i].isExtracted) continue;
        ss.str("");
        //ss << graphTitle << extractedValuesVec[i].toFileName("B");
        ss << graphTitle << extractedValuesVec[i].toFileName("KbT");
        
        /*
        graphFile.clear();
        graphFile.open(ss.str().c_str());
        //cout<<ss.str()<<endl;
        if (!graphFile.is_open() || !graphFile.good()){
            cout << "Can't open graphFile" << endl;
            exit(1);
        }
        graphFile.close();
         */
        
        graphFile.clear();
        graphFile.open(ss.str().c_str(), ios::app);
        if (!graphFile.is_open() || !graphFile.good()){
            cout << "Can't open graphFile" << endl;
            exit(1);
        }
        ss.str("");
        
        //cout<<"Df"<<endl;
        ss << std::scientific << std::setprecision(14) << extractedValuesVec[i].sp.KbT;
        //ss << std::scientific << std::setprecision(14) << extractedValuesVec[i].sp.bMag;
        
        
        for (int j = 0; j < cellCombinations.size(); j++) {
            for (int k = 0; k < extractedValuesVec.size(); k++) {
                simParameters sp = extractedValuesVec[i].sp;
                sp.celA = cellCombinations[j];
                sp.celB = cellCombinations[j];
                if(isMatchParameters(sp, extractedValuesVec[k].sp) && extractedValuesVec[i].paramName == extractedValuesVec[k].paramName){
                    ss << " " << std::scientific << std::setprecision(14) << cellCombinations[j];
                    ss << " " << std::scientific << std::setprecision(14) << extractedValuesVec[k].valueDbl;
                    extractedValuesVec[k].isExtracted = true;
                    break;
                }
            }
        }
        graphFile << ss.str() << endl;
        graphFile.close();
    }
    
    
    cout << "Done writing new graph files" << endl;
    
    system("echo -e '\a'");
    cout<<"done!"<<endl;
    return 0;
    
}


