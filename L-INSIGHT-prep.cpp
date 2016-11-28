#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <set>
#include <cstdint>
#include "text_tools.hpp"

// using namespace fitreg;
using namespace std;
using namespace hui;

vector<double> computeWattersons_a(unsigned long numValues);

int main(int argc, char *argv[]){
    string inFile;
    string outDir;

    if (argc == 1 || (argc == 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h"))){
        cout << "usage: LINSIGHT-prep [-h] -i INSIGHT_INPUT_FILE -o INSIGHT_DATA_DIRECTORY\n";
        exit(0);
    }

    if (argc % 2 != 1){
        cerr << "Error: odd number of arguments provided.\n";
        exit(1);
    }

    for (int i = 1; i < argc; i +=2){
        string label(argv[i]);
        string value(argv[i + 1]);
        if (label == "--input-file" || label == "-i"){
            inFile = value;
        }
        else if (label == "--output-dir" || label == "-o"){
            outDir = value;
        }
        else{
            cerr << "ERROR: Unknown command line argument: " << label << ".\n";
            exit(1);
        }
    }

    if (inFile == ""){
        cerr << "ERROR: Input INSIGHT file is not provided.\n";
        exit(1);
    }

    if (outDir == ""){
        cerr << "ERROR: Output directory is not provided.\n";
        exit(1);
    }

    // A hack to prevent very small lambda and theta
    const double minThetaAndLambda = 1.e-5;

    ifstream inHandle;
    inHandle.open(inFile.c_str());

    if (!inHandle.is_open()){
        cerr << "Error: cannot open file: " << inFile << endl;
        exit(1);
    }

    string entry;
    getline(inHandle, entry);
    vector<string> items = TextTools::splitByWhiteSpace(entry);
    unsigned long nbSamples;
    if (items[0] == "samples"){
        nbSamples = TextTools::string2long(items[1]);
    }
    else{
        throw runtime_error("Cannot get sample size!\n");
    }

    vector<double> watternsonsMatrix = computeWattersons_a(nbSamples);

    double thetaTimesA = 0;
    double lambda = 0;
    char * outBuffer = (char *) malloc(sizeof(uint64_t) + 5 * sizeof(float));
    set<string> chromSet;
    ofstream outHandle;
    while (getline(inHandle, entry)){
        vector<string> items = TextTools::splitByWhiteSpace(entry);
        if (items[0] == "block"){
            double theta = TextTools::string2double(items[3]);
            lambda = TextTools::string2double(items[5]);
            thetaTimesA = theta * watternsonsMatrix[nbSamples];
            if (thetaTimesA < minThetaAndLambda){
                thetaTimesA = minThetaAndLambda;
            }
            if (lambda < minThetaAndLambda){
                lambda = minThetaAndLambda;
            }
        }
        else if (items[0] == "site"){
            vector<string> tmp = TextTools::split(items[1], ':');
            string chr = tmp[0];
            uint64_t pos = TextTools::string2long(tmp[1]);
            double pZeqXmaj = 0;
            double pZeqXmin = 0;
            double alleleFreq = 0;
            if (items[2] == "M"){
                pZeqXmaj = TextTools::string2double(items[3]);
            }
            else if (items[2] == "L"){
                pZeqXmaj = TextTools::string2double(items[3]);
                pZeqXmin = TextTools::string2double(items[4]);
                // dummy allele frequency. We essentially use the three category representation of SFS. 
                alleleFreq = 0.05;
            }
            else if (items[2] == "H"){
                pZeqXmaj = TextTools::string2double(items[3]);
                pZeqXmin = TextTools::string2double(items[4]);
                // dummy allele frequency. We essentially use the three category representation of SFS. 
                alleleFreq = 0.5;
            }
            else{
                throw runtime_error("Wrong site type!\n");
            }

            *((uint64_t *) outBuffer) = pos;
            float * numPara = (float *) (outBuffer + sizeof(uint64_t));
            numPara[0] = (float) alleleFreq;
            numPara[1] = (float) pZeqXmaj;
            numPara[2] = (float) pZeqXmin;
            numPara[3] = (float) thetaTimesA;
            numPara[4] = (float) lambda;
            if (chromSet.find(chr) == chromSet.end()){
                outHandle.close();
                string outFile = outDir + "/" + chr + ".fit";
                outHandle.open(outFile.c_str(), ios::binary | ios::out);

                if (!outHandle.is_open()){
                    cerr << "Error: cannot open file: " << outFile << endl;
                    exit(1);
                }

                chromSet.insert(chr);
            }
            outHandle.write(outBuffer, sizeof(uint64_t) + 5 * sizeof(float));
        }
    }
    outHandle.close();

}

vector<double> computeWattersons_a(unsigned long numValues){
    if (numValues < 1){
        throw runtime_error("Illegel number of samples!\n");
    }
    vector<double> valueArray(numValues + 1);
    valueArray[0] = 0.0;
    valueArray[1] = 0.0;

    for(unsigned long k=1; k<numValues;k++) {
        valueArray[k+1] = valueArray[k] + (1.0 / k);
    }
    return valueArray;
}
