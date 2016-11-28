#include <cmath>
#include <cstdlib>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <ctime>
#include "text_tools.hpp"
#include <fenv.h>
#include <iomanip>
#include <dirent.h>
#include <map>
#include <stdint.h>
#include "fitreg_util.hpp"

#define GOMPERTZ_SIGMOID

using namespace std;
using namespace hui;

vector<double> readParameters(string inFile);

inline vector<double> calculateSelectionVariables(const vector<double> &para, const vector<double> &features, size_t featureSize, double rhoShift){
    if (para.size() % 2 == 1){
        cerr << "Error: The number of parameters is odd number!\n";
        cerr << "       Please check the parameter file.\n";
        exit(1);
    }
    if (features.size() * 2 + 2 != para.size()){
        cerr << "Error: The number of features does not match the number of weights!\n";
        cerr << "       Please check the input files.\n";
        exit(1);
    }
    double rhoLatent = para[0];
    unsigned long halfSize = para.size() / 2;
    double gammaPrimeLatent = para[halfSize];
    for (size_t i = 0; i < features.size(); i ++){
        double f = features[i];
        rhoLatent += para[i + 1] * f;
        gammaPrimeLatent += para[i + 1 + halfSize] * f;
    }

#ifdef GOMPERTZ_SIGMOID
    // a = 1
    // b = -3
    // c = -1
    // f(x) = a * exp(b * exp(c * x))
    // rhoLatent = -1;
    double rho = exp(rhoShift * exp(-rhoLatent));

    // using logistic sigmoid for gamma
    double gammaPrime = 1. / (1. + exp(-gammaPrimeLatent));
#else
    double rho = 1. / (1. + exp(-rhoLatent));
    double gammaPrime = 1. / (1. + exp(-gammaPrimeLatent));
#endif
    vector<double> selection(2);
    selection[0] = rho;
    selection[1] = gammaPrime;
    return selection;
}

int main(int argc, char *argv[]){

    string featureFile;
    string paraFile;
    string outFile;

    if (argc == 1 || (argc == 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h"))){
        cout << "usage: LINSIGHT-score [-h] -f FEATURE_FILE -p PARAMETER_FILE -o OUTPUT_FILE\n";
        exit(0);
    }

    if (argc % 2 != 1){
        cerr << "Error: odd number of arguments provided.\n";
        exit(1);
    }

    for (int i = 1; i < argc; i +=2){
        string label(argv[i]);
        string value(argv[i + 1]);
        if (label == "--feature-file" || label == "-f"){
            featureFile = value;
        }
        else if (label == "--parameter-file" || label == "-p"){
            paraFile = value;
        }
        else if (label == "--output-file" || label == "-o"){
            outFile = value;
        }
        else{
            cerr << "ERROR: Unknown command line arugment: " << label << ".\n";
            exit(1);
        }
    }

    if (featureFile == ""){
        cerr << "ERROR: Feature file is not provided.\n";
        exit(1);
    }

    if (paraFile == ""){
        cerr << "ERROR: Parameter file is not provided.\n";
        exit(1);
    }

    if (outFile == ""){
        cerr << "ERROR: Output file is not provided.\n";
        exit(1);
    }

    // parameters for sigmoid neuron
    double rhoShift = -3.;

    vector<double> para = readParameters(paraFile);
    // the number of features
    size_t featureSize = 0;

    string featureChr;
    uint64_t featureStart;
    uint64_t featureEnd;
    vector<double> features;

    // read features
    ifstream featureHandle;
    featureHandle.open(featureFile.c_str(), ios::in);
    if (!featureHandle.is_open()){
        cerr << "Error: cannot open feature file!\n";
        exit(1);
    }
    string line;

    ofstream out;

    out.open(outFile.c_str(), ios::out);
    if (!out.is_open()){
        cerr << "Error: cannot open output file!\n";
        exit(1);
    }

    // out << "#chr\tstart\tend\trho\tgamma\n";
    // if the first time of reading file
    while (readNextFeatureEntry(featureHandle, featureChr, featureStart, featureEnd, features, featureSize)){
        if (featureSize == 0){
            featureSize = features.size();
        }
        vector<double> selection = calculateSelectionVariables(para, features, featureSize, rhoShift);
        // out << featureChr << "\t" << featureStart << "\t" << featureEnd << "\t" << selection[0] << "\t" << selection[1] << endl;
        out << featureChr << "\t" << featureStart << "\t" << featureEnd << "\t" << selection[0] << endl;
    }
}

vector<double> readParameters(string inFile){
    ifstream file;
    file.open(inFile.c_str());
    if(!file.is_open()){
        cerr << "Error: cannot open parameter file!\n";
        exit(1);
    }

    file.seekg (0, file.end);
    long length = file.tellg();// get file size

    vector<double> para;
    long i;
    // loop backward over the file
    for(i = length - 2; i > 0; i-- )
    {   
        file.seekg(i);
        char c = file.get();
        if( c == '\r' || c == '\n' )//new line?
             break;
    }
    // check whether we should start from the begining of the file
    if (i == 0){
        file.seekg(0);
    }


    string line;
    getline(file, line);//read last line

    vector<string> tmp = TextTools::split(line, '\t');
    for (vector<string>::const_iterator it = tmp.begin(); it != tmp.end(); it ++){
        para.push_back(TextTools::string2double(*it));
    }

    file.close();
    return para;
}
