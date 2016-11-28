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
#include <cstdint>
#include "fitreg_util.hpp"

#define GOMPERTZ_SIGMOID
// #define DEBUG

using namespace std;
using namespace hui;

inline void updateGradient(vector<double> &gradient, const float *evolPara, const vector<double> &para, const vector<double> &features, const vector<double> &beta, double alleleFreqCutoff, double rhoShift, double gammaPrimeShift){

        double valA, valB, valC, valD, valE, valF;
        float alleleFreq = evolPara[0];
        float pZeqXmaj = evolPara[1];
        float pZeqXmin = evolPara[2];
        float thetaTimesA = evolPara[3];
        float lambda = evolPara[4];

        // add bias terms
        double rhoLatent = para[0];
        unsigned long halfSize = para.size() / 2;
        double gammaPrimeLatent = para[halfSize];

        // calculate latent variables
        for (size_t i = 0; i < features.size(); i ++){
            double f = features[i];
            rhoLatent += para[i + 1] * f;
            gammaPrimeLatent += para[i + 1 + halfSize] * f;
        }

#ifdef GOMPERTZ_SIGMOID
        // a = 1
        // b = rhoShift
        // c = -1
        // f(x) = a * exp(b * exp(c * x))
        // rhoLatent = -1;
        double rho = exp(rhoShift * exp(-rhoLatent));

        // Using logistic sigmoid function for gamma
        double gammaPrime = 1. / (1. + exp(-gammaPrimeLatent));
#else
        double rho = 1. / (1. + exp(-rhoLatent));
        double gammaPrime = 1. / (1. + exp(-gammaPrimeLatent));
#endif
        // eta is fixed to zero (ignore positive selection)
        double eta = 0;
       
        // the definition of gamma in INSIGHT
        double gamma = gammaPrime * beta[0];

        // for numerical stability
        if (lambda < 0.00001)
            lambda = 0.00001;

        if (thetaTimesA < 0.00001)
            thetaTimesA = 0.00001;

        // modified from INSIGHT
        // monomorphic sites
        // if (type == 'M'){
        if (alleleFreq == 0){
            valA = 1.0;
            valB = pZeqXmaj;
            valC = lambda*((1-pZeqXmaj)/3.0 - pZeqXmaj);
            valD = -pZeqXmaj*thetaTimesA;
            valE = pZeqXmaj*thetaTimesA*lambda;
            valF = (1-thetaTimesA) * ( (1-lambda) * pZeqXmaj + lambda/3 * (1-pZeqXmaj) );
        }
        // else if (type == 'L'){
        else if (alleleFreq < alleleFreqCutoff && alleleFreq > 0){
            valA = thetaTimesA / 3.0;
            valB = 0.0;
            valC = 0.0;
            valD = 1.0;
            valE = -lambda;
            valA *= pZeqXmaj;
            valF =
                 lambda/3  * (beta[0] * (1-pZeqXmaj) + beta[2] * (1-pZeqXmin) ) * thetaTimesA / 3 +
                (1-lambda) * (beta[0] *   pZeqXmaj   + beta[2] *   pZeqXmin   ) * thetaTimesA / 3;
        }
        // else if (type == 'H'){
        else if (alleleFreq >= alleleFreqCutoff && alleleFreq <= 1){
            valA = 0.0;
            valB = 0.0;
            valC = 0.0;
            valD = 0.0;
            valE = 0.0;

            valF = ((1.0-2.0*lambda/3.0) * (pZeqXmaj+pZeqXmin) + (2.0*lambda/3.0) * (1-pZeqXmaj-pZeqXmin) )
                    * beta[1]*thetaTimesA/3.0;
        }
        else{
            cout << "### " << alleleFreq << endl;
            throw runtime_error("Unknown site type in input data!\n");
        }

        double derivRho     = valA * (valB + valC*eta + valD*gamma + valE*eta*gamma) - valF;
        double derivGamma   = valA * (valD + valE*eta);

        // likelihood
        double valP = derivRho*rho + valF;

        // the derivatives of loss function
        derivRho = -derivRho / valP;
        derivGamma = -derivGamma / valP;
        double derivGammaPrime = derivGamma / beta[0];

        // the intermediate values used in the gradient computation
        // use of the intermediate results and significantly reduce computational cost
#ifdef GOMPERTZ_SIGMOID
        double rhoIntermediate = derivRho * -1 * rhoShift * exp(rhoShift * exp(-rhoLatent)) * exp(-rhoLatent);

        // Using the original logistic sigmoid function for gamma
        double gammaPrimeIntermediate = derivGammaPrime * gammaPrime * gammaPrime * exp(-gammaPrimeLatent);
#else
        double rhoIntermediate = derivRho * rho * rho * exp(-rhoLatent);
        double gammaPrimeIntermediate = derivGammaPrime * gammaPrime * gammaPrime * exp(-gammaPrimeLatent);
#endif

        gradient[0] += rhoIntermediate;
        gradient[halfSize] += gammaPrimeIntermediate;
        for (size_t i = 0; i < features.size(); i ++){
            double f = features[i];
            gradient[i + 1] += rhoIntermediate * f;
            gradient[halfSize + i + 1] += gammaPrimeIntermediate * f;
        }
}

inline void updateParameters(vector<double> &para, vector<double> &gradient, vector<double> &histGrad, vector<double> histGrad2, double learningRateRho, double learningRateGammaPrime, double epsilon, unsigned long miniBatchSize){
    for (size_t i = 0; i < gradient.size(); i++){
        // normalize gradient
        double deriv = gradient[i];
        histGrad[i] += deriv * deriv;
        if (i < para.size() / 2){
            para[i] -= deriv / (sqrt(histGrad[i]) + epsilon) * learningRateRho;
        }
        else{
            para[i] -= deriv / (sqrt(histGrad[i]) + epsilon) * learningRateGammaPrime;
        }
        gradient[i] = 0;
    }
}

int main(int argc, char *argv[]){
    string inDir;
    string featureFile;
    string outputFile;
    unsigned long nbEpoch = 0;
    double learningRateRho = 0;
    double learningRateGammaPrime = 0;
    vector<double> beta(3);

    if (argc == 1 || (argc == 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h"))){
        cout << "usage: LINSIGHT-fit [-h] -d INSIGHT_DATA_DIRECTORY -f FEATURE_FILE\n";
        cout << "                    -o OUTPUT_FILE -n EPOCH_NUMBER -b1 BETA_1 -b2 BETA_2\n";
        cout << "                    [-r1 RHO_LEARNING_RATE] [-r2 GAMMA_LEARNING_RATE]\n";
        exit(0);
    }

    if (argc % 2 != 1){
        cerr << "Error: odd number of arguments provided.\n";
        exit(1);
    }

    for (int i = 1; i < argc; i +=2){
        string label(argv[i]);
        string value(argv[i + 1]);
        if (label == "--data-dir" || label == "-d"){
            inDir = value;
        }
        else if (label == "--feature-file" || label == "-f"){
            featureFile = value;
        }
        else if (label == "--output-file" || label == "-o"){
            outputFile = value;
        }
        else if (label == "--epoch-number" || label == "-n"){
            nbEpoch = TextTools::string2long(value);
        }
        else if (label == "--rho-rate" || label == "-r1"){
            learningRateRho = TextTools::string2double(value);
        }
        else if (label == "--gamma-rate" || label == "-r2"){
            learningRateGammaPrime = TextTools::string2double(value);
        }
        else if (label == "--beta-one" || label == "-b1"){
            beta[0] = TextTools::string2double(value);
        }
        else if (label == "--beta-two" || label == "-b2"){
            beta[1] = TextTools::string2double(value);
        }
        else{
            cerr << "ERROR: Unknown command line arugment: " << label << ".\n";
            exit(1);
        }
    }

    beta[2] = 1. - beta[0] - beta[1];

    if (inDir == ""){
        cerr << "Error: Input directory of phylogenetic/polymorphic data is not provided." << endl;
        exit(1);
    }

    if (featureFile == ""){
        cerr << "Error: Feature file is not provided." << endl;
        exit(1);
    }

    if (outputFile == ""){
        cerr << "Error: Out file is not provided." << endl;
        exit(1);
    }

    if (nbEpoch == 0){
        cerr << "Error: Epoch number is not provided." << endl;
        exit(1);
    }

    if (learningRateRho == 0){
        learningRateRho = 0.001;
    }

    if (learningRateGammaPrime == 0){
        learningRateGammaPrime = learningRateRho * 10;
    }

    if (beta[0] <= 0 || beta[0] >= 1){
        cerr << "Error: Invalid Beta_1: " << beta[0] << endl;
        exit(1);
    }

    if (beta[1] <= 0 || beta[1] >= 1){
        cerr << "Error: Invalid Beta_2: " << beta[1] << endl;
        exit(1);
    }

    if (beta[2] <= 0 || beta[2] >= 1){
        cerr << "Error: Invalid Beta_3: " << beta[2] << endl;
        exit(1);
    }

    double alleleFreqCutoff = 0.15;
    // parameters for sigmoid neuron
    double rhoShift = -3.;
    double gammaPrimeShift = -3.;

    cout << "Original CMD:";
    for (int i = 0; i < argc; i ++){
        cout << " " << argv[i];
    }
    cout << endl;
    cout << "Number of epoch = " << nbEpoch << endl;
    cout << "Learning rate of rho = " << learningRateRho << endl;
    cout << "Learning rate of gamma = " << learningRateGammaPrime << endl;
    cout << "Beta1 = " << beta[0] << endl;
    cout << "Beta2 = " << beta[1] << endl;
    cout << "Beta3 = " << beta[2] << endl;

    const unsigned long chunkSize = 1e6;
    const unsigned long miniBatchSize = 100;
    // const unsigned long miniBatchSize = 1;
    const double epsilon = 1.e-6;

    map<string, string> chrFiles = readFitregDirectory(inDir);

    size_t singleEntrySize = sizeof(uint64_t) + sizeof(float) * 5;
    size_t bufferSize = chunkSize * singleEntrySize;
    char *buffer = (char *) malloc(bufferSize);
    if (buffer == NULL){
        throw runtime_error("Cannot allocate buffer!\n");
    }



    // the vector of feature sizes
    // vector<unsigned long> featureSizes;
    size_t featureSize = 0;
    // the number of weights (except bias)
    unsigned long nbWeights = 0;

    // the gradient of minibatch
    vector<double> gradient;

    // the weights
    vector<double> para;

    // the historical gradient
    vector<double> histGrad;
    vector<double> histGrad2;

    unsigned long counter = 0;
    ofstream out;
    out.open(outputFile.c_str(), ios::out);
    if (!out.is_open()){
        cerr << "Error: cannot open output file!\n";
        exit(1);
    }

    for (unsigned long epoch = 0; epoch < nbEpoch; epoch ++){
        // read features
        ifstream featureHandle;
        featureHandle.open(featureFile.c_str(), ios::in);
        if (!featureHandle.is_open()){
            cerr << "Error: cannot open feature file!\n";
            exit(1);
        }
        string line;

        // reset counter for mini-batch
        counter = 0;

        string featureChr;
        uint64_t featureStart;
        uint64_t featureEnd;
        vector<double> features;

        // read the first feature
        readNextFeatureEntry(featureHandle, featureChr, featureStart, featureEnd, features, featureSize);

        // if the first time of reading file
        if (featureSize == 0){
            featureSize = features.size();
            nbWeights = featureSize;

            cout << "Number of weights = " << nbWeights << endl;
            gradient = vector<double>(2 * (1 + nbWeights));
            histGrad = vector<double>(2 * (1 + nbWeights));
            histGrad2 = vector<double>(2 * (1 + nbWeights));
            para  = vector<double>(2 * (1 + nbWeights));
        }

        // read data
        uint64_t totalSite = 0;
        uint64_t totalGoodSite = 0;
        for (map<string, string>::const_iterator it = chrFiles.begin(); it != chrFiles.end(); it ++){
#ifdef DEBUG
            cout << "### Chrom = " << it->first << "\t" << it->second << endl;
#endif
            string chr = it->first;
            string fileDir = it->second;
            ifstream inHandle;
            inHandle.open(fileDir.c_str(), ios::binary | ios::in);
            if (!inHandle.is_open()){
                throw runtime_error("Cannot open data file!\n");
            }
            size_t nbExtractedEntry = 0;

            do {
                inHandle.read(buffer, bufferSize);
                nbExtractedEntry = inHandle.gcount();
                for (size_t i = 0; i < nbExtractedEntry; i += singleEntrySize){
                    char *p = buffer + i;
                    uint64_t pos = *((uint64_t *) p);
                    bool matchFlag = false;
                    if (chr == featureChr && pos >= featureStart && pos < featureEnd){
                        // find corresponding features
                        matchFlag = true;
                    }
                    else if (chr.compare(featureChr) > 0 || (chr == featureChr && pos >= featureEnd)){
                        while (readNextFeatureEntry(featureHandle, featureChr, featureStart, featureEnd, features, featureSize)){
                            if (chr.compare(featureChr) < 0 || pos < featureStart){
                                // cannot find corresponding features
                                break;
                            }
                            else if (chr == featureChr && pos >= featureStart && pos < featureEnd){
                                // find corresponding features
                                matchFlag = true;
                                break;
                            }
                        }
                    }

                    totalGoodSite += matchFlag;
                    totalSite ++;
                    
                    if (matchFlag){
                        // increase the number of count
                        counter ++;
                    
                        updateGradient(gradient, (float *) (p + sizeof(uint64_t)), para, features, beta, alleleFreqCutoff, rhoShift, gammaPrimeShift);
                        
                        // update parameters
                        if (counter == miniBatchSize){
                            updateParameters(para, gradient, histGrad, histGrad2, learningRateRho, learningRateGammaPrime, epsilon, miniBatchSize);

                            counter = 0;
                        }
                    }
                }
#ifdef DEBUG
                cout << nbExtractedEntry << "\t" << nbExtractedEntry / singleEntrySize << "\t" << nbExtractedEntry % singleEntrySize << endl;
                unsigned long kk = 0;
                for (vector<double>::const_iterator it = para.begin(); it != para.end(); it ++){
                    cout << "\t" << *it;
                    kk ++;
                    if (kk == para.size() / 2){
                        cout << endl;
                    }
                }
                cout << endl << endl;
#endif
            } while (nbExtractedEntry);
        }

        featureHandle.close();

#ifdef DEBUG
        cout << "Total site = " << totalSite << "; total good site = " << totalGoodSite << endl;
#endif

        // write parameters after each epoch
        for (size_t m = 0; m < para.size(); m ++){
            if (m == 0){
                out << para[m];
            }
            else{
                out << "\t" << para[m];
            }
        }
        out << endl;
    }
    out.close();

    free(buffer);
}
