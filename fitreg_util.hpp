#ifndef FITREG_UTIL_HPP
#define FITREG_UTIL_HPP

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

using namespace std;
using namespace hui;

// the size of the line buffer
// It must be large enough!!
const size_t lineBufferSize = 1000000;

map<string, string> readFitregDirectory(string inDir){
    // read list of fitreg files
    DIR *dir;
    struct dirent *ent;
    map<string, string> chrFiles;
    string ending(".fit");
    if ((dir = opendir (inDir.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            string fileName(ent->d_name);
            if (fileName.length() >= ending.length()) {
                if (fileName.compare (fileName.length() - ending.length(), ending.length(), ending) == 0){
                    vector<string> tokens = TextTools::split(fileName, '.');
                    string chr = tokens[0];
                    string dir = inDir + "/" + fileName;
                    chrFiles[chr] = dir;
                }
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        //throw runtime_error("Could not open directory!\n");
        cerr << "Error: cannot open INSIGHT data directory!\n";
        exit(1);
    }
    return chrFiles;
}

bool readNextFeatureEntry(ifstream &handle, string &chr, uint64_t &start, uint64_t &end, vector<double> &features, size_t featureSize){
    string oldChr = chr;
    uint64_t oldEnd = end;

    string line;
    // skip command lines and null lines here
    while (line == "" || line[0] == '#'){
        istream &state = getline(handle, line);
        // check end of file
        if (!state){
            return false;
        }
    }
    
    char buffer[lineBufferSize];
    if (lineBufferSize < line.length() + 10){
        throw runtime_error("Line is longer than buffer size!\n");
    }

    strcpy(buffer, line.c_str());
    char *token;
    // chromosome
    if ((token = strtok(buffer, "\t")) != NULL){
        chr = token;
    }
    else{
        throw runtime_error("Cannot obtain genomic chr number!\n");
    }

    // start
    if ((token = strtok(NULL, "\t")) != NULL){
        start = strtoull(token, NULL, 0);
    }
    else{
        throw runtime_error("Cannot obtain genomic chr number!\n");
    }

    // end
    if ((token = strtok(NULL, "\t")) != NULL){
        end = strtoull(token, NULL, 0);
    }
    else{
        throw runtime_error("Cannot obtain genomic chr number!\n");
    }

    features.clear();
    while ((token = strtok(NULL, "\t")) != NULL){
        features.push_back(atof(token));
    }

    // check validation of features
    if (featureSize != 0 && features.size() != featureSize){
        cerr << "Error:" << endl;
        cerr << line << endl;
        cerr << features.size() << " != " << featureSize << endl;
        throw runtime_error("The number of features do not mathch the title line!\n");
    }

    // check whether the features are sorted and proper
    if (chr == oldChr){
        if (start < oldEnd){
            throw runtime_error("Feature file is not sorted properly or has overlapping intervals. Please check the input file!!\n");
        }

        if (end < start){
            throw runtime_error("Start of a feature is greater than its end. Please check the input file!!\n");
        }
    }

    return true;
}

#endif
