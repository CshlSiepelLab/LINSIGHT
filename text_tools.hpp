#ifndef TEXT_TOOLS_HPP
#define TEXT_TOOLS_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdint>

using namespace std;

namespace hui{
    class TextTools{
        public:
            static double string2double(string str){
                return string2double(str.c_str());
            }
            
            static double string2double(const char str[]){
                stringstream ss(str);
                double y;
                ss >> y;
                if (ss.fail()){
                    throw runtime_error("Conversion to double type failed!\n");
                }
                return y;
            }

            static long string2long(string str){
                return string2long(str.c_str());
            }

            static long string2long(const char str[]){
                stringstream ss(str);
                long y;
                ss >> y;
                if (ss.fail()){
                    throw runtime_error("Conversion to double type failed!\n");
                }
                return y;
            }

            static long string2uint64_t(const char str[]){
                stringstream ss(str);
                uint64_t y;
                ss >> y;
                if (ss.fail()){
                    throw runtime_error("Conversion to double type failed!\n");
                }
                return y;
            }

            static long string2uint64_t(string str){
                return string2uint64_t(str.c_str());
            }

            static string long2string(long x){
                ostringstream ss;
                ss << x;
                return ss.str();
            }

            static vector<string> split(const string &str, char delim){
                stringstream ss(str);
                string item;
                vector<string> tokens;
                while(std::getline(ss, item, delim)) {
                    tokens.push_back(item);
                }
                return tokens;
            }

            static vector<string> splitByWhiteSpace(const string &str){
                stringstream ss(str);
                string item;
                string token;
                vector<string> tokens;
                while(std::getline(ss, item, '\t')) {
                    stringstream iss(item);
                    while(std::getline(iss, token, ' ')) {
                        tokens.push_back(token);
                    }
                }
                return tokens;
            }
    };
}

#endif
