#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <future>
#include "tinyxml2.h"
#include "H5Cpp.h"


template<typename Out>
void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);

std::map<std::string, std::vector<float>> prepareMotorConfig(const std::vector <std::string>&motorLog, const std::vector <std::string>&motorVal);

void prepareMonitor(const std::vector<std::string>& monitorLogName, const std::vector <std::vector <std::string>>& monitorLog, std::map<std::string, std::vector<float>>& motor, std::map<std::string, std::vector<float>>& monitor);

void verbose(const int lvl, const std::string& out);


void appendUnits(H5::DataSet dset, const std::string& units);

H5::DataSet writeH5Data(const int& RANK, int& dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, float* data);
H5::DataSet writeH5Data(const int& RANK, int& dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, std::string data);
