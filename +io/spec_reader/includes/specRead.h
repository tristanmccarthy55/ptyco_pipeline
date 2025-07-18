#include <getopt.h>
#include <omp.h>
#include <sys/stat.h>
#include <chrono>
#include <fstream>
#include <future>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>

#include "help.h"
#include "spec_reader_utils.h"
#include "tinyxml2.h"


#ifndef SPECREAD_H
#define SPECREAD_H
class specRead {
    uint lineidx;
    int start;
    std::string line;
    std::string scanNr_str;
    std::string startString;  //("#S 5398");
    std::string endString;    //("#X 5398");
    std::string startScanString;
    std::string motorString;
    std::vector<std::vector<std::string>> scanContent;
    std::vector<std::string> scanContentName;
    std::vector<std::string> motorConfig;
    std::vector<std::string> motorUpdate;
    std::vector<std::string> metadata;
    std::future<std::map<std::string, std::vector<float>>> motorRes;


   public:
    int scan;
    std::string jsonOut;

    std::vector<std::string> header;
    std::map<std::string, std::vector<float>> motor;
    std::map<std::string, std::vector<float>> monitor;
    std::map<std::string, std::vector<std::string>> metadataContainer;

    specRead(int& scanNr, std::ifstream& datFile);
    void read_from_file(const std::vector<std::string>& inputFile);
    void read_from_orchestra(const std::string& orchestra_path);
    void output_stream(const bool hdf5, const std::string output, const uint& scanNrIdx, const uint scanNrSz, const std::string& xmlLayoutFile);

};
#endif
