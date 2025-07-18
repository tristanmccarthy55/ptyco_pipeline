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
#include "specRead.h"
// #include "xmlParser.cpp"

using namespace std;
using namespace tinyxml2;
using namespace H5;

namespace fs = std::experimental::filesystem;

int GLOBAL_verbosity = 0;
fs::path binaryName;


int main(int argc, char** argv) {
    binaryName = fs::path(std::string(argv[0]));
    globalArgs_t globalArgs = ProcessArgs(argc, argv);

    if (globalArgs.scanNr.size() == 0) {
        std::cout << "The scan number has to be specified!" << std::endl;
        display_usage();
        return 1;
    }

    std::string specFile = globalArgs.source;
    std::ifstream datFile(specFile);  //("specES1_started_2016_06_28_0954.dat");
    std::vector<std::string> inputFile;
    std::string line;

    // std::string scanNr_str;
    // std::string startString;
    // std::string endString;
    // std::string motorString;
    // bool readStatus = false;
    // bool lineRead = false;

    datFile.seekg(0);
    while (datFile.peek() != EOF) {
        getline(datFile, line);
        // if (line.length()>0){
        inputFile.push_back(line);
        // }

        // lineRead = false;
        // // inputFile.push_back(line);
        // for (uint scansID=0; (scansID < globalArgs.scanNr.size() && !lineRead); scansID++){
        // 	scanNr_str = std::to_string(globalArgs.scanNr[scansID]);
        // 	startString = "#S " + scanNr_str;
        // 	endString = "#X " + scanNr_str;
        // 	motorString = "#O";
        //
        // 	if (line.compare(0,motorString.size(),motorString)==0 && !lineRead){
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==true && line.compare(0,endString.size(),endString)==0 && !lineRead){
        // 		inputFile.push_back(line);
        // 		readStatus=false;
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==true && !lineRead){
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==false && line.compare(0,startString.size(),startString)==0 && !lineRead){
        // 		readStatus=true;
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // }
    }
    // cout << line << endl;

    std::vector<std::unique_ptr<specRead>> storage;

    storage.reserve(globalArgs.scanNr.size() + 1);
    // std::cout << storage.max_size() << std::endl;

    if (datFile.is_open()) {
#pragma omp critical
        for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
            storage.push_back(std::make_unique<specRead>(globalArgs.scanNr[scanNrIdx], datFile));
            // std::cout << globalArgs.scanNr[scanNrIdx] << std::endl;
        }
        if (globalArgs.orchestra.length() > 0) {
#pragma omp parallel for
            for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
                storage[scanNrIdx]->read_from_file(inputFile);
                storage[scanNrIdx]->read_from_orchestra(globalArgs.orchestra);
                storage[scanNrIdx]->output_stream(globalArgs.hdf5, globalArgs.output, scanNrIdx, globalArgs.scanNr.size(), globalArgs.xmlLayoutFile);
            }
        } else {
#pragma omp parallel for
            for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
                storage[scanNrIdx]->read_from_file(inputFile);
                storage[scanNrIdx]->output_stream(globalArgs.hdf5, globalArgs.output, scanNrIdx, globalArgs.scanNr.size(), globalArgs.xmlLayoutFile);
            }
        }
        if (!globalArgs.hdf5 && globalArgs.output.length() == 0) {
            // send to terminal
            // #pragma omp critical
            for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
                std::cout << storage[scanNrIdx]->jsonOut;
            }
        } 

    } else {
        cout << "Unable to open file" << endl;
    }
    datFile.close();
}
