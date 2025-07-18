#include "spec_reader_utils.h"

using namespace std;
using namespace H5;

extern int GLOBAL_verbosity;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      if (!item.empty()){
        *(result++) = item;
      }

    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



std::map<std::string, vector<float>> prepareMotorConfig(const std::vector <std::string>&motorLog, const std::vector <std::string>&motorVal){
    std::map<std::string, vector<float>> motor;
    std::vector <std::string> motorConfigList;
    std::vector <std::string> motorValList;
    // std::cout << motorLog.size() << std::endl;
    for (uint vecID=0; vecID<motorLog.size(); vecID++){
        // std::cout << motorLog[vecID] << std::endl;
        istringstream iss(motorLog[vecID]);
        motorConfigList = split(motorLog[vecID], ' ');
        iss.str(std::string());
        iss.clear();
        istringstream issVal(motorVal[vecID]);
        motorValList = split(motorVal[vecID], ' ');
        issVal.str(std::string());
        issVal.clear();

        for (uint motorID=1; motorID<motorConfigList.size(); motorID++){
            motor[motorConfigList[motorID]].push_back(std::stof(motorValList[motorID]));
        }
    }

    return motor;
}

void prepareMonitor(const std::vector<std::string>& monitorLogName, const std::vector <std::vector <std::string>>& monitorLog, std::map<std::string, std::vector<float>>& motor, std::map<std::string, vector<float>>& monitor){
    std::vector <float> monitorEntry;

    for (uint nameID=1; nameID<monitorLogName.size(); nameID++){

        monitorEntry.clear();
        // monitorEntry.reserve(monitorLog.size());
        for (uint vecID=0; vecID<monitorLog.size(); vecID++){
            // monitorEntryArray[vecID] = std::stof(monitorLog[vecID][nameID-1]);
            try{
                monitorEntry.push_back(std::stof(monitorLog[vecID][nameID-1]));
            } catch(...){
                throw std::runtime_error("Failed to convert monitorLog");
            }
            
        }
        // cout << monitorLogName[nameID] << endl;

        if ( motor.find(monitorLogName[nameID]) == motor.end()){
            monitor[monitorLogName[nameID]] = monitorEntry;
        }
        else {

            motor[monitorLogName[nameID]] = monitorEntry;
        }


    }

return;
}


void verbose(const int lvl, const std::string& out){
    if (lvl <= GLOBAL_verbosity){
			std::cout << out << std::endl;
    }
}

void appendUnits(H5::DataSet dset, const std::string& units) {
    H5std_string ATTR_NAME("units");
    const char* s[1] = {units.c_str()};
    hsize_t dims1[] = {1};
    DataSpace sid1(1, dims1);
    StrType tid1(0, H5T_VARIABLE);
    Attribute attr = dset.createAttribute(ATTR_NAME, tid1, sid1);
    attr.write(tid1, (void*)s);
}

DataSet writeH5Data(const int& RANK, int & dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, float* data) {
    H5std_string DATASET_NAME(dsetName);

    hsize_t dimsf[1];
    dimsf[0] = dimsfSz;
    
    DataSpace dataspace(RANK, dimsf);
    IntType datatype(PredType::NATIVE_FLOAT);
    datatype.setOrder(H5T_ORDER_LE);
    DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
    dataset.write(data, PredType::NATIVE_FLOAT);
    return dataset;
}

DataSet writeH5Data(const int& RANK, int & dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, std::string data) {
    H5std_string DATASET_NAME(dsetName);

    hsize_t dimsf[] = {(hsize_t)dimsfSz};
    const char* s[1] = {data.c_str()};
    DataSpace dataspace(1, dimsf);
    StrType datatype(0, H5T_VARIABLE);
    DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
    dataset.write((void*)s, datatype);
    return dataset;
}