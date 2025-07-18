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

// #include "help.h"
// #include "specRead.h"
// #include "spec_reader_utils.h"
// #include "tinyxml2.h"
#include "xmlParser.h"
#include "exprtk.cpp"

using namespace std;
using namespace tinyxml2;
using namespace H5;

xmlParser::xmlParser(){};

template <typename h5GroupData>
void xmlParser::writeAttributes(tinyxml2::XMLElement* element, std::string indent, h5GroupData group) {
    for (tinyxml2::XMLElement* child = element->FirstChildElement("attribute"); child != 0; child = child->NextSiblingElement("attribute")) {
        if (child->ToElement()->Attribute("name")) {
            std::string attrName = child->ToElement()->Attribute("name");
            verbose(3, (indent + "Writing attribute " + attrName + "."));
            if (child->ToElement()->Attribute("source")) {
                // std::cout << ", source = " << child->ToElement()->Attribute("source");

                // std::cout << std::endl;
                try {
                    H5std_string ATTR_NAME(attrName);
                    if (child->ToElement()->Attribute("source")) {
                        std::string source = child->ToElement()->Attribute("source");
                        if (source.size() == 4 && source.compare(0, 4, "spec") == 0) {
                        } else if (source.size() == 8 && source.compare(0, 8, "constant") == 0) {
                            const char* s[1] = {child->ToElement()->Attribute("value")};
                            hsize_t dims1[] = {1};
                            DataSpace sid1(1, dims1);
                            StrType tid1(0, H5T_VARIABLE);
                            Attribute attr = group.createAttribute(ATTR_NAME, tid1, sid1);
                            attr.write(tid1, (void*)s);
                        }
                    }
                } catch (...) {
                    verbose(1, (indent + "Failed to write attribute " + attrName + "."));
                }
            }
        }
    }
}


void xmlParser::parseXML(tinyxml2::XMLElement* element, std::string indent, H5::Group group, std::map<std::string, std::vector<float>>& motor, std::map<std::string, std::vector<float>>& monitor, std::map<std::string, std::vector<std::string>>& metadataContainer) {
    indent += "\t";
    // write attributes
    writeAttributes(element, indent, group);

    // write datasets
    for (tinyxml2::XMLElement* child = element->FirstChildElement("dataset"); child != 0; child = child->NextSiblingElement("dataset")) {
        if (child->ToElement()->Attribute("name")) {
            std::string dsetName = child->ToElement()->Attribute("name");
            // std::cout << indent + "dataset name = " << dsetName;
            verbose(3, (indent + "Writing dataset " + dsetName + "."));
            H5std_string DATASET_NAME(dsetName);
            const int RANK = 1;

            if (child->ToElement()->Attribute("source")) {
                std::string source = child->ToElement()->Attribute("source");
                if (source.size() == 4 && source.compare(0, 4, "spec") == 0) {
                    try{
                    std::string specVal;
                    bool writeData = true;
                    if (child->ToElement()->Attribute("source")) {
                        specVal = child->ToElement()->Attribute("specval");
                    } else {
                        throw std::runtime_error("Expected attribute 'specval'.");
                    }
                    verbose(3, (indent + "Loading from spec"));
                    DataSet dataset;
                    std::string dtype = "float";

                    if (motor.find(specVal) != motor.end()) {
                        // found entry in motor
                        // std::cout << "Motor Value: " << motor[specVal][0] << std::endl;
                        verbose(3, (indent + "Found motor entry in spec for " + dsetName + " (" + specVal + ")."));
                        float* data;
                        data = &(motor[specVal][0]);
                        int dimsfSz = motor[specVal].size();
                        dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                    } else if (monitor.find(specVal) != monitor.end()) {
                        // found entry in monitor
                        verbose(3, (indent + "Found monitor entry in spec for " + dsetName + " (" + specVal + ")."));
                        // std::cout << "Monitor Value: " << monitor[specVal].size() << std::endl;
                        float* data;
                        data = &(monitor[specVal][0]);
                        int dimsfSz = monitor[specVal].size();
                        dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                    } else if (metadataContainer.find(specVal) != metadataContainer.end()) {
                        verbose(3, (indent + "Found metadata entry in spec for " + dsetName + " (" + specVal + ")."));
                        int dimsfSz;
                        dtype = metadataContainer[specVal][1];
                        if (dtype.compare(0,6,"string")==0){
                            dimsfSz = 1;
                            dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, metadataContainer[specVal][3]);
                        } else if ((dtype.compare(0,5,"float")==0) || ((dtype.compare(0,3,"int"))==0)){
                            float* data;
                            std::vector<float> metadataEntry;
                            uint numEntries = std::stoi(metadataContainer[specVal][2]);
                            dimsfSz = numEntries;
                            for(uint ii=0; ii<numEntries; ii++){
                                metadataEntry.push_back(std::stof(metadataContainer[specVal][3+ii]));
                            }
                            data = &(metadataEntry[0]);
                            // std::cout << metadataEntry[0] << std::endl;
                            dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                        }

                    } else {
                        writeData = false;
                        verbose(1, (indent + "Did not find any entry for " + dsetName + " (" + specVal + ")."));
                    }
                    if (writeData) {
                        xmlParser::writeAttributes(child, indent, dataset);
                        if (child->ToElement()->Attribute("units")) {
                            verbose(3, (indent + "Appending units as attribute."));
                            appendUnits(dataset, child->ToElement()->Attribute("units"));
                        }
                    }
                    } catch (...) {

                    }

                } else if (source.size() == 8 && source.compare(0, 8, "constant") == 0) {
                    verbose(3, (indent + "Parsing static value"));
                    try {
                    std::string dtype;
                    if (child->ToElement()->Attribute("type")) {
                        dtype = child->ToElement()->Attribute("type");
                    } else {
                        throw std::runtime_error("Expected attribute 'type'.");
                    }

                    if (dtype.compare(0, 6, "string") == 0) {
                        const char* s[1] = {child->ToElement()->Attribute("value")};
                        hsize_t dims1[] = {1};
                        DataSpace sid1(1, dims1);
                        StrType tid1(0, H5T_VARIABLE);
                        DataSet dataset = group.createDataSet(DATASET_NAME, tid1, sid1);
                        dataset.write((void*)s, tid1);
                        xmlParser::writeAttributes(child, indent, dataset);
                        if (child->ToElement()->Attribute("units")) {
                            verbose(3, (indent + "Appending units as attribute."));
                            appendUnits(dataset, child->ToElement()->Attribute("units"));
                        }
                    } else if (dtype.compare(0, 5, "float") == 0) {
                        float val = std::stof(child->ToElement()->Attribute("value"));
                        float* data = &val;
                        hsize_t dimsf[1];
                        dimsf[0] = 1;
                        DataSpace dataspace(RANK, dimsf);
                        IntType datatype(PredType::NATIVE_FLOAT);
                        datatype.setOrder(H5T_ORDER_LE);
                        DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
                        dataset.write(data, PredType::NATIVE_FLOAT);
                        xmlParser::writeAttributes(child, indent, dataset);
                        if (child->ToElement()->Attribute("units")) {
                            verbose(3, (indent + "Appending units as attribute."));
                            appendUnits(dataset, child->ToElement()->Attribute("units"));
                        }
                    } else {
                        throw std::runtime_error("Expected value of type 'string' or 'float'.");
                    }
                    } catch (...) {
                    // throw std::runtime_error("Failed to parse static value");
                    }
                } else if (source.size() == 8 && source.compare(0, 8, "relative") == 0) {
                    std::string expressionString = child->ToElement()->Attribute("expression");
                    typedef exprtk::symbol_table<float> symbol_table_t;
                    typedef exprtk::expression<float> expression_t;
                    typedef exprtk::parser<float> parser_t;
                    std::string specVal;
                    bool useDefault = false;
                    float* data;

                    if (child->ToElement()->Attribute("specval")) {
                        specVal = child->ToElement()->Attribute("specval");
                    } else {
                        throw std::runtime_error("Expected attribute 'specval'.");
                    }
                    
                    symbol_table_t symbol_table;
                    if (motor.find(specVal) != motor.end()) {
                        symbol_table.add_variable(specVal, motor[specVal][0]);
                    } else if (monitor.find(specVal) != monitor.end()) {
                        symbol_table.add_variable(specVal, monitor[specVal][0]);
                    } else if (metadataContainer.find(specVal) != metadataContainer.end()) {
                        float val = std::stof(metadataContainer[specVal][3]);
                        symbol_table.add_variable(specVal, val);
                    } else {
                        useDefault = true;
                    }

                    if (!useDefault){
                        symbol_table.add_constants();

                        expression_t expression;
                        expression.register_symbol_table(symbol_table);

                        parser_t parser;
                        parser.compile(expressionString, expression);
                        float val = expression.value();
                        data = &val;
                    } else {
                        float val = std::stof(child->ToElement()->Attribute("default"));
                        data = &val;
                    }
                    try{
                    hsize_t dimsf[1];
                    dimsf[0] = 1;
                    DataSpace dataspace(RANK, dimsf);
                    IntType datatype(PredType::NATIVE_FLOAT);
                    datatype.setOrder(H5T_ORDER_LE);
                    DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
                    dataset.write(data, PredType::NATIVE_FLOAT);
                    xmlParser::writeAttributes(child, indent, dataset);
                    if (child->ToElement()->Attribute("units")) {
                        verbose(3, (indent + "Appending units as attribute."));
                        appendUnits(dataset, child->ToElement()->Attribute("units"));
                    }
                    } catch (...){

                    }
                } else {
                    verbose(1, ("Unknown source type " + source));
                }
                // std::cout << ", source = " << child->ToElement()->Attribute("source");
            } else {
                throw std::runtime_error("Expected attribute 'source'.");
            }
        }
    }

    // write groups
    for (tinyxml2::XMLElement* child = element->FirstChildElement("group"); child != 0; child = child->NextSiblingElement("group")) {
        if (child->ToElement()->Attribute("name")) {
            std::string groupName = child->ToElement()->Attribute("name");
            if (child->ToElement()->Attribute("dependency")){
                if (H5Lexists(group.getId(), child->ToElement()->Attribute("dependency"), H5P_DEFAULT) > 0) {
                    verbose(2, (indent + "Writing group " + groupName + "."));
                    Group subgroup;
                    if (H5Lexists(group.getId(), child->ToElement()->Attribute("name"), H5P_DEFAULT) > 0) {
                        // std::cout << "group exists" << std::endl;
                        Group tmpsubgroup(group.openGroup(child->ToElement()->Attribute("name")));
                        subgroup = tmpsubgroup;
                    } else {
                        Group tmpsubgroup(group.createGroup(child->ToElement()->Attribute("name")));
                        subgroup = tmpsubgroup;
                    }
                    // Group subgroup(group.createGroup(child->ToElement()->Attribute("name")));
                    parseXML(child, indent, subgroup, motor, monitor, metadataContainer);
                } else {
                    continue;
                }
            } else {
                verbose(2, (indent + "Writing group " + groupName + "."));
                Group subgroup;
                if (H5Lexists(group.getId(), child->ToElement()->Attribute("name"), H5P_DEFAULT) > 0) {
                    // std::cout << "group exists" << std::endl;
                    Group tmpsubgroup(group.openGroup(child->ToElement()->Attribute("name")));
                    subgroup = tmpsubgroup;
                } else {
                    Group tmpsubgroup(group.createGroup(child->ToElement()->Attribute("name")));
                    subgroup = tmpsubgroup;
                }
                // Group subgroup(group.createGroup(child->ToElement()->Attribute("name")));
                parseXML(child, indent, subgroup, motor, monitor, metadataContainer);
            }
        }

        if (child->ToElement()->Attribute("spec_default")) {
            // spec_default defines the location where we dump the rest, that is any undefined data
            // std::cout << ", spec_default = " << child->ToElement()->Attribute("spec_default");
            tinyxml2::XMLElement* currentChild = child;
            std::string tmpString;
            tmpString = currentChild->ToElement()->Attribute("name");
            specDefaultPath = "/" + tmpString + specDefaultPath;
            while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
                tmpString = currentChild->Parent()->ToElement()->Attribute("name");
                specDefaultPath = "/" + tmpString + specDefaultPath;
                currentChild = currentChild->Parent()->ToElement();
            }
        }
        // std::cout << std::endl;
    }

    // create hardlink 
    for (tinyxml2::XMLElement* child = element->FirstChildElement("hardlink"); child != 0; child = child->NextSiblingElement("hardlink")) {
        std::string linkName;
        if (child->ToElement()->Attribute("name")) {
            linkName = child->ToElement()->Attribute("name");
            verbose(2, (indent + "Creating hardlink for " + linkName + "."));
        }
        std::string targetString;
        if (child->ToElement()->Attribute("target")) {
            targetString = child->ToElement()->Attribute("target");
        } else {
            throw std::runtime_error("Expected attribute 'target'.");
        }
        
        string linkPath;
        string tmpString;
        linkPath = "/" + linkName;
        tinyxml2::XMLElement* currentChild = child;
        while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
            tmpString = currentChild->Parent()->ToElement()->Attribute("name");
            linkPath = "/" + tmpString + linkPath;
            currentChild = currentChild->Parent()->ToElement();
        }
        try{
            group.link(H5L_TYPE_HARD, targetString, linkPath);
        } catch (...) {
            verbose(1, (indent + "Failed to create hardlink " + linkName + "."));
        }
    }
    // create softlink 
    for (tinyxml2::XMLElement* child = element->FirstChildElement("softlink"); child != 0; child = child->NextSiblingElement("softlink")) {
        std::string linkName;
        if (child->ToElement()->Attribute("name")) {
            linkName = child->ToElement()->Attribute("name");
            verbose(2, (indent + "Creating softlink for " + linkName + "."));
        }
        std::string targetString;
        if (child->ToElement()->Attribute("target")) {
            targetString = child->ToElement()->Attribute("target");
        } else {
            throw std::runtime_error("Expected attribute 'target'.");
        }
        
        string linkPath;
        string tmpString;
        linkPath = "/" + linkName;
        tinyxml2::XMLElement* currentChild = child;
        while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
            tmpString = currentChild->Parent()->ToElement()->Attribute("name");
            linkPath = "/" + tmpString + linkPath;
            currentChild = currentChild->Parent()->ToElement();
        }
        try{
            group.link(H5L_TYPE_HARD, targetString, linkPath);
        } catch (...) {
            verbose(1, (indent + "Failed to create softlink " + linkName + "."));
        }
    }
}