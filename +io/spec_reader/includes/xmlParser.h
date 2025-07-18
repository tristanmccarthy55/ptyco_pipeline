
#ifndef XMLPARSER_H
#define XMLPARSER_H
#include "specRead.h"
#include "tinyxml2.h"
// #include "Cpp.h"

class xmlParser {
   public:
    std::string specDefaultPath;

    xmlParser();
    template <typename h5GroupData>
    void writeAttributes(tinyxml2::XMLElement* element, std::string indent, h5GroupData group);

    void parseXML(tinyxml2::XMLElement* element, std::string indent, H5::Group group, std::map<std::string, std::vector<float>>& motor, std::map<std::string, std::vector<float>>& monitor, std::map<std::string, std::vector<std::string>>& metadataContainer);
};
#endif