#include "help.h"

namespace fs = std::experimental::filesystem;

extern int GLOBAL_verbosity;

void display_usage(void) {
    std::cout << std::endl
              << std::endl
              << binaryName.filename() << std::endl
              << std::endl;

    std::cout << "USAGE: " << std::endl;
    std::cout << "\t " + binaryName.string() + " [options]" << std::endl;
    std::cout << "ARGUMENTS: " << std::endl;
    std::cout << "\t -s | --source          SPEC source file" << std::endl;
    std::cout << "\t -o | --output          output file" << std::endl;
    std::cout << "\t -v | --verbosity       set verbosity (integer)" << std::endl;
    std::cout << "\t -h | --help            this usage info" << std::endl
              << std::endl;
    std::cout << "\t --scanNr               specify a single scan number or a range (e.g. 247-320)" << std::endl;
    std::cout << "\t --hdf5                 save output as HDF5" << std::endl;
    std::cout << "\t --orchestra			path to orchestra directory" << std::endl;
    std::cout << "\t --xmlLayout			path to the XML layout file for writing NEXUS files" << std::endl;

    exit(EXIT_FAILURE);
}

globalArgs_t ProcessArgs(int argc, char** argv) {
    globalArgs_t globalArgs;

    const char* const short_opts = "s:o:v:0123:h";
    const option long_opts[] = {
        {"source", required_argument, nullptr, 's'},
        {"output", required_argument, nullptr, 'o'},
        {"verbosity", required_argument, nullptr, 'v'},
        {"scanNr", required_argument, nullptr, '0'},
        {"hdf5", no_argument, nullptr, '1'},
        {"orchestra", required_argument, nullptr, '2'},
        {"xmlLayout", required_argument, nullptr, '3'},
        {"help", no_argument, nullptr, 'h'},
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
            case 's':
                try {
                    globalArgs.source = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'o':
                try {
                    globalArgs.output = fs::path(std::string(optarg));
                    fs::path outdir = globalArgs.output;

                    break;
                } catch (...) {
                    display_usage();
                }

            case 'v':
                try {
                    GLOBAL_verbosity = std::stoi(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case '0':
                try {
                    auto scanRange = split(optarg, ',');
                    for (auto scanID : scanRange) {
                        auto scanValRange = split(scanID, '-');
                        if (scanValRange.size() > 1) {
                            for (int ii = 0; std::stoi(scanValRange[0]) + ii <= std::stoi(scanValRange[scanValRange.size() - 1]); ii++) {
                                globalArgs.scanNr.push_back(std::stoi(scanValRange[0]) + ii);
                            }
                        } else {
                            globalArgs.scanNr.push_back(std::stoi(scanID));
                        }
                    }

                    break;
                } catch (...) {
                    display_usage();
                }
            case '1':
                try {
                    globalArgs.hdf5 = true;
                    break;
                } catch (...) {
                    display_usage();
                }
            case '2':
                try {
                    globalArgs.orchestra = optarg;
                    break;
                } catch (...) {
                    display_usage();
                }
            case '3':
                try {
                    globalArgs.xmlLayoutFile = optarg;
                    break;
                } catch (...) {
                    display_usage();
                }
            case 'h':  // -h or --help
            case '?':  // Unrecognized option
                display_usage();
                break;
            default:
                display_usage();
                break;
        }
    }
    return globalArgs;
}
