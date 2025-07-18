#include <iostream>
#include <experimental/filesystem>
#include <getopt.h>
#include <vector>
#include "spec_reader_utils.h"

extern std::experimental::filesystem::path binaryName;

#ifndef HELP_H
#define HELP_H
struct globalArgs_t {
	std::string source;     /* -s option */
	std::string output;	/* -o option */
	int verbosity;		/* -v option */
	std::vector<int> scanNr;		/* --scanNr option */
	bool hdf5=false;	/* --hdf5 option */
	std::string orchestra;	/* --hdf5 option */
	std::string xmlLayoutFile;	/* --xmlLayout option */
};
#endif /* HELP_H */

void display_usage( void);

globalArgs_t ProcessArgs(int argc, char** argv);
