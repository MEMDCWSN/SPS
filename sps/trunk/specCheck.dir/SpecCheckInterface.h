///////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_INFO_INTERFACE_H__
#define __CLUSTER_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "SpecSet.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class SpecCheckInterface {

  // input files filenames (pklbin)
  vector<string> inputFiles;

  // input files filenames
  string inputFilesNames;

  // the result vector
  map<int, vector<int> > invalid;

  // output file and directory
  string outFileName;

  // Input data (specsets)
  vector<SpecSet> inputData;



  // load string vector file
  int  loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles);
  // load pklbin file
  int loadSpecsetPklbin(string &inputFilesName, SpecSet &specset);
  // load mzxml file
  int loadSpecsetMzxml(string &inputFilesName, SpecSet &specset);
  // load mgf file
  int loadSpecsetMgf(string &inputFilesName, SpecSet &specset);
  // load specset file
  int loadSpecset(string &inputFilesName, SpecSet &specset);

 public:


  // Constructors and destructor
  SpecCheckInterface();
  ~SpecCheckInterface();

  // Option parsing
  int processOptions(int argc, char **argv);

  // Output help
  int help(ostream &);

  // Output version information
  int version(ostream &);

  // Output error messages
  int error(const string &);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
