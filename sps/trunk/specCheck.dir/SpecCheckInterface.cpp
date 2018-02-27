///////////////////////////////////////////////////////////////////////////////
#include "SpecCheckInterface.h"


#include "CommandLineParser.h"
#include "ParameterList.h"
#include "utils.h"
#include "Logger.h"
#include "mzxml.h"


#include "copyright.h"

#define DEFAULT_INPUT_FILES_LIST      "spectra/pklbin_files.txt"
#define DEFAULT_ORIGINAL_FILES_LIST   "spectra/input_index.txt"
#define DEFAULT_OUTFILE_NAME          "check.txt"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecCheckInterface::SpecCheckInterface() :
  outFileName(DEFAULT_OUTFILE_NAME),
  inputFilesNames(DEFAULT_INPUT_FILES_LIST)
{
}

SpecCheckInterface::~SpecCheckInterface()
{
}

///////////////////////////////////////////////////////////////////////////////
int  SpecCheckInterface::loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles)
{
  std::vector<std::string> file;
  std::string line;
  file.clear();
  std::ifstream infile (inputFilesNames.c_str(), std::ios_base::in | std::ios_base::binary);

  if(!infile) return 0;

  while (getline(infile, line, '\n')) {
    inputFiles.push_back (line);
  }

  //cout << "for: " << inputFilesNames << " got:" << endl;
  //for(int i = 0 ; i < inputFiles.size() ; i++)
  //  cout << inputFiles[i] << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::loadSpecsetPklbin(string &inputFilesName, SpecSet &specset)
{
  if(specset.LoadSpecSet_pklbin(inputFilesName.c_str()) <= 0) {
    cerr << "Unable to load file: " << inputFilesName << endl;
    return 0;
  }
  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::loadSpecsetMzxml(string &inputFilesName, SpecSet &specset)
{
  // auxilizary vector needed
  vector<short> msLevel;
  int ret = LoadMzxml( (char * const)(inputFilesName.c_str()), specset, & msLevel, 2);
  if(! ret) {
    cerr << "Error loading mzxml file " << inputFilesName << endl;
    return 0;
  }
  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::loadSpecsetMgf(string &inputFilesName, SpecSet &specset)
{
  int fileLoaded = specset.LoadSpecSet_mgf(inputFilesName.c_str());
  if(!fileLoaded) {
    cout << "Error loading mgf file: " << inputFilesName;
    return 0;
  }
  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::loadSpecset(string &inputFilesName, SpecSet &specset)
{
  //cout << "Loading file " << inputFilesName << endl;

  // get extension
  string aux = inputFilesName.substr(inputFilesName.find_last_of(".") + 1);
  // change to lower case
  std::transform(aux.begin(), aux.end(), aux.begin(), ::tolower);

  if(aux == "mzxml")
    return loadSpecsetMzxml(inputFilesName, specset);

  if(aux == "pklbin")
    return loadSpecsetPklbin(inputFilesName, specset);

  if(aux == "mgf")
    return loadSpecsetMgf(inputFilesName, specset);

  cout << "   unknown file extension " << aux << endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::processOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------


  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));

  listOptions.push_back(CommandLineParser::Option("-input-list",      "INPUT_LIST",       true));
  listOptions.push_back(CommandLineParser::Option("-inputspectra",    "INPUT_SPECS_MS",   true));

  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",               "PARAMETER_FILE",   true));


  ////////////////////////////////////////////////////////////////////////////////
  // Execute the command line parser
  CommandLineParser clp(argc, argv, 0, listOptions);

  ////////////////////////////////////////////////////////////////////////////////
  // Checking for errors
  string parser_error;
  if (!clp.validate(parser_error)) {
    ERROR_MSG(parser_error);
    return -1;
  }


  ////////////////////////////////////////////////////////////////////////////////
  // The parameters' values
  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (commandLineParams.exists("PARAMETER_FILE")) {

    string parameterFilename = commandLineParams.getValue("PARAMETER_FILE");

    ParameterList ip;
    ip.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag not set)
    commandLineParams.addList(ip, false);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);



  ////////////////////////////////////////////////////////////////////////////////
  // File load section

  // input spectra file to load
  if (commandLineParams.exists("INPUT_LIST")) {
    inputFilesNames = commandLineParams.getValue("INPUT_LIST").c_str();

    // load input files names
    if(!loadStringVector(inputFilesNames, inputFiles)) {
      cerr << "Unable to load file: " << inputFilesNames << endl;
      exit(0);
    }
  }


  if (commandLineParams.exists("INPUT_SPECS_MS")) {
    string aux = commandLineParams.getValue("INPUT_SPECS_MS").c_str();
    splitText(aux.c_str(), inputFiles, ";|");
  }


  // Keep record of loaded files
  int filesLoaded = 0;

  cout  << "Loading " << inputFiles.size() << " files." << endl;

  for(int i = 0 ; i < inputFiles.size() ; i++) {
    SpecSet spec;
    if(loadSpecset(inputFiles[i], spec)) {
      inputData.push_back(spec);
      filesLoaded++;
    }
  }

  cout <<  filesLoaded << " files loaded." << endl;

  // exit if no file was loaded
  if(!filesLoaded) {
    cerr << "No file loaded." << endl;
    return 0;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  Processing section

  // check for zero PM

  cout << "Checking " << inputData.size() << " specsets." << endl;
  for(int i = 0 ; i < inputData.size() ; i++) {
    cout << "  #" << i << " Checking " << inputData[i].size() << " spectra." << endl;
    vector<int> aux;
    for(int j = 0 ; j < inputData[i].size() ; j++) {
      if(inputData[i][j].parentMass < 2.0) {
        aux.push_back(j);
        cout << "found: " << i << " ; " << j << " -> " << inputData[i][j].parentMass << endl;
      } //else
      //  cout << "value: " << i << " ; " << j << " -> " << inputData[i][j].parentMass << " - " << inputData[i][j].parentMZ << endl;
    }
    if(aux.size())
      invalid.insert(pair<int, vector<int> >(i, aux));
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  output section

  map<int, vector<int> >::iterator it;
  for(it = invalid.begin() ; it != invalid.end() ; it++) {
    int fileIdx = it->first;
    string fn = inputFiles[fileIdx];

    cout << "------- " << fn << " -------" << endl;

    for(int i = 0 ; i > it->second.size() ; i++) {
      if(!i) cout << " ; ";
      cout << it->second[i];
    }
    cout << endl;
  }


   // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: specCheck [OPTION]\n";
  sout << "Options:\n";
  sout << "  --input-list   FILE         Input spectra list file (txt format)\n";
  sout << "  --inputspectra FILE         Input spectra list, separated by |\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::version(ostream & sout)
{
  sout << "specCheck 1.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecCheckInterface::error(const string & a)
{
  cerr << "specCheck: invalid or missing option: " << a << endl;

  cerr << "Type 'specCheck --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
