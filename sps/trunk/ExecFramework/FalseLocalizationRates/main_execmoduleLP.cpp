//
//  SpsModuleExecutor - stand alone execution driver for SpSModules
//
#include "CommandLineParser.h"
#include "ExecModuleFactoryLP.h"
#include "Logger.h"
#include "ParameterList.h"
#include "ParallelThreadedExecution.h"

#include "Logger.h"

#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <time.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
string getCurrentTimeString(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc < 3)
  {
    cerr << "Usage: main_execmodule moduleName parametersFileName [options]" << endl;
    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("t", "NUM_THREADS", true));
  listOptions.push_back(CommandLineParser::Option("lf", "LOG_FILE_NAME", true));
  listOptions.push_back(CommandLineParser::Option("ll", "LOG_LEVEL", true));

  // Command line parameters for ExecSpectraExtraction for CCMS Workflows
  listOptions.push_back(CommandLineParser::Option("ccms_output_library", "OUTPUT_MGF", true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_library", "EXISTING_LIBRARY_MGF", true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_spectradir", "SPECTRA_DIR", true));
  listOptions.push_back(CommandLineParser::Option("ccms_results_dir", "RESULTS_DIR", true));
  listOptions.push_back(CommandLineParser::Option("ccms_newresults_dir", "NEWLIBRARYRESULTS_DIR", true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_spectradatafile", "INPUTSPECTRA_DATAFILE", true));
  listOptions.push_back(CommandLineParser::Option("ccms", "CCMS", true));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parser_error;
  if (!clp.validate(parser_error))
  {
    ERROR_MSG(parser_error);
    return -2;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  // Load the parameter file
  ParameterList ip;
  std::string paramsFileName = argv[2];
  std::string::size_type idx = paramsFileName.rfind(".");
  if( idx != std::string::npos &&
      paramsFileName.substr(idx + 1) == "xml"){
    if (!ip.readFromProteosafeXMLFile(argv[2]))
    {
        ERROR_MSG("Can not read parameter file [" << argv[2] << "].");
        return -3;
    }
  }
  else{
    if (!ip.readFromFile(argv[2]))
    {
        ERROR_MSG("Can not read parameter file [" << argv[2] << "].");
        return -3;
    }
  }

  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  ip.addList(commandLineParams, true);

  int logLevel = ip.getValueInt("LOG_LEVEL", 0);
  if (ip.exists("LOG_FILE_NAME"))
  {
    string logFileName = ip.getValue("LOG_FILE_NAME");
    Logger::setDefaultLogger(Logger::getLogger(logFileName, logLevel));
  }
  else
  {
    Logger::setDefaultLogger(Logger::getLogger(logLevel));
  }

  DEBUG_TRACE;
  // Fill the factory with all known module types
  ExecModuleFactoryLP::RegisterAllModules();

  DEBUG_TRACE;
  // Get the exemplar from the factory
  ExecBase * moduleExemplar = ExecModuleFactoryLP::getModule(argv[1]);

  if (moduleExemplar == 0)
  {
    ERROR_MSG("Module name [" << argv[1] << "] not found in factory.");
    ExecModuleFactoryLP::cleanup();
    return -4;
  }

  DEBUG_TRACE;
  // Clone the exemplar and create a real module with the desired parameters
  ExecBase * moduleReal = moduleExemplar->clone(ip);

  string jobDoneFile(moduleReal->getName());
  jobDoneFile += "_";
  string numNode = moduleReal->m_params.getValue("NUM_SPLIT");
  jobDoneFile += numNode;
  jobDoneFile += "_results.param";

  DEBUG_VAR(jobDoneFile);

  // Validate all the parameters before invoking
  std::string errorString;
  if (!moduleReal->validateParams(errorString))
  {
    ERROR_MSG(errorString);
    string errorString2("Invalid parameters [");
    errorString2 += argv[2];
    errorString2 += "] to module of type [";
    errorString2 += argv[1];
    errorString2 += "]";
    ERROR_MSG("Invalid parameters [" << argv[2] << "] to module of type ["
        << argv[1] << "].");
    delete moduleReal;
    ExecModuleFactoryLP::cleanup();
    return -3;
  }

  DEBUG_TRACE;
  // Load all the input data from the specified files
  if (!moduleReal->loadInputData())
  {
    ERROR_MSG("Loading input data for module of type [" << argv[1] << "]");
    delete moduleReal;
    ExecModuleFactoryLP::cleanup();
    return -4;
  }

  DEBUG_TRACE;
  
  bool returnStatus;
  if(!ip.exists("NUM_THREADS") || ip.getValueInt("NUM_THREADS") <= 1)
  {
    DEBUG_TRACE;
    returnStatus = moduleReal->invoke();
    DEBUG_VAR(returnStatus);
  }
  else
  {
    int numThreads = ip.getValueInt("NUM_THREADS");
    DEBUG_VAR(numThreads);
    //returnStatus = moduleReal->invoke();
    ParallelThreadedExecution exec(moduleReal);
    returnStatus = exec.invoke(numThreads);
    DEBUG_VAR(returnStatus);
  }
    
  if (!returnStatus)
  {
    ERROR_MSG("Invoking module of type [" << argv[1] << "]");
    delete moduleReal;
    ExecModuleFactoryLP::cleanup();
    return -5;
  }

  DEBUG_TRACE;
  // Save any result data
  if (!moduleReal->saveOutputData())
  {
    ERROR_MSG("Saving output data for module of type [" << argv[1] << "]");
    delete moduleReal;
    ExecModuleFactoryLP::cleanup();
    return -6;
  }

  DEBUG_TRACE;

  delete moduleReal;
  ExecModuleFactoryLP::cleanup();
  return 0;
}

