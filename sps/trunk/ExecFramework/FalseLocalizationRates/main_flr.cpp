#include "CommandLineParser.h"
#include "ExecModuleFactoryLP.h"
#include "ExecMergeConvert.h"
#include "Logger.h"
#include "ParameterList.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "copyright.h"

#include "Logger.h"
#include "ExecFlr.h"
#include "ExecLpSolver.h"
#include "StatusFile.h"

#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <time.h>

using namespace specnets;
using namespace std;

const string DEFAULT_INPUT_FILE_BASE = "spectra/specs_ms";
const string DEFAULT_INPUT_FILES_LIST = "spectra/pklbin_files.txt";

enum Stage
{
  STAGE_BEGIN = 0,
  STAGE_PAIRSPECTRA = 1,
  STAGE_LPSOLVER = 2,
  STAGE_FLR = 3,
};

// System Defines
#define TEST_VALID {          \
  DEBUG_VAR(isValid);         \
  if (!isValid) {             \
    ERROR_MSG(errorString);   \
    return false;             \
  }                           \
}

#define TEST_RETURN(__module, __var) {                                                    \
  if(__var.size() == 0) {                                                               \
    ERROR_MSG("invoking " << __module << ": empty return set \"" << #__var << "\".");     \
  /*  return false; */                                                                     \
  }                                                                                   \
}

#define TEST_RETURN_STATUS(__module) {                            \
  DEBUG_VAR( returnStatus);                                     \
  if (!returnStatus) {                                          \
    ERROR_MSG("invoking " << __module << " exited in error.");    \
    return false;                                              \
  }                                                             \
}

#define TEST_SAVE_OUPUT_DATA(__module) {                                      \
  DEBUG_VAR(returnStatus);                                                  \
  if (!returnStatus) {                                                      \
    ERROR_MSG("Saving output data for " << __module << " exited in error.");  \
    return false;                                                           \
  }                                                                         \
}

//-----------------------------------------------------------------------------
string getCurrentDirectory(void)
{
  string currentWorkindDir;
  currentWorkindDir = getcwd(NULL, 1024);
  return currentWorkindDir;

//  int dummy;
//  dummy = system("pwd > pwd.tmp");
//  ifstream ifs("pwd.tmp");
//  char buf[1024];
//  ifs.getline(buf, 1024);
//  return buf;
}

//-----------------------------------------------------------------------------
string getProjPath(const ParameterList & pl, const string & addPath)
{
  string projDir = pl.getValue("PROJECT_DIR", "");
  bool relDir = pl.getValueBool("RELATIVE_DIR", true);

  return getPath(projDir, addPath, relDir);
}

//-----------------------------------------------------------------------------
string getCurrentTimeString(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);
}

//-----------------------------------------------------------------------------

void addDefaultParameterValues(ParameterList &p)
{
  // Basic parameters
  p.addIfDoesntExist("TOLERANCE_PEAK", "0.45");
  p.addIfDoesntExist("MAX_MOD_COUNT", "2");
  p.addIfDoesntExist("NORMALIZE_SPECTRUM", "1");
  p.addIfDoesntExist("UNIQUE_MODIFIED_PEPTIDES", "0");
  p.addIfDoesntExist("REINDEX_SCANS", "0");
  p.addIfDoesntExist("NUM_GROUPING_THRESHOLDS",30);
  p.addIfDoesntExist("GROUPING_THRESHOLD_STEPS",.015);


  // Grid parameters

  p.addIfDoesntExist("GRID_TYPE", "sge");
  p.addIfDoesntExist("GRID_NUMNODES", "0");
  p.addIfDoesntExist("GRID_NUMCPUS", "4");
  p.addIfDoesntExist("GRID_SGE_EXE_DIR", "");
  p.addIfDoesntExist("GRID_PARAMS", "-l h_vmem=1G");

  string currDir = getCurrentDirectory();
  p.addIfDoesntExist("PROJECT_DIR", currDir);
  p.addIfDoesntExist("RELATIVE_DIR", "1");
}
// -------------------------------------------------------------------------
/**
 * Loads list of filenames into an ordered collection of SpecSets
 * @param fileList list of filenames
 * @param outputSpectra output collection of spectra
 * @return true if successful, false otherwise
 */
bool loadSpecsMS(const string& exeDir,
                 const string& fileList,
                 SpecSet &outputSpectra)
{
  vector<string> names;
  if (!readFilesFromFile(fileList, names))
  {
    ERROR_MSG("readFilesFromFile() failed for " << fileList);
    return false;
  }

  if (names.size() > 1)
  {
    ERROR_MSG("There should only be a single filtered spectrum file!");
    return false;
  }
  if (!ExecMergeConvert::loadSpecset(exeDir, names[0], &outputSpectra))
  {
    ERROR_MSG("ExecMergeConvert::loadSpecset() failed for " << names[0]);
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
bool loadResultSpectrum(ParameterList &ip,
                        const char * specSetListFile,
                        PeptideSpectrumMatchSet &psms,
                        SpecSet &outputResults,
                        bool resetScans)
{
  vector<vector<string> > spectraFileList;
  map<string, unsigned int> spectraFileListHeader;
  vector<string> requiredHeader;
  requiredHeader.push_back("Path");
  vector<int> requiredHeaderIndex;

  if (!DelimitedTextReader::loadDelimitedFile(specSetListFile,
                                              "\t",
                                              "#",
                                              spectraFileListHeader,
                                              spectraFileList,
                                              requiredHeader,
                                              requiredHeaderIndex))
  {
    ERROR_MSG("Unable to load spectrum files!" << specSetListFile);
    return false;
  }

  vector<int> outputResultIndex(psms.m_psmSet.size(), -1);
  std::tr1::unordered_set<string> passedSpectra;
  outputResults.resize(psms.m_psmSet.size());

  for (int i = 0; i < spectraFileList.size(); i++)
  {
    DEBUG_TRACE;
    //get file path
    string path = spectraFileList[i][requiredHeaderIndex[0]];
    DEBUG_VAR(path);

    FilenameManager spectrumFm(path.c_str());

    DEBUG_VAR(spectrumFm.extension);
    DEBUG_VAR(spectrumFm.filename);

    ParameterList loadParams(ip);
    loadParams.setValue("INPUT_SPECTRA", path);

    DEBUG_VAR(path);

    SpecSet currSet(1);

    ExecMergeConvert* loader = new ExecMergeConvert(loadParams,
        &currSet);

    // let ExecMergeConvert take care of loading spectra
    if (!loader->loadInputData())
    {
      return false;
    }

    // let ExecMergeConvert pre-process spectra (set tolerances, activation, etc)
    if (!loader->invoke())
    {
      ERROR_MSG("Unable to load file! " << path);
      return false;
    }

    currSet.setFilename(spectrumFm.filename);
    delete loader;

    if (resetScans)
    {
      //ONLY BECAUSE INSPECT DOES NOT DEAL WITH SCAN NUMBERS!
      //REMOVE ME!
      DEBUG_MSG("Reindexing scans");
      for (int k = 0; k < currSet.size(); k++)
      {
        currSet[k].scan = k+1;
        currSet[k].msLevel = 2;
      }
    }

    currSet.index();

    for (int j = 0; j < psms.m_psmSet.size(); j++)
    {
      psmPtr currPsm = psms.m_psmSet[j];
      string psmFilename;

      FilenameManager psmFm(currPsm->m_spectrumFile.c_str());

      stringstream key;
      key.str("");
      key.clear();

      key << psms.m_psmSet[j]->m_scanNum << "$" << psmFm.filename;

      if (psmFm.filename.compare(spectrumFm.filename) == 0)
      {
        Spectrum * currScan = currSet.getIndex(key.str().c_str());

        if (currScan != NULL)
        {
          outputResults[j] = *currScan;
          outputResultIndex[j] = outputResults.size() - 1;
        }
        else
        {
          DEBUG_VAR(currPsm->m_scanNum);
          DEBUG_TRACE;
        }
      }
    }
  }

  for (int i = 0; i < outputResultIndex.size(); i++)
  {
    if (outputResultIndex[i] > -1)
    {
      psms.m_psmSet[i]->m_spectrum = &(outputResults[i]);
    }
  }

  ostringstream cmd;
  cmd.str("");
// generate sequential name, 1 based index
  cmd << DEFAULT_INPUT_FILE_BASE << ".pklbin";

  string pklbinPath = getProjPath(ip, cmd.str());

// List of generated data files
  vector<string> pklbinFileList;

  DEBUG_VAR(pklbinPath);
  pklbinFileList.push_back(pklbinPath);

// save it with the name we want (pklbin)
  DEBUG_MSG("Saving \'" << pklbinPath << "\' ...");
  outputResults.savePklBin(pklbinPath.c_str());

  // save input file list
  string aux;
  aux = getProjPath(ip, DEFAULT_INPUT_FILES_LIST);
  if (!writeFileIndex(aux.c_str(), pklbinFileList))
  {
    return false;
  }

  return true;
}

// -------------------------------------------------------------------------
bool performPairSpectra(ParameterList & ip,
                        PeptideSpectrumMatchSet &peptideResults,
                        SpecSet &specSet,
                        bool gridExecutionFlag)
{
  DEBUG_TRACE;
  DEBUG_VAR(peptideResults.size());
  DEBUG_VAR(specSet.size());
  ExecPairSpectra modulePairSpectra(ip,
      &peptideResults,
      &specSet);
  DEBUG_TRACE;

  string errorString;
  bool isValid = modulePairSpectra.validateParams(errorString);
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  bool returnStatus;

  if (!modulePairSpectra.loadInputData())
  {
    return false;
  }

  if (!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") <= 0)
  {
    returnStatus = modulePairSpectra.invoke();
  }
  else
  {
    DEBUG_MSG("Grid execution not supported for ExecPairSpectra");
    returnStatus = modulePairSpectra.invoke();
  }
// Test for return status
  DEBUG_VAR( returnStatus);
  if (!returnStatus)
  {
    ERROR_MSG("invoking modulePairSpectra exited in error.");
    return false;
  }

  DEBUG_TRACE;

  modulePairSpectra.saveOutputData();
  return true;
}

// -------------------------------------------------------------------------
bool performLpSolver(ParameterList & ip, bool gridExecutionFlag)
{
  DEBUG_TRACE;
  ExecLpSolver moduleLpSolver(ip);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleLpSolver.validateParams(errorString);
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  isValid = moduleLpSolver.loadInputData();
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  bool returnStatus;

  if (!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") <= 0)
  {
    returnStatus = moduleLpSolver.invoke();
  }
  else
  {
    DEBUG_TRACE;
    int numNodes = ip.getValueInt("GRID_NUMNODES");

    string gridType = ip.getValue("GRID_TYPE");
    if (gridType == "pbs")
    {
      ParallelPbsExecution exec(&moduleLpSolver,
          gridExecutionFlag,
          !gridExecutionFlag,
          false);
      returnStatus = exec.invoke(numNodes);
    }
    else if (gridType == "sge")
    {
      ParallelSgeExecution exec(&moduleLpSolver,
          gridExecutionFlag,
          !gridExecutionFlag,
          false);
      returnStatus = exec.invoke(numNodes);
    }
  }
  // Test for return status
  DEBUG_VAR( returnStatus);
  if (!returnStatus)
  {
    ERROR_MSG("invoking moduleFlr exited in error.");
    return false;
  }

  DEBUG_TRACE;

  return true;
}

bool runFlrSingle(ParameterList & ip, bool gridExecutionFlag)
{
  DEBUG_TRACE;
  ExecFlr moduleFlr(ip);
  DEBUG_TRACE;
}

// -------------------------------------------------------------------------
bool performFlr(ParameterList & ip, bool gridExecutionFlag)
{
  DEBUG_TRACE;
  ExecFlr moduleFlr(ip);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleFlr.validateParams(errorString);
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  isValid = moduleFlr.loadInputData();
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  bool returnStatus;

  if (!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") <= 0)
  {
    returnStatus = moduleFlr.invoke();
  }
  else
  {
    DEBUG_TRACE;
    int numNodes = ip.getValueInt("GRID_NUMNODES");

    string gridType = ip.getValue("GRID_TYPE");
    if (gridType == "pbs")
    {
      ParallelPbsExecution exec(&moduleFlr,
          gridExecutionFlag,
          !gridExecutionFlag,
          false);
      returnStatus = exec.invoke(numNodes);
    }
    else if (gridType == "sge")
    {
      ParallelSgeExecution exec(&moduleFlr,
          gridExecutionFlag,
          !gridExecutionFlag,
          false);
      returnStatus = exec.invoke(numNodes);
    }
  }
  // Test for return status
  DEBUG_VAR( returnStatus);
  if (!returnStatus)
  {
    ERROR_MSG("invoking moduleFlr exited in error.");
    return false;
  }

  DEBUG_TRACE;

  moduleFlr.saveOutputData();
  return true;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  int logLevel = 5;
  Logger::setDefaultLogger(Logger::getLogger(logLevel));

  string initialStageString("");
  string finalStageString("");
  string statusFileName("status.txt");

  bool resumeFlag = false;
  bool gridExecutionFlag = false;
  bool runGenoMSFlag = false;
  bool runMergedFlag = false;

  bool showHelp = false;
  for (size_t i = 0; i < argc; i++)
  {
    string arg(argv[i]);
    if (arg.compare("--help") == 0)
    {
      showHelp = true;
    }
  }

  if (argc < 2 || showHelp)
  {
    if (showHelp)
    {
      cout << "Usage: main_flr [PARAM FILE] [OPTION]..." << endl << endl;
      cout << "Optional arguments are listed below " << endl;
      cout << "  -i  <intialstage>   begin processing at specified stage:"
          << endl;
      cout << "                         begin,pairspectra,lpsolver,flr" << endl;
      cout
          << "  -f  <finalstage>   end processing after completing specified stage:"
          << endl;
      cout << "                        lpsolver,flr" << endl;

      cout << "  -g                  execution is on a grid" << endl;
      cout << "  -lf <filename>      name of log file for output" << endl;
      cout << "  -ll <loglevel>      log level for debug/warn/error output:"
          << endl;
      cout << "                         9 for errors only" << endl;
      cout << "                         5 for warnings and errors" << endl;
      cout << "                         0 for all debug output" << endl;
      cout << "  -s                   execute a single step then exit" << endl;
    }
    else
    {
      cerr << "main_specnets: insufficient arguments" << endl;
      cerr << "Try \'main_specnets --help\' for more information." << endl
          << endl;
    }

    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;

    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;

    return -1;
  }
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("i", "INITIAL_STAGE", 1));
  listOptions.push_back(CommandLineParser::Option("f", "FINAL_STAGE", 1));
  listOptions.push_back(CommandLineParser::Option("g", "GRID_EXECUTION", 0));
  listOptions.push_back(CommandLineParser::Option("lf", "LOG_FILE_NAME", 1));
  listOptions.push_back(CommandLineParser::Option("ll", "LOG_LEVEL", 1));
  listOptions.push_back(CommandLineParser::Option("s", "SINGLE_STEP", 0));

  CommandLineParser clp(argc, argv, 1, listOptions);
  string parser_error;
  if (!clp.validate(parser_error))
  {
    ERROR_MSG(parser_error);
    return -2;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  ParameterList ip;
  ip.readFromFile(argv[1]);
  ip.writeToFile("debug_sps.params");

  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  ip.addList(commandLineParams, true);
  ip.writeToFile("debug_wcommand.params");

  logLevel = ip.getValueInt("LOG_LEVEL", 5);
  if (ip.exists("LOG_FILE_NAME"))
  {
    string logFileName = ip.getValue("LOG_FILE_NAME");
    Logger::setDefaultLogger(Logger::getLogger(logFileName, logLevel));
  }
  else
  {
    Logger::setDefaultLogger(Logger::getLogger(logLevel));
  }

  // add divert for segfault
  addSegFaultDivert();

  DEBUG_TRACE;

  if (!ip.exists("EXE_DIR"))
  {

    // extract EXE_DIR from command line
    string exeDir(argv[0]);

    // find last /, and remove from that point on
    size_t found = exeDir.find_last_of("/\\");
    string aux = exeDir.substr(0, found);

    //string mainSpecnetsStr = "/main_specnets";
    //exeDir.erase(exeDir.length() - mainSpecnetsStr.length(), mainSpecnetsStr.length());

    // remove /ExecFramework, if it exists
    ip.setValue("EXE_DIR", aux);
  }

  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = commandLineParams.getValue("INITIAL_STAGE");
  }
  DEBUG_VAR(initialStageString);

  if (ip.exists("FINAL_STAGE"))
  {
    finalStageString = commandLineParams.getValue("FINAL_STAGE");
  }
  DEBUG_VAR(finalStageString);

  addDefaultParameterValues(ip);
  ip.writeToFile("debug_default.params");

  if (ip.exists("EXE_DIR"))
  {
    string exeDir = ip.getValue("EXE_DIR");

    // if path begins with '~', exit program.
    if (exeDir[0] == '~')
    {
      cout
          << "EXE_DIR path begins with tilde (~). Paths beginning with tilde are not supported."
          << endl;
      exit(0);
    }

    // In case there is a "/" at the end of EXE_DIR.. remove it
    if (exeDir.length() > 2 && exeDir[exeDir.length() - 1] == '/')
    {
      exeDir = exeDir.substr(0, exeDir.length() - 1);
      ip.setValue("EXE_DIR", exeDir);
    }
  }

  if (ip.exists("GRID_EXE_DIR"))
  {
    // In case there is a "/" at the end of EXE_DIR.. remove it
    string gridExeDir = ip.getValue("GRID_EXE_DIR");
    if (gridExeDir.length() > 2 && gridExeDir[gridExeDir.length() - 1] == '/')
    {
      gridExeDir = gridExeDir.substr(0, gridExeDir.length() - 1);
      ip.setValue("GRID_EXE_DIR", gridExeDir);
    }
  }

  if (ip.exists("GRID_SGE_EXE_DIR"))
  {
    // In case there is a "/" at the end of GRID_SGE_EXE_DIR.. remove it
    string gridSgeExeDir = ip.getValue("GRID_SGE_EXE_DIR");
    if (gridSgeExeDir.length() > 2
        && gridSgeExeDir[gridSgeExeDir.length() - 1] == '/')
    {
      gridSgeExeDir = gridSgeExeDir.substr(0, gridSgeExeDir.length() - 1);
      ip.setValue("GRID_SGE_EXE_DIR", gridSgeExeDir);
    }
  }

  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = ip.getValue("INITIAL_STAGE");
  }
  DEBUG_VAR(initialStageString);
  if (initialStageString.empty())
  {
    initialStageString = "begin";
  }
  DEBUG_VAR(initialStageString);

  if (ip.exists("FINAL_STAGE"))
  {
    finalStageString = ip.getValue("FINAL_STAGE");
  }
  DEBUG_VAR(finalStageString);
  if (finalStageString.empty())
  {
    finalStageString = "flr";
  }
  DEBUG_VAR(finalStageString);

  map<string, Stage> map_stage;
  map_stage["begin"] = STAGE_BEGIN;
  map_stage["pairspectra"] = STAGE_PAIRSPECTRA;
  map_stage["lpsolver"] = STAGE_LPSOLVER;
  map_stage["flr"] = STAGE_FLR;

  if (map_stage.find(initialStageString) == map_stage.end())
  {
    ERROR_MSG("Unknown starting stage [" << initialStageString << "]");
    return -1;
  }

  if (map_stage.find(finalStageString) == map_stage.end())
  {
    ERROR_MSG("Unknown final stage [" << finalStageString << "]");
    return -1;
  }

  //Start the status as "running" and write itout
  writeStatusFile(statusFileName, "Running");

  int initialStage = map_stage[initialStageString];
  DEBUG_VAR(initialStage);
  int finalStage = map_stage[finalStageString];
  DEBUG_VAR(finalStage);

  bool res;
  res = mkdir_if_not_exist("spectra");
  if (res)
  {
    DEBUG_MSG("Made directory \'spectra\'");
  }

  // get LD_LIBRARY_PATH from system
  char *curLibPath = getenv("LD_LIBRARY_PATH");

  // Build the needed library path
  string libPath;
  libPath = ip.getValue("EXE_DIR");
  // set LD_LIBRARY_PATH to EXE_DIR + /libs.
  libPath += "/libs";

  string fullLibPath;
  // check if LD_LIBRARY_PATH is already defined.
  if (curLibPath)
  {
    // if it is, check if it contains the path we want.
    fullLibPath = curLibPath;
    // if the library path IS NOT contained in the path variable, add it, and set the environment variable.
    if (fullLibPath.find(libPath) == string::npos)
    {
      fullLibPath += ':';
      fullLibPath += libPath;
      mysetenv("LD_LIBRARY_PATH", fullLibPath.c_str());
    }
  }
  else
  {
    // if LD_LIBRARY_PATH is not defined,, define it with what we want.
    mysetenv("LD_LIBRARY_PATH", libPath.c_str());
  }

  if (commandLineParams.exists("GRID_EXECUTION"))
  {
    gridExecutionFlag = true;
  }
  DEBUG_VAR(gridExecutionFlag);

  //---------------------------------
  // Load amino acid masses
  //---------------------------------
  AAJumps jumps(1); // Amino acid masses
  if (ip.exists("AMINO_ACID_MASSES"))
  {
    DEBUG_MSG("Loading amino acid masses from [" << ip.getValue("AMINO_ACID_MASSES") << "]");
    if (!jumps.loadJumps(ip.getValue("AMINO_ACID_MASSES").c_str(), true))
    {
      ERROR_MSG("Unable to load amino acid jumps");
      ERROR_MSG("Aborting!");
      exit(-1);
    }
  }
  else
  {
    DEBUG_MSG("No amino acid masses loaded. Using defaults");
  }

  string exeDir = ip.getValue("EXE_DIR");
  string convertCmd = exeDir + "/convert ";

  PeptideSpectrumMatchSet peptideResults;
  SpecSet filteredSpectra;

  if (initialStage == STAGE_BEGIN)
  {
    DEBUG_TRACE;

    //load results
    if (!loadPsmResults(ip, peptideResults))
    {
      return false;
    }

    if (!loadResultSpectrum(ip,
            ip.getValue("INPUT_SPECTRA_FILE_LIST").c_str(),
            peptideResults,
            filteredSpectra,
            ip.getValueBool("REINDEX_SCANS")))
    {
      return false;
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else
  {
    if (initialStage == STAGE_PAIRSPECTRA)
    {
      //load results
      if (!loadPsmResults(ip, peptideResults))
      {
        return false;
      }

      if (!loadSpecsMS(ip.getValue("EXE_DIR"), DEFAULT_INPUT_FILES_LIST, filteredSpectra))
      {
        ERROR_MSG("loadSpecsMS() exited in error. Exiting program");
        exit(0);
      }
      filteredSpectra.index();
      peptideResults.addSpectraByFilename(&filteredSpectra,false);
    }
  }

  DEBUG_VAR(peptideResults.size());
  DEBUG_VAR(filteredSpectra.size());

  if (finalStage == STAGE_BEGIN)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  if (initialStage <= STAGE_PAIRSPECTRA)
  {
    DEBUG_TRACE;

    if (!performPairSpectra(ip,
            peptideResults,
            filteredSpectra,
            ip.getValueBool("GRID_EXECUTION")))
    {
      exit(0);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }

  if (finalStage == STAGE_PAIRSPECTRA)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  if (initialStage <= STAGE_LPSOLVER)
  {
    DEBUG_TRACE;
    if (!performLpSolver(ip,
            ip.getValueBool("GRID_EXECUTION")))
    {
      exit(0);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }

  if (finalStage == STAGE_LPSOLVER)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  if (initialStage <= STAGE_FLR)
  {
    DEBUG_TRACE;
    if (!performLpSolver(ip,
            ip.getValueBool("GRID_EXECUTION")))
    {
      exit(0);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }

  if (finalStage == STAGE_FLR)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

}

