//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_performs.h"
#include "main_specnets_perform_assembly.h"
#include "main_specnets_perform_filteralign.h"
#include "main_specnets_perform_genoms.h"
#include "main_specnets_perform_specprotalign.h"

// Module Includes
#include "AlignmentPenaltyBased.h"
#include "CommandLineParser.h"
#include "FdrPeptide.h"
#include "Logger.h"

//#include "ExecReportProteinCoverage.h"
#include "db_fasta.h"
#include "ExecMergeConvert.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
#include "utils.h"  // stringSplit
#include "clusters.h"
#include "ClusterData.h"
#include "abruijn.h"
#include "copyright.h"
#include "SpecSet.h"
#include "Specific.h"
#include "StatusFile.h"

// System Includes
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

using namespace specnets;
using namespace std;


//-----------------------------------------------------------------------------

char *curLibPath = NULL;


//-----------------------------------------------------------------------------
void addDefaultParameterValues(ParameterList &p)
{
  // Basic parameters
  p.addIfDoesntExist("TOLERANCE_PEAK", "0.4");
  p.addIfDoesntExist("TOLERANCE_PM", "1.5");
  p.addIfDoesntExist("RESOLUTION", "0.1");

  // Preprocessing parameters
  p.addIfDoesntExist("CLUSTER_MIN_SIZE", "1");
  p.addIfDoesntExist("CLUSTER_MODEL", "LTQ_TRYP");
  p.addIfDoesntExist("MSCLUSTER_MIX_PROB ", "0.05");
  p.addIfDoesntExist("INSTRUMENT_TYPE", "IT");
  p.addIfDoesntExist("MIN_SPECTRUM_QUALITY", "0");
  p.addIfDoesntExist("CORRECT_PM", "no");
  p.addIfDoesntExist("GUESS_CHARGE", "no");
  p.addIfDoesntExist("PEPNOVO_OUTDIR", "spectra");

  // Alignment parameters
  p.addIfDoesntExist("AA_DIFF_COUNT", "2");
  p.addIfDoesntExist("MIN_SHIFT", "0");
  p.addIfDoesntExist("MIN_MOD_MASS", "-100");
  p.addIfDoesntExist("MAX_MOD_MASS", "100");
  p.addIfDoesntExist("MAX_NUM_MODS", "1");
  p.addIfDoesntExist("MIN_RATIO", "0.35");

  p.addIfDoesntExist("MAX_PVALUE", "0.045");
  p.addIfDoesntExist("MIN_MATCHED_PEAKS", "6"); // Minimum number of matched peaks to consider a spectral alignment
  p.addIfDoesntExist("MAX_AA_JUMP", "2");
  p.addIfDoesntExist("MIN_OVERLAP_AREA", "0.45");
  p.addIfDoesntExist("PENALTY_PTM", "-200"); // Set to "0" for SpecNets
  p.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
  p.addIfDoesntExist("FILTER_TRIGS", "no"); // Set to "no" for SpecNets
  p.addIfDoesntExist("PARTIAL_OVERLAPS", "1"); // Set to "0" for SpecNets
  p.addIfDoesntExist("TAGS_FILTER", "");
  p.addIfDoesntExist("TAGS_MATCH_FLANK", "1");
  p.addIfDoesntExist("TAGS_MATCH_COUNT", "2");

  // Comparative Shotgun Protein Sequencing (CSPS) parameters
  //p.addIfDoesntExist("CLUSTALW_EXE_DIR",      "");
  // p.addIfDoesntExist("CLUSTALW_MINSCORE", "250");
  p.addIfDoesntExist("CLUSTALW_MINSCORE", "1000000"); // Disabled by default - activate explicitly (i.e., set to ~250) for cSPS projects
  //  p.addIfDoesntExist("FORCE_REFERENCE",       "-1");

  // De novo sequencing parameters
  p.addIfDoesntExist("SPSPATH_MIN_NUM_PEAKS", "5");
  p.addIfDoesntExist("ADD_ENDPOINTS", "0");
  p.addIfDoesntExist("SPSPATH_MIN_NUM_SPECS", "2");
  p.addIfDoesntExist("PARALLEL_PATHS", "0");
  p.addIfDoesntExist("SPS_MIN_EDGES_TO_COMPONENT", "1");
  p.addIfDoesntExist("MIN_METACONTIG_SCORE", "3.35");
  p.addIfDoesntExist("MIN_METACONTIG_SIZE", "0");

  // tagsearch/matchma parameters
  p.addIfDoesntExist("TAG_LEN", "6");
  p.addIfDoesntExist("DOUBLE_AA_JUMPS", "1");
  p.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0"); // Set to 2 for SpecNets
  p.addIfDoesntExist("MAX_NUM_TAGS", "0");
  p.addIfDoesntExist("MAX_NUM_MODS", "2");
  p.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7"); // Minimum number of matched peaks between spectrum/database to accept PSM
  p.addIfDoesntExist("TAG_MATCH_TOP_SCORING_ONLY", "1");

  // Networks parameters (pathproj)
  //  p.addIfDoesntExist("MIN_PERC_EXPINT",   "0.01");
  //  p.addIfDoesntExist("MIN_PERC_TP",       "0.01");

  // Grid parameters
  p.addIfDoesntExist("GRID_TYPE", "pbs");
  p.addIfDoesntExist("GRID_NUMNODES", "-1");
  p.addIfDoesntExist("GRID_NUMCPUS", "1");
  p.addIfDoesntExist("GRID_EXE_DIR", "");
  p.addIfDoesntExist("GRID_SGE_EXE_DIR", "");

  DEBUG_VAR(p.getValue("GRID_PARAMS", "XXX"));
  p.addIfDoesntExist("GRID_PARAMS", "-l vmem=8G");

  string currDir = getCurrentDirectory();
  p.addIfDoesntExist("PROJECT_DIR", currDir);
  p.addIfDoesntExist("RELATIVE_DIR", "0");

  // Reporting parameters
  p.addIfDoesntExist("REPORT_DIR", "report");
  p.addIfDoesntExist("REPORT_DYNAMIC", "1");
  p.addIfDoesntExist("SPS_PROJECTS", "./sps_projects.txt");

}

//-----------------------------------------------------------------------------
/**
 * Loads raw MS/MS spectra from user, converts them to pklbin, then saves pklbin
 *   spectra along with lists of filenames to the project spectra directory. ExecMergeConvert
 *   is used to do any pre-processing of input spectra (set tolerances, set activation, etc.)
 * @param ip list of parameters from user
 * @param inputFileNames input file names delimeted by ";"
 * @param loadedSpectra output collection of SpecSets representing the converted spectra
 * @return true if successful, false if not
 */
//-----------------------------------------------------------------------------
bool loadInitialData(ParameterList & ip,
                     string inputFileNames,
                     vector<string>& loadedSpectraFiles)
{
  DEBUG_TRACE;
  //---------------------------------
  // Load spectra
  //---------------------------------

  string pklbinPath;

  // List of generated data files
  vector<string> pklbinFileList;
  //
  vector<string> filenames;
  ostringstream cmd;
  stringSplit(inputFileNames, filenames, ";");

  vector<string> set_filenames;

  if (ip.exists("SET_FILENAMES"))
  {
    stringSplit(ip.getValue("SET_FILENAMES"), set_filenames, ";");
    if (set_filenames.size() != filenames.size())
    {
      ERROR_MSG("SET_FILENAMES length (" << set_filenames.size() << ") must match INPUT_SPECS_MS length (" << filenames.size() << ")");
      return false;
    }
  }
  vector<pair<int, int> > loadedIndices;

  loadedSpectraFiles.resize(filenames.size());

  DEBUG_TRACE;
  unsigned int baseIdx, totalSpectrumCount = 0;

  bool merge_convert_error = false;

//#pragma omp parallel for
  for (unsigned int i = 0; i < filenames.size(); i++)
  {
    DEBUG_MSG("OPENING FILE");
    DEBUG_VAR (filenames[i]);

    ParameterList loadParams(ip);
    loadParams.setValue("INPUT_SPECTRA", filenames[i]);
    if (set_filenames.size() > 0)
    {
      loadParams.setValue("SET_FILENAMES", set_filenames[i]);
    }
    SpecSet tempSpecs;

    ExecMergeConvert* loader = new ExecMergeConvert(loadParams, &tempSpecs);

    // let ExecMergeConvert take care of loading spectra
    if (!loader->loadInputData())
    {
      delete loader;
      merge_convert_error = true;
      continue;
      //return false;
    }

    // let ExecMergeConvert pre-process spectra (set tolerances, activation, etc)
    if (!loader->invoke())
    {
      delete loader;
      merge_convert_error = true;
      continue;
    }

    cmd.str("");
    // generate sequential name, 1 based index
    cmd << SPECTRA_DIR << "/" << DEFAULT_INPUT_FILE_BASE << i + 1 << ".pklbin";
    pklbinPath = getProjPath(ip, cmd.str());
    DEBUG_VAR(pklbinPath);
    #pragma omp critical
    {
      pklbinFileList.push_back(pklbinPath);

      loadedSpectraFiles[i] = pklbinPath;
    }

    delete loader;

    if (!ExecMergeConvert::saveSpecset(pklbinPath, &tempSpecs))
    {
      merge_convert_error = true;
      continue;
    }

    DEBUG_TRACE;

    if(tempSpecs.size() == 0)
    {
      ERROR_MSG("Input spectra set is empty for file " << filenames[i]);
      //return false;
    }

    #pragma omp critical
    {
      totalSpectrumCount += tempSpecs.size();
    }

  }

  if(merge_convert_error){
      return false;
  }

  // Check if there are no spectra, in which case it fails
  if(totalSpectrumCount == 0)
  {
    ERROR_MSG("Input spectra set is empty");
    ERROR_MSG("Please make sure input files are 32-bit uncompressed mzXML files that contain MS/MS spectra");
    return false;
  }

  if (!ip.exists("SET_FILENAMES"))
  {
    set_filenames = filenames;
  }

  // save input file list
  string aux;
  aux = getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_USER_INPUT_FILES_LIST;
  if (!writeFileIndex(aux.c_str(), set_filenames))
  {
    return false;
  }

  DEBUG_TRACE;

  // save pklbin file list
  aux = getProjPath(ip, DEFAULT_INPUT_FILES_LIST);
  if (!writeFileIndex(aux.c_str(), pklbinFileList))
  {
    return false;
  }

  DEBUG_TRACE;
  return true;
}

#if 0
//-----------------------------------------------------------------------------
/**
 * Does the merging of separate sets of spectra into one. This sets scan #s in the merged set so
 *  they are unique and saves a mapping to keep track of them.
 * @param ip input list of user parameters
 * @param loadedSpectra input set of SpecSets, represents current working set. This will have size 0 after this function returns
 * @param ms2spectra output set of merged spectra w/ modified scan #s
 * @return true if successful, false if not
 */
//-----------------------------------------------------------------------------
bool mergeSeparateSpecs(ParameterList& ip,
                        vector<SpecSet>& loadedSpectra,
                        SpecSet& ms2spectra,
                        vector<vector<int> >& specMapping)
{
  DEBUG_VAR(loadedSpectra.size());

  if (loadedSpectra.size() == 1)
  {
    ms2spectra.clear();
    ms2spectra.swap(loadedSpectra[0]);
    loadedSpectra.resize(0);

    for (int i = 0; i < ms2spectra.size(); i++)
    {
      vector<int> aux;
      aux.push_back(0);
      aux.push_back(i);
      specMapping.push_back(aux);
    }
  }
  else
  {
    // count # of spectra
    unsigned int totalSpectrumCount = 0;
    for (unsigned int i = 0; i < loadedSpectra.size(); i++)
    {
      totalSpectrumCount += loadedSpectra[i].size();
    }
    totalSpectrumCount = 0;

    for (unsigned int i = 0; i < loadedSpectra.size(); i++)
    {

      ms2spectra.swapAppendSpecSet(loadedSpectra[i], false);

      for (unsigned int j = totalSpectrumCount; j < ms2spectra.size(); j++)
      {
        // make sure scan #s are unique
        ms2spectra[j].scan = j + 1;

        vector<int> aux;
        aux.push_back(i);
        aux.push_back(j - totalSpectrumCount);
        specMapping.push_back(aux);
      }
      totalSpectrumCount = ms2spectra.size();
      DEBUG_VAR(totalSpectrumCount);
    }
  }

  // save input mapping between input spectra files and combined spectra files
  string aux3 = getProjPath(ip, DEFAULT_INPUT_MAPPING);
  if (!Save_binArray(aux3.c_str(), specMapping))
  {
    ERROR_MSG("Failed to save to \'" << DEFAULT_INPUT_MAPPING << "\'");
    return false;
  }
  loadedSpectra.resize(0);
  return true;
}
#endif


#if 0
#include <sys/resource.h>
int setStackSize(int newSize)
{
  const rlim_t kStackSize = newSize * 1024 * 1024; // min stack size in MB
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  DEBUG_MSG("Stack soft limit: " << (int)(rl.rlim_cur) );
  DEBUG_MSG("Stack hard limit: " << (int)(rl.rlim_max) );

  if (result == 0)
  {
    if (rl.rlim_cur < kStackSize)
    {
      rl.rlim_cur = kStackSize;
      if(rl.rlim_cur > rl.rlim_max)
      rl.rlim_cur = rl.rlim_max;
      DEBUG_MSG("Changing stack soft limit to : " << rl.rlim_cur );
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)
      {
        int aa = errno;
        WARN_MSG("setrlimit returned result = " << result);
        WARN_MSG("errno = " << aa);
        string aux = strerror(aa);
        WARN_MSG("error = " << aux);
        return 0;
      }
      DEBUG_MSG("Stack soft limit changed to : " << kStackSize );
    }
  }

  // ...

  return 1;
}
#endif

//-----------------------------------------------------------------------------
void makeDirectoryIfNotExist(ParameterList & ip,
                             const char * directoryName)
{
  string fullName = getProjPath(ip, directoryName);
  bool res;
  res = mkdir_if_not_exist(fullName.c_str());
  if (res)
  {
    DEBUG_MSG("Made directory \'"<< fullName << "\'");
  }
}

//-----------------------------------------------------------------------------
int showVersion(void)
{
  cout << PROGRAM_NAME << endl;
  cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
  cout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
  cout << "SH1: " << XSTR(GIT_SH1) << endl;

  cout << endl;

  cout << COPYRIGHT1 << endl;
  cout << COPYRIGHT2 << endl;
  cout << endl;

  return 0;
}


//-----------------------------------------------------------------------------
void loadScanSpecificPenalties(string & filename,
                               map<int, map<int, float> > & penalties)
{
  vector<string> header;
  vector<vector<string> > lines;
  vector<string> requiredHeader;
  vector<int> requiredHeaderIndex;
  DelimitedTextReader::loadDelimitedFile(filename.c_str(),
                            "\t",
                            "",
                            header,
                            lines,
                            requiredHeader,
                            requiredHeaderIndex);
  for (int i = 0; i < lines.size(); i++) {
    //DEBUG_VAR(lines[i][0]);
    string scanString = lines[i][0];
    string massString = lines[i][1];
    string scoreString = lines[i][2];
    int scan = atoi(scanString.c_str());
    int mass = atoi(massString.c_str());
    float score = atof(scoreString.c_str());
    penalties[scan][mass] = score;
  }
  return;
}

//-----------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  int logLevel = 5;
  Logger::setDefaultLogger(Logger::getLogger(logLevel));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

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
      showHelp = true;

    if (arg.compare("--version") == 0)
      return showVersion();
  }

  if (argc < 2 || showHelp)
  {
    if (showHelp)
    {
      cout << "Usage: main_specnets [PARAM FILE] [OPTION]..." << endl << endl;
      cout << "Optional arguments are listed below " << endl;
      cout << "  -i  <intialstage>   begin processing at specified stage:"
          << endl;
      cout
          << "                         ms2deconv, mscluster, scoring, prmclustering"
          << endl;
      cout
          << "                         filterpairs, filteraligns, alignment, filterstarpairs"
          << endl;
      cout
          << "                         penaltygen, specnets, assembly, metaassembly"
          << endl;
      cout
          << "                         contigprotalign, spectaggen, specprotalign, protprotalign"
          << endl;
      cout
          << "                         homologyassembly, genoms, statprotseqs, report"
          << endl;
      cout
          << "  -f  <finalstage>   end processing after completing specified stage:"
          << endl;
      cout
          << "                         ms2deconv, mscluster, scoring, prmclustering"
          << endl;
      cout
          << "                         filterpairs, filteraligns, alignment, filterstarpairs"
          << endl;
      cout
          << "                         penaltygen, specnets, assembly, metaassembly"
          << endl;
      cout
          << "                         contigprotalign, spectaggen, specprotalign, protprotalign"
          << endl;
      cout
          << "                         homologyassembly, genoms, statprotseqs, report"
          << endl;
      cout << "  -proteosafe         input params are passed by ProteoSAFe"
          << endl;
      cout << "  -g                  execution is on a grid" << endl;
      cout << "  -lf <filename>      name of log file for output" << endl;
      cout << "  -ll <loglevel>      log level for debug/warn/error output:"
          << endl;
      cout << "                         9 for errors only" << endl;
      cout << "                         5 for warnings and errors" << endl;
      cout << "                         0 for all debug output" << endl;
      cout << "  -s                   execute a single step then exit" << endl;
      cout << "  -z                   resume an earlier run and prepare"
          << endl;
      cout
          << "                       specprotalign for running on a remote grid"
          << endl;
      cout << "       NOTE: this option is only used when intermediate results"
          << endl;
      cout << "             were saved during a previous run where FilterPairs"
          << endl;
      cout << "             was executed on a remote grid" << endl;
      cout << "  -zno                 resume an earlier run and execute"
          << endl;
      cout << "                       remaining steps on a single thread"
          << endl;
      cout << "       NOTE: this option is only used when intermediate results"
          << endl;
      cout << "             were saved during a previous run where FilterPairs"
          << endl;
      cout << "             was executed on a remote grid" << endl;
      cout << "  -q                   Run in GenoMS mode" << endl;
      cout << "       NOTE: In this mode, only 4 stages are valid: mscluster, "
          << endl;
      cout << "             scoring, genoms, and report" << endl;
      cout << "  -m                   Run in integrated mode (GenoMS + CSPS)"
          << endl;
      cout
          << "  -x  <execmode>        execution mode (sps, specnets, or signatures):"
          << endl << endl;
    }
    else
    {
      cerr << "main_specnets: insufficient arguments" << endl;
      cerr << "Type \'main_specnets --help\' for more information." << endl
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
  listOptions.push_back(CommandLineParser::Option("z", "RESUME_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("zno",
                                                  "RESUME_SINGLE_FLAG",
                                                  0));
  listOptions.push_back(CommandLineParser::Option("q", "GENOMS_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("m", "MERGE_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("x", "EXECUTION_MODE", 1));
  listOptions.push_back(CommandLineParser::Option("ccms_numnodes",
                                                  "GRID_NUMNODES",
                                                  1));

  CommandLineParser clp(argc, argv, 1, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "main_specnets: " << parserError << endl;
    cerr << "Type \'main_specnets --help\' for more information." << endl
        << endl;
    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;
    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  if (commandLineParams.exists("GENOMS_FLAG"))
    commandLineParams.setValue("GENOMS_FLAG", "1");

  if (commandLineParams.exists("MERGE_FLAG"))
    commandLineParams.setValue("MERGE_FLAG", "1");

  ParameterList ip;
  ip.readFromFile(argv[1]);
  ip.writeToFile("debug_sps.params");

  //---------------------------------------------------------------------------
  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  //---------------------------------------------------------------------------
  ip.addList(commandLineParams, true);
  ip.writeToFile("debug_wcommand.params");

  string execMode = ip.getValue("EXECUTION_MODE", "sps");

  if (ip.getValueInt("GENOMS_FLAG", 0) == 1)
    runGenoMSFlag = true;

  if (ip.getValueInt("MERGE_FLAG", 0) == 1)
    runMergedFlag = true;

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

  DEBUG_VAR(execMode);
  DEBUG_MSG("GenoMS_FLAG=" + ip.getValue("GENOMS_FLAG", "X"));
  DEBUG_MSG("Merge_FLAG=" + ip.getValue("MERGE_FLAG", "X"));

  DEBUG_TRACE;

  DEBUG_MSG("SH1:" << XSTR(GIT_SH1));

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
    string mainSpecnetsStr = "/ExecFramework";
    if (aux.rfind(mainSpecnetsStr) == aux.length() - mainSpecnetsStr.length())
      aux.erase(aux.length() - mainSpecnetsStr.length(),
                mainSpecnetsStr.length());
    else
    {
      mainSpecnetsStr = "\\ExecFramework";
      if (aux.rfind(mainSpecnetsStr) == aux.length() - mainSpecnetsStr.length())
        aux.erase(aux.length() - mainSpecnetsStr.length(),
                  mainSpecnetsStr.length());
    }

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

  DEBUG_VAR(commandLineParams.getValue("SINGLE_STEP"));

  if (ip.exists("MIN_PVALUE"))
  {
    ERROR_MSG("MIN_PVALUE is a deprecated variable, must use MAX_PVALUE instead");
    return -1;
  }

  DEBUG_VAR(ip.getValue("GRID_PARAMS", "XXX"));

  addDefaultParameterValues(ip);
  ip.writeToFile("debug_default.params");

  ip.addIfDoesntExist("OUTPUT_SPECTRA_PATH", getProjPath(ip, "./spectra"));
  if (ip.exists("EXE_DIR"))
  {
    string exeDir = ip.getValue("EXE_DIR");

    // if path begins with '~', exit program.
    if (exeDir[0] == '~')
    {
      cout
          << "EXE_DIR path begins with tilde (~). Paths beginning with tilde are not supported."
          << endl;
      return (0);
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
    finalStageString = "report";
  }
  DEBUG_VAR(finalStageString);

  map<string, Stage> map_stage;
  map_stage["begin"] = STAGE_BEGIN;
  map_stage["ms2deconv"] = STAGE_MS2DECONV;
  map_stage["mscluster"] = STAGE_MSCLUSTER;
  map_stage["scoring"] = STAGE_SCORING;
  map_stage["prmclustering"] = STAGE_PRMCLUSTERING;
  map_stage["filterpairs"] = STAGE_FILTERPAIRS;
  map_stage["filteraligns"] = STAGE_FILTERALIGNS;
  map_stage["alignment"] = STAGE_ALIGNMENT;
  map_stage["filterstarpairs"] = STAGE_FILTERSTARPAIRS;
  map_stage["penaltygen"] = STAGE_PENALTYGENERATE;
  map_stage["assembly"] = STAGE_ASSEMBLY;
  map_stage["metaassembly"] = STAGE_METAASSEMBLY;
  map_stage["specnets"] = STAGE_SPECNETS;
  map_stage["contigprotalign"] = STAGE_CONTIGPROTALIGN;
  map_stage["spectaggen"] = STAGE_SPECTAGGENERATE;
  map_stage["specprotalign"] = STAGE_SPECPROTALIGN;
  map_stage["protprotalign"] = STAGE_PROTPROTALIGN;
  map_stage["homologyassembly"] = STAGE_HOMOLOGYASSEMBLY;
  map_stage["genoms"] = STAGE_GENOMS;
  map_stage["merge"] = STAGE_MERGE;
  map_stage["statprotseqs"] = STAGE_STATPROTSEQS;
  map_stage["report"] = STAGE_REPORT;

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

  if (commandLineParams.exists("RESUME_FLAG")
      || commandLineParams.exists("RESUME_SINGLE_FLAG"))
  {
    resumeFlag = true;
  }
  DEBUG_VAR(resumeFlag);

  if (commandLineParams.exists("GRID_EXECUTION"))
  {
    gridExecutionFlag = true;
  }
  DEBUG_VAR(gridExecutionFlag);

  //---------------------------------------------------------------------------
  // Initialize environment
  //---------------------------------------------------------------------------
  DEBUG_MSG("Spectral Networks 2.0.0: session started on: " << getCurrentTimeString());
  DEBUG_MSG("Starting stage: [" << initialStageString << "] on: " << getCurrentTimeString());

  ip.writeToFile(getProjPath(ip, "debug_initial.params"));

  makeDirectoryIfNotExist(ip, SPECTRA_DIR.c_str());
  makeDirectoryIfNotExist(ip, SPECTRA_GRID_SCORING_DIR.c_str());
  makeDirectoryIfNotExist(ip, ALIGNS_DIR.c_str());
  if (commandLineParams.exists("PROTEOSAFE_FLAG")) {
    makeDirectoryIfNotExist(ip, "aligns/params");
  }
  makeDirectoryIfNotExist(ip, "specnets");
  makeDirectoryIfNotExist(ip, "homology");

  makeDirectoryIfNotExist(ip, "homology/grid_t");
  makeDirectoryIfNotExist(ip, "homology/grid_d");
  makeDirectoryIfNotExist(ip, "homology/grid2_t");
  makeDirectoryIfNotExist(ip, "homology/grid2_d");
 
  makeDirectoryIfNotExist(ip, "report");
  makeDirectoryIfNotExist(ip, "ReportData");
  makeDirectoryIfNotExist(ip, "assembly");

  // get LD_LIBRARY_PATH from system
  curLibPath = getenv("LD_LIBRARY_PATH");

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
  jumps.saveJumps(getProjPath(ip, "homology/specprotalign_aa_masses.txt").c_str());

  string exeDir = ip.getValue("EXE_DIR");
  string convertCmd = exeDir + "/convert ";

  vector<string> ms2Spectra_separate; // MS/MS spectra in separate files
  string inputPklbinFilesList = getProjPath(ip, DEFAULT_INPUT_FILES_LIST);

  if (initialStage == STAGE_BEGIN)
  {
    DEBUG_TRACE;

    // Build the needed library path
    string libPath = exeDir + "/libs";
    addEnvironmentVariable(convertCmd, "LD_LIBRARY_PATH", libPath);

    bool ret = loadInitialData(ip, ip.getValue("INPUT_SPECS_MS"), ms2Spectra_separate);

    if(!ret)
    {
      ERROR_MSG("Failed loading initial data.");
      writeStatusFile(statusFileName, "Error");
      return(1);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return(0);
    }
  }
  else
  {
    if (!readFilesFromFile(inputPklbinFilesList, ms2Spectra_separate))
    {
      ERROR_MSG("readFilesFromFile() failed for " << inputPklbinFilesList);
      return(0);
    }
  }

  if (finalStage == STAGE_BEGIN)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  //TEST_RETURN("moduleBegin", *ms2Spectra_separate);
  DEBUG_VAR(ms2Spectra_separate.size());

  if (ms2Spectra_separate.size() == 0)
  {
    ERROR_MSG("ms2Spectra_separate size is 0!");
    writeStatusFile(statusFileName, "Error");
    exit(-2);
  }

  vector<string> deconvMs2Spectra_separate; // MS/MS spectra in separate files
  string deconvPklbinFilesList = 
          getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_DECONV_FILES_LIST;

  if (initialStage <= STAGE_MS2DECONV && ip.getValueInt("DECONV_MS2", 0) > 0)
  {

    if (!performMS2Deconv(ip, ms2Spectra_separate, deconvMs2Spectra_separate))
    {
      return (0);
    }

    DEBUG_VAR(ms2Spectra_separate.size());

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else if (ip.getValueInt("DECONV_MS2", 0) > 0)
  {
    if (!readFilesFromFile(deconvPklbinFilesList, deconvMs2Spectra_separate))
    {
      ERROR_MSG("readFilesFromFile() failed for " << inputPklbinFilesList);
      return (0);
    }
  }
  else
  {
    deconvMs2Spectra_separate = ms2Spectra_separate;
    deconvPklbinFilesList = inputPklbinFilesList;
  }

  if (finalStage == STAGE_MS2DECONV)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  string clusterTool = ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST);

  const bool performClustering = ip.getValueInt("CLUSTER_MIN_SIZE", 0) > 0;
  SpecSet ms2Spectra;

  if (initialStage <= STAGE_MSCLUSTER)
  {

    DEBUG_MSG("Starting stage McCluster on: " << getCurrentTimeString());
    if (!performMsCluster(ip, deconvMs2Spectra_separate, ms2Spectra))
    {
      ERROR_MSG("Problem encountered during MsCluster stage");
      writeStatusFile(statusFileName, "Error");
      exit(-1);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else
  {
    // load clustered spectra
    DEBUG_MSG("Bypassing MsCluster stage");
    string ms2SpectraFile = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_MS_FILE;
    if (ms2Spectra.loadPklBin(ms2SpectraFile.c_str()) <= 0) {
      ERROR_MSG("Could not load input MS/MS spectra! [" << ms2SpectraFile << "]");
    }

    DEBUG_VAR(ms2Spectra.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleMsCluster", ms2Spectra);

  DEBUG_VAR(ms2Spectra.size());

  if (ms2Spectra.size() == 0)
  {
    ERROR_MSG("ms2Spectra size is 0!");
    writeStatusFile(statusFileName, "Error");
    exit(-2);
  }

  if (finalStage == STAGE_MSCLUSTER)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  //---------------------------------------------------------------------------
  // PEPNOVO/SCORING STAGE
  //---------------------------------------------------------------------------
  SpecSet prmSpectra;
  if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
  {
    prmSpectra = ms2Spectra;
  }
  else
  {
    if (initialStage <= STAGE_SCORING)
    {
      DEBUG_MSG("Starting stage Scoring on: " << getCurrentTimeString());
      if (!performScoring(ip,
                          ms2Spectra,
                          prmSpectra,
                          gridExecutionFlag,
                          resumeFlag)) {
        ERROR_MSG("Problem encountered during Scoring stage");
        writeStatusFile(statusFileName, "Error");
        exit(-3);
      }

      if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0) && 
          !ip.exists("NO_GRID_SCORING") && (!resumeFlag && !gridExecutionFlag)) {
        // If we are doing a grid execution (and are not actually on the grid)
        //    and we aren't resuming... then exit (we'll resume execution later)
        DEBUG_MSG("Files for grid execution have been saved.");
        DEBUG_MSG("Restart with -z option when grid execution has been completed.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }

      if (runGenoMSFlag || runMergedFlag) {
        // need to save a copy of the scored spectra as MGF so GenoMS works
        if (!ExecMergeConvert::saveSpecset(
            getProjPath(ip, SPECTRA_DIR) + "/" +SPECS_SCORED_MGF_FILE,
                                           &prmSpectra)) {
          ERROR_MSG("Failed to save to " << 
              getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_SCORED_MGF_FILE);
          exit(-3);
        }
      }

      if (commandLineParams.exists("SINGLE_STEP")) {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Scoring stage");
      string prmFile = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_SCORED_FILE;
      if (prmSpectra.loadPklBin(prmFile.c_str()) <= 0) {
        ERROR_MSG("Could not load PRM spectra! [" << prmFile << "]");
      }
      DEBUG_VAR(prmSpectra.size());
    }
  }
  
  DEBUG_VAR(prmSpectra.size());

  // Test for empty return data structures
  TEST_RETURN("moduleScoring", prmSpectra);

  if (finalStage == STAGE_SCORING)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  if (performClustering && clusterTool == CLUSTER_PRMS
      && initialStage <= STAGE_PRMCLUSTERING)
  {
    DEBUG_MSG("Starting stage PRM Clustering on: " << getCurrentTimeString());
    if (!performPrmClustering(ip, prmSpectra, prmSpectra))
    {
      ERROR_MSG("Problem encountered during PRM Clustering stage");
      writeStatusFile(statusFileName, "Error");
      exit(-3);
    }

    if (runGenoMSFlag || runMergedFlag)
    {
      // need to save a copy of the scored spectra as MGF so GenoMS works
      string mgfFile = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_SCORED_MGF_FILE;
      if (!ExecMergeConvert::saveSpecset(mgfFile, &prmSpectra))
      {
        ERROR_MSG("Failed to save to " << mgfFile);
        exit(-3);
      }
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing PRM Clustering stage");
    DEBUG_VAR(prmSpectra.size());
  }

  for (int i = 0; i < prmSpectra.size(); i++)
  {
    prmSpectra[i].msFragType = Spectrum::FragType_PRM;
  }

  if (finalStage == STAGE_PRMCLUSTERING)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }
  
  DEBUG_VAR(prmSpectra.size());

  //--------------------------------------------------------------------------
  // GenoMS STAGE
  //--------------------------------------------------------------------------
  if (runGenoMSFlag || runMergedFlag)
  {
    if (initialStage <= STAGE_GENOMS)
    {

#if INCLUDE_GENOMS == 0
      ERROR_MSG("GenoMS is not available. Returning in error.");
      exit(-STAGE_GENOMS);
      return false;
#endif

      DEBUG_MSG("Starting GenoMS");

      if (!performGenoMS(ip))
      {
        ERROR_MSG("Problem encountered during GenoMS");
        writeStatusFile(statusFileName, "Error");
        exit(-1);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing GenoMS stage");
    }
    /*
     SpecSet testSpectra;

     testSpectra.loadPklBin( getProjPath(ip, "./spectra/stars.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/stars.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./spectra/specs_ms.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/spec_ms2.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./assembly/sps_seqs.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./assembly/sps_seqs.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./homology/contigs_midx_all.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/contigs_midx_all.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./homology/contigs.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/contigs.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./homology/homglue_matches.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/homglue_matches.mgf").c_str());
     DEBUG_TRACE;
     abinfo_t abinfo;
     Load_abinfo("./assembly/component_info.bin", abinfo);
     dumpAbInfo("./assembly/genoms_component_info.txt", abinfo);
     Save_abinfo_v1_0("./assembly/component_info.sps.bin",abinfo);
     */
  }

  if (runMergedFlag && initialStage <= STAGE_GENOMS)
  {

    DEBUG_MSG("Stage Genoms: " << STAGE_GENOMS);
    
    performGenoMSRename(statusFileName);

  }
  else if (runGenoMSFlag)
  {

    //If we are running GenoMS, then we need to update FASTA_DATABASE to point to our output

    ip.setValue("FASTA_DATABASE", "./protid.fasta");

    if (initialStage <= STAGE_REPORT)
    {
      DEBUG_MSG("Starting Report stage");
      ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                       ios_base::out | ios_base::binary);
      spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
          << ip.getValue("TOLERANCE_PM") << "\n";
      spsProj.close();

      if (!performReport(ip))
      {
        ERROR_MSG("Problem encountered during Report stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_REPORT);
      }
      //--------------------------------------------------------------------------
      //Set up for relaunching
      //--------------------------------------------------------------------------
      if (!generateRelaunchScript(ip))
      {
        ERROR_MSG("Problem encountered during relaunch script creation");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_REPORT);
      }

    }
    else
    {
      DEBUG_MSG("Bypassing Report stage");
    }
    writeStatusFile(statusFileName, "Finished");
    return 0;
  }

  SpectrumPairSet filteredPairs;

  vector<TwoValues<float> > ratios;
  vector<TwoValues<float> > means;
  vector<float> varTerms;
  list<vector<float> > alignStats;
  vector<vector<float> > specStats;
  vector<unsigned int> idxKept;
  vector<TwoValues<float> > pvalues;

  DEBUG_VAR(prmSpectra.size());
  //---------------------------------------------------------------------------
  // FLITERPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERPAIRS)
  {

	if( ip.getValueBool("TAG_FILTER", false) ) {

		prmSpectra.saveTags(getProjPath(ip, "./aligns/sequence.tags"),
				            3,//tag length
				            1,//gap
				            ip.getValueInt("MAX_TAG_SIZE", 50)*2,
				            ip.getValueDouble("TOLERANCE_PEAK", 0.5),
				            ip.getValue("PAIRS_MATCH_MODE", "") != "cosine");

		ip.setValue("INPUT_TAGS", getProjPath(ip, "./aligns/sequence.tags"));
	}

    DEBUG_MSG("Starting stage FilterPairs on: " << getCurrentTimeString());
    if (!performFilterPairs(ip,
                            prmSpectra,
                            ms2Spectra,
                            filteredPairs,
                            ratios,
                            means,
                            varTerms,
                            alignStats,
                            specStats,
                            idxKept,
                            pvalues,
                            gridExecutionFlag,
                            resumeFlag))
    {
      ERROR_MSG("Problem encountered during Filter Pairs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-4);
    }

    if (commandLineParams.exists("PROTEOSAFE_FLAG"))
    {
      bool res = spsSystem("cp aligns/*.params aligns/params/");
      if (!res)
      {
        ERROR_MSG("Failed to copy params files to aligns/params !!!");
        exit(-4);
      }
    }

    if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
        && (!resumeFlag && !gridExecutionFlag))
    {
      // If we are doing a grid execution (and are not actually on the grid)
      //    and we aren't resuming... then exit (we'll resume execution later)
      DEBUG_MSG("Files for grid execution have been saved.");
      DEBUG_MSG("Restart with -z option when grid execution has been completed.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

    if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
    {
      filteredPairs.saveToBinaryFile(getProjPath(ip,
                                                 "./aligns/pairs_cosine.bin"));
    }
    else
    {
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs_raw.bin"));
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

  }
  else
  {
    DEBUG_MSG("Bypassing Filter Pairs stage");
    filteredPairs.loadFromBinaryFile(getProjPath(ip, "./aligns/pairs_raw.bin"));
    DEBUG_VAR(filteredPairs.size());
    Load_binArray(getProjPath(ip, "aligns/ratios.bin").c_str(), ratios);
    DEBUG_VAR(ratios.size());
    Load_binArray(getProjPath(ip, "aligns/means.bin").c_str(), means);
    DEBUG_VAR(means.size());
    Load_binArray(getProjPath(ip, "aligns/vars.bin").c_str(), varTerms);
    DEBUG_VAR(varTerms.size());
  }

  if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
  {
    DEBUG_MSG("Parameters include PAIRS_MATCH_MODE=cosine, exiting after ExecFilterPairs");
    return 0;
  }

  // Test for empty return data structures
  //  TEST_RETURN("moduleFilterPairs", ratios);
  //  TEST_RETURN("moduleFilterPairs", means);
  //  TEST_RETURN("moduleFilterPairs", varTerms);
  TEST_RETURN("moduleFilterPairs", filteredPairs);
  //  TEST_RETURN("moduleFilterPairs", pvalues);
  //  TEST_RETURN("moduleFilterPairs", idxKept);
  //  TEST_RETURN("moduleFilterPairs", alignStats);
  //  TEST_RETURN("moduleFilterPairs", specStats);

  if (finalStage == STAGE_FILTERPAIRS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  PeptideSpectrumMatchSet filterPsmSet;
  string psmFilename = ip.getValue("INPUT_FILTER_PSMS");
  if (!psmFilename.empty()) {
    DEBUG_MSG("Loading filter PSM file.. [" << psmFilename << "]");
    if (!filterPsmSet.loadFromFile(ip.getValue("INPUT_FILTER_PSMS").c_str())) {
      ERROR_MSG("Error reading filter PSM file.");
      return false;
    }
  } 
  DEBUG_VAR(filterPsmSet.size());

  psmFilename = ip.getValue("INPUT_FILTER_MSGFDB_PSMS");
  if (!psmFilename.empty()) {
    DEBUG_MSG("Loading filter MSGFDB PSM file.. [" << psmFilename << "]");
    if (!filterPsmSet.loadFromFile(ip.getValue("INPUT_FILTER_MSGFDB_PSMS").c_str())) {
      ERROR_MSG("Error reading filter MSGFDB PSM file.");
      return false;
    }
  } 
  DEBUG_VAR(filterPsmSet.size());

  //---------------------------------------------------------------------------
  // FLITERALIGNS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERALIGNS)
  {
    DEBUG_MSG("Starting stage FilterAligns on: " << getCurrentTimeString());
    if (!performFilterAligns(ip,
                             filteredPairs,
                             prmSpectra,
                             filterPsmSet,
                             ratios,
                             means,
                             varTerms,
                             idxKept,
                             pvalues))
    {
      ERROR_MSG("Problem encountered during Filter Aligns stage");
      writeStatusFile(statusFileName, "Error");
      exit(-4);
    }
    DEBUG_VAR(filteredPairs.size());

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
    else
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Aligns stage");
    filteredPairs.loadFromBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
    DEBUG_VAR(filteredPairs.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleFilterPairs", filteredPairs);

  if (finalStage == STAGE_FILTERALIGNS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  if (commandLineParams.exists("RESUME_SINGLE_FLAG"))
  {
    ip.setValue("GRID_NUMNODES", "-1");
  }

  SpecSet pairAlignments;
  SpecSet starSpectraOnly;
  SpecSet starSpectra;
  vector<unsigned int> alignedSpectra;

  //---------------------------------------------------------------------------
  // ALIGNMENT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ALIGNMENT)
  {
    DEBUG_MSG("Starting stage Alignment on: " << getCurrentTimeString());
    if (!performAlignment(ip,
                          jumps,
                          prmSpectra,
                          filteredPairs,
                          pairAlignments,
                          starSpectraOnly,
                          starSpectra,
                          alignedSpectra))
    {
      ERROR_MSG("Problem encountered during Alignment stage");
      writeStatusFile(statusFileName, "Error");
      exit(-5);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Alignment stage");
    string starSpectraFile = getProjPath(ip, SPECTRA_DIR) + "/" + "stars.pklbin";
    if (starSpectra.loadPklBin(starSpectraFile.c_str()) <= 0) {
      ERROR_MSG("Could not load star spectra! [" << starSpectraFile << "]");
    }
    DEBUG_VAR(starSpectra.size());
  }

  // Test for empty return data structures
  //TEST_RETURN("ExecAlignment", pairAlignments);
  //TEST_RETURN("ExecAlignment", starSpectraOnly);
  TEST_RETURN("ExecAlignment", starSpectra);
  //TEST_RETURN("ExecAlignment", alignedSpectra);

  if (finalStage == STAGE_ALIGNMENT)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  vector<vector<float> > starRatios;
  SpecSet matchedPeaks;

  SpectrumPairSet filteredStarPairs;
  filteredStarPairs.copy(filteredPairs);
  
  //---------------------------------------------------------------------------
  // FILTERSTARPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERSTARPAIRS)
  {
    if (!performFilterStarPairs(ip,
                                filteredStarPairs,
                                starSpectra,
                                starRatios,
                                matchedPeaks))
    {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Star Pairs stage");
    //    DEBUG_MSG("Not implemented yet" );
    //    Load_binArray("./aligns/ratios.bin", starRatios);
    //    DEBUG_VAR(starRatios.size());
    //    matchedPeaks->loadPklBin("");
    //    DEBUG_VAR(matchedPeaks.size());

    if (!filteredStarPairs.loadFromBinaryFile("aligns/pairs_stars.bin"))
    {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage (loading filteredStarPairs)");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }
  }

  if (finalStage == STAGE_FILTERSTARPAIRS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  DEBUG_MSG("Reloading Star specturm");
  // Need to reload the star spectra because they were altered by FilterStarPairs
  if (starSpectra.loadPklBin(getProjPath(ip, "./spectra/stars.pklbin").c_str())
      <= 0)
  {
    ERROR_MSG("Problem encountered while reloading star specturm");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_ASSEMBLY);
  }
  DEBUG_VAR(starSpectra.size());

  if (execMode == string("signatures"))
  {
    return 0;
  }

  // Test for empty return data structures
  //TEST_RETURN("ExecFilterStarPairs", starRatios);
  //TEST_RETURN("ExecFilterStarPairs", matchedPeaks);

  //---------------------------------------------------------------------------
  // Spectral Networks STAGE
  //---------------------------------------------------------------------------
  DB_fasta dbAll;
  string dbFileName;
  DEBUG_MSG(ip.getValue("FASTA_DATABASE"));
  if (ip.exists("FASTA_DATABASE"))
  {
    dbFileName = ip.getValue("FASTA_DATABASE");

    //If we are integrating, then add GenoMS protein sequences to the DB
    if (runMergedFlag && initialStage <= STAGE_FILTERSTARPAIRS)
    {
      DEBUG_MSG("Concatenating " << dbFileName << " and genoms.protid.fasta");
      if (!concatenateFiles("./genoms.protid.fasta",
                            dbFileName,
                            "./protid.fasta"))
      {
        ERROR_MSG("Problem encountered concatenating fasta files");
        writeStatusFile(statusFileName, "Error");
        exit(-100);
      }
      dbFileName = "./protid.fasta";
      ip.setValue("FASTA_DATABASE", "./protid.fasta");
    }
    if (dbAll.Load(dbFileName.c_str()) <= 0)
    {
      ERROR_MSG("Problem encountered during Spectral Networks stage (loading " << dbFileName << ")");
      writeStatusFile(statusFileName, "Error");
      exit(-100);
    }
  }

  //---------------------------------------------------------------------------
  // Create penalty matrices for new alignment
  //---------------------------------------------------------------------------
  float resolution = ip.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
  float maxPeakEquivalents = ip.getValueDouble("MAX_PEAK_EQUIVALENTS", 1.5);
  float minPeakEquivalents = ip.getValueDouble("MIN_PEAK_EQUIVALENTS", 1.0);
  float minFrequency = ip.getValueDouble("MIN_PENALTY_FREQUENCY", 0.005);
  float unknownPenalty = ip.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY",
                                           1.0);
  float knownModPenalty = ip.getValueDouble("PENALTY_ALIGNMENT_KNOWN_PENALTY",
                                            0.01);
  float minModMass = ip.getValueDouble("MIN_MOD_MASS", -100.0);
  float maxModMass = ip.getValueDouble("MAX_MOD_MASS", 100.0);
  float unknownMultiplier =
      ip.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 1.5);
      
  PenaltyMatrix penaltyMatrixBlosum(jumps,
                                    resolution,
                                    knownModPenalty,
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);

  PenaltyMatrix penaltyMatrixMods(jumps,
                                  resolution,
                                  knownModPenalty,
                                  unknownPenalty,
                                  unknownMultiplier,
                                  minModMass,
                                  maxModMass);
  map<int, map<int, float> > scanSpecificPenalties;

  bool penaltyAlign = ip.getValueBool("PENALTY_ALIGNMENT");
  DEBUG_VAR(penaltyAlign);
  if (penaltyAlign) {
  
    if (ip.exists("BLOSUM_PENALTY_FILE")) {
      string blosumFilename = ip.getValue("BLOSUM_PENALTY_FILE");
      DEBUG_VAR(blosumFilename);
      if (!penaltyMatrixBlosum.loadFromBlosum(blosumFilename, maxPeakEquivalents)) {
        ERROR_MSG("Unable to load BLOSUM_PENALTY_FILE [" << blosumFilename << "]");
        exit(-1);
      }
    }

    string blosumMatFilename(getProjPath(ip, "homology/specprotalign_blosum_penal.txt"));
    DEBUG_VAR(blosumMatFilename);
    penaltyMatrixBlosum.saveMatrix(blosumMatFilename);
  
    if (ip.exists("KNOWN_MODS_FILE")) {
      string knowmModsFileName = ip.getValue("KNOWN_MODS_FILE");
      DEBUG_VAR(knowmModsFileName);
      if (!penaltyMatrixMods.loadKnownModifications(knowmModsFileName)) {
        ERROR_MSG("Unable to load known mods [" << knowmModsFileName << "]");
        exit(-1);
      }
    }
    if (ip.exists("CLEAVAGE_PENALTY_FILE")) {
      string cleavagePenaltiesFileName = ip.getValue("CLEAVAGE_PENALTY_FILE");
      DEBUG_VAR(cleavagePenaltiesFileName);
      if (!penaltyMatrixMods.loadCleavagePenalties(cleavagePenaltiesFileName)) {
        ERROR_MSG("Unable to load cleavage penalties [" << cleavagePenaltiesFileName << "]");
        exit(-1);
      }
    }
  
    //---------------------------------------------------------------------------
    // Penalty Matrix Generation STAGE
    //---------------------------------------------------------------------------
    if (initialStage <= STAGE_PENALTYGENERATE) {
      if (!performPenaltyGen(ip,
                             filteredPairs,
                             prmSpectra,
                             penaltyMatrixMods,
                             scanSpecificPenalties)) {
        ERROR_MSG("Problem encountered during PenaltyGen stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_PENALTYGENERATE);
      }

      if (commandLineParams.exists("SINGLE_STEP")) {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }
    } else {
      string modFileName = getProjPath(ip, "homology/specprotalign_mod_penal.txt");
      string knowmModsFileName = getProjPath(ip, "homology/specprotalign_known_mods.txt");
      string cleavagePenaltiesFileName = getProjPath(ip, "homology/specprotalign_cleave_pen.txt");
      DEBUG_MSG("Loading: " << modFileName);
      DEBUG_MSG("Loading: " << knowmModsFileName);
      DEBUG_MSG("Loading: " << cleavagePenaltiesFileName);
      if (!penaltyMatrixMods.load(modFileName,
                                  knowmModsFileName,
                                  cleavagePenaltiesFileName)) {
        ERROR_MSG("Error loading mod penalties from 3 files.");
        return false;
      }
      string scanPenaltiesFileName = getProjPath(ip, "homology/specprotalign_scan_specific_pen.txt");
      loadScanSpecificPenalties(scanPenaltiesFileName, scanSpecificPenalties);
    }
  } // if (penaltyAlign)

  if (finalStage == STAGE_PENALTYGENERATE)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  //---------------------------------------------------------------------------
  // SPECNETS EXEC MODE (as opposed to SPS)
  //   The name "main_specnets" is a bit misleading.. prob should be "main_sps"
  //---------------------------------------------------------------------------
  if (execMode == string("specnets"))
  {
    // Entering Spectral Networks mode
    PeptideSpectrumMatchSet * input_psms = new PeptideSpectrumMatchSet;
    PeptideSpectrumMatchSet * output_psms = new PeptideSpectrumMatchSet;
    SpecSet * psms_spectra = new SpecSet;
    SpecSet * psms_midx = new SpecSet;
    vector<vector<int> > *psms_mp = new vector<vector<int> >;
    SpecSet * snets_contigs = new SpecSet;
    SpecSet * snets_midx = new SpecSet;
    float resolution = ip.getValueFloat("@", .1);

    vector<vector<int> > *snets_mp = new vector<vector<int> >;

    if (ip.exists("INSPECT_PSMS"))
    {
      //load in inspect seed PSMS
      if (!input_psms->loadInspectResultsFile(ip.getValue("INSPECT_PSMS").c_str(),
                                              ip.getValueBool("SCAN_ZERO_INDEX",
                                                              1)))
      {
        ERROR_MSG("Unable to load inspect file! " << ip.getValue("INSPECT_PSMS"));
        writeStatusFile(statusFileName, "Error");
        exit(-100);
      }
    }

    if (not performExecMainSpecnets(ip,
                                    &ms2Spectra,
                                    &prmSpectra,
                                    &starSpectra,
                                    &filteredStarPairs,
                                    &dbAll,
                                    &penaltyMatrixMods,
                                    &penaltyMatrixBlosum,
                                    output_psms,
                                    input_psms,
                                    psms_spectra,
                                    psms_midx,
                                    psms_mp,
                                    snets_contigs,
                                    snets_midx,
                                    snets_mp))
    {
      ERROR_MSG("Problem encountered during ExecMainSpecnets stage");
      writeStatusFile(statusFileName, "Error");
      exit(-100);
    }

    delete input_psms;
    delete output_psms;
    delete psms_spectra;
    delete psms_midx;
    delete psms_mp;
    delete snets_contigs;
    delete snets_midx;
    delete snets_mp;

    return 0;
  }
  //---------------------------------------------------------------------------
  // END OF SPECNETS EXEC MODE
  //---------------------------------------------------------------------------

  Clusters contigShifts;
  abinfo_t contigAbinfo;

  //---------------------------------------------------------------------------
  // ASSEMBLY STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ASSEMBLY and execMode == string("sps"))
  {
    if (!performAssembly(ip,
                         starSpectra,
                         filteredStarPairs,
                         contigShifts,
                         contigAbinfo))
    {
      ERROR_MSG("Problem encountered during Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
    contigShifts.Save(getProjPath(ip, "assembly/path_spectra_as_cluster.txt").c_str());

    bool res = FileCopy(getProjPath(ip, "assembly/sps_seqs.pklbin").c_str(),
                        getProjPath(ip, "assembly/old_sps_seqs.pklbin").c_str());
    //spsSystem("cp assembly/sps_seqs.pklbin assembly/old_sps_seqs.pklbin");
    res = FileCopy(getProjPath(ip, "assembly/path_spectra_as_cluster.txt").c_str(),
                   getProjPath(ip, "assembly/old_path_spectra_as_cluster.txt").c_str());
    //= spsSystem("cp assembly/path_spectra_as_cluster.txt assembly/old_path_spectra_as_cluster.txt");
    res = FileCopy(getProjPath(ip, "assembly/component_info.bin").c_str(),
                   getProjPath(ip, "assembly/old_component_info.bin").c_str());
    //= spsSystem("cp assembly/component_info.bin assembly/old_component_info.bin");
  }
  else
  {
    DEBUG_MSG("Bypassing Assembly stage");
    if (contigShifts.Load(getProjPath(ip, "assembly/old_path_spectra_as_cluster.txt").c_str()) <= 0
        and contigShifts.Load(getProjPath(ip, "assembly/path_spectra_as_cluster.txt").c_str()) <= 0)
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
    if (contigShifts.consensus.loadPklBin(getProjPath(ip, "assembly/old_sps_seqs.pklbin").c_str()) <= 0
        and contigShifts.consensus.loadPklBin(getProjPath(ip, "assembly/sps_seqs.pklbin").c_str()) <= 0)
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
    DEBUG_MSG("Loading contigAbinfo...");
    if (!Load_abinfo(getProjPath(ip, "assembly/component_info.bin").c_str(), contigAbinfo))
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
  }

  // Test for empty return data structures
  TEST_RETURN("ExecAssembly", contigShifts);

  if (finalStage == STAGE_ASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  DEBUG_MSG("Reloading Star specturm");
  // Need to reload the star spectra because they were altered by Assembly
  if (starSpectra.loadPklBin(getProjPath(ip, "./spectra/stars.pklbin").c_str())
      <= 0)
  {
    ERROR_MSG("Problem encountered while reloading star specturm");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_ASSEMBLY);
  }
  DEBUG_VAR(starSpectra.size());

  if (initialStage <= STAGE_METAASSEMBLY
      and ip.getValueInt("MIN_METACONTIG_SIZE", 0) > 0)
  {
    SpectrumPairSet contigPairs;
    Clusters metaContigs;
    abinfo_t metaContigAbinfo;

    if (!performContigAlignment(ip, contigShifts, contigPairs))
    {
      ERROR_MSG("Problem encountered during MetaAssembly alignment stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }

    if (!performMetaAssembly(ip,
                             contigShifts,
                             contigPairs,
                             contigAbinfo,
                             starSpectra,
                             metaContigs,
                             metaContigAbinfo))
    {
      ERROR_MSG("Problem encountered during MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }

    if (commandLineParams.exists("SINGLE_STEP")) {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

    contigShifts = metaContigs;
    contigAbinfo = metaContigAbinfo;
  }
  else if (ip.getValueInt("MIN_METACONTIG_SIZE", 0) > 0)
  {
    DEBUG_MSG("Bypassing MetaAssembly stage");
    if (contigShifts.Load(getProjPath(ip, "assembly/path_spectra_as_cluster.txt").c_str()) <= 0)
    {
      ERROR_MSG("Problem encountered while skipping MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }
    if (!Load_abinfo(getProjPath(ip, "assembly/component_info.bin").c_str(), contigAbinfo))
    {
      ERROR_MSG("Problem encountered while skipping MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }
  } else if (initialStage == STAGE_METAASSEMBLY &&
             ip.getValueInt("MIN_METACONTIG_SIZE", 0) <= 0) {
      ERROR_MSG("Initial stage set to META-ASSEMBLY but MIN_METACONTIG_SIZE <= 0");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
  }

  if (finalStage == STAGE_METAASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  DEBUG_VAR(contigShifts.size());

  TEST_RETURN("ExecMetaAssembly", contigShifts);

  DEBUG_MSG("Reloading Star specturm");
  // Need to reload the star spectra because they were altered by Meta Assembly
  if (starSpectra.loadPklBin(getProjPath(ip, "./spectra/stars.pklbin").c_str())
      <= 0)
  {
    ERROR_MSG("Problem encountered while reloading star specturm");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_ASSEMBLY);
  }
  DEBUG_VAR(starSpectra.size());


  DB_fasta dbDecoys;
  if (!ip.exists("FASTA_DATABASE"))
  {
    ERROR_MSG("FASTA_DATABASE not set");
    return (0);
  }
  if (dbDecoys.Load(dbFileName.c_str()) <= 0)
  {
    ERROR_MSG("Unable to load FASTA_DATABASE [" << dbFileName << "]");
    return (0);
  }
  dbDecoys.replaceDecoyShuffled();
  
  if (ip.exists("FASTA_DATABASE"))
  {
    // Test for empty return data structures
    //TEST_RETURN("ExecTagSearch", alignedSpectra);
    DEBUG_VAR(contigShifts.size());

    //---------------------------------------------------------------------------
    // Contig/Protein Alignment STAGE
    //---------------------------------------------------------------------------
    SpecSet matchedContigsAll; // Matched versions of spectra significantly matched to a protein (All)
    SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein

    PeptideSpectrumMatchSet psmSetContig;
    PeptideSpectrumMatchSet psmSetContigDecoy;
    PeptideSpectrumMatchSet psmSetContigFdr;

    if (initialStage <= STAGE_CONTIGPROTALIGN)
    {
      DEBUG_VAR(contigShifts.size());
      if (!performContigProtAlign(ip,
                                  contigShifts.consensus,
                                  dbAll,
                                  dbDecoys,
                                  penaltyMatrixBlosum,
                                  penaltyMatrixMods,
                                  scanSpecificPenalties,
                                  contigAbinfo,
                                  filterPsmSet,
                                  matchedContigsAll,
                                  matchedContigs,
                                  psmSetContig,
                                  psmSetContigDecoy,
                                  psmSetContigFdr,
                                  gridExecutionFlag,
                                  resumeFlag))
      {
        ERROR_MSG("Problem encountered during ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
          && (!resumeFlag && !gridExecutionFlag))
      {
        // If we are doing a grid execution (and are not actually on the grid)
        //    and we aren't resuming... then exit (we'll resume execution later)
        DEBUG_MSG("Files for grid execution have been saved.");
        DEBUG_MSG("Restart with -z option when grid execution has been completed.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Contig/Protein Alignment stage");
      if (matchedContigs.loadPklBin(getProjPath(ip, "homology/contigs.pklbin").c_str(),
                                    getProjPath(ip, "homology/contigs_psm.txt").c_str(),
                                    getProjPath(ip,
                                                "homology/contigs_midx.pklbin").c_str())
          <= 0)
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      psmSetContig.getPSMSet(&matchedContigs);

      DEBUG_MSG("Loading " << getProjPath(ip, "homology/contigs.pklbin").c_str());
      if (matchedContigsAll.loadPklBin(getProjPath(ip,
                                                   "homology/contigs.pklbin").c_str())
          <= 0)
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      DEBUG_MSG("Loading " << getProjPath(ip, "homology/contigs_psm.txt").c_str());
      if (!psmSetContig.loadFromFile(getProjPath(ip, "homology/contigs_psm.txt").c_str()))
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      DEBUG_MSG("Loading " << getProjPath(ip, "homology/contigs_psm_dec.txt").c_str());
      if (!psmSetContigDecoy.loadFromFile(getProjPath(ip,
                                                      "homology/contigs_psm_dec.txt").c_str()))
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }
    }

    //SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein
    //SpecSet matchedMassIndices; // Per-spectrum indices of matched spectrum/protein masses

    // Test for empty return data structures
    //TEST_RETURN("ExecContigProtAlign", matchedContigs);
    DEBUG_VAR(contigShifts.size());
    DEBUG_VAR(matchedContigs.size());
    DEBUG_VAR(matchedContigsAll.size());

    if (finalStage == STAGE_CONTIGPROTALIGN)
    {
      DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

    //---------------------------------------------------------------------------
    // Spectrum Tag Generation STAGE
    //---------------------------------------------------------------------------
    if (initialStage <= STAGE_SPECTAGGENERATE
        && ip.getValueInt("ALIGN_STARS", 0) > 0)
    {
    
      if (!performSpecTagGen(ip,
                             contigAbinfo,
                             contigShifts.consensus,
                             starSpectra,
                             matchedContigsAll,
                             psmSetContig,
                             psmSetContigDecoy))
      {
        ERROR_MSG("Problem encountered during SpecTagGen stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_SPECTAGGENERATE);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }

    } else {
      // Nothing to do here - output is used from file always
    }

    //---------------------------------------------------------------------------
    // Spectrum/Protein Alignment STAGE
    //---------------------------------------------------------------------------
    if (initialStage <= STAGE_SPECPROTALIGN
        && ip.getValueInt("ALIGN_STARS", 0) > 0)
    {
      PeptideSpectrumMatchSet psmSetContigTag;
      PeptideSpectrumMatchSet psmSetContigTagDecoy;
      SpecSet starPrmAll; // Matched versions of spectra significantly matched to a protein (All)
      SpecSet starPrm; // Matched versions of spectra significantly matched to a protein
      PeptideSpectrumMatchSet psmSetSpectra;
      PeptideSpectrumMatchSet psmSetSpectraDecoy;
      PeptideSpectrumMatchSet psmSetSpectraFdr;

      DEBUG_TRACE;
      if (!performSpecProtAlign(ip, starSpectra, prmSpectra, // PRM spectra for "sprinkling" into star spectra
                                dbAll,
                                dbDecoys,
                                psmSetContigTag,
                                psmSetContigTagDecoy,
                                penaltyMatrixBlosum,
                                penaltyMatrixMods,
                                scanSpecificPenalties,
                                contigAbinfo,
                                filterPsmSet,
                                starPrmAll,
                                starPrm,
                                psmSetSpectra,
                                psmSetSpectraDecoy,
                                psmSetSpectraFdr,
                                gridExecutionFlag,
                                resumeFlag))
      {
        ERROR_MSG("Problem encountered during SpecProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_SPECPROTALIGN);
      }

      if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
          && (!resumeFlag && !gridExecutionFlag))
      {
        // If we are doing a grid execution (and are not actually on the grid)
        //    and we aren't resuming... then exit (we'll resume execution later)
        DEBUG_MSG("Files for grid execution have been saved.");
        DEBUG_MSG("Restart with -z option when grid execution has been completed.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return (0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Spectrum/Protein Alignment stage");
    }

    if (finalStage == STAGE_SPECPROTALIGN)
    {
      DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

  //---------------------------------------------------------------------------
  // Protein/Protein Alignment STAGE
  //---------------------------------------------------------------------------
  DEBUG_TRACE;
  string indexFileName; // File names of clustalw .aln output files
  if (initialStage <= STAGE_PROTPROTALIGN)
  {

    //---------------------------------------------------------------------------
    // Find reference protein for protein/protein sequence alignments
    //---------------------------------------------------------------------------
    unsigned int refProtIdx = 0, refProtMatchCount = 0;
    set<unsigned int> dbIndexes;

    //---------------------------------------------------------------------------
    // Assign the reference protein for protein/protein sequence alignments
    //---------------------------------------------------------------------------
    if (ip.exists("FORCE_REFERENCE"))
    {
      DEBUG_MSG("POINT A-prime");
      // Force use of the supplied reference
      refProtIdx = ip.getValueInt("FORCE_REFERENCE");
      DEBUG_MSG("Point B-prime");
      for (unsigned int i = 0; i < matchedContigs.size(); i++)
      {
        DEBUG_MSG("Point C-prime");
        dbIndexes.insert(matchedContigs[i].psmList.front()->m_dbIndex);
      }
      DEBUG_VAR(refProtIdx);

    }
    else
    {
      // Use the most common protein match as the reference
      DEBUG_MSG("Point A.a");
      vector<unsigned int> numMatches(dbAll.size());
      DEBUG_MSG("Point A.b");
      for (unsigned int i = 0; i < dbAll.size(); i++)
        numMatches[i] = 0;
      DEBUG_MSG("Point A");
      for (unsigned int i = 0; i < matchedContigs.size(); i++)
      {
        DEBUG_VAR(matchedContigs.size());
        DEBUG_VAR(i);
        DEBUG_VAR(matchedContigs[i].psmList.size());
        if (matchedContigs[i].psmList.size() > 0)
        {
          DEBUG_VAR(matchedContigs[i].psmList.front());
          int index = matchedContigs[i].psmList.front()->m_dbIndex;
          DEBUG_MSG("Point B");
          dbIndexes.insert(index);
          DEBUG_MSG("Point C");
          if (index >= 0 and index < dbAll.size())
          {
            DEBUG_MSG("Point D");
            //    numMatches[index]++;  // Select by highest number of matched contigs
            DEBUG_VAR(numMatches.size());
            DEBUG_VAR(index);

            DEBUG_VAR(matchedContigs[i].psmList.front()->m_matchedPeaks.size());
            numMatches[index] +=
                matchedContigs[i].psmList.front()->m_matchedPeaks.size(); // Select by highest number of matched masses
            if (numMatches[index] > refProtMatchCount)
            {
              DEBUG_MSG("Point E");
              refProtIdx = index;
              refProtMatchCount = numMatches[index];
            }
          }

          DEBUG_VAR(refProtIdx);
          DEBUG_VAR(refProtMatchCount);
        }
        //Since we didn't write out the PMS matches for these contigs (they are GenoMS contigs)
        else if (matchedContigs[i].psmList.size() == 0
            && (runGenoMSFlag || runMergedFlag))
        {

          int index = 0;
          DEBUG_MSG("Point B_NEC");
          dbIndexes.insert(index);
          DEBUG_MSG("Point C_NEC");
          if (index >= 0 and index < dbAll.size())
          {
            DEBUG_MSG("Point D_NEC");
            //    numMatches[index]++;  // Select by highest number of matched contigs
            DEBUG_VAR(numMatches.size());
            DEBUG_VAR(index);

            //DEBUG_VAR(matchedContigs[i].psmList.front()->m_matchedPeaks.size());
            numMatches[index] += matchedContigs[i].size(); // Select by highest number of matched masses
            if (numMatches[index] > refProtMatchCount)
            {
              DEBUG_MSG("Point E_NEC");
              refProtIdx = index;
              refProtMatchCount = numMatches[index];
            }
          }

          DEBUG_VAR(refProtIdx);
          DEBUG_VAR(refProtMatchCount);

        }
      }
    } // if (ip.exists("FORCE_REFERENCE")) {

    if (!performProtProtAlign(ip, refProtIdx, dbIndexes, dbAll))
    {
      ERROR_MSG("Problem encountered during ProtProtAlign stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_PROTPROTALIGN);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Protein/Protein Alignment stage");
  }

    DEBUG_VAR(dbAll.size());
    DEBUG_VAR(matchedContigs.size());

    //---------------------------------------------------------------------------
    // Homology-based Assembly STAGE
    //---------------------------------------------------------------------------
    SpecSet cspsContigs;
    SpecSet cspsMatchedIndices;
    vector<vector<int> > cspsMatchedProts;
    if (initialStage <= STAGE_HOMOLOGYASSEMBLY)
    {
      DEBUG_TRACE;
      if (!performHomologyAssembly(ip,
              matchedContigs,
              dbAll,
              cspsContigs,
              cspsMatchedIndices,
              cspsMatchedProts))
      {
        ERROR_MSG("Problem encountered during cSPS HomologyAssembly stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_HOMOLOGYASSEMBLY);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        return(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Homology Assembly stage");
      if (cspsContigs.loadPklBin(getProjPath(ip,
                  "homology/homglue_matches.pklbin").c_str())
          <= 0
          or cspsMatchedIndices.loadPklBin(getProjPath(ip,
                  "homology/homglue_matches_midx.pklbin").c_str())
          <= 0
          or Load_binArray<int> (getProjPath(ip,
                  "homology/homglue_matches_mp.bin").c_str(),
              cspsMatchedProts) <= 0)
      {
        ERROR_MSG("Problem encountered while skipping SpecProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_HOMOLOGYASSEMBLY);
      }
    }

    // Test for empty return data structures
    TEST_RETURN("ExecHomologyAssembly", cspsContigs);
    TEST_RETURN("ExecHomologyAssembly", cspsMatchedIndices);
    TEST_RETURN("ExecHomologyAssembly", cspsMatchedProts);

    //---------------------------------------------------------------------------
    // MERGE RESULTS IF RUNNING IN INTEGRATIVE MODE
    //---------------------------------------------------------------------------
    if (runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)
    {
      //We need to rename all of our results files
      if (rename("./spectra/stars.pklbin", "./spectra/csps.stars.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./assembly/sps_seqs.pklbin",
                 "./assembly/csps.sps_seqs.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./assembly/component_info.bin",
                 "./assembly/csps.component_info.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_mp.bin", "./homology/csps.contigs_mp.bin")
          < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_midx.pklbin",
                 "./homology/csps.contigs_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_mp_all.bin",
                 "./homology/csps.contigs_mp_all.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      /*if (rename("./homology/contigs_midx_all.pklbin",
       "./homology/csps.contigs_midx_all.pklbin") < 0)
       {
       ERROR_MSG("Problem encountered renaming CSPS files");
       writeStatusFile(statusFileName, "Error");
       exit(-STAGE_MERGE);
       }*/

      if (!FileCopy("homology/contigs.pklbin",
                    "spectra/csps.contigs.pklbin"))
      {
        ERROR_MSG("Problem encountered copying contig.pklbin");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./spectra/contigs_indices.bin",
                 "./spectra/csps.contigs_indices.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_ref_mp.bin",
                 "./homology/csps.homglue_ref_mp.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_ref_midx.pklbin",
                 "./homology/csps.homglue_ref_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches.pklbin",
                 "./homology/csps.homglue_matches.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches_mp.bin",
                 "./homology/csps.homglue_matches_mp.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches_midx.pklbin",
                 "./homology/csps.homglue_matches_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      ifstream ifile("homology/ref_sps_names.txt");
      if (ifile)
      {

        if (rename("./homology/ref_sps_names.txt",
                   "./homology/csps.ref_sps_names.txt") < 0)
        {
          ERROR_MSG("Problem encountered renaming CSPS files");
          writeStatusFile(statusFileName, "Error");
          exit(-STAGE_MERGE);
        }
      }

    } // if(runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)

    if (runMergedFlag && initialStage <= STAGE_MERGE)
    {
      if (!performMergeOfCSPSAndGenoMS(ip))
      {
        ERROR_MSG("Problem encounctered during merge of CSPS and GenoMS");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
    }
  }
  else // if(ip.exists("FASTA_DATABASE"))
  {
    //
    //  --- need to fix performReports to allow reporting de novo results only
    //
    //  WARN_MSG("No FASTA_DATABASE specified, reporting only de novo results.\n");
    WARN_MSG("No FASTA_DATABASE specified, returning without database matches.\n");
    return 0;
  }

  if (finalStage == STAGE_HOMOLOGYASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  bool generatedStats = false;
  if (initialStage <= STAGE_STATPROTSEQS
      && (ip.exists("INPUT_SPEC_IDS") || ip.exists("INPUT_CLUST_SPEC_IDS")))
  {
    if (!performStatProtSeqs(ip))
    {
      ERROR_MSG("Problem encountered during StatProtSeqs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_STATPROTSEQS);
    }
    generatedStats = true;
    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }
  }

  if (generatedStats && finalStage == STAGE_STATPROTSEQS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    return (0);
  }

  //---------------------------------------------------------------------------
  // REPORT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_REPORT)
  {
    ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                     ios_base::out | ios_base::binary);
    spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
        << ip.getValue("TOLERANCE_PM") << "\n";
    spsProj.close();

    if (!performReport(ip))
    {
      ERROR_MSG("Problem encountered during Report stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_REPORT);
    }
    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      return (0);
    }

  }
  else
  {
    DEBUG_MSG("Bypassing reports stage");
  }
  //--------------------------------------------------------------------------
  //Set up for relaunching
  //--------------------------------------------------------------------------
  if (!generateRelaunchScript(ip))
  {
    ERROR_MSG("Problem encountered during relaunch script creation");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_REPORT);
  }

  //---------------------------------------------------------------------------
  // END
  //---------------------------------------------------------------------------
  writeStatusFile(statusFileName, "Finished");
  return 0;
} // END

