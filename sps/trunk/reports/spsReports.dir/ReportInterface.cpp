///////////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include "ReportInterface.h"

#include "ParameterList.h"
#include "CommandLineParser.h"
#include "Specific.h"

#include "Logger.h"


#include "ReportRendererHtml.h"
#include "ReportRendererHtmlClient.h"
#include "ReportRendererHtmlSingle.h"
#include "ReportRendererHtmlDynamic.h"
#include "ReportRendererPdf.h"

#include "ReportConfig.h"

#include "Timer.h"
#include "Tokenizer.h"

#include "copyright.h"

#include "hpdf.h"


#define TSV_SUFFIX  ".tsv"

#define TABLES_DIR_DEFAULT "ReportData"

///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
using namespace std;

namespace spsReports {



//void error_handler  (HPDF_STATUS   error_no,
//                HPDF_STATUS   detail_no,
//                void         *user_data)
//{
//    printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT)error_no,
//                (HPDF_UINT)detail_no);
//}



///////////////////////////////////////////////////////////////////////////////
ReportInterface::ReportInterface() :
   m_tableNameHeader("tableHeader.txt"),
   m_tableNameProtein("tableProtein.txt"),
   m_tableNameProteinCoverage("tableProteinCoverage.txt"),
   m_tableNameContig("tableContig.txt"),
   m_tableNameCluster("tableCluster.txt"),
   m_tableNameSpectra("tableSpectra.txt"),
   m_tableConfig("tableConfig.txt"),
   m_verbose(false)
{
}
///////////////////////////////////////////////////////////////////////////////
ReportInterface::~ReportInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
string ReportInterface::composeFileName(const string &projectDir, const string &fileName)
{
  // Compose output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  // add filename
  aux += fileName;
  // return composed filename
  return aux;
}

/////////////////////
int ReportInterface::buildDirectoryPath(const string &dir)
{
  if (! dir.empty())
    for (string::size_type i = 0, j = min(dir.find('/', i), dir.size()); i < dir.size(); i = j + 1, j = min(dir.find('/', i), dir.size()))
      if (! dir.substr(0, j).empty())
        if (MKDIR(dir.substr(0, j).c_str()) == -1)
        //if (mkdir(dir.substr(0, j).c_str()) == -1)
          if (errno != EEXIST) {
            cerr << "error: cannot create directory '" << dir.substr(0, j) << "': " << strerror(errno) << endl;
            return ERROR;
          }
   return OK;
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::dump_abruijn(ReportGeneratorData &data)
{
  // load needed data files
  SpsFiles spsFiles;
  spsFiles.loadData(data);
  // Dump
  spsFiles.dump_abruijn(cout);
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::dump_binArray(ReportGeneratorData &data)
{
  // load needed data files
  SpsFiles spsFiles;
  spsFiles.loadData(data);
  // generate dump
  spsFiles.dump(cout);
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::buildTables(ReportGeneratorData &data)
{
  Timer_c timer;

  // verbose output
  DEBUG_MSG("++++ Generating report tables ++++");

  // build directory path for tables
  if(buildDirectoryPath(data.tablesDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report tables.");
    return ERROR;
  }

  // set the table filenames
  ReportTableHeader           tHeader(data.tablesDir,           m_tableNameHeader);
  ReportTableProtein          tProtein(data.tablesDir,          m_tableNameProtein);
  ReportTableProteinCoverage  tProteinCoverage(data.tablesDir,  m_tableNameProteinCoverage);
  ReportTableContig           tContig(data.tablesDir,           m_tableNameContig);
  ReportTableClusterConsensus tCluster(data.tablesDir,          m_tableNameCluster);
  ReportTableInputSpectra     tInputSpectra(data.tablesDir,     m_tableNameSpectra);

  // set the configuration file filename
  ReportConfig config(data.tablesDir, m_tableConfig);

  // load needed data files
  SpsFiles spsFiles;
  if(spsFiles.loadData(data) == ERROR) {
    ERROR_MSG("Needed files were not loaded");
    return ERROR;
  }

  // table generator object
  ReportTableGenerator rtg;
  // set data files
  rtg.setSpsFiles(&spsFiles);
  // init
  rtg.init(data);
  // generate the tables
  int ret = rtg.buildTables(&tHeader, &tProtein, &tProteinCoverage, &tContig, &tCluster, &tInputSpectra);
  if(ret == ERROR) {
    ERROR_MSG("Error generating tables.");
    return ERROR;
  }

  // save the tables
  DEBUG_MSG("Saving the tables");
  tHeader.saveTable();
  tProtein.saveTable();
  tProteinCoverage.saveTable();
  tContig.saveTable();
  tCluster.saveTable();
  tInputSpectra.saveTable();

  // save the tables is tsv format
  DEBUG_MSG("Saving the tables is TSV format");
  string aux = m_tableNameProtein;
  aux += TSV_SUFFIX;
  tProtein.setTableFilename(aux.c_str());

  aux = m_tableNameProteinCoverage;
  aux += TSV_SUFFIX;
  tProteinCoverage.setTableFilename(aux);

  aux = m_tableNameContig;
  aux += TSV_SUFFIX;
  tContig.setTableFilename(aux);

  aux = m_tableNameCluster;
  aux += TSV_SUFFIX;
  tCluster.setTableFilename(aux);

  aux = m_tableNameSpectra;
  aux += TSV_SUFFIX;
  tInputSpectra.setTableFilename(aux);

  tProtein.saveTable(FILE_SEPARATOR_TSV, true);
  tProteinCoverage.saveTable(FILE_SEPARATOR_TSV, true);
  tContig.saveTable(FILE_SEPARATOR_TSV, true);
  tCluster.saveTable(FILE_SEPARATOR_TSV, true);
  tInputSpectra.saveTable(FILE_SEPARATOR_TSV, true);

  // save config data
  config.setPwd(data.pwd);
  config.writeFile();

  // done
  DEBUG_MSG("Tables saved.");
  // verbose output
  DEBUG_MSG("---- Report tables generation took " << timer.stop());

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportPdf(ReportGeneratorData &data)
{
  Timer_c timer;

  DEBUG_MSG("++++ Building PDF report ++++");

  if(buildDirectoryPath(data.outDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report.");
    return ERROR;
  }

  // static report renderer for contigs page
  ReportRendererPdf rrPDF;
  // init report
  int ret = rrPDF.generateReport(data);
  // verbose output
  DEBUG_MSG("---- Full report generation took " << timer.stop());

  return ret;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtmlDynamic(ReportGeneratorData &data)
{
  // start timers
  Timer_c timer;

  DEBUG_MSG("++++ Building dynamic report ++++");

  if(buildDirectoryPath(data.outDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report.");
    return ERROR;
  }

  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  // static report renderer for contigs page
  ReportRendererHtmlDynamic rrHtmlDyn;
  // Specify project directory (needed for input datat)
  int ret =  rrHtmlDyn.generateReport(data);
  // verbose output
  DEBUG_MSG("---- Full report generation took " << timer.stop());

  return ret;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtmlSingle(ReportGeneratorData &data)
{
  // start timers
  Timer_c timer;

  DEBUG_MSG("++++ Building static HTML report in a single file ++++");

  if(buildDirectoryPath(data.outDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report.");
    return ERROR;
  }

  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  // load needed data files
  SpsFiles spsFiles;
  if(spsFiles.loadData(data) == ERROR) {
    ERROR_MSG("Needed files were not loaded");
    return ERROR;
  }

  // static report renderer for contigs page
  ReportRendererHtmlSingle rrHtmlSingle;
  // set data files
  rrHtmlSingle.setSpsFiles(&spsFiles);
  // Specify project directory (needed for input datat)
  int ret =  rrHtmlSingle.generateReport(data);
  // verbose output
  DEBUG_MSG("---- Full report generation took " << timer.stop());

  return ret;
}///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtmlClient(ReportGeneratorData &data)
{
  // start timers
  Timer_c timer;

  DEBUG_MSG("++++ Building static HTML report with pagination on client ++++");

  if(buildDirectoryPath(data.outDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report.");
    return ERROR;
  }

  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  // load needed data files
  SpsFiles spsFiles;
  if(spsFiles.loadData(data) == ERROR) {
    ERROR_MSG("Needed files were not loaded");
    return ERROR;
  }

  // static report renderer for contigs page
  ReportRendererHtmlClient rrHtmlClient;
  // set data files
  rrHtmlClient.setSpsFiles(&spsFiles);
  // Specify project directory (needed for input datat)
  int ret =  rrHtmlClient.generateReport(data);
  // verbose output
  DEBUG_MSG("---- Full report generation took " << timer.stop());

  return ret;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtml(ReportGeneratorData &data)
{

  DEBUG_MSG("++++ Building HTML report ++++");

  if(buildDirectoryPath(data.outDir) == ERROR) {
    ERROR_MSG("Could not build directory path for report.");
    return ERROR;
  }

  // start timers
  Timer_c timer, timer2;

  // load needed data files
  SpsFiles spsFiles;
  if(spsFiles.loadData(data) == ERROR) {
    ERROR_MSG("Needed files were not loaded");
    return ERROR;
  }
  DEBUG_MSG("---- File loading section took " << timer2.stop());

  // static report renderer for contigs page
  ReportRendererHtml rrHtml;
  // set data files
  rrHtml.setSpsFiles(&spsFiles);
  // render the report
  int ret = rrHtml.generateReport(data);
  // verbose output
  if(data.cpu == 1) DEBUG_MSG("---- Full report generation took " << timer.stop());

  return ret;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::parseOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------
  // Get the execultable directory.
  string str = argv[0];
  size_t found;
  found = str.find_last_of("/\\");
  string exeDir = ".";
  if(found != string::npos)
    exeDir = str.substr(0, found);

  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
  listOptions.push_back(CommandLineParser::Option("-debug",           "debug",            false));
  listOptions.push_back(CommandLineParser::Option("-verbose",         "verbose",          false));

  // Data tables names and project location
  listOptions.push_back(CommandLineParser::Option("-exe-dir",         "EXE_DIR",          true));
  listOptions.push_back(CommandLineParser::Option("-project-dir",     "PROJECT_DIR",      true));
  listOptions.push_back(CommandLineParser::Option("-tables-dir",      "TABLES_DIR",       true));
  listOptions.push_back(CommandLineParser::Option("-table-header",    "TABLE_HEADER",     true));
  listOptions.push_back(CommandLineParser::Option("-table-protein",   "TABLE_PROTEIN",    true));
  listOptions.push_back(CommandLineParser::Option("-table-contig",    "TABLE_CONTIG",     true));
  listOptions.push_back(CommandLineParser::Option("-table-cluster",   "TABLE_CLUSTER",    true));
  listOptions.push_back(CommandLineParser::Option("-table-spectra",   "TABLE_SPECTRA",    true));

  listOptions.push_back(CommandLineParser::Option("-project-dir-server", "PROJECT_DIR_SERVER", true));


  // input data files

  listOptions.push_back(CommandLineParser::Option("-abruijn",           "FILE_COMP",        true));
  listOptions.push_back(CommandLineParser::Option("-sps-seqs",          "FILE_SEQS",        true));
  listOptions.push_back(CommandLineParser::Option("-contigs-mp",        "FILE_MP2",         true));
  listOptions.push_back(CommandLineParser::Option("-contigs-midx",      "FILE_MIDX2",       true));
  listOptions.push_back(CommandLineParser::Option("-homglue-ref-mp",    "FILE_REFMP",       true));
  listOptions.push_back(CommandLineParser::Option("-homglue-ref-midx",  "FILE_REFMIDX",     true));
  listOptions.push_back(CommandLineParser::Option("-stars",             "FILE_STARS",       true));
  listOptions.push_back(CommandLineParser::Option("-consensus-spectra", "FILE_MS",          true));
  listOptions.push_back(CommandLineParser::Option("-input-spectra-list","FILE_INDEX",       true));
  listOptions.push_back(CommandLineParser::Option("-cluster-file",      "FILE_CLUSTER",   true));
  listOptions.push_back(CommandLineParser::Option("-cluster-ms",        "FILE_CLUSTERMS",   true));
  listOptions.push_back(CommandLineParser::Option("-cluster-scan",      "FILE_CLUSTERSCAN", true));
  listOptions.push_back(CommandLineParser::Option("-proteins",          "FILE_FASTA",       true));


  // variables
  listOptions.push_back(CommandLineParser::Option("-title",             "TITLE",            true));

  listOptions.push_back(CommandLineParser::Option("-job",               "JOB",              true));
  listOptions.push_back(CommandLineParser::Option("-user",              "USER",             true));
  listOptions.push_back(CommandLineParser::Option("-password",          "REPORT_PWD",       true));
  listOptions.push_back(CommandLineParser::Option("-server",            "SERVER",           true));
  listOptions.push_back(CommandLineParser::Option("-cells-per-line",    "CELLS_PER_LINE",   true));

  listOptions.push_back(CommandLineParser::Option("-peakmasstol",       "TOLERANCE_PEAK",   true));
  listOptions.push_back(CommandLineParser::Option("-parentmasstol",     "TOLERANCE_PM",     true));
  listOptions.push_back(CommandLineParser::Option("-resolution",        "RESOLUTION",       true));

  listOptions.push_back(CommandLineParser::Option("-outdir",            "REPORT_DIR",       true));
  listOptions.push_back(CommandLineParser::Option("-prefix",            "PREFIX",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",           "OUTFILE",          true));
  listOptions.push_back(CommandLineParser::Option("-format",            "FORMAT",           true));
  listOptions.push_back(CommandLineParser::Option("-target",            "TARGET",           true));

  listOptions.push_back(CommandLineParser::Option("-no-msms-images",    "NO_MSMS_IMAGES",   true));
  listOptions.push_back(CommandLineParser::Option("-tool",              "TOOL",             true));
  listOptions.push_back(CommandLineParser::Option("-no-clusters",       "NO_CLUSTERS",      false));
  listOptions.push_back(CommandLineParser::Option("-cpu-count",         "CPUS",             true));

  // models and mass shifts
  listOptions.push_back(CommandLineParser::Option("-aa",                "AMINO_ACID_MASSES",true));
  listOptions.push_back(CommandLineParser::Option("-annotation-model",  "ANNOTATION_MODEL", true));
  listOptions.push_back(CommandLineParser::Option("-shift-value",       "SHIFT_VALUE",      true));

  listOptions.push_back(CommandLineParser::Option("-annotation-model-prm",  "ANNOTATION_MODEL_PRM", true));
  listOptions.push_back(CommandLineParser::Option("-shift-value-prm",       "SHIFT_VALUE_PRM",      true));


  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",                 "PARAMETER_FILE",   true));

  // Build tables command
  listOptions.push_back(CommandLineParser::Option("-build-tables",      "buildTables",      false));

  // Build report command
  listOptions.push_back(CommandLineParser::Option("-build-html",        "REPORT_HTML_TYPE",  true));

  // Build report command
  listOptions.push_back(CommandLineParser::Option("-build-pdf",         "REPORT_PDF",       false));


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


  // add EXE_DIR parameter
  commandLineParams.addIfDoesntExist("EXE_DIR", exeDir);


  ////////////////////////////////////////////////////////////////////////////////
  // Process the commands
  return processOptions(commandLineParams);

}
//////////////////////////////////////////////////////////////////////////////////
// Process options
//////////////////////////////////////////////////////////////////////////////////
int ReportInterface::processOptions(ParameterList &ip)
{
  ParameterList commandLineParams;

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (ip.exists("PARAMETER_FILE")) {

    string parameterFilename = ip.getValue("PARAMETER_FILE");

    commandLineParams.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag set)
  }

  commandLineParams.addList(ip, true);

  ////////////////////////////////////////////////////////////////////////////////
  // help message control
  bool showHelp = true;

  // Data holder
  ReportGeneratorData data;

  ////////////////////////////////////////////////////////////////////////////////
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // sets the "verbose" flag
  if (commandLineParams.exists("verbose"))
    m_verbose = true;

  ////////////////////////////////////////////////////////////////////////////////
  // Boolean parameters


  ////////////////////////////////////////////////////////////////////////////////
  // File locations
  ////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  // File load section

  // Executables directory
  if (commandLineParams.exists("EXE_DIR"))
    data.exeDir = commandLineParams.getValue("EXE_DIR");

  // Project directory
  if (commandLineParams.exists("PROJECT_DIR")) {
    data.projectDir  = commandLineParams.getValue("PROJECT_DIR");
    data.outDir      = composeFileName(data.projectDir, data.outDir);
  }

  // Project directory on server
  if (commandLineParams.exists("PROJECT_DIR_SERVER"))
    data.targetProjectDir = commandLineParams.getValue("PROJECT_DIR_SERVER");
  else
    data.targetProjectDir = data.projectDir;

  // Server cgi-bin location
  if (commandLineParams.exists("CELLS_PER_LINE"))
    data.cellsPerLine  = getInt(commandLineParams.getValue("CELLS_PER_LINE").c_str());

  // Server cgi-bin location
  if (commandLineParams.exists("SERVER"))
    data.server  = commandLineParams.getValue("SERVER");


  // Tables directory
  if (commandLineParams.exists("TABLES_DIR"))
    data.tablesDir = commandLineParams.getValue("TABLES_DIR");
  else
    data.tablesDir   = TABLES_DIR_DEFAULT;


  // table file names
  if (commandLineParams.exists("TABLE_HEADER"))
    m_tableNameHeader = commandLineParams.getValue("TABLE_HEADER");

  if (commandLineParams.exists("TABLE_PROTEIN"))
    m_tableNameProtein = commandLineParams.getValue("TABLE_PROTEIN");

  if (commandLineParams.exists("TABLE_CONTIG"))
    m_tableNameContig = commandLineParams.getValue("TABLE_CONTIG");

  if (commandLineParams.exists("TABLE_CLUSTER"))
    m_tableNameCluster = commandLineParams.getValue("TABLE_CLUSTER");

  if (commandLineParams.exists("TABLE_SPECTRA"))
    m_tableNameSpectra = commandLineParams.getValue("TABLE_SPECTRA");

  // configuration file name
  if (commandLineParams.exists("TABLE_CONFIG"))
    m_tableConfig = commandLineParams.getValue("TABLE_CONFIG");


  // Font file location
  if (commandLineParams.exists("FONT_LOCATION"))
    ;


  data.aasFileDir = data.exeDir;
  data.annotationModelDir = data.exeDir;


  ////////////////////////////////////////////////////////////////////////////////
  // Annotation models & mass shifts

  // Amino-acid masses file
  if (commandLineParams.exists("AMINO_ACID_MASSES")) {
    string aux = commandLineParams.getValue("AMINO_ACID_MASSES");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      data.aasFile = aux.substr(found + 1);
      data.aasFileDir = aux.substr(0, found);
    } else {
      data.aasFile = aux;
      data.aasFileDir = "." ;
    }
  }

  // default

  if(commandLineParams.exists("ANNOTATION_MODEL")) {
    string aux = commandLineParams.getValue("ANNOTATION_MODEL");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      data.annotationModel = aux.substr(found + 1);
      data.annotationModelDir = aux.substr(0, found);
    } else {
      data.annotationModelDir = ".";
      data.annotationModel = aux;
    }
  }

  // Annotation model for PRM

  if(commandLineParams.exists("ANNOTATION_MODEL_PRM")) {
    string aux = commandLineParams.getValue("ANNOTATION_MODEL_PRM");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      data.annotationModelPrm = aux.substr(found + 1);
      data.annotationModelDirPrm = aux.substr(0, found);
    } else {
      data.annotationModelDirPrm = ".";
      data.annotationModelPrm = aux;
    }
  }


  if (commandLineParams.exists("SHIFT_VALUE"))
    data.massShift += commandLineParams.exists("SHIFT_VALUE");

  if (commandLineParams.exists("SHIFT_VALUE_PRM"))
    data.massShiftPrm += commandLineParams.exists("SHIFT_VALUE_PRM");


  ////////////////////////////////////////////////////////////////////////////////
  //  Input data section - files needed to generate the tables

  // abruijn graph
  if(commandLineParams.exists("FILE_COMP")) {
    data.filenameAbruijn = commandLineParams.getValue("FILE_COMP");
    //data.absoluteFilenameAbruijn = true;
  }

  if(commandLineParams.exists("FILE_SEQS")) {
    data.filenameSpsSeqs = commandLineParams.getValue("FILE_SEQS");
    //data.absoluteFilenameSpsSeqs = true;
  }

  if(commandLineParams.exists("FILE_MP2")) {
    data.filenameContigsMp = commandLineParams.getValue("FILE_MP2");
    //data.absoluteFilenameContigsMp = true;
  }

  if(commandLineParams.exists("FILE_MIDX2")) {
    data.filenameContigsMidx = commandLineParams.getValue("FILE_MIDX2");
    //data.absoluteFilenameContigsMidx = true;
  }

  if(commandLineParams.exists("FILE_REFMP")) {
    data.filenameHomglueRefMp = commandLineParams.getValue("FILE_REFMP");
    //data.absoluteFilenameHomglueRefMp = true;
  }

  if(commandLineParams.exists("FILE_REFMIDX")) {
    data.filenameHomglueRefMidx = commandLineParams.getValue("FILE_REFMIDX");
    //data.absoluteFilenameHomglueRefMidx = true;
  }

  if(commandLineParams.exists("FILE_STARS")) {
    data.filenameStarSpectra = commandLineParams.getValue("FILE_STARS");
    //data.absoluteFilenameStarSpectra = true;
  }

  if(commandLineParams.exists("FILE_MS")) {
    data.filenameConsensusSpectra = commandLineParams.getValue("FILE_MS");
    //data.absoluteFilenameConsensusSpectra = true;
  }

  if(commandLineParams.exists("FILE_INDEX")) {
    data.filenameInputSpectraList = commandLineParams.getValue("FILE_INDEX");
    //data.absoluteFilenameInputSpectraList = true;
  }

  if(commandLineParams.exists("FILE_CLUSTERMS")) {
    data.filenameClusterMS = commandLineParams.getValue("FILE_CLUSTERMS");
    //data.absoluteFilenameClusterMS = true;
  }

  if(commandLineParams.exists("FILE_CLUSTER")) {
    data.filenameCluster = commandLineParams.getValue("FILE_CLUSTER");
    //data.absoluteFilenameClusterMS = true;
  }
  if(commandLineParams.exists("FILE_CLUSTERSCAN")) {

    data.filenameScanFiles = commandLineParams.getValue("FILE_CLUSTERSCAN");
    //data.absoluteFilenameScanFiles = true;
  }

  if(commandLineParams.exists("FILE_FASTA")) {
    data.filenameProteins = commandLineParams.getValue("FILE_FASTA");
    //data.absoluteFilenameProteins = true;
  }

  if(commandLineParams.exists("MATCHED_CONTIGS")) {
    data.filenameContigSpectra = commandLineParams.getValue("MATCHED_CONTIGS");
    //data.absoluteFilenameContigSpectra = true;
  }

  if(commandLineParams.exists("MATCHED_CONTIGS_IDX")) {
    data.filenameContigIndices = commandLineParams.getValue("MATCHED_CONTIGS_IDX");
    //data.absoluteFilenameContigIndices = true;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("REPORT_DIR")) {
    data.outDir = commandLineParams.getValue("REPORT_DIR");
    //data.tablesDir = data.outDir;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters

  // spectrumInfo
  if (commandLineParams.exists("spectruminfo"))
    ;
    //ReportInterfacetrum.setSpectrumInfo(commandLineParams.getValue("spectruminfo"));

  // peak mass tolerance
  if (commandLineParams.exists("TOLERANCE_PEAK"))
    data.peakMassTol = commandLineParams.getValueFloat("TOLERANCE_PEAK");

  // parent mass tolerance
  if (commandLineParams.exists("TOLERANCE_PM"))
    data.peakMassTol = commandLineParams.getValueFloat("TOLERANCE_PM");

  // parent mass tolerance
  if (commandLineParams.exists("RESOLUTION"))
    data.resolution = commandLineParams.getValueFloat("RESOLUTION");


  ////////////////////////////////////////////////////////////////////////////////
  //  Report related parameters

  // Job name
  if (commandLineParams.exists("JOB"))
    data.job = commandLineParams.getValue("JOB");

  // User
  if (commandLineParams.exists("USER"))
    data.user = commandLineParams.getValue("USER");

  // Report password
  if (commandLineParams.exists("REPORT_PWD"))
    data.pwd = commandLineParams.getValue("REPORT_PWD");

  // tool used
  if (commandLineParams.exists("TOOL"))
    data.tool = commandLineParams.getValueInt("TOOL");


  ////////////////////////////////////////////////////////////////////////////////
  //  Image related parameters

  // Title
  if (commandLineParams.exists("TITLE"))
    ;
    //ReportInterfacetrum.setTitle(commandLineParams.getValue("TITLE"));

  if (commandLineParams.exists("NO_MSMS_IMAGES"))
    data.displayLevel = 3;


  ////////////////////////////////////////////////////////////////////////////////
  //  other parameters

  // Executables directory
  if (commandLineParams.exists("CPUS"))
    data.cpu = commandLineParams.getValueInt("CPUS");

  if (commandLineParams.exists("NO_CLUSTERS"))
    data.noClusters = true;

  if (commandLineParams.exists("ALLOW_REALIGN"))
    data.allowRealign = (bool)commandLineParams.getValueInt("ALLOW_REALIGN");

  if (commandLineParams.exists("DYNAMIC"))
    data.dynamic = (bool)commandLineParams.getValueInt("DYNAMIC");

  ////////////////////////////////////////////////////////////////////////////////
  //  Action commands

  // the debug flag
  if (commandLineParams.exists("debug")) {
    dump_abruijn(data);
    dump_binArray(data);
    return 0;
  }

  // Build text tables
  if (commandLineParams.exists("buildTables")) {
    showHelp = false;
    if(buildTables(data) == ERROR) {
      ERROR_MSG("Failed to build report tables.");
      return ERROR;
    }
  }


  // Build HTML reports
  if (commandLineParams.exists("REPORT_HTML_TYPE")) {
    showHelp = false;

    string type = commandLineParams.getValue("REPORT_HTML_TYPE");

    if(type == "dynamic") {
      if(generateReportHtmlDynamic(data) == ERROR) {
        ERROR_MSG("Failed to build html 'dynamic' report.");
        return ERROR;
      }
    }

    if(type == "client") {
      if(generateReportHtmlClient(data) == ERROR) {
        ERROR_MSG("Failed to build html 'client' report.");
        return ERROR;
      }
    }

    if(type == "single") {
      if(generateReportHtmlSingle(data) == ERROR) {
        ERROR_MSG("Failed to build html 'single' report.");
        return ERROR;
      }
    }

    if(type == "static") {
      if(generateReportHtml(data) == ERROR) {
        ERROR_MSG("Failed to build html 'static' report.");
        return ERROR;
      }
    }

  }

  // Build PDF reports
  if (commandLineParams.exists("REPORT_PDF")) {
    showHelp = false;
    if(generateReportPdf(data) == ERROR) {
      ERROR_MSG("Failed to build PDF report.");
      return ERROR;
    }
  }


  if(showHelp)
    help(cout);


  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: spsReports [OPTION]\n";
  sout << "Options:\n";
  sout << '\n';
  sout << "  --build-tables              Build report data tables\n";
  sout << "  --build-html TYPE           Build html report. TYPE = static|client|single|dynamic.\n";
  sout << "  --build-pdf                 Build PDF report\n";
  sout << '\n';
  sout << "  --exe-dir PATH              Executables directory (defaults to spsReports exec dir)\n";
  sout << "  --project-dir PATH          Project directory (defaults to current directory)\n";
  sout << "  --tables-dir PATH           Table files directory (defaults to project directory)\n";
  sout << "  --project-dir-server        Directory for project after relocation. Used by dynamic reports. Defaults to --project-dir";
  sout << '\n';
  sout << "  --tables-header FILENAME    Header table filename\n";
  sout << "  --tables-protein FILENAME   Proteins table filename\n";
  sout << "  --tables-contig FILENAME    Contigs table filename\n";
  sout << "  --tables-cluster FILENAME   Cluster consensus table filename\n";
  sout << "  --tables-spectra FILENAME   Input spectra table filename\n";
  sout << '\n';
  sout << "  --aa FILE                   Amino acids file (txt format)\n";
  sout << "  --annotation-model FILE     Annotation model file (defaults to './model_cid.txt')\n";
  sout << '\n';
  sout << "  --peakmasstol NUMBER        Amino acid peak mass tolerance (defaults to 0.45)\n";
  sout << "  --parentmasstol NUMBER      Parent mass tolerance\n";
  sout << "  --resolution NUMBER         Resolution\n";
  sout << "  --shift-value               Specify and apply mass shift\n";
  sout << '\n';
  sout << "  --abruijn FILE              Abruijn graph\n";
  sout << "  --sps-seqs FILE             \n";
  sout << "  --contigs-mp FILE           \n";
  sout << "  --contigs-midx FILE         \n";
  sout << "  --homglue-ref-midx FILE     \n";
  sout << "  --stars FILE                Star spectra file\n";
  sout << "  --consensus-spectra FILE    Consensus spectra file\n";
  sout << "  --input-spectra-list FILE   Input spectra list file\n";
  sout << "  --cluster-ms FILE           \n";
  sout << "  --cluster-scan FILE         \n";
  sout << "  --cluster-file FILE         Cluster file name (either from MsCluster or PrmClust\n";
  sout << "  --proteins FILE             Fasta file containg protein information\n";
  sout << '\n';
  sout << "  --tool NUMBER               Indicates tool(s) used. 1 for SPS, 2 for GenoMS, 3 for both\n";
  sout << "  --no-clusters               Generate reports without the cluster layer\n";
  sout << '\n';
  sout << "  --server URL                Path to server cgi-bin directory\n";
  sout << "  --user STRING               User name\n";
  sout << "  --job STRING                job name\n";
  sout << "  --title STRING              Specifies title\n";
  sout << "  --cells-per-line NUMBER     Specifies number of aa cells per line in protein coverage view\n";
  sout << "  --no-msms-images            Specifies that no MS/MS images are generated\n";
  sout << '\n';
  sout << "  --outdir PATH               Output directory (defaults to <--project-dir>/report)\n";
  sout << '\n';
  sout << "  --verbose                   Display progress information\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::version(ostream & sout)
{
  sout << "spsReports 3.0." << XSTR(SPS_VERSION) << '\n';
  sout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
  sout << "SH1: " <<  XSTR(GIT_SH1) << endl;
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::error(const string & a)
{
  cerr << "spsReports: invalid or missing option: " << a << endl;

  cerr << "Type 'spsReports --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
