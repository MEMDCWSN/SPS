///////////////////////////////////////////////////////////////////////////////
#include "ClusterInfoInterface.h"

#include "CommandLineParser.h"
#include "ParameterList.h"
#include "utils.h"
#include "Logger.h"
#include "OutputTable.h"

#include "copyright.h"

#include "main_specnets_defs.h"

#define DEFAULT_BIN_FILES_FILENAME    "spectra/bin_files.txt"
#define DEFAULT_INPUT_FILES_LIST      "spectra/pklbin_files.txt"
#define DEFAULT_ORIGINAL_FILES_LIST   "spectra/input_index.txt"
#define DEFAULT_SUMMARY_FILEMAME      "cluster_summary.tsv"
#define DEFAULT_AUTFILE_NAME          "ClusterInfo.tsv"

#define DEFAULT_SPECSMS_FILENAME      "spectra/specs_ms.pklbin"

#define CSV_SEP '\t'
#define CSV_SEP2 '\t'

///////////////////////////////////////////////////////////////////////////////
namespace specnets
{
///////////////////////////////////////////////////////////////////////////////
  ClusterInfoInterface::ClusterInfoInterface() :
      outdir("."), outFileName(DEFAULT_AUTFILE_NAME),
          binFilesNames(DEFAULT_BIN_FILES_FILENAME),
          inputFilesNames(DEFAULT_INPUT_FILES_LIST),
          originalFilesNames(DEFAULT_ORIGINAL_FILES_LIST),
          outSummaryFilename(DEFAULT_SUMMARY_FILEMAME), projectdir(".")
  {
  }

  ClusterInfoInterface::~ClusterInfoInterface()
  {
  }

///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::loadStringVector(string &inputFilesNames,
                                             vector<std::string> &inputFiles)
  {
    std::vector<std::string> file;
    std::string line;
    file.clear();
    std::ifstream infile(inputFilesNames.c_str(),
                         std::ios_base::in | std::ios_base::binary);

    if (!infile)
      return 0;

    while (getline(infile, line, '\n'))
    {
      inputFiles.push_back(line);
    }

    cout << "for: " << inputFilesNames << " got:" << endl;
    for (int i = 0; i < inputFiles.size(); i++)
      cout << inputFiles[i] << endl;

    return 1;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::writeCsvFile(ParameterList & commandLineParams,
                                         vector<string> *inputFiles,
                                         vector<vector<vector<int> > > *scanNumbers,
                                         string &outFileName,
                                         vector<SpecSet> &specs)
  {
    // data holder
    OutputTable data;
    // get the data
    //bool fileLoaded = clusterData.load(projectdir.c_str());
    //cout << "# of loaded clusters: " << clusterData.data.size() << endl;
    //bool fileLoaded = clusterData.loadMsClusterData(projectdir.c_str());

    string aux = projectdir;
    if (aux.length())
      if (aux[aux.length() - 1] != '/')
        aux += '/';

    string scoredClusterFile;
    if (commandLineParams.exists("SCORED_CLUSTER_FILE")) {
      scoredClusterFile =
          commandLineParams.getValue("SCORED_CLUSTER_FILE").c_str();
    } else {
      scoredClusterFile = aux;
      scoredClusterFile += SPECTRA_DIR + "/" + SPECS_SCORED_CLUST_FILE;
    }
DEBUG_VAR(scoredClusterFile);

    bool fileLoaded = clusterData.loadBinaryFile(scoredClusterFile.c_str());

    if (!fileLoaded)
    {
      string msClusterFile;
      if (commandLineParams.exists("MS_CLUSTER_FILE")) {
        msClusterFile =
            commandLineParams.getValue("MS_CLUSTER_FILE").c_str();
      } else {
        msClusterFile = aux;
        msClusterFile += SPECTRA_DIR + "/" + SPECS_MS_CLUST_FILE;
      }
DEBUG_VAR(msClusterFile);

      fileLoaded = clusterData.loadBinaryFile(msClusterFile.c_str());
    }

    if (!fileLoaded)
    {
      cout << "Unable to load cluster file." << endl;
      return 0;
    }

    cout << "# of loaded clusters: " << clusterData.getClusterCount() << endl;
    clusterData.getCsvData(inputFiles, scanNumbers, specs, data);


    // if error, say so
    if (!data.printToCSV(outFileName.c_str(), "\t"))
    {
      cerr << "ERROR: ClusterInfo: could not open file for writing: "
          << outFileName << endl;
      return -1;
    }

    // exit OK
    return 1;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::WriteClusterSummary(std::map<string, string> &groups,
                                                string groupName,
                                                vector<string> *inputFiles)
  {
    // 1) cluster index, 1 to numClusters
    // 2) number of spectra in cluster
    // 3) cluster Spectrum::parentMass
    // 4) cluster Spectrum::precursorCharge
    // 5) cluster (Spectrum::parentMass + Spectrum::precursorCharge -1) / Spectrum::precursorCharge
    // 6) summed Spectrum::precursor_intensity for all spectra in the cluster

    OutputTable data;
    // add headers
    vector<pair<string, bool> > row;
    row.push_back(pair<string, bool>("cluster index", true));
    row.push_back(pair<string, bool>("number of spectra", true));
    row.push_back(pair<string, bool>("parent mass", true));
    row.push_back(pair<string, bool>("precursor charge", true));
    row.push_back(pair<string, bool>("precursor mass", true));
    row.push_back(pair<string, bool>("sum(precursor intensity)", true));
    // add group names
    int groupNameSize = groupName.size();
    for (map<std::string, std::string>::iterator it = groups.begin();
        it != groups.end(); it++){
      cout<<"Getting summary for cluster: "<<it->first<<"\t"<<it->second<<endl;
      row.push_back(pair<string, bool>(it->first.substr(groupNameSize), true));
    }
    // store headers
    data.values.push_back(row);

    // get summary data
    clusterData.getSummary(data, inputData, specs_ms, inputFiles, groups);


    // if error, say so
    if (!data.printToCSV(outSummaryFilename.c_str(), "\t"))
    {
      cerr << "ERROR: ClusterInfo: could not open summary file for writing: "
          << outSummaryFilename << endl;
      return -1;
    }

    // exit OK
    return 1;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::processOptions(int argc, char **argv)
  {
    //--------------------------------------------------------------------------------------------
    // Initialize directories used
    //--------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------
    // Parse the command line parameters
    //--------------------------------------------------------------------------------------------
    vector<CommandLineParser::Option> listOptions;

    //  c++ instruction                            cmd line option    parameter name    takes value?

    listOptions.push_back(CommandLineParser::Option("-help", "help", false));
    listOptions.push_back(CommandLineParser::Option("-version",
                                                    "VERSION",
                                                    false));

    listOptions.push_back(CommandLineParser::Option("-project-dir",
                                                    "PROJECT_DIR",
                                                    true));
    listOptions.push_back(CommandLineParser::Option("-inputspectra",
                                                    "INPUT_SPECTRA",
                                                    true));
    listOptions.push_back(CommandLineParser::Option("-originalspectra",
                                                    "ORIGINAL_SPECTRA",
                                                    true));

    listOptions.push_back(CommandLineParser::Option("-inputbin",
                                                    "INPUT_BIN",
                                                    true));

    listOptions.push_back(CommandLineParser::Option("-cluster-spectra",
                                                    "CLUSTER_SPECTRA",
                                                    true));

    listOptions.push_back(CommandLineParser::Option("-scored-cluster",
                                                    "SCORED_CLUSTER_FILE",
                                                    true));

    listOptions.push_back(CommandLineParser::Option("-ms-cluster",
                                                    "MS_CLUSTER_FILE",
                                                    true));

    listOptions.push_back(CommandLineParser::Option("-outdir", "OUTDIR", true));
    listOptions.push_back(CommandLineParser::Option("-outfile",
                                                    "OUTFILE",
                                                    true));
    listOptions.push_back(CommandLineParser::Option("-out-summary-file",
                                                    "OUT_SUMMARY_FILE",
                                                    true));

    // parameter file
    listOptions.push_back(CommandLineParser::Option("-p",

                                                    "PARAMETER_FILE",
                                                    true));

    ////////////////////////////////////////////////////////////////////////////////
    // Execute the command line parser
    CommandLineParser clp(argc, argv, 0, listOptions);

    ////////////////////////////////////////////////////////////////////////////////
    // Checking for errors
    string parser_error;
    if (!clp.validate(parser_error))
    {
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

    if (commandLineParams.exists("PARAMETER_FILE"))
    {

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

    // A file must be loaded
    bool fileLoaded = false;

    // input spectra file to load
    if (commandLineParams.exists("INPUT_SPECTRA"))
    {
      inputFilesNames = commandLineParams.getValue("INPUT_SPECTRA").c_str();
DEBUG_VAR(inputFilesNames);
    }

    // input spectra file to load
    if (commandLineParams.exists("ORIGINAL_SPECTRA"))
    {
      originalFilesNames =
          commandLineParams.getValue("ORIGINAL_SPECTRA").c_str();
DEBUG_VAR(originalFilesNames);
    }

    // input spectra file to load
    if (commandLineParams.exists("INPUT_BIN"))
    {
      binFilesNames = commandLineParams.getValue("INPUT_BIN").c_str();
DEBUG_VAR(binFilesNames);
    }

    // cluster binary file to load
    if (commandLineParams.exists("PROJECT_DIR"))
    {
      projectdir = commandLineParams.getValue("PROJECT_DIR").c_str();
DEBUG_VAR(projectdir);
    }

    // load input files names
    string aux = projectdir;
    if (aux[aux.length() - 1] != '/')
      aux += '/';
    aux += inputFilesNames;
DEBUG_VAR(aux);

    if (!loadStringVector(aux, inputFiles))
    {
      cerr << "Unable to load file: " << aux << endl;
      exit(0);
    }

    // load original files names
    aux = projectdir;
    if (aux[aux.length() - 1] != '/')
      aux += '/';
    aux += originalFilesNames;
DEBUG_VAR(aux);

    if (!loadStringVector(aux, originalFiles))
    {
      cerr << "Unable to load file: " << aux << endl;
      exit(0);
    }

    // read input data into memory structure
    for (int i = 0; i < inputFiles.size(); i++)
    {
      SpecSet spec;
      //cout << inputFiles[i] << " loaded" << endl;
      aux = projectdir;
      if (aux[aux.length() - 1] != '/')
        aux += '/';
      aux += inputFiles[i];
DEBUG_VAR(aux);
      if (spec.loadPklBin(aux.c_str()) <= 0)
      {
        cerr << "Unable to load file: " << aux << endl;
        //return 0;
      }
      
      cout<<"Loaded: " <<inputFiles[i]<<endl;
      for(int j = 0; j < spec.size(); j++){
          spec[j].resize(0);
      }
      
      inputData.push_back(spec);
    }

    // states if scan numbers should be read from pklbin (default)
    bool readFromPklbin = true;

    if (commandLineParams.exists("INPUT_BIN"))
    {

      aux = projectdir;
      if (aux[aux.length() - 1] != '/')
        aux += '/';
      aux += binFilesNames;
DEBUG_VAR(aux);

      //cout << "chekcing for " << aux << endl;
      if (loadStringVector(aux, scanFileNames))
      {
        //cout << "reading scan numbers from .bin" << endl;
        for (int i = 0; i < scanFileNames.size(); i++)
        {
          vector<vector<int> > aux;
          Load_binArray(scanFileNames[i].c_str(), aux);
          //cout << aux.size() << " ; " << aux[0].size() << " ; " << aux[0][0] << " ; " << aux[0][1] << endl;
          scanNumbers.push_back(aux);
        }
        readFromPklbin = false;
      }
    }

    if (readFromPklbin)
    {
      //cout << "reading scan numbers from .pklbin" << endl;
      for (int i = 0; i < inputData.size(); i++)
      {
        vector<vector<int> > aux;
        //cout << "Has " << spec.size() << " spectra" << endl;
        for (int j = 0; j < inputData[i].size(); j++)
        {
          vector<int> aux2;
          aux2.push_back(inputData[i][j].scan);
          aux2.push_back(0);
          aux.push_back(aux2);
        }
        scanNumbers.push_back(aux);
      }
    }

    // load specs_ms
    string specs_msFilename = DEFAULT_SPECSMS_FILENAME;
    if (commandLineParams.exists("CLUSTER_SPECTRA"))
      string specs_msFilename = commandLineParams.getValue("CLUSTER_SPECTRA");

    aux = projectdir;
    if (aux[aux.length() - 1] != '/')
      aux += '/';
    aux += specs_msFilename;
DEBUG_VAR(aux);
    if (specs_ms.loadPklBin(aux.c_str()) <= 0)
      cerr << "Unable to load file: " << aux << endl;
    
    cout<<"Clearing Specs_ms peaks"<<endl;
    for(int j = 0; j < specs_ms.size(); j++){
          specs_ms[j].resize(0);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //  File save section

    std::map<string, string> groups;
    std::string groupName("GROUP_");
    commandLineParams.getGroups(groups, groupName);

    ////////////////////////////////////////////////////////////////////////////////
    //  File save section

    // Output directory
    if (commandLineParams.exists("OUTDIR"))
      outdir = commandLineParams.getValue("OUTDIR");
DEBUG_VAR(outdir);

    // Output file name.
    if (commandLineParams.exists("OUTFILE"))
      outFileName = commandLineParams.getValue("OUTFILE");
DEBUG_VAR(outFileName);

    if (commandLineParams.exists("OUT_SUMMARY_FILE"))
      outSummaryFilename = commandLineParams.getValue("OUT_SUMMARY_FILE");
DEBUG_VAR(outSummaryFilename);

    // compose output filename
    aux = outdir;
    if (aux[aux.length() - 1] != '/')
      aux += '/';
    aux += outFileName;
DEBUG_VAR(aux);

    // write the clusterInfo csv file
    writeCsvFile(commandLineParams, &originalFiles, &scanNumbers, aux, inputData);
    
    cout<<"Done Writing Cluster Info"<<endl;

    // write the summary file
    WriteClusterSummary(groups, groupName, &originalFiles);

    // return status ok
    return 0;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::help(ostream & sout)
  {
    version(sout);

    sout << "Usage: clusterinfo [OPTION]\n";
    sout << "Options:\n";
    sout << "  --project-dir DIRECTORY     Project directory\n";
    sout
         << "  --cluster-spectra FILE      Cluster spectra file (defaults to 'spectra/specs_ms.pklbin')\n";
    sout << "  --clusterbin FILE           Clustering binary data file\n";
    sout
         << "  --inputspectra FILE         Input spectra file list (txt format)\n";
    sout
         << "  --scored-cluster FILE       Scored Cluster spectra file (defaults to 'spectra/specs_scored.clust')\n";
    sout
         << "  --ms-cluster FILE           MS Cluster spectra file (defaults to 'spectra/specs_ms.clust')\n";
    sout
         << "  --originalspectra FILE      Original spectra file list (txt format)\n";
    sout << '\n';
    sout << "  --outdir PATH               Output directory (defaults to .)\n";
    sout
         << "  --outfile FILE              Output filename (defaults to 'ClusterInfo.tsv')\n";
    sout
         << "  --out-summary-file FILE     Output filename for the summary file (defaults to 'cluster_summary.tsv')\n";
    sout << '\n';
    sout << "  --p FILE                    Read parameters from file\n";
    sout << '\n';
    sout << "  --help                      Display this help and exit\n";
    sout
        << "  --version                    Output version information and exit\n";
    sout << endl;

    return 1;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::version(ostream & sout)
  {
    sout << "clusterinfo 1.0." << XSTR(SPS_VERSION) << '\n';
    sout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
    sout << "SH1: " << XSTR(GIT_SH1) << endl;
    sout << COPYRIGHT1 << '\n';
    sout << COPYRIGHT2;
    sout << "\n\n";

    return 1;
  }
///////////////////////////////////////////////////////////////////////////////
  int ClusterInfoInterface::error(const string & a)
  {
    cerr << "clusterinfo: invalid or missing option: " << a << endl;

    cerr << "Type 'clusterinfo --help' for more information." << endl << endl;

    return -1;
  }
///////////////////////////////////////////////////////////////////////////////
}
;
// namespace specnets
///////////////////////////////////////////////////////////////////////////////
