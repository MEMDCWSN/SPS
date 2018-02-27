///////////////////////////////////////////////////////////////////////////////
#include "PairsInfoInterface.h"


#include "CommandLineParser.h"
#include "ParameterList.h"
#include "utils.h"
#include "Logger.h"
#include "SpecnetsGraph.h"
#include <set>

#include "copyright.h"

#define DEFAULT_INPUT_FILE            "aligns/pairs.bin"
#define DEFAULT_OUTFILE_NAME          "pairs.txt"

#define CSV_SEP '\t'


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
PairsInfoInterface::PairsInfoInterface() :
  outdir("."),
  outFileName(DEFAULT_OUTFILE_NAME),
  inputFilename(DEFAULT_INPUT_FILE),
  projectdir("."),
  m_minScore1(-1.0),
  m_minScore2(-1.0),
  m_edgeTopKBoth(-1.0),
  m_edgeTopKOne(-1.0),
  m_exists_minScore1(false),
  m_exists_minScore2(false),
  m_exists_edgeTopKBoth(false),
  m_exists_edgeTopKOne(false),
  m_exists_component_size(false)
{
}

PairsInfoInterface::~PairsInfoInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
int PairsInfoInterface::writeOutFile(void)
{
  // 1) index of the first spectrum
  // 2) index of the second spectrum
  // 3) shift 1
  // 4) shift 2
  // 5) score 1
  // 6) score 2

  vector<vector<string> > data;

  vector<vector< float > > ranks;

  vector<float> rank_cutoffs;

  float mzdiff_cutoff = 0.0;

  //cout<<"EDGES " <<pairSet.size()<<endl;

  int max_index = 0;
  for(int i = 0; i < pairSet.size(); i++){
		int index1 = pairSet[i].spec1;
		int index2 = pairSet[i].spec2;
		max_index = max(max(index1, max_index), index2);
  }

  ranks.resize(max_index+1);
  rank_cutoffs.resize(max_index+1);
  for(int i = 0; i < pairSet.size(); i++){
  	if(abs(pairSet[i].shift1) < mzdiff_cutoff) continue;
  	ranks[pairSet[i].spec1].push_back(pairSet[i].score1);
  	ranks[pairSet[i].spec2].push_back(pairSet[i].score1);
  }


  for(int i = 0; i < ranks.size(); i++){
  	//sort each guy by score1
  	sort(ranks[i].begin(), ranks[i].end());
  }


  int k_top = 0;
  if(m_exists_edgeTopKBoth)
  	k_top = m_edgeTopKBoth;
  if(m_exists_edgeTopKOne)
  	k_top = m_edgeTopKOne;
  for(int i = 0; i < ranks.size(); i++){
  	int cutoff_idx = min(k_top - 1, (int)ranks[i].size() - 1);
  	float cutoff = 1.f;
  	if(cutoff_idx < 0)
  		cutoff = 1.f;
  	else
  		cutoff = ranks[i][ranks[i].size() - cutoff_idx - 1];
  	rank_cutoffs[i] = cutoff;
  }
  
  SpectrumPairSet filtered_pairset;
  for(int i = 0 ; i < pairSet.size() ; i++) {

    if(m_exists_minScore1)
      if(m_minScore1 > pairSet[i].score1)
        continue;

    if(m_exists_minScore2)
      if(m_minScore2 > pairSet[i].score2)
        continue;

    //EDGE_TOPK_BOTH=K --> only outputs edges where the rank of the edge (by score1) is <=K for both adjacent nodes
    if(m_exists_edgeTopKBoth)
      if(rank_cutoffs[pairSet[i].spec1] > pairSet[i].score1 || rank_cutoffs[pairSet[i].spec2] > pairSet[i].score1)
        continue;

    //EDGE_TOPK_ONE=K --> only outputs edges where the rank of the edge (by score1) is <=K for at least one adjacent node
    if(m_exists_edgeTopKOne)
      if(rank_cutoffs[pairSet[i].spec1] > pairSet[i].score1 && rank_cutoffs[pairSet[i].spec2] > pairSet[i].score1)
        continue;


    //cout<<pairSet[i].shift1 <<" "<< mzdiff_cutoff<<endl;
    if(abs(pairSet[i].shift1) < mzdiff_cutoff) continue;
    
    filtered_pairset.push_back(pairSet[i]);
    
  }
  
  
  //Filter minimium component sizes
  if(m_exists_component_size){
        //Don't do anything if component sizes are really small
        if(m_component_size > 0){
            //Creating graph and filtering
            SpecnetsGraph spec_graph(filtered_pairset);
            
            spec_graph.filter_graph_component_size(m_component_size);
            
            pairSet.resize(0);
            
            std::vector<unsigned int> deleted_edges = spec_graph.get_pairs_deleted();
            
            std::set<unsigned int> deleted_edgs_set;
            
            for(int i = 0; i < deleted_edges.size(); i++){
                deleted_edgs_set.insert(deleted_edges[i]);
            }
            
            for(int i = 0; i < filtered_pairset.size(); i++){
                if(deleted_edgs_set.find(i) == deleted_edgs_set.end()){
                    pairSet.push_back(filtered_pairset[i]);
                }
            }
            
            filtered_pairset.resize(0);
        }
        else{
            cout<<"No Component Filtering"<<endl;
        }
  }
  else{
      pairSet.resize(0);
      pairSet = filtered_pairset;
  }
  
  map<unsigned int, unsigned int> node_to_connected_component;
  SpecnetsGraph spec_graph(pairSet);
  spec_graph.get_node_to_connected_component_mapping(node_to_connected_component);
  
  for(int i = 0 ; i < pairSet.size() ; i++) {

    vector<string> row;
    stringstream aux;

    // spectrum 1 index
    int aaa = pairSet[i].spec1 + 1;
    aux << aaa;
    row.push_back(aux.str());

    // spectrum 2 index
    aux.str(std::string());
    aaa = pairSet[i].spec2 + 1;
    aux << aaa;
    row.push_back(aux.str());

    // shift 1 value
    aux.str(std::string());
    aux << pairSet[i].shift1;
    row.push_back(aux.str());

    // shift 2 value
    aux.str(std::string());
    aux << pairSet[i].shift2;
    row.push_back(aux.str());

    aux.str(std::string());
    aux << pairSet[i].score1;
    row.push_back(aux.str());

    aux.str(std::string());
    aux << pairSet[i].score2;
    row.push_back(aux.str());
    
    
    // connected component
    aux.str(std::string());
    aux << node_to_connected_component[pairSet[i].spec1];
    row.push_back(aux.str());

    data.push_back(row);
  }

  // open file to write to
  ofstream file(outFileName.c_str(), ios::out | ios::binary);
    // if error, say so
  if(!file.is_open()) {
    cerr << "ERROR: PairsInfo: could not open file for writing: " << outFileName << endl;
    return -1;
  }

  // output file

  // Write records
  for(int i = 0 ; i < data.size() ; i++ ) {
    // Get record
    vector<string> & row = data[i];
    // cycle thru the row
    for(int j = 0 ; j < row.size() ; j++) {
      // write the separator
      if(j)
        file << CSV_SEP;
      // write the element
      file << row[j];
    }
    // write the end of line
    file << endl;
  }

  // close the file
  file.close();

  // exit OK
  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int PairsInfoInterface::processOptions(int argc, char **argv)
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


  listOptions.push_back(CommandLineParser::Option("-project-dir",     "PROJECT_DIR",      true));
  listOptions.push_back(CommandLineParser::Option("-input-file",      "INPUT_FILE",       true));

  listOptions.push_back(CommandLineParser::Option("-outdir",          "OUTDIR",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",         "OUTFILE",          true));

  listOptions.push_back(CommandLineParser::Option("-min-score1",      "MIN_SCORE1",       true));
  listOptions.push_back(CommandLineParser::Option("-min-score2",      "MIN_SCORE2",       true));

  listOptions.push_back(CommandLineParser::Option("-edge-topk-both",  "EDGE_TOPK_BOTH",   true));
  listOptions.push_back(CommandLineParser::Option("-edge-topk-one",   "EDGE_TOPK_ONE",    true));
  
  // filter component sizes
  listOptions.push_back(CommandLineParser::Option("-max-component-size",   "MAXCOMPONENTSIZE",    true));

  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",               "PARAMETER_FILE",   true));


  ////////////////////////////////////////////////////////////////////////////////
  // Execute the command line parser
  CommandLineParser clp(argc, argv, 0, listOptions);

  ////////////////////////////////////////////////////////////////////////////////
  // Checking for errors
  string parser_error;
  if (!clp.validate(parser_error)) {
    cerr << parser_error << endl;
    //ERROR_MSG(parser_error);
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
  if (commandLineParams.exists("INPUT_FILE")) {
    inputFilename = commandLineParams.getValue("INPUT_FILE").c_str();
  }

  // cluster binary file to load
  if (commandLineParams.exists("PROJECT_DIR")) {
    projectdir = commandLineParams.getValue("PROJECT_DIR").c_str();
  }

  // load the paris file
  if(pairSet.loadFromBinaryFile(inputFilename)  < 0) {
    cerr << "Unable to load file: " << inputFilename << endl;
    return 0;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //  Data parameters section

  // min score 1
  if (commandLineParams.exists("MIN_SCORE1")) {
    m_minScore1 = commandLineParams.getValueFloat("MIN_SCORE1");
    m_exists_minScore1 = true;
  }

  // min score 2
  if (commandLineParams.exists("MIN_SCORE2")) {
    m_minScore2 = commandLineParams.getValueFloat("MIN_SCORE2");
    m_exists_minScore2 = true;
  }

  // EDGE_TOPK_BOTH
  if (commandLineParams.exists("EDGE_TOPK_BOTH")) {
    m_edgeTopKBoth = commandLineParams.getValueFloat("EDGE_TOPK_BOTH");
    m_exists_edgeTopKBoth = true;
  }

  // EDGE_TOPK_ONE
  if (commandLineParams.exists("EDGE_TOPK_ONE")) {
    m_edgeTopKOne = commandLineParams.getValueFloat("EDGE_TOPK_ONE");
    m_exists_edgeTopKOne = true;
  }
  
  // EDGE_TOPK_ONE
  if (commandLineParams.exists("MAXCOMPONENTSIZE")) {
    m_component_size = commandLineParams.getValueInt("MAXCOMPONENTSIZE");
    m_exists_component_size = true;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("OUTDIR"))
    outdir = commandLineParams.getValue("OUTDIR");

  // Output file name.
  if (commandLineParams.exists("OUTFILE"))
    outFileName = commandLineParams.getValue("OUTFILE");

  // compose output filename
  string aux = outdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += outFileName;

  // write the clusterInfo csv file
  writeOutFile();

  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int PairsInfoInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: pairsinfo [OPTION]\n";
  sout << "Options:\n";
  sout << "  --project-dir DIRECTORY     Project directory\n";
  sout << "  --input-file FILE           Input pairs file\n";
  sout << '\n';
  sout << "  --outdir PATH               Output directory (defaults to .)\n";
  sout << "  --outfile FILE              Output filename (defaults to 'pairs.txt')\n";
  sout << '\n';
  sout << "  --min-score1 VALUE          Only output pairs with score >= VALUE for the 1st node\n";
  sout << "  --min-score2 VALUE          Only output pairs with score >= VALUE for the 2nd node\n";
  sout << "  --edge-topk-both VALUE      Only outputs edges where the rank of the edge (by score) is <= VALUE for both adjacent nodes\n";
  sout << "  --edge-topk-one VALUE       Only outputs edges where the rank of the edge (by score) is <= VALUE for at least one adjacent node\n";
  sout << '\n';
  sout << "  --max-component-size VALUE  Maximum Connected Component Size in Network. Filters edges in component to meet criterion.\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int PairsInfoInterface::version(ostream & sout)
{
  sout << "pairsinfo 1.0." << XSTR(SPS_VERSION) << '\n';
  sout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
  sout << "SH1: " <<  XSTR(GIT_SH1) << endl;
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int PairsInfoInterface::error(const string & a)
{
  cerr << "pairsinfo: invalid or missing option: " << a << endl;
  cerr << "Type 'pairsinfo --help' for more information." << endl << endl;
  return -1;
}


}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
