//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "denovo.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "ClusterSet.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define USING_CANCER 0

#define DEBUG_ 0


using namespace std;
using namespace sps;
using namespace specnets;


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 5) {
    cerr << "Usage: main_filter_psms_by_network <psm_file> <contig_info> <star_spectra> <out_psms>" << endl;
    return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("outputdir", "RESULTS_OUTPUT_DIR", 1));

  CommandLineParser clp(argc, argv, 4, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_filter_psms_by_network <psm_file> <contig_info> <out_psms>" << endl;
    cerr << "  " << parserError << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  string specfdrfile(argv[1]);
  string compfile(argv[2]);
  string starspecfile(argv[3]);

  DEBUG_MSG("Loading Spectra PSM file [" << specfdrfile << "]");
  PeptideSpectrumMatchSet	psmSet;
  if (!psmSet.loadFromFile(specfdrfile.c_str())) {
    ERROR_MSG("Loading Spectra PSM file [" << specfdrfile << "]");
    return -6;
  }
  DEBUG_VAR(psmSet.size());

  DEBUG_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
  abinfo_t contigAbinfo;
  if (!Load_abinfo(compfile.c_str(), contigAbinfo)) {
    ERROR_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
    return -9;
  }
  DEBUG_VAR(contigAbinfo.size());

  DEBUG_MSG("Loading Star Spectra file [" << starspecfile << "]");
  SpecSet starSpectra;
  if (starSpectra.loadPklBin(starspecfile.c_str()) <= 0) {
    ERROR_MSG("Loading Star Spectra file [" << starspecfile << "]");
    return -11;
  }
  //DEBUG_VAR(starSpectra.size());


  DEBUG_MSG("Iterating through all contigs/spectra...");


  //----------------------------------------------------------
  // Iterate through the abinfo and get all contigs and spectra
  //----------------------------------------------------------
  map<int, int> specToContig;  
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr = contigAbinfo.begin();
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr_end = contigAbinfo.end();
  for ( ; itr != itr_end; itr++) {
    int contigIndex = itr->first + 1;  // We prefer 1-based contigs

    //DEBUG_MSG("C = " << contigIndex);
    vector<int> specs = itr->second.first.first;
//    totalSpecsInContigs += specs.size();
    for (int i = 0; i < specs.size(); i++) {

      if (specs[i] >= starSpectra.size()) {
        WARN_MSG(specs[i] << " out of bounds!");
        continue;
      }

      int specIndex;
#if 0
      specIndex = specs[i] + 1;  // We prefer 1-based spectra
#else
      int verIndex = specs[i];
      //DEBUG_MSG("  V = " << verIndex);
      specIndex = starSpectra[verIndex].scan;
#endif
      //DEBUG_MSG("  S = " << specIndex);
      if (specToContig.find(specIndex) != specToContig.end()) {
        //WARN_MSG(specIndex << " was found twice!"<< "  " << contigIndex << " and " << specToContig[specIndex ]);
      } else {
        specToContig[specIndex] = contigIndex;
        //contigToSpecSet[contigIndex].insert(specIndex);
        //totalSpecsInContigs++;
        //contigResults[contigIndex].numSpectra++;
      }
    }
  }
  DEBUG_VAR(specToContig.size());

  // Filter out PSMs that are not in contigs  
  PeptideSpectrumMatchSet	psmSetOutput;
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    int scan = psmSet[iPsm]->m_scanNum;
    if (specToContig.find(scan) != specToContig.end()) {
      psmSetOutput.push_back(psmSet[iPsm]);
    }
  }
  DEBUG_VAR(psmSetOutput.size());

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[4];
        
  DEBUG_VAR(psmSetOutput.size());
  DEBUG_MSG("Saving file [" << outputFile << "]");
  if (!psmSetOutput.saveToFile(outputFile.c_str(), true, false)) {
    ERROR_MSG("Error saving to PSM file [ " << outputFile << "]");
    return -1;
  }
  
  
  DEBUG_MSG("Done.");
  return 0;
}


