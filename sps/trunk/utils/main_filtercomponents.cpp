//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "SpecnetsGraph.h"
#include "tuple.h"

#include <stdlib.h>

#define DEBUG_ 0
#define MAX_COMPONENT_SIZE 20

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 5) {
    cerr << "Usage: main_filtercomponents in_pairset_file fullpairs_file ratios_file out_pairset_file" << endl;
    return -1;
  }

  SpectrumPairSet fullPairSet;
  SpectrumPairSet inPairSet;
  SpectrumPairSet outPairSet;

  inPairSet.loadFromBinaryFile(argv[1]);
  DEBUG_VAR(inPairSet.size());

  if (!fullPairSet.loadFromBinaryFile(argv[2])) {
    ERROR_MSG("Could not load: " << argv[2]);
    return -2;
  }
  DEBUG_VAR(fullPairSet.size());

  vector<TwoValues<float> > ratios;
  if (!Load_binArray(argv[3], ratios)) {
    ERROR_MSG("Could not load: " << argv[3]);
    return -1;
  }
  DEBUG_VAR(ratios.size());

  // Reset the scores to ratio scores for filtering
  int nIndex = 0;
  for (int i = 0; i < inPairSet.size(); i++) {
    SpectrumPair & thePair = inPairSet[i];
    while(thePair.spec1 != fullPairSet[nIndex].spec1) nIndex++;
    while(thePair.spec2 != fullPairSet[nIndex].spec2) nIndex++;
    thePair.score1 = min(ratios[nIndex][0], ratios[nIndex][1]);
    //DEBUG_MSG(thePair.spec1 << "  " << thePair.spec2 << "  " << nIndex << "  " << thePair.score1);
  }

  SpecnetsGraph spec_graph(inPairSet);
            
  spec_graph.filter_graph_component_size(MAX_COMPONENT_SIZE);
            
  outPairSet.resize(0);
            
  std::vector<unsigned int> deleted_edges = spec_graph.get_pairs_deleted();

  std::set<unsigned int> deleted_edgs_set;
  for(int i = 0; i < deleted_edges.size(); i++){
    deleted_edgs_set.insert(deleted_edges[i]);
  }
            
  for(int i = 0; i < inPairSet.size(); i++){
    if(deleted_edgs_set.find(i) == deleted_edgs_set.end()){
      outPairSet.push_back(inPairSet[i]);
    }
  }
            
  DEBUG_VAR(outPairSet.size());
  outPairSet.saveToBinaryFile(argv[4]);

  return 0;
}
