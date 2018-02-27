#include "abruijn.h"
#include "SpecSet.h"

#include <list>
#include <string>
#include <vector>

namespace specnets
{
  class AAJumps;
  class ParameterList;
  class SpectrumPairSet;
  class SpecSet;
  class Clusters;
  class PeptideSpectrumMatchSet;
  class PenaltyMatrix;
  class DB_fasta;

  bool performMS2Deconv(ParameterList & ip,
                        vector<string>& inputSpecs,
                        vector<string>& outputSpecs);

  bool performMsCluster(ParameterList & ip,
                        vector<string>& inputFilesList,
                        SpecSet &clustSpectra);

  bool performScoring(ParameterList & ip,
                      SpecSet & inputSpectra,
                      SpecSet & outputSpectra,
                      bool gridExecutionFlag,
                      bool resume);

  bool performPrmClustering(ParameterList & ip,
                            SpecSet & inputSpectra,
                            SpecSet & outputSpectra);

  bool performStatProtSeqs(ParameterList & ip);

  bool performProtProtAlign(ParameterList & ip,
                            unsigned int refProtIdx,
                            set<unsigned int> dbIndexes,
                            DB_fasta & db);

  bool performHomologyAssembly(ParameterList & ip,
                               SpecSet & spectra,
                               DB_fasta & db,
                               SpecSet & contigShifts,
                               SpecSet & contigMatchedIndices,
                               vector<vector<int> > & contigMatchedProts);

  bool performReport(ParameterList & ip);

  bool performExecMainSpecnets(ParameterList & ip,
                               SpecSet * msSpectra,
                               SpecSet * scoredSpectra,
                               SpecSet * starSpectra,
                               SpectrumPairSet * pairs,
                               DB_fasta * db,
                               PenaltyMatrix * penaltyMatrixBlosum,
                               PenaltyMatrix * penaltyMatrixMods,
                               PeptideSpectrumMatchSet * psms,
                               PeptideSpectrumMatchSet * origPsms,
                               SpecSet * psms_spectra,
                               SpecSet * psms_midx,
                               vector<vector<int> > * psms_mp,
                               SpecSet * snets_contigs,
                               SpecSet * snets_midx,
                               vector<vector<int> > * snets_mp);
}

