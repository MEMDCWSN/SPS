#include "spectrum.h"

#include <list>
#include <string>
#include <vector>

namespace specnets
{
  class AAJumps;
  class ParameterList;
  class PeptideSpectrumMatchSet;
  class SpectrumPairSet;
  class SpecSet;

  bool performFilterPairs(ParameterList & ip,
                          SpecSet & inputSpectra,
                          SpecSet & inputSpectraMS2,
                          SpectrumPairSet & filteredPairs,
                          vector<TwoValues<float> > & ratios,
                          vector<TwoValues<float> > & means,
                          vector<float> & varTerms,
                          list<vector<float> > & alignStats,
                          vector<vector<float> > & specStats,
                          std::vector<unsigned int> & idxKept,
                          std::vector<TwoValues<float> > & pvalues,
                          bool gridExecutionFlag,
                          bool resume);

  bool performFilterAligns(ParameterList & ip,
                           SpectrumPairSet & filteredPairs,
                           SpecSet & prmSpectra,
                           PeptideSpectrumMatchSet & filterPsmSet,
                           vector<TwoValues<float> > & ratios,
                           vector<TwoValues<float> > & means,
                           vector<float> & varTerms,
                           std::vector<unsigned int> & idxKept,
                           std::vector<TwoValues<float> > & pvalues);

  bool performAlignment(ParameterList & ip,
                        AAJumps & jumps,
                        SpecSet &inputSpectra,
                        SpectrumPairSet &inputPairs,
                        SpecSet &pairAlignments,
                        SpecSet &starSpectraOnly,
                        SpecSet &starSpectra,
                        vector<unsigned int> &alignedSpectra);

  bool performFilterStarPairs(ParameterList & ip,
                              SpectrumPairSet & inputPairs,
                              SpecSet & starSpectra,
                              vector<vector<float> > & ratios,
                              SpecSet & matchedPeaks);
}

