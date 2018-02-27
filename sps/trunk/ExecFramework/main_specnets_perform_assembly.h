#include "abruijn.h"

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

  bool performAssembly(ParameterList & ip,
                       SpecSet & starSpectra,
                       SpectrumPairSet & starPairs,
                       Clusters & contigShifts,
                       abinfo_t & contigAbinfo);

  bool performContigAlignment(ParameterList & ip,
                              Clusters & contigs,
                              SpectrumPairSet & contigPairs);

  bool performMetaAssembly(ParameterList & ip,
                           Clusters & contigs,
                           SpectrumPairSet & contigAligns,
                           abinfo_t & contigAbruijn,
                           SpecSet & starSpectra,
                           Clusters & outputMetaContigs,
                           abinfo_t & metaContigAbruijn);
}

