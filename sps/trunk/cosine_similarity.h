#ifndef COSINE_SIMILARITY
#define COSINE_SIMILARITY

#include <string>
#include <vector>
#include "utils.h"
#include "spectrum.h"
#include "spectrum_scoring.h"
#include "PeptideSpectrumMatch.h"

using namespace std;

namespace specnets
{
  class CosineSimilarity
  {
    vector<string> m_ions_to_extract;

  public:

    CosineSimilarity(vector<string> ions_to_extract)
    {
      this->m_ions_to_extract = ions_to_extract;
    }

    float similarity(const Spectrum &spec1,
                     const Spectrum &spec2,
                     const string &annotation1,
                     const string &annotation2,
                     const MS2ScoringModel &model,
                     int do_sqrt_filter);
  };
}

#endif
