#include "cosine_similarity.h"

#include <string>
#include <vector>


using namespace std;

namespace specnets
{
  // Computes to cosine between two normalized vectors of the same size
  inline float cosine(vector<float> &u, vector<float> &v)
  {
    if (u.size() != v.size())
      return -1000;
    float c = 0;
    for (unsigned int i = 0; i < u.size(); i++)
      c += u[i] * v[i];
    return c;
  }

  inline void normalize_spectrum(Spectrum &s)
  {

    float n = 0.f;
    for (int peakIdx = 0; peakIdx < s.size(); peakIdx++)
    {
      n += s[peakIdx][1] * s[peakIdx][1];
    }
    n = sqrt(n);
    for (int peakIdx = 0; peakIdx < s.size(); peakIdx++)
    {
      s[peakIdx][1] /= n;
    }
  }

  inline void sqrt_filter(Spectrum &s)
  {
    //Experimental Step to Sqrt intesities
    for (int peakIdx = 0; peakIdx < s.size(); peakIdx++)
    {
      s[peakIdx][1] = sqrt(s[peakIdx][1]);
    }

  }

  inline void normalize(vector<float> &v)
  {

    float n = 0;
    for (unsigned int i = 0; i < v.size(); i++)
      n += v[i] * v[i];
    if (n > 0)
    {
      n = sqrt(n);
      for (unsigned int i = 0; i < v.size(); i++)
        v[i] /= n;
    }
  }

  inline void extractIons(const PeptideSpectrumMatch &psm,
                          const int peptideLength,
                          const MS2ScoringModel &model,
                          const vector<string> &ionsToExtract,
                          vector<float> &ions)
  {
    ions.resize(peptideLength * ionsToExtract.size());
    for (unsigned int i = 0; i < ions.size(); i++)
      ions[i] = 0;

    // Select target ions
    unsigned int baseIdx = 0;
    for (unsigned int ionType = 0; ionType < ionsToExtract.size(); ionType++)
    {
      for (unsigned int peakIdx = 0; peakIdx < psm.m_peakAnnotations.size(); peakIdx++)
      {
        if (psm.m_peakAnnotations[peakIdx].first and psm.m_peakAnnotations[peakIdx].first->name.compare(ionsToExtract[ionType]) == 0)
        {
          ions[baseIdx + psm.m_peakAnnotations[peakIdx].second - 1] = (*psm.m_spectrum)[peakIdx][1];
        }
      }
      baseIdx += peptideLength;
    }
  }

  float CosineSimilarity::similarity(const Spectrum &spec1,
                                     const Spectrum &spec2,
                                     const string &annotation1,
                                     const string &annotation2,
                                     const MS2ScoringModel &model,
                                     int do_sqrt_filter)
  {
    vector<float> temp_ion_intensity1;
    vector<float> temp_ion_intensity2;
    vector<float> temp_ion_mass1;
    vector<float> temp_ion_mass2;

    string annotation1_clean = annotation1;
    string annotation2_clean = annotation2;

    string allIons("all");

    PeptideSpectrumMatch psm1;
    PeptideSpectrumMatch psm2;

    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;

    psm1.m_spectrum = &spec1_temp;
    psm2.m_spectrum = &spec2_temp;

    AAJumps jumps(1);
    int psm1Length = jumps.getPeptideLength(annotation1);
    int psm2Length = jumps.getPeptideLength(annotation2);

    psm1.annotate(annotation1, allIons, model, 0, 0, .45);
    psm2.annotate(annotation2, allIons, model, 0, 0, .45);

    normalize_spectrum(spec1_temp);
    normalize_spectrum(spec2_temp);

    if (do_sqrt_filter)
    {
      sqrt_filter(spec1_temp);
      sqrt_filter(spec2_temp);
    }

    extractIons(psm1,
                psm1Length,
                model,
                m_ions_to_extract,
                temp_ion_intensity1);

    extractIons(psm2,
                psm2Length,
                model,
                m_ions_to_extract,
                temp_ion_intensity2);

    //normalize(temp_ion_intensity1);
    //normalize(temp_ion_intensity2);

    float cosine_val = cosine(temp_ion_intensity1, temp_ion_intensity2);

    return cosine_val;
  }
}
