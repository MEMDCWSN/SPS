/**
 * Annotation of ion statistics
 */

/**
 * Calculate explained intensity for input spectrum
 * @param spectrum - annotated Spectrum (@see spectrum)
 * @param ionNames - comma delimited string indicating which ion names to consider (@see spectrum_scoring for
 * information on ftIonFragment names)
 * @param jumps - amino acid masses to use. @see aminoacid.h.
 */
#include "Logger.h"
#include "SpectrumAnnotStatistics.h"

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentExplainedIntensity(const PeptideSpectrumMatch &psm,
                                                           SpectrumAnnotParameter &ionParam)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    AAJumps jumps(1);

    vector<string> ion_types;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    vector<float> masses;

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

    float total_intensity = 0;
    float identified_intensity = 0;

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        total_intensity += (*spectrum)[i][1]; // add to total intensity regardless of annotation
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        { // make sure peak is annotated
          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            identified_intensity += (*spectrum)[i][1];
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              { // include only listed peaks
                identified_intensity += (*spectrum)[i][1];
              }
            }
          }
        }
      }
    }
    if (total_intensity > 0)
    {
      return (identified_intensity / total_intensity) * 100;
    }
    else
    {
      return 0;
    }
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentExplainedIntensity(const PeptideSpectrumMatchSet &psmSet,
                                                           SpectrumAnnotParameter &ionParam)
  {

    float total_intensity = 0;
    float identified_intensity = 0;

    AAJumps jumps(1);

    for (int psmIdx = 0; psmIdx < psmSet.size(); psmIdx++)
    {

      PeptideSpectrumMatch& psm = *(psmSet[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return -1.0;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

      vector<string> ion_types;

      splitText(ionNames.c_str(), ion_types, (const char*)",");

      vector<float> masses;

      double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

      for (unsigned i = 0; i < spectrum->size(); i++)
      {
        if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
            <= (peptide_mass - 37))
        { // ignore small masses or masses larger than we expect.
          total_intensity += (*spectrum)[i][1]; // add to total intensity regardless of annotation
          const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
          if (curr_frag != NULL)
          { // make sure peak is annotated
            if (ionNames.compare("all") == 0)
            { // include all identified peaks
              identified_intensity += (*spectrum)[i][1];
            }
            else
            {
              for (int j = 0; j < ion_types.size(); j++)
              {
                if (ion_types[j].compare(curr_frag->name) == 0)
                { // include only listed peaks
                  identified_intensity += (*spectrum)[i][1];
                }
              }
            }
          }
        }
      }
    }
    if (total_intensity > 0)
    {
      return (identified_intensity / total_intensity) * 100;
    }
    else
    {
      return 0;
    }
  }
  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentExplainedPeaks(const PeptideSpectrumMatch &psm,
                                                       SpectrumAnnotParameter &ionParam)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    AAJumps jumps(1);

    vector<string> ion_types;

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

    float total_peaks = 0;
    float identified_peaks = 0;

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        total_peaks++;// add to total regardless of annotation
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        { // make sure peak is annotated
          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            identified_peaks++;
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              { // include only listed peaks
                identified_peaks++;
              }
            }
          }
        }
      }
    }
    if (total_peaks > 0)
    {
      return ((float)identified_peaks / (float)total_peaks) * 100;
    }
    else
    {
      return 0;
    }
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentExplainedPeaks(const PeptideSpectrumMatchSet &psmSet,
                                                       SpectrumAnnotParameter &ionParam)
  {

    float total_peaks = 0;
    float identified_peaks = 0;
    AAJumps jumps(1);

    for (int psmIdx = 0; psmIdx < psmSet.size(); psmIdx++)
    {
      PeptideSpectrumMatch& psm = *(psmSet[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return -1.0;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

      vector<string> ion_types;

      splitText(ionNames.c_str(), ion_types, (const char*)",");

      double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

      for (unsigned i = 0; i < spectrum->size(); i++)
      {
        if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
            <= (peptide_mass - 37))
        { // ignore small masses or masses larger than we expect.
          total_peaks++;// add to total regardless of annotation
          const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
          if (curr_frag != NULL)
          { // make sure peak is annotated
            if (ionNames.compare("all") == 0)
            { // include all identified peaks
              identified_peaks++;
            }
            else
            {
              for (int j = 0; j < ion_types.size(); j++)
              {
                if (ion_types[j].compare(curr_frag->name) == 0)
                { // include only listed peaks
                  identified_peaks++;
                }
              }
            }
          }
        }
      }
    }
    if (total_peaks > 0)
    {
      return ((float)identified_peaks / (float)total_peaks) * 100;
    }
    else
    {
      return 0;
    }
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentObservedIons(const PeptideSpectrumMatch &psm,
                                                     SpectrumAnnotParameter &ionParam)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    AAJumps jumps(1);

    vector<string> ion_types;

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

    int peptide_length = jumps.getPeptideLength(psm.m_annotation);

    int num_ions; //number of different ion types
    int identified_ions = 0; //number of identified ion peaks

    if (ionNames.compare("all") == 0)
    {
      num_ions = psm.m_ionTypes.size();
    }
    else
    {
      num_ions = ion_types.size();
    }

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        { // make sure peak is annotated
          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            identified_ions++;
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              { // include only listed peaks
                identified_ions++;
              }
            }
          }
        }
      }
    }

    if (num_ions > 0 && peptide_length - 1 > 0)
    {
      return ((float)identified_ions / ((float)num_ions
          * ((float)peptide_length - 1))) * 100;
    }
    else
    {
      return 0;
    }
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::percentObservedIons(const PeptideSpectrumMatchSet &psmSet,
                                                     SpectrumAnnotParameter &ionParam)
  {
    AAJumps jumps(1);
    long peptide_length = 0;
    long num_ions; //number of different ion types
    long identified_ions = 0; //number of identified ion peaks

    for (int psmIdx = 0; psmIdx < psmSet.size(); psmIdx++)
    {
      PeptideSpectrumMatch& psm = *(psmSet[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return -1.0;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

      vector<string> ion_types;

      splitText(ionNames.c_str(), ion_types, (const char*)",");

      double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

      peptide_length += jumps.getPeptideLength(psm.m_annotation);

      if (ionNames.compare("all") == 0)
      {
        num_ions += psm.m_ionTypes.size();
      }
      else
      {
        num_ions += ion_types.size();
      }

      for (unsigned i = 0; i < spectrum->size(); i++)
      {
        if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
            <= (peptide_mass - 37))
        { // ignore small masses or masses larger than we expect.
          const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
          if (curr_frag != NULL)
          { // make sure peak is annotated
            if (ionNames.compare("all") == 0)
            { // include all identified peaks
              identified_ions++;
            }
            else
            {
              for (int j = 0; j < ion_types.size(); j++)
              {
                if (ion_types[j].compare(curr_frag->name) == 0)
                { // include only listed peaks
                  identified_ions++;
                }
              }
            }
          }
        }
      }
    }

    if (num_ions > 0 && peptide_length - 1 > 0)
    {
      return ((float)identified_ions / ((float)num_ions
          * ((float)peptide_length - 1))) * 100;
    }
    else
    {
      return 0;
    }
  }

  unsigned int SpectrumAnnotStatistics::numNonEmptySpectra(const PeptideSpectrumMatchSet &psmSet)
  {
    unsigned int count = 0;
    for (int psmIdx = 0; psmIdx < psmSet.size(); psmIdx++)
    {
      PeptideSpectrumMatch& psm = *(psmSet[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return -1;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() > 5)
      {
        count++;
      }
    }
    return count;
  }

  // -------------------------------------------------------------------------
  int SpectrumAnnotStatistics::totalPeaks(const PeptideSpectrumMatch &psm)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    int total_peaks = 0;

    AAJumps jumps(1);

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        total_peaks++;// add to total regardless of annotation
      }
    }
    return total_peaks;
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::parentMassErrorPPM(const PeptideSpectrumMatch &psm,
                                                    int charge)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    AAJumps jumps(1);

    if (charge <= 0)
    { //charge is not defined
      charge = spectrum->parentCharge;
    }

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);
    double monoisotopic_mass = peptide_mass;
    double oneisotopic_mass = peptide_mass + AAJumps::massHion; /* one C13 value */

    if (charge <= 0)
    {
      WARN_MSG("Warning: Charge not defined! Calculating from parentMass.");
      charge = peptide_mass / spectrum->parentMass;
      if (charge <= 0)
      {
        WARN_MSG("Charge not able to be calculated from ParentMass, using charge 1");
        charge = 1;
      }
    }

    /* set to parent mass/charge */
    monoisotopic_mass = (monoisotopic_mass + ((double)charge
        * AAJumps::massHion) + AAJumps::massH2O) / (double)charge;
    oneisotopic_mass = (oneisotopic_mass + ((double)charge * AAJumps::massHion)
        + AAJumps::massH2O) / (double)charge;

#ifdef DEBUG
    std::cout << "monoisotopic_mass " << monoisotopic_mass << endl;
    std::cout << "oneisotopic_mass " << oneisotopic_mass << endl;
    std::cout << "charge " << charge << endl;
#endif

#ifdef DEBUG
    std::cout << "ParentMZ " << spectrum->parentMZ << endl;
#endif

    float monoisotopic_ppm = 0;
    float oneisotopic_ppm = 0;

    if (spectrum->parentMZ > 0)
    {
      //if parentMZ is correctly set.
      /* get PPM */
      monoisotopic_ppm = 1000000.0 * ((spectrum->parentMZ - monoisotopic_mass)
          / monoisotopic_mass);
      oneisotopic_ppm = 1000000.0 * ((spectrum->parentMZ - oneisotopic_mass)
          / monoisotopic_mass);
    }
    else
    {
      //if parentMass is the only thing set.
      if (spectrum->parentCharge > 0)
      {
        //if charge is set, parentMass is scaled by charge.
        float parentMZ = (spectrum->parentMass + (AAJumps::massHion
            * (spectrum->parentCharge - 1))) / spectrum->parentCharge;
        monoisotopic_ppm = 1000000.0 * ((parentMZ - monoisotopic_mass)
            / monoisotopic_mass);
        oneisotopic_ppm = 1000000.0 * ((parentMZ - oneisotopic_mass)
            / monoisotopic_mass);
      }
      else
      {
        //if charge isn't set, parentMass is equivalent to mz.
        monoisotopic_ppm = 1000000.0 * ((spectrum->parentMass
            - monoisotopic_mass) / monoisotopic_mass);
        oneisotopic_ppm = 1000000.0
            * ((spectrum->parentMass - oneisotopic_mass) / monoisotopic_mass);
      }
    }

#ifdef DEBUG
    std::cout << "monoisotopic_ppm " << monoisotopic_ppm << endl;
    std::cout << "oneisotopic_ppm " << oneisotopic_ppm << endl;
#endif

    if (abs(monoisotopic_ppm) > abs(oneisotopic_ppm))
    {
      return oneisotopic_ppm;
    }
    else
    {
      return monoisotopic_ppm;
    }
  }

  // -------------------------------------------------------------------------
  void SpectrumAnnotStatistics::fragmentMassErrorDaHistogram(const PeptideSpectrumMatch &psm,
                                                             SpectrumAnnotParameter &ionParam,
                                                             float maxError,
                                                             vector<TwoValues<
                                                                 float> > &bins)
  {
    const unsigned int defautNumBins = 60;
    if (bins.size() == 0)
    {
      bins.resize(defautNumBins);
      float binWidth = (maxError * 2.0) / ((float)defautNumBins);
      bins[0][0] = 0.0 - maxError;
      for (unsigned int i = 1; i < bins.size(); i++)
      {
        bins[i][0] = bins[i - 1][0] + binWidth;
      }
    }

    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return;
    }

    Spectrum * spectrum = psm.m_spectrum;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    AAJumps jumps(1);

    vector<string> ion_types;

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    vector<float> masses;

    float peptide_mass;

    //generate srm and prm masses, no offsets.
    vector<float> prm_masses;
    vector<float> srm_masses;

    vector<float> massErrors;

    jumps.getPRMandSRMMasses(psm.m_annotation,
                             prm_masses,
                             srm_masses,
                             peptide_mass);

    Spectrum prmSpec;
    prmSpec.resize(prm_masses.size());
    for (unsigned int i = 0; i < prm_masses.size(); i++)
    {
      prmSpec[i][0] = prm_masses[i];
    }
    prmSpec.sortPeaks();

    Spectrum srmSpec;
    srmSpec.resize(srm_masses.size());
    for (unsigned int i = 0; i < srm_masses.size(); i++)
    {
      srmSpec[i][0] = srm_masses[i];
    }
    srmSpec.sortPeaks();

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        { // make sure peak is annotated
          float baseMass = ((*spectrum)[i][0] * ((float)curr_frag->charge))
              - curr_frag->massOffset;
          float trueBaseMass = (curr_frag->isNTerm)
              ? prmSpec[prmSpec.findClosest(baseMass)][0]
              : srmSpec[srmSpec.findClosest(baseMass)][0];
          float trueIonMass = (trueBaseMass + curr_frag->massOffset)
              / ((float)curr_frag->charge);
          float DaError = (*spectrum)[i][0] - trueIonMass;

          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            massErrors.push_back(DaError);
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              { // include only listed peaks
                massErrors.push_back(DaError);
              }
            }
          }
        }
      }
    }
    getHistogramInfo(massErrors, bins);
  }

  void SpectrumAnnotStatistics::fragmentMassErrorDaHistogram(const PeptideSpectrumMatchSet &psms,
                                                             SpectrumAnnotParameter &ionParam,
                                                             float maxError,
                                                             vector<TwoValues<
                                                                 float> > &bins)
  {
    const unsigned int defautNumBins = 60;
    if (bins.size() == 0)
    {
      bins.resize(defautNumBins);
      float binWidth = (maxError * 2.0) / ((float)defautNumBins);
      bins[0][0] = 0.0 - maxError;
      for (unsigned int i = 1; i < bins.size(); i++)
      {
        bins[i][0] = bins[i - 1][0] + binWidth;
      }
    }

    for (int psmIdx = 0; psmIdx < psms.size(); psmIdx++)
    {

      PeptideSpectrumMatch& psm = *(psms[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      fragmentMassErrorDaHistogram(psm, ionParam, maxError, bins);
    }
  }

  // -------------------------------------------------------------------------
  void SpectrumAnnotStatistics::fragmentMassErrorPPMHistogram(const PeptideSpectrumMatch &psm,
                                                              SpectrumAnnotParameter &ionParam,
                                                              float maxErrorPPM,
                                                              vector<TwoValues<
                                                                  float> > &bins)
  {
    const unsigned int defautNumBins = 60;
    if (bins.size() == 0)
    {
      bins.resize(defautNumBins);
      float binWidth = (maxErrorPPM * 2.0) / ((float)defautNumBins);
      bins[0][0] = 0.0 - maxErrorPPM;
      for (unsigned int i = 1; i < bins.size(); i++)
      {
        bins[i][0] = bins[i - 1][0] + binWidth;
      }
    }

    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return;
    }

    Spectrum * spectrum = psm.m_spectrum;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    AAJumps jumps(1);

    vector<string> ion_types;

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    vector<float> masses;

    float peptide_mass;

    //generate srm and prm masses, no offsets.
    vector<float> prm_masses;
    vector<float> srm_masses;

    vector<float> massErrors;

    jumps.getPRMandSRMMasses(psm.m_annotation,
                             prm_masses,
                             srm_masses,
                             peptide_mass);

    Spectrum prmSpec;
    prmSpec.resize(prm_masses.size());
    for (unsigned int i = 0; i < prm_masses.size(); i++)
    {
      prmSpec[i][0] = prm_masses[i];
    }
    prmSpec.sortPeaks();

    Spectrum srmSpec;
    srmSpec.resize(srm_masses.size());
    for (unsigned int i = 0; i < srm_masses.size(); i++)
    {
      srmSpec[i][0] = srm_masses[i];
    }
    srmSpec.sortPeaks();

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        { // make sure peak is annotated
          float baseMass = ((*spectrum)[i][0] * ((float)curr_frag->charge))
              - curr_frag->massOffset;
          float trueBaseMass = (curr_frag->isNTerm)
              ? prmSpec[prmSpec.findClosest(baseMass)][0]
              : srmSpec[srmSpec.findClosest(baseMass)][0];
          float trueIonMass = (trueBaseMass + curr_frag->massOffset)
              / ((float)curr_frag->charge);
          float ppmError = 1000000.0 * ((*spectrum)[i][0] - trueIonMass)
              / trueIonMass;

          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            massErrors.push_back(ppmError);
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              { // include only listed peaks
                massErrors.push_back(ppmError);
              }
            }
          }
        }
      }
    }
    getHistogramInfo(massErrors, bins);
  }

  void SpectrumAnnotStatistics::fragmentMassErrorPPMHistogram(const PeptideSpectrumMatchSet &psms,
                                                              SpectrumAnnotParameter &ionParam,
                                                              float maxErrorPPM,
                                                              vector<TwoValues<
                                                                  float> > &bins)
  {
    const unsigned int defautNumBins = 60;
    if (bins.size() == 0)
    {
      bins.resize(defautNumBins);
      float binWidth = (maxErrorPPM * 2.0) / ((float)defautNumBins);
      bins[0][0] = 0.0 - maxErrorPPM;
      for (unsigned int i = 1; i < bins.size(); i++)
      {
        bins[i][0] = bins[i - 1][0] + binWidth;
      }
    }

    for (int psmIdx = 0; psmIdx < psms.size(); psmIdx++)
    {
      PeptideSpectrumMatch& psm = *(psms[psmIdx].get());

      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      fragmentMassErrorPPMHistogram(psm, ionParam, maxErrorPPM, bins);
    }
  }

  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::parentMassErrorDa(const PeptideSpectrumMatch &psm,
                                                   int charge)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      ERROR_MSG("Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    if (charge <= 0)
    { //charge is not defined
      charge = spectrum->parentCharge;
    }
    AAJumps jumps(1);

    double peptide_mass = jumps.getPeptideMass(psm.m_annotation);
    double monoisotopic_mass = peptide_mass;
    double oneisotopic_mass = peptide_mass + AAJumps::massHion; /* one C13 value */

    if (charge <= 0)
    {
      WARN_MSG("Warning: Charge not defined! Calculating from parentMass.");
      charge = peptide_mass / spectrum->parentMass;
      if (charge <= 0)
      {
        WARN_MSG("Charge not able to be calculated from ParentMass, using charge 1");
        charge = 1;
      }
    }

    /* set to parent mass/charge */
    monoisotopic_mass = (monoisotopic_mass + (charge * AAJumps::massHion)
        + AAJumps::massH2O) / charge;
    oneisotopic_mass = (oneisotopic_mass + (charge * AAJumps::massHion)
        + AAJumps::massH2O) / charge;

    float monoisotopic_da;
    float oneisotopic_da;

    if (spectrum->parentMZ > 0)
    {
      /* get error in Da */
      monoisotopic_da = spectrum->parentMZ - monoisotopic_mass;
      oneisotopic_da = spectrum->parentMZ - oneisotopic_mass;
    }
    else
    {
      if (spectrum->parentCharge > 0)
      {
        float parentMZ = (spectrum->parentMass + (AAJumps::massHion
            * (spectrum->parentCharge - 1))) / spectrum->parentCharge;
        monoisotopic_da = parentMZ - monoisotopic_mass;
        oneisotopic_da = parentMZ - oneisotopic_mass;
      }
      else
      {
        monoisotopic_da = spectrum->parentMass - monoisotopic_mass;
        oneisotopic_da = spectrum->parentMass - oneisotopic_mass;
      }
    }

    if (abs(monoisotopic_da) > abs(oneisotopic_da))
    {
      return oneisotopic_da;
    }
    else
    {
      return monoisotopic_da;
    }
  }
  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::observedBreaks(const PeptideSpectrumMatch &psm,
                                                SpectrumAnnotParameter &ionParam)
  {
    if (psm.m_spectrum == (Spectrum *)NULL)
    {
      WARN_MSG("Warning: Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = psm.m_spectrum;

    string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

    AAJumps jumps(1);

    int peptide_length = jumps.getPeptideLength(psm.m_annotation);
    float peptide_mass = jumps.getPeptideMass(psm.m_annotation);

    vector<string> ion_types;

    splitText(ionNames.c_str(), ion_types, (const char*)",");

    set<short> breaks; //keep track of breaks for which we've seen an ion
    //set<short>::iterator breaks_it; //breaks iterator

    //breaks_it = breaks.begin();

    for (unsigned i = 0; i < spectrum->size(); i++)
    {
      if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
          <= (peptide_mass - 37))
      { // ignore small masses or masses larger than we expect.
        const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
        if (curr_frag != NULL)
        {
          if (ionNames.compare("all") == 0)
          { // include all identified peaks
            breaks.insert(psm.m_peakAnnotations[i].second);
          }
          else
          {
            for (int j = 0; j < ion_types.size(); j++)
            {
              if (ion_types[j].compare(curr_frag->name) == 0)
              {
                breaks.insert(psm.m_peakAnnotations[i].second);
              }
            }
          }
        }
      }
    }
    if (peptide_length - 1 > 0)
    {
      return ((float)breaks.size() / ((float)peptide_length - 1)) * 100;
    }
    else
    {
      return 0;
    }
  }
  // -------------------------------------------------------------------------
  float SpectrumAnnotStatistics::observedBreaks(const PeptideSpectrumMatchSet &psmSet,
                                                SpectrumAnnotParameter &ionParam)
  {

    AAJumps jumps(1);

    float totalBreaks = 0;
    float totalPepLength = 0;

    for (int psmIdx = 0; psmIdx < psmSet.size(); psmIdx++)
    {
      PeptideSpectrumMatch& psm = *(psmSet[psmIdx].get());
      if (psm.m_spectrum == (Spectrum *)NULL)
      {
        ERROR_MSG("Spectrum " << psmIdx << " is not defined!");
        return -1.0;
      }

      Spectrum * spectrum = psm.m_spectrum;

      if (spectrum->size() == 0)
      {
        continue;
      }

      string ionNames = ionParam.getFragSpecificIonNames(Spectrum::activationToString(spectrum->msFragType));

      int peptide_length = jumps.getPeptideLength(psm.m_annotation);
      float peptide_mass = jumps.getPeptideMass(psm.m_annotation);

      vector<string> ion_types;

      splitText(ionNames.c_str(), ion_types, (const char*)",");

      set<short> breaks; //keep track of breaks for which we've seen an ion
      //set<short>::iterator breaks_it; //breaks iterator

      //breaks_it = breaks.begin();

      for (unsigned i = 0; i < spectrum->size(); i++)
      {
        if ((*spectrum)[i][0] >= AAJumps::minAAmass && (*spectrum)[i][0]
            <= (peptide_mass - 37))
        { // ignore small masses or masses larger than we expect.
          const ftIonFragment* curr_frag = psm.m_peakAnnotations[i].first;
          if (curr_frag != NULL)
          {
            if (ionNames.compare("all") == 0)
            { // include all identified peaks
              breaks.insert(psm.m_peakAnnotations[i].second);
            }
            else
            {
              for (int j = 0; j < ion_types.size(); j++)
              {
                if (ion_types[j].compare(curr_frag->name) == 0)
                {
                  breaks.insert(psm.m_peakAnnotations[i].second);
                }
              }
            }
          }
        }
      }
      totalBreaks += (float)breaks.size();
      totalPepLength += (float)(peptide_length - 1);
    }
    if (totalPepLength > 0)
    {
      return (totalBreaks / totalPepLength) * 100;
    }
    else
    {
      return 0;
    }
  }
}
