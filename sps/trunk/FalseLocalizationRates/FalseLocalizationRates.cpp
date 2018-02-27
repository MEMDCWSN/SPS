/*
 * FalseLocalizationRates.cpp
 *
 *  Created on: May 14, 2011
 *      Author: jsnedecor
 */

#include "FalseLocalizationRates.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  FalseLocalizationRates::FalseLocalizationRates()
  {
    //void
  }

  // -------------------------------------------------------------------------
  FalseLocalizationRates::FalseLocalizationRates(AAJumps jumps,
                                                 MS2ScoringModel model,
                                                 float peakTol)
  {
    m_jumps = jumps;
    m_model = model;
    m_peakTol = peakTol;
  }
  // -------------------------------------------------------------------------
  FalseLocalizationRates::~FalseLocalizationRates(void)
  {
    //void
  }
  // -------------------------------------------------------------------------
  // helper function for groupVariants. Ignores precursor masses
  float FalseLocalizationRates::getTotalIonCurrent(PeptideSpectrumMatch &psm1,
                                                   PeptideSpectrumMatch &psm2)
  {
    Spectrum * spectrum = psm1.m_spectrum;

    if (spectrum == NULL)
    {
      WARN_MSG("Spectrum not defined!");
    }

    float total = 0;
    for (int i = 0; i < spectrum->size(); i++)
    {

      if (psm1.m_peakAnnotations[i].first == NULL
          && psm2.m_peakAnnotations[i].first == NULL)
      {
        // total += (*spectrum)[i][1]; // there is no associated annotation
      }
      else
      {
        if (psm1.m_peakAnnotations[i].first != NULL
            && psm1.m_peakAnnotations[i].first->isIF
            || psm2.m_peakAnnotations[i].first != NULL
                && psm2.m_peakAnnotations[i].first->isIF)
        {
          total += (*spectrum)[i][1];
        }
      }
    }
    return total;
  }
  // -------------------------------------------------------------------------
  float FalseLocalizationRates::getTotalIonCurrent(Spectrum &modifiedSpectrum,
                                                   Spectrum &allMassesSpectrum)
  {
    float totalIonCurrent = 0;
    for (int i = 0; i < allMassesSpectrum.size(); i++)
    {
      float currMass = allMassesSpectrum[i][0];
      vector<int> matches;
      modifiedSpectrum.findMatches(currMass, m_peakTol, matches);

      for (int j = 0; j < matches.size(); j++)
      {
        totalIonCurrent += modifiedSpectrum[matches[j]][1];
      }
    }
    return totalIonCurrent;

  }
  // -------------------------------------------------------------------------
  bool FalseLocalizationRates::groupHasValidVariant(vector<string> &validSites,
                                                    const vector<float> &modifications,
                                                    vector<string> &variants,
                                                    vector<string> &validVariants)
  {
    if (validSites.size() != modifications.size())
    {
      ERROR_MSG("Valid sites and modification vectors different sizes!");
      return false;
    }

    std::tr1::unordered_map<float, unsigned int> modToValidSites;

    for (unsigned int i = 0; i < modifications.size(); i++)
    {
      modToValidSites[modifications[i]] = i;
    }

    bool validVariant = false;

    for (int varIndex = 0; varIndex < variants.size(); varIndex++)
    {
      size_t aaStart;
      aaStart = variants[varIndex].find_first_of("(");

      vector<bool> foundVariant(modifications.size());
      int modCount = 0;

      while (aaStart != string::npos)
      {
        if (aaStart != string::npos)
        {
          size_t aaEnd = variants[varIndex].find_first_of(",", aaStart);
          string modAA = variants[varIndex].substr(aaStart + 1, aaEnd - aaStart
              - 1);
          size_t modEnd = variants[varIndex].find_first_of(")", aaEnd + 1);
          stringstream ss;

          float modification;
          ss << variants[varIndex].substr(aaEnd + 1, modEnd - aaEnd - 1);
          ss >> modification;

          if (modToValidSites.find(modification) == modToValidSites.end())
          {
            WARN_MSG("Unable to find modification in variant: " << modification);
            return false;
          }
          else
          {
            string currValidSites = validSites[modToValidSites[modification]];
            vector < string > currValidSitesVec;
            stringSplit2(currValidSites, currValidSitesVec, ":");

            for (int k = 0; k < currValidSitesVec.size(); k++)
            {
#ifdef DEBUG
              DEBUG_VAR(currValidSitesVec[k]);
              DEBUG_VAR(aaStart);
              DEBUG_VAR(variants[varIndex]);
              DEBUG_VAR(modCount);
#endif
              if (currValidSitesVec[k].compare("nterm") == 0)
              {
                if (aaStart == 0)
                {
                  foundVariant[modCount] = true;
                }
              }
              else if (currValidSitesVec[k].compare(modAA) == 0)
              {
                foundVariant[modCount] = true;
              }
            }
          }
        }
        modCount++;
        aaStart = variants[varIndex].find_first_of("(", aaStart + 1);
      }

      bool allValid = true;

      for (int i = 0; i < foundVariant.size(); i++)
      {
        allValid = allValid and foundVariant[i];
      }
      if (allValid)
      {
        validVariants.push_back(variants[varIndex]);
        validVariant = true;
      }
    }

    return validVariant;
  }
  // -------------------------------------------------------------------------
  bool FalseLocalizationRates::generateTheoreticalSpectrum(string &unmodifiedPeptide,
                                                           Spectrum &unmodifiedSpectrum,
                                                           vector<
                                                               vector<string> > &variants,
                                                           vector<float> &variantQuantities,
                                                           Spectrum &allPossiblePeaks,
                                                           Spectrum &outputSpectrum)
  {

    Timer tm;
    //Generate expected Intensities annotations;
    tm.start();
    PeptideSpectrumMatch expectedIntensitiesPSM;
    expectedIntensitiesPSM.m_spectrum = &unmodifiedSpectrum;
    string includeIons = "all";

    DEBUG_VAR(unmodifiedPeptide);

    if (!expectedIntensitiesPSM.annotate(unmodifiedPeptide,
                                         includeIons,
                                         m_model,
                                         0,
                                         0,
                                         m_peakTol,
                                         false))
    {
      return false;
    }

    std::tr1::unordered_map<string, int> unmodifiedSpectrumMap;
    vector < vector<int> > unmodifiedSpectrumIons;
    expectedIntensitiesPSM.mapIons(unmodifiedSpectrumIons,
                                   unmodifiedSpectrumMap);

    vector<float> peaksPercentIntensity(allPossiblePeaks.size(), 0);
    vector < std::tr1::unordered_set<string>
        > peakIons(allPossiblePeaks.size());

    for (int varIndex = 0; varIndex < variants.size(); varIndex++)
    {
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allPossiblePeaks;
      string includeIons = "all";
      for (int j = 0; j < variants[varIndex].size(); j++)
      {
        psm.annotate(variants[varIndex][j],
                     includeIons,
                     m_model,
                     0,
                     0,
                     m_peakTol,
                     false,
                     false,
                     true);
      }

      for (int peakIndex = 0; peakIndex < psm.m_spectrum->size(); peakIndex++)
      {
        if (psm.m_peakAnnotations[peakIndex].first != NULL)
        {
          peaksPercentIntensity[peakIndex] += variantQuantities[varIndex];
          stringstream ss;
          ss << psm.m_peakAnnotations[peakIndex].first->name
              << psm.m_peakAnnotations[peakIndex].second;
          string test = ss.str();
          peakIons[peakIndex].insert(ss.str());
        }
      }
    }

    int count = 0;
    Spectrum newSpectrum;
    newSpectrum.resize(peaksPercentIntensity.size());

    for (int peakIndex = 0; peakIndex < peaksPercentIntensity.size(); peakIndex++)
    {
      if (peaksPercentIntensity[peakIndex] > 0)
      {
        TwoValues<float> currValue;
        currValue[0] = 0;
        currValue[1] = 0;

        std::tr1::unordered_set<string>::const_iterator it;

        for (it = peakIons[peakIndex].begin(); it != peakIons[peakIndex].end(); it++)
        {
          string ionName = *it;

          if (unmodifiedSpectrumMap.find(ionName)
              != unmodifiedSpectrumMap.end())
          {
            float expectedIntensity = 0;

            for (int k = 0; k
                < unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]].size(); k++)
            {
              int currPeakIndex =
                  unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]][k];
              expectedIntensity += unmodifiedSpectrum[currPeakIndex][1];
            }

            currValue[1] += expectedIntensity
                * peaksPercentIntensity[peakIndex];
          }
          else
          {
            currValue[1] += 1 * peaksPercentIntensity[peakIndex];
          }
          currValue[0] = allPossiblePeaks[peakIndex][0];
        }
        newSpectrum[count] = currValue;
        count++;
      }
    }
    newSpectrum.resize(count);
    outputSpectrum = newSpectrum;

    tm.stop();
    DEBUG_MSG("Check theoretical time " << tm.duration());

    return true;
  }
  // -------------------------------------------------------------------------
  bool FalseLocalizationRates::generateTheoreticalSpectrum(string &unmodifiedPeptide,
                                                           Spectrum &unmodifiedSpectrum,
                                                           vector<string> &variants,
                                                           vector<float> &variantQuantities,
                                                           Spectrum &allPossiblePeaks,
                                                           Spectrum &outputSpectrum)
  {

    Timer tm;
    //Generate expected Intensities annotations;
    tm.start();
    PeptideSpectrumMatch expectedIntensitiesPSM;
    expectedIntensitiesPSM.m_spectrum = &unmodifiedSpectrum;
    string includeIons = "all";

    DEBUG_VAR(unmodifiedPeptide);

    if (!expectedIntensitiesPSM.annotate(unmodifiedPeptide,
                                         includeIons,
                                         m_model,
                                         0,
                                         0,
                                         m_peakTol,
                                         false))
    {
      return false;
    }

    std::tr1::unordered_map<string, int> unmodifiedSpectrumMap;
    vector < vector<int> > unmodifiedSpectrumIons;
    expectedIntensitiesPSM.mapIons(unmodifiedSpectrumIons,
                                   unmodifiedSpectrumMap);

    vector<float> peaksPercentIntensity(allPossiblePeaks.size(), 0);
    vector < std::tr1::unordered_set<string>
        > peakIons(allPossiblePeaks.size());

    for (int varIndex = 0; varIndex < variants.size(); varIndex++)
    {
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allPossiblePeaks;
      string includeIons = "all";
      psm.annotate(variants[varIndex],
                   includeIons,
                   m_model,
                   0,
                   0,
                   m_peakTol,
                   false);

      for (int peakIndex = 0; peakIndex < psm.m_spectrum->size(); peakIndex++)
      {
        if (psm.m_peakAnnotations[peakIndex].first != NULL)
        {
          peaksPercentIntensity[peakIndex] += variantQuantities[varIndex];
          stringstream ss;
          ss << psm.m_peakAnnotations[peakIndex].first->name
              << psm.m_peakAnnotations[peakIndex].second;
          string test = ss.str();
          peakIons[peakIndex].insert(ss.str());
        }
      }
    }

    int count = 0;
    Spectrum newSpectrum;
    newSpectrum.resize(peaksPercentIntensity.size());

    for (int peakIndex = 0; peakIndex < peaksPercentIntensity.size(); peakIndex++)
    {
      if (peaksPercentIntensity[peakIndex] > 0)
      {
        TwoValues<float> currValue;
        currValue[0] = 0;
        currValue[1] = 0;

        std::tr1::unordered_set<string>::const_iterator it;

        for (it = peakIons[peakIndex].begin(); it != peakIons[peakIndex].end(); it++)
        {
          string ionName = *it;

          if (unmodifiedSpectrumMap.find(ionName)
              != unmodifiedSpectrumMap.end())
          {
            float expectedIntensity = 0;

            for (int k = 0; k
                < unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]].size(); k++)
            {
              int currPeakIndex =
                  unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]][k];
              expectedIntensity += unmodifiedSpectrum[currPeakIndex][1];
            }

            currValue[1] += expectedIntensity
                * peaksPercentIntensity[peakIndex];
          }
          else
          {
            currValue[1] += 1 * peaksPercentIntensity[peakIndex];
          }
          currValue[0] = allPossiblePeaks[peakIndex][0];
        }
        newSpectrum[count] = currValue;
        count++;
      }
    }
    newSpectrum.resize(count);
    outputSpectrum = newSpectrum;

    tm.stop();
    DEBUG_MSG("Check theoretical time " << tm.duration());

    return true;
  }
  // -------------------------------------------------------------------------
  //helper function for generateVariantSequences
  void generateVariantFromPosString(string &unmodifiedPeptide,
                                    string &positionString,
                                    const string &modIndices,
                                    const vector<float> &modifications,
                                    string &outputPeptide)
  {
    vector<unsigned int> positions; //current positions
    vector<float> newMods;
    size_t position = 0;

    while (position != string::npos)
    {
      position = positionString.find_first_of(modIndices, position);

      if (position != string::npos)
      {
        stringstream ss;
        positions.push_back(position);
        ss << positionString[position];
        int modIndex;
        ss >> modIndex;
        newMods.push_back(modifications[modIndex - 1]);
        position = position + 1;
      }
    }

    PeptideSpectrumMatch::insertModifications(unmodifiedPeptide,
                                              newMods,
                                              positions,
                                              outputPeptide);
  }

  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateVariantSequences(string &unmodifiedPeptide,
                                                        const vector<float> &modifications,
                                                        vector<string> &outputVariants)
  {
    stringstream ss;
    stringstream modIndices;
    ss << std::setw(unmodifiedPeptide.length() - modifications.size() + 1)
        << std::setfill('0');
    for (int i = 0; i < modifications.size(); i++)
    {
      ss << i + 1;
      modIndices << i + 1;
    }
    string input = ss.str();
    DEBUG_VAR(input);

    string variant;
    generateVariantFromPosString(unmodifiedPeptide,
                                 input,
                                 modIndices.str(),
                                 modifications,
                                 variant);
    outputVariants.push_back(variant);

    while (next_permutation(input.begin(), input.end()))
    {
      generateVariantFromPosString(unmodifiedPeptide,
                                   input,
                                   modIndices.str(),
                                   modifications,
                                   variant);
      outputVariants.push_back(variant);
    }
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateVariantSequences(string &unmodifiedPeptide,
                                                        const vector<vector<
                                                            float> > &modifications,
                                                        vector<vector<string> > &outputVariants)
  {
    stringstream ss;
    stringstream modIndices;
    ss << std::setw(unmodifiedPeptide.length() - modifications.size() + 1)
        << std::setfill('0');
    for (int i = 0; i < modifications.size(); i++)
    {
      ss << i + 1;
      modIndices << i + 1;
    }
    string input = ss.str();
    DEBUG_VAR(input);

    vector < vector<float> > tempModifications;
    vector<float> currTempMods(modifications.size());

    for (int i = 0; i < modifications.size(); i++)
    {
      vector < vector<float> > newTempModifications;
      for (int j = 0; j < modifications[i].size(); j++)
      {
        vector < vector<float> > temp = tempModifications;

        if (temp.size() > 0)
        {
          for (int k = 0; k < temp.size(); k++)
          {
            temp[k].push_back(modifications[i][j]);
            newTempModifications.push_back(temp[k]);
          }
        }
        else
        {
          vector<float> newTemp;
          newTemp.push_back(modifications[i][j]);
          newTempModifications.push_back(newTemp);
        }
      }
      tempModifications = newTempModifications;
    }

    string variant;

    vector < string > variantString;
    variantString.push_back(input);

    while (next_permutation(input.begin(), input.end()))
    {
      variantString.push_back(input);
    }
    outputVariants.resize(variantString.size());
    for (int i = 0; i < variantString.size(); i++)
    {
      for (int j = 0; j < tempModifications.size(); j++)
      {
        string variant;
        generateVariantFromPosString(unmodifiedPeptide,
                                     variantString[i],
                                     modIndices.str(),
                                     tempModifications[j],
                                     variant);
        outputVariants[i].push_back(variant);
      }
    }
  }
  // -------------------------------------------------------------------------
  bool getIndices(const int &v1, const int &v2, int &i, int &j)
  {
    if (v1 == v2) //we don't get similarity for the same indices
      return false;

    if (v1 < v2) //v1 is row, v2 is column
    {
      i = v1;
      j = v2 - v1 - 1;
    }
    else
    {
      i = v2;
      j = v1 - v2 - 1;
    }
  }
  // -------------------------------------------------------------------------
  bool getVariantFromIndex(const int &i, const int &j, int &v1, int &v2)
  {
    v1 = i;
    v2 = i + j + 1;
  }
  // -------------------------------------------------------------------------
  void getVariantIndices(const vector<string> &variantGroup, const vector<
      string> &variantSeq, vector<int> &variantSeqIndices)
  {
    variantSeqIndices.resize(variantGroup.size());
    for (int i = 0; i < variantSeq.size(); i++)
    {
      for (int j = 0; j < variantGroup.size(); j++)
      {
        if (variantSeq[i].compare(variantGroup[j]) == 0)
        {
          variantSeqIndices[j] = i;
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void getVariantIndices(const vector<string> &variantGroup, const vector<
      vector<string> > &variantSeq, vector<int> &variantSeqIndices)
  {
    variantSeqIndices.resize(variantGroup.size());
    for (int i = 0; i < variantSeq.size(); i++)
    {
      for (int j = 0; j < variantGroup.size(); j++)
      {
        if (variantSeq[i][0].compare(variantGroup[j]) == 0)
        {
          variantSeqIndices[j] = i;
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  float getIntensity(Spectrum &spectrum,
                     const float mass,
                     const float peakTol,
                     int &startIndex)
  {
    vector<int> matches;

    spectrum.findMatches(mass, peakTol, matches, startIndex);

    float intensity = 0;
    for (int i = 0; i < matches.size(); i++)
    {
      intensity += spectrum[matches[i]][1];
    }
    if (matches.size() > 0)
      startIndex = matches[matches.size() - 1];
    return intensity;
  }
  // -------------------------------------------------------------------------
  bool getMinimumFromMatrix(vector<vector<float> > &variantMatrix,
                            float totalIntensity,
                            int &minRow,
                            int &minColumn,
                            float &cutoff)
  {
    float minimumValue = 1000;
    for (int row = 0; row < variantMatrix.size(); row++)
    {
      for (int column = 0; column < variantMatrix[row].size(); column++)
      {
        if (minimumValue > variantMatrix[row][column] / totalIntensity)
        {
          minimumValue = variantMatrix[row][column] / totalIntensity;
          minRow = row;
          minColumn = column;
        }
      }
    }
#ifdef DEBUG
    DEBUG_VAR(minimumValue);
#endif
    return cutoff > minimumValue;
  }
  // -------------------------------------------------------------------------
  float getAverageDistance(vector<string> &group1,
                           vector<string> &group2,
                           vector<string> &variantSeq,
                           vector<vector<float> > &variantMatrix)
  {
    vector<int> group1Indices;
    vector<int> group2Indices;
    getVariantIndices(group1, variantSeq, group1Indices);
    getVariantIndices(group2, variantSeq, group2Indices);

    int count = 0;
    float total = 0;
    for (int i = 0; i < group1Indices.size(); i++)
    {
      for (int j = 0; j < group2Indices.size(); j++)
      {
        count++;
        int row, column = 0;
        if (getIndices(group1Indices[i], group2Indices[j], row, column))
        {
          total += variantMatrix[row][column];
        }
      }
    }
    return (float)total / (float)count;
  }
  // -------------------------------------------------------------------------
  float getAverageDistance(vector<string> &group1,
                           vector<string> &group2,
                           vector<vector<string> > &variantSeq,
                           vector<vector<float> > &variantMatrix)
  {
    vector<int> group1Indices;
    vector<int> group2Indices;
    getVariantIndices(group1, variantSeq, group1Indices);
    getVariantIndices(group2, variantSeq, group2Indices);

    int count = 0;
    float total = 0;
    for (int i = 0; i < group1Indices.size(); i++)
    {
      for (int j = 0; j < group2Indices.size(); j++)
      {
        count++;
        int row, column = 0;
        if (getIndices(group1Indices[i], group2Indices[j], row, column))
        {
          total += variantMatrix[row][column];
        }
      }
    }
    return (float)total / (float)count;
  }
  // -------------------------------------------------------------------------
  bool getMinimumFromGroup(vector<vector<string> > &variantGroups,
                           vector<vector<string> > &variantSeq,
                           vector<vector<float> > &variantMatrix,
                           float totalIntensity,
                           int &group1Index,
                           int &group2Index,
                           float &cutoff)
  {
    float minimum = 1000;
    for (int i = 0; i < variantGroups.size(); i++)
    {
      for (int j = i + 1; j < variantGroups.size(); j++)
      {
        vector < string > *group1 = &(variantGroups[i]);
        vector < string > *group2 = &(variantGroups[j]);
        if (group1->size() > 0 && group2->size() > 0)
        {
          float averageDistance = getAverageDistance(*group1,
                                                     *group2,
                                                     variantSeq,
                                                     variantMatrix)
              / totalIntensity;

          if (minimum > averageDistance)
          {
            minimum = averageDistance;
            group1Index = i;
            group2Index = j;
          }
        }
      }
    }
    return minimum < cutoff;
  }
  // -------------------------------------------------------------------------
  bool getMinimumFromGroup(vector<vector<string> > &variantGroups,
                           vector<string> &variantSeq,
                           vector<vector<float> > &variantMatrix,
                           float totalIntensity,
                           int &group1Index,
                           int &group2Index,
                           float &cutoff)
  {
    float minimum = 1000;
    for (int i = 0; i < variantGroups.size(); i++)
    {
      for (int j = i + 1; j < variantGroups.size(); j++)
      {
        vector < string > *group1 = &(variantGroups[i]);
        vector < string > *group2 = &(variantGroups[j]);
        if (group1->size() > 0 && group2->size() > 0)
        {
          float averageDistance = getAverageDistance(*group1,
                                                     *group2,
                                                     variantSeq,
                                                     variantMatrix)
              / totalIntensity;

          if (minimum > averageDistance)
          {
            minimum = averageDistance;
            group1Index = i;
            group2Index = j;
          }
        }
      }
    }
    return minimum < cutoff;
  }
  // -------------------------------------------------------------------------
  void initializeVariantGroups(vector<string> &variantSeq,
                               vector<vector<string> > &variantGroups,
                               std::tr1::unordered_map<string, unsigned int> &variantMap,
                               const vector<float> &quantities,
                               vector<float> &variantQuantities)
  {
    variantQuantities.resize(quantities.size());
    for (int i = 0; i < variantSeq.size(); i++)
    {
      vector < string > currGroup;
      currGroup.push_back(variantSeq[i]);
      variantGroups.push_back(currGroup);
      variantMap[variantSeq[i]] = i;
      variantQuantities[i] = quantities[i];
    }
  }
  // -------------------------------------------------------------------------
  void initializeVariantGroups(vector<vector<string> > &variantSeq,
                               vector<vector<string> > &variantGroups,
                               std::tr1::unordered_map<string, unsigned int> &variantMap,
                               const vector<float> &quantities,
                               vector<float> &variantQuantities)
  {
    variantQuantities.resize(quantities.size());
    for (int i = 0; i < variantSeq.size(); i++)
    {
      vector < string > currGroup;
      currGroup.push_back(variantSeq[i][0]);
      variantGroups.push_back(currGroup);
      variantMap[variantSeq[i][0]] = i;
      variantQuantities[i] = quantities[i];
    }
  }
  // -------------------------------------------------------------------------
  void getInverse(const vector<int> &input, vector<int> &output, int size)
  {
    vector<int> temp(size, 1);
    for (int i = 0; i < input.size(); i++)
    {
      temp[input[i]] = 0;
    }

    for (int i = 0; i < temp.size(); i++)
    {
      if (temp[i] == 1)
      {
        output.push_back(i);
      }
    }
  }
  // -------------------------------------------------------------------------
  void mergeGroups(vector<vector<float> > &variantMatrix,
                   vector<vector<string> > &variantGroups,
                   std::tr1::unordered_map<string, unsigned int> &variantMap,
                   vector<float> &variantQuantities,
                   int &variantGroup1Index,
                   int &variantGroup2Index)
  {
    //merge variants in group
    vector < string > newGroup;
    for (int i = 0; i < variantGroups[variantGroup1Index].size(); i++)
    {
      string currVariant = variantGroups[variantGroup1Index][i];
      newGroup.push_back(currVariant);
      variantMap[currVariant] = variantGroup1Index;
    }
    for (int i = 0; i < variantGroups[variantGroup2Index].size(); i++)
    {
      string currVariant = variantGroups[variantGroup2Index][i];
      newGroup.push_back(currVariant);
      variantMap[currVariant] = variantGroup1Index;
    }
    //set to new group
    variantGroups[variantGroup1Index] = newGroup;
    variantGroups[variantGroup2Index].resize(0);
    variantQuantities[variantGroup1Index]
        = variantQuantities[variantGroup1Index]
            + variantQuantities[variantGroup2Index];
    variantQuantities[variantGroup2Index] = -1;

    //zero out group1
    for (int i = 0; i < variantGroups.size(); i++)
    {
      if (variantQuantities[i] >= 0)
      {
        int row;
        int column;
#ifdef DEBUG
        DEBUG_VAR(variantGroup1Index);
#endif
        if (getIndices(i, variantGroup1Index, row, column))
        {
#ifdef DEBUG
          DEBUG_VAR(i);
          DEBUG_VAR(row);
          DEBUG_VAR(column);
#endif
          variantMatrix[row][column] = 0;
        }
      }
    }

    //zero out group2
    for (int i = 0; i < variantGroups.size(); i++)
    {
      if (variantQuantities[i] >= 0)
      {
        int row;
        int column;
        if (getIndices(i, variantGroup2Index, row, column))
        {
          variantMatrix[row][column] = INT_MAX;
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  void mergeGroups2(vector<vector<string> > &variantGroups,
                    std::tr1::unordered_map<string, unsigned int> &variantMap,
                    vector<float> &variantQuantities,
                    int &variantGroup1Index,
                    int &variantGroup2Index)
  {
    //merge variants in group
    vector < string > newGroup;
    for (int i = 0; i < variantGroups[variantGroup1Index].size(); i++)
    {
      string currVariant = variantGroups[variantGroup1Index][i];
      newGroup.push_back(currVariant);
      variantMap[currVariant] = variantGroup1Index;
    }

    for (int i = 0; i < variantGroups[variantGroup2Index].size(); i++)
    {
      string currVariant = variantGroups[variantGroup2Index][i];
      newGroup.push_back(currVariant);
      variantMap[currVariant] = variantGroup1Index;
    }

    //set to new group
    variantGroups[variantGroup1Index] = newGroup;
    variantGroups[variantGroup2Index].resize(0);
    variantQuantities[variantGroup1Index]
        = variantQuantities[variantGroup1Index]
            + variantQuantities[variantGroup2Index];
    variantQuantities[variantGroup2Index] = -1;
  }
  // -------------------------------------------------------------------------
  bool contains(const vector<int> &vector1, const vector<int> &vector2) //does any element of vector1 show up in vector2
  {
    bool passedVector = false;
    for (int i = 0; i < vector1.size(); i++)
    {
      for (int j = 0; j < vector2.size(); j++)
      {
        if (vector1[i] == vector2[j])
        {
          passedVector = true;
        }
      }
    }
    return passedVector;
  }
  // -------------------------------------------------------------------------
  // -- use max distance, all distinguishing intensity --
  void recalculateDistances(vector<vector<float> > &variantMatrix,
                            vector<
                                std::tr1::unordered_map<string, vector<int> > > &peaksMappedToVariant,
                            vector<
                                std::tr1::unordered_map<string, vector<int> > > &variantMappedToPeaks,
                            const vector<vector<string> > &variantGroups,
                            const vector<string> &variantSeq,
                            const int &variantGroupIndex,
                            const Spectrum &allMassesSpectrum,
                            Spectrum &modifiedSpectrum,
                            const float &peakTol)
  {
    const vector<string> * currVariantGroup =
        &(variantGroups[variantGroupIndex]);
    vector<int> variantGroupIndices;
    getVariantIndices(*currVariantGroup, variantSeq, variantGroupIndices);

    int startIndex = -1;
    for (int peakIdx = 0; peakIdx < allMassesSpectrum.size(); peakIdx++)
    {
      std::tr1::unordered_map < string, vector<int> > *currMap
          = &(peaksMappedToVariant[peakIdx]);

      float intensity = getIntensity(modifiedSpectrum,
                                     allMassesSpectrum[peakIdx][0],
                                     peakTol,
                                     startIndex);
#ifdef DEBUG
      DEBUG_VAR(allMassesSpectrum[peakIdx][0]);
      DEBUG_VAR(intensity);
#endif

      std::tr1::unordered_map<string, vector<int> >::const_iterator ionIter;

      int count = 0;

      std::tr1::unordered_set < string > passedMatches; //if we have a peak with two assignments, we only count it once.

      for (ionIter = currMap->begin(); ionIter != currMap->end(); ionIter++)
      {
        count++;
#ifdef DEBUG
        DEBUG_VAR(count);
#endif
        const vector<int> * currVariants = &(ionIter->second);

#ifdef DEBUG
        DEBUG_VAR(peakIdx);
        DEBUG_VAR(ionIter->first);
        for (int var = 0; var < currVariants->size(); var++)
        {
          DEBUG_VAR((*currVariants)[var]);
        }
#endif

        vector<int> inverseVariants;
        getInverse(*currVariants, inverseVariants, variantSeq.size());

#ifdef DEBUG
        for (int var = 0; var < inverseVariants.size(); var++)
        {
          DEBUG_VAR(inverseVariants[var]);
        }
#endif
        if (contains(variantGroupIndices, *currVariants))
        {
          for (int i = 0; i < inverseVariants.size(); i++)
          {
            int row, column;

            stringstream ss;
            if (getIndices(variantGroupIndex, inverseVariants[i], row, column))
            {

              ss << row << "_" << column;
              if (passedMatches.count(ss.str()) == 0)
              {
#ifdef DEBUG
                DEBUG_VAR(row);
                DEBUG_VAR(column);
                DEBUG_VAR(variantMatrix[row][column]);
#endif
                variantMatrix[row][column] += intensity;
                passedMatches.insert(ss.str());
              }
            }
          }
        }

        if (contains(variantGroupIndices, inverseVariants))
        {
          for (int i = 0; i < currVariants->size(); i++)
          {
            int row, column;

            stringstream ss;
            if (getIndices(variantGroupIndex, (*currVariants)[i], row, column))
            {
              ss << row << "_" << column;
              if (passedMatches.count(ss.str()) == 0)
              {
#ifdef DEBUG
                DEBUG_VAR(row);
                DEBUG_VAR(column);
                DEBUG_VAR(variantMatrix[row][column]);
#endif
                variantMatrix[row][column] += intensity;
                passedMatches.insert(ss.str());
              }
            }
          }
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void initializePeakMap(vector<string> &variantSeq,
                         Spectrum &allMassesSpectrum,
                         vector<std::tr1::unordered_map<string, vector<int> > > &peaksMappedToVariant,
                         vector<std::tr1::unordered_map<string, vector<int> > > &variantMappedToPeaks,
                         MS2ScoringModel &model)
  {
    string includeIons = "all";

    for (int i = 0; i < variantSeq.size(); i++)
    {
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allMassesSpectrum;

      psm.annotate(variantSeq[i], includeIons, model, 0, 0, 0.001, false);

      for (int j = 0; j < allMassesSpectrum.size(); j++)
      {
        const ftIonFragment * currFrag = psm.m_peakAnnotations[j].first;
        if (currFrag != NULL && currFrag->isIF)
        {
          stringstream ss;
          ss << currFrag->name << psm.m_peakAnnotations[j].second;

          peaksMappedToVariant[j][ss.str()].push_back(i);
          variantMappedToPeaks[i][ss.str()].push_back(j);
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void initializePeakMap(vector<vector<string> > &variantSeq,
                         Spectrum &allMassesSpectrum,
                         vector<std::tr1::unordered_map<string, vector<int> > > &peaksMappedToVariant,
                         vector<std::tr1::unordered_map<string, vector<int> > > &variantMappedToPeaks,
                         MS2ScoringModel &model)
  {
    string includeIons = "all";

    for (int i = 0; i < variantSeq.size(); i++)
    {
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allMassesSpectrum;

      for (int j = 0; j < variantSeq[i].size(); j++)
      {
        psm.annotate(variantSeq[i][j],
                     includeIons,
                     model,
                     0,
                     0,
                     0.001,
                     false,
                     false,
                     true);
      }

      for (int j = 0; j < allMassesSpectrum.size(); j++)
      {
        const ftIonFragment * currFrag = psm.m_peakAnnotations[j].first;
        if (currFrag != NULL && currFrag->isIF)
        {
          stringstream ss;
          ss << currFrag->name << psm.m_peakAnnotations[j].second;

          peaksMappedToVariant[j][ss.str()].push_back(i);
          variantMappedToPeaks[i][ss.str()].push_back(j);
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void printMatrix(vector<string> &variantSeq,
                   vector<vector<float> > &variantMatrix)
  {
    for (int varIdx = 0; varIdx < variantSeq.size(); varIdx++)
    {
      cout << variantSeq[varIdx] << "\t";
    }
    cout << endl;
    for (int row = 0; row < variantMatrix.size(); row++)
    {
      cout << variantSeq[row] << "\t";
      int prefixSize = variantSeq.size() - variantMatrix[row].size();
      for (int prefixIdx = 0; prefixIdx < prefixSize - 1; prefixIdx++)
      {
        cout << "\t";
      }
      for (int column = 0; column < variantMatrix[row].size(); column++)
      {
        cout << variantMatrix[row][column] << "\t";
      }
      cout << endl;
    }
  }
  // -------------------------------------------------------------------------
  void initializeDistanceMatrix(vector<vector<string> > &variantSeq,
                                vector<vector<float> > &variantMatrix,
                                Spectrum &allMassesSpectrum,
                                Spectrum &modifiedSpectrum,
                                vector<std::tr1::unordered_map<string, vector<
                                    int> > > &peaksMappedToVariant,
                                float &peakTol)
  {
    //initialize similarity matrix
    for (int i = 0; i < variantMatrix.size(); i++)
    {
      int rowSize = variantMatrix.size() - i - 1;
      if (rowSize > 0)
      {
        variantMatrix[i].resize(rowSize, 0);
      }
      else
      {
        variantMatrix[i].resize(0);
      }
    }

    int startIndex = -1;

    for (int peakIdx = 0; peakIdx < allMassesSpectrum.size(); peakIdx++)
    {
      std::tr1::unordered_map < string, vector<int> > *currMap
          = &(peaksMappedToVariant[peakIdx]);

      float intensity = getIntensity(modifiedSpectrum,
                                     allMassesSpectrum[peakIdx][0],
                                     peakTol,
                                     startIndex);
#ifdef DEBUG
      DEBUG_VAR(allMassesSpectrum[peakIdx][0]);
      DEBUG_VAR(intensity);
#endif

      std::tr1::unordered_map<string, vector<int> >::const_iterator ionIter;

      int count = 0;

      std::tr1::unordered_set < string > passedMatches; //if we have a peak with two assignments, we only count it once.

      for (ionIter = currMap->begin(); ionIter != currMap->end(); ionIter++)
      {
        count++;
#ifdef DEBUG
        DEBUG_VAR(count);
#endif
        const vector<int> * currVariants = &(ionIter->second);
#ifdef DEBUG
        DEBUG_VAR(peakIdx);
        DEBUG_VAR(ionIter->first);

        for (int var = 0; var < currVariants->size(); var++)
        {
          DEBUG_VAR((*currVariants)[var]);
        }
#endif

        vector<int> inverseVariants;
        getInverse(*currVariants, inverseVariants, variantSeq.size());

#ifdef DEBUG
        for (int var = 0; var < inverseVariants.size(); var++)
        {
          DEBUG_VAR(inverseVariants[var]);
        }
#endif

        for (int i = 0; i < currVariants->size(); i++) //row
        {
          for (int j = 0; j < inverseVariants.size(); j++) //column
          {
            int row, column = -1;

            stringstream ss;
            if (getIndices((*currVariants)[i], inverseVariants[j], row, column))
            {
              ss << row << "_" << column;
              if (passedMatches.count(ss.str()) == 0)
              {
#ifdef DEBUG
                DEBUG_VAR(row);
                DEBUG_VAR(column);
                DEBUG_VAR(variantMatrix[row][column]);
#endif
                variantMatrix[row][column] += intensity;
                passedMatches.insert(ss.str());
              }
            }
          }
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void initializeDistanceMatrix(vector<string> &variantSeq,
                                vector<vector<float> > &variantMatrix,
                                Spectrum &allMassesSpectrum,
                                Spectrum &modifiedSpectrum,
                                vector<std::tr1::unordered_map<string, vector<
                                    int> > > &peaksMappedToVariant,
                                float &peakTol)
  {
    //initialize similarity matrix
    for (int i = 0; i < variantMatrix.size(); i++)
    {
      int rowSize = variantMatrix.size() - i - 1;
      if (rowSize > 0)
      {
        variantMatrix[i].resize(rowSize, 0);
      }
      else
      {
        variantMatrix[i].resize(0);
      }
    }

    int startIndex = -1;

    for (int peakIdx = 0; peakIdx < allMassesSpectrum.size(); peakIdx++)
    {
      std::tr1::unordered_map < string, vector<int> > *currMap
          = &(peaksMappedToVariant[peakIdx]);

      float intensity = getIntensity(modifiedSpectrum,
                                     allMassesSpectrum[peakIdx][0],
                                     peakTol,
                                     startIndex);
#ifdef DEBUG
      DEBUG_VAR(allMassesSpectrum[peakIdx][0]);
      DEBUG_VAR(intensity);
#endif

      std::tr1::unordered_map<string, vector<int> >::const_iterator ionIter;

      int count = 0;

      std::tr1::unordered_set < string > passedMatches; //if we have a peak with two assignments, we only count it once.

      for (ionIter = currMap->begin(); ionIter != currMap->end(); ionIter++)
      {
        count++;
#ifdef DEBUG
        DEBUG_VAR(count);
#endif
        const vector<int> * currVariants = &(ionIter->second);
#ifdef DEBUG
        DEBUG_VAR(peakIdx);
        DEBUG_VAR(ionIter->first);

        for (int var = 0; var < currVariants->size(); var++)
        {
          DEBUG_VAR((*currVariants)[var]);
        }
#endif

        vector<int> inverseVariants;
        getInverse(*currVariants, inverseVariants, variantSeq.size());

#ifdef DEBUG
        for (int var = 0; var < inverseVariants.size(); var++)
        {
          DEBUG_VAR(inverseVariants[var]);
        }
#endif

        for (int i = 0; i < currVariants->size(); i++) //row
        {
          for (int j = 0; j < inverseVariants.size(); j++) //column
          {
            int row, column = -1;

            stringstream ss;
            if (getIndices((*currVariants)[i], inverseVariants[j], row, column))
            {
              ss << row << "_" << column;
              if (passedMatches.count(ss.str()) == 0)
              {
#ifdef DEBUG
                DEBUG_VAR(row);
                DEBUG_VAR(column);
                DEBUG_VAR(variantMatrix[row][column]);
#endif
                variantMatrix[row][column] += intensity;
                passedMatches.insert(ss.str());
              }
            }
          }
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::groupVariants4(Spectrum &modifiedSpectrum,
                                              Spectrum &allMassesSpectrum,
                                              vector<vector<string> > &variantSeq,
                                              const vector<float> &quantities,
                                              const vector<float> &massShifts,
                                              vector<vector<string> > &variantGroups,
                                              vector<float> &variantQuantities,
                                              vector<float> &minDistinguishingIntensity,
                                              float minimumGroupingPercentIntensity)
  {
    DEBUG_VAR(minimumGroupingPercentIntensity);
    Timer tm;
    tm.start();
    vector < std::tr1::unordered_map<string, vector<int> >
        > peaksMappedToVariant(allMassesSpectrum.size()); //peak indices associated with peak annotation and variant indices
    //variantIndexMappedToPeaks[peakIndex][peakAnnotation][variantIndex]

    vector < std::tr1::unordered_map<string, vector<int> >
        > variantMappedToPeaks(variantSeq.size()); //variantIndices associated with peak indices. variantMappedToPeaks[variant][peakIndex]

    tm.start();

    initializePeakMap(variantSeq,
                      allMassesSpectrum,
                      peaksMappedToVariant,
                      variantMappedToPeaks,
                      m_model);

    tm.stop();
    DEBUG_MSG("STEP 1 " << tm.duration());

    vector < vector<float> > variantMatrix(variantSeq.size()); //similarity matrix

    initializeDistanceMatrix(variantSeq,
                             variantMatrix,
                             allMassesSpectrum,
                             modifiedSpectrum,
                             peaksMappedToVariant,
                             m_peakTol);

    float totalIntensity = getTotalIonCurrent(modifiedSpectrum,
                                              allMassesSpectrum);

    tm.stop();
    DEBUG_MSG("STEP 2 " << tm.duration());

    std::tr1::unordered_map<string, unsigned int> variantMap;
#ifdef DEBUG
    printMatrix(variantSeq,variantMatrix);
#endif
    tm.start();
    //create variant groups
    initializeVariantGroups(variantSeq,
                            variantGroups,
                            variantMap,
                            quantities,
                            variantQuantities);
    tm.stop();
    DEBUG_MSG("initialization " << tm.duration());

    //find minimum from matrix
    int variantGroup1, variantGroup2 = 0;
    while (getMinimumFromGroup(variantGroups,
                               variantSeq,
                               variantMatrix,
                               totalIntensity,
                               variantGroup1,
                               variantGroup2,
                               minimumGroupingPercentIntensity))
    {
      tm.start();
#ifdef DEBUG
      DEBUG_VAR(variantGroup1);
      DEBUG_VAR(variantGroup2);
#endif
      //merge groups
      mergeGroups2(variantGroups,
                   variantMap,
                   variantQuantities,
                   variantGroup1,
                   variantGroup2);

      tm.stop();
      DEBUG_MSG("ITERATION " << tm.duration());
    }

  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::groupVariants4(Spectrum &modifiedSpectrum,
                                              Spectrum &allMassesSpectrum,
                                              vector<string> &variantSeq,
                                              const vector<float> &quantities,
                                              const vector<float> &massShifts,
                                              vector<vector<string> > &variantGroups,
                                              vector<float> &variantQuantities,
                                              vector<float> &minDistinguishingIntensity,
                                              float minimumGroupingPercentIntensity)
  {
    DEBUG_VAR(minimumGroupingPercentIntensity);
    Timer tm;
    tm.start();
    vector < std::tr1::unordered_map<string, vector<int> >
        > peaksMappedToVariant(allMassesSpectrum.size()); //peak indices associated with peak annotation and variant indices
    //variantIndexMappedToPeaks[peakIndex][peakAnnotation][variantIndex]

    vector < std::tr1::unordered_map<string, vector<int> >
        > variantMappedToPeaks(variantSeq.size()); //variantIndices associated with peak indices. variantMappedToPeaks[variant][peakIndex]

    tm.start();

    initializePeakMap(variantSeq,
                      allMassesSpectrum,
                      peaksMappedToVariant,
                      variantMappedToPeaks,
                      m_model);

    tm.stop();
    DEBUG_MSG("STEP 1 " << tm.duration());

    vector < vector<float> > variantMatrix(variantSeq.size()); //similarity matrix

    initializeDistanceMatrix(variantSeq,
                             variantMatrix,
                             allMassesSpectrum,
                             modifiedSpectrum,
                             peaksMappedToVariant,
                             m_peakTol);

    float totalIntensity = getTotalIonCurrent(modifiedSpectrum,
                                              allMassesSpectrum);

    tm.stop();
    DEBUG_MSG("STEP 2 " << tm.duration());

    std::tr1::unordered_map<string, unsigned int> variantMap;
#ifdef DEBUG
    printMatrix(variantSeq,variantMatrix);
#endif
    tm.start();
    //create variant groups
    initializeVariantGroups(variantSeq,
                            variantGroups,
                            variantMap,
                            quantities,
                            variantQuantities);
    tm.stop();
    DEBUG_MSG("initialization " << tm.duration());

    //find minimum from matrix
    int variantGroup1, variantGroup2 = 0;
    while (getMinimumFromGroup(variantGroups,
                               variantSeq,
                               variantMatrix,
                               totalIntensity,
                               variantGroup1,
                               variantGroup2,
                               minimumGroupingPercentIntensity))
    {
      tm.start();
#ifdef DEBUG
      DEBUG_VAR(variantGroup1);
      DEBUG_VAR(variantGroup2);
#endif
      //merge groups
      mergeGroups2(variantGroups,
                   variantMap,
                   variantQuantities,
                   variantGroup1,
                   variantGroup2);

      tm.stop();
      DEBUG_MSG("ITERATION " << tm.duration());
    }

  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::groupVariants3(Spectrum &modifiedSpectrum,
                                              Spectrum &allMassesSpectrum,
                                              vector<string> &variantSeq,
                                              const vector<float> &quantities,
                                              const vector<float> &massShifts,
                                              vector<vector<string> > &variantGroups,
                                              vector<float> &variantQuantities,
                                              vector<float> &minDistinguishingIntensity,
                                              float minimumGroupingPercentIntensity)
  {
    DEBUG_VAR(minimumGroupingPercentIntensity);
    Timer tm;
    tm.start();
    vector < std::tr1::unordered_map<string, vector<int> >
        > peaksMappedToVariant(allMassesSpectrum.size()); //peak indices associated with peak annotation and variant indices
    //variantIndexMappedToPeaks[peakIndex][peakAnnotation][variantIndex]

    vector < std::tr1::unordered_map<string, vector<int> >
        > variantMappedToPeaks(variantSeq.size()); //variantIndices associated with peak indices. variantMappedToPeaks[variant][peakIndex]

    tm.start();

    initializePeakMap(variantSeq,
                      allMassesSpectrum,
                      peaksMappedToVariant,
                      variantMappedToPeaks,
                      m_model);

    tm.stop();
    DEBUG_MSG("STEP 1 " << tm.duration());

    vector < vector<float> > variantMatrix(variantSeq.size()); //similarity matrix

    initializeDistanceMatrix(variantSeq,
                             variantMatrix,
                             allMassesSpectrum,
                             modifiedSpectrum,
                             peaksMappedToVariant,
                             m_peakTol);

    float totalIntensity = getTotalIonCurrent(modifiedSpectrum,
                                              allMassesSpectrum);

    tm.stop();
    DEBUG_MSG("STEP 2 " << tm.duration());

    std::tr1::unordered_map<string, unsigned int> variantMap;

#ifdef DEBUG
    printMatrix(variantSeq,variantMatrix);
#endif
    tm.start();
    //create variant groups
    initializeVariantGroups(variantSeq,
                            variantGroups,
                            variantMap,
                            quantities,
                            variantQuantities);
    tm.stop();
    DEBUG_MSG("initialization " << tm.duration());

    //find minimum from matrix
    int minRow, minColumn;

    while (getMinimumFromMatrix(variantMatrix,
                                totalIntensity,
                                minRow,
                                minColumn,
                                minimumGroupingPercentIntensity))
    {
      tm.start();
#ifdef DEBUG
      DEBUG_VAR(minRow);
      DEBUG_VAR(minColumn);
#endif
      int variantGroup1;
      int variantGroup2;
      getVariantFromIndex(minRow, minColumn, variantGroup1, variantGroup2);
#ifdef DEBUG
      DEBUG_VAR(variantGroup1);
      DEBUG_VAR(variantGroup2);
#endif
      //merge groups
      mergeGroups(variantMatrix,
                  variantGroups,
                  variantMap,
                  variantQuantities,
                  variantGroup1,
                  variantGroup2);

#ifdef DEBUG
      printMatrix(variantSeq,variantMatrix);
#endif
      recalculateDistances(variantMatrix,
                           peaksMappedToVariant,
                           variantMappedToPeaks,
                           variantGroups,
                           variantSeq,
                           variantGroup1,
                           allMassesSpectrum,
                           modifiedSpectrum,
                           m_peakTol);

#ifdef DEBUG
      printMatrix(variantSeq,variantMatrix);
#endif
      tm.stop();
      DEBUG_MSG("ITERATION " << tm.duration());
    }
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::groupVariants2(Spectrum &modifiedSpectrum,
                                              vector<string> &variantSeq,
                                              const vector<float> &quantities,
                                              const vector<float> &massShifts,
                                              vector<vector<string> > &variantGroups,
                                              vector<float> &variantQuantities,
                                              vector<float> &minDistinguishingIntensity,
                                              float minimumGroupingPercentIntensity)
  {
    Timer tm;
    std::tr1::unordered_map<string, unsigned int> variantMap; //maps between variants and group indices

    string includeIons = "all";

    //set up initial groups
    for (int i = 0; i < variantSeq.size(); i++)
    {
      vector < string > currGroup;
      currGroup.push_back(variantSeq[i]);
      variantMap[variantSeq[i]] = i;
      variantQuantities.push_back(quantities[i]);
      minDistinguishingIntensity.push_back(1.0);
      variantGroups.push_back(currGroup);
    }

    int lastIndex = variantGroups.size();

    for (int i = 0; i < lastIndex; i++)
    {
      //check to see whether group has disappeared.
      if (variantGroups[i].size() == 0)
      {
        continue;
      }

      string kAnnotation = variantGroups[i][0];

      PeptideSpectrumMatch kMatch;
      kMatch.m_spectrum = &modifiedSpectrum;
      kMatch.annotate(kAnnotation, includeIons, m_model, 0, 0, m_peakTol, false);
      tm.start();

      for (int j = i + 1; j < lastIndex; j++)
      {
        //check to see whether group has disappeared.
        if (variantGroups[j].size() == 0)
        {
          continue;
        }

        string jAnnotation = variantGroups[j][0];

        float distinguishingIntensity = 0;
        PeptideSpectrumMatch jMatch;
        jMatch.m_spectrum = &modifiedSpectrum;
        jMatch.annotate(variantGroups[j][0],
                        includeIons,
                        m_model,
                        0,
                        0,
                        m_peakTol,
                        false);

        float totalIonCurrent = getTotalIonCurrent(jMatch, kMatch);
#ifdef DEBUG
        cout << kAnnotation << endl;
        cout << jAnnotation << endl;
#endif

        for (int k = 0; k < modifiedSpectrum.size(); k++)
        {
          const ftIonFragment * jFrag = jMatch.m_peakAnnotations[k].first;
          const ftIonFragment * kFrag = kMatch.m_peakAnnotations[k].first;

          if (jFrag != NULL && kFrag != NULL)
          {
            //allow for peaks which have two possible ion interpretations to be considered in total intensity
            if (jFrag->name.compare(kFrag->name) != 0
                || jMatch.m_peakAnnotations[k].second
                    != kMatch.m_peakAnnotations[k].second)
            {
#ifdef DEBUG
              if (jFrag != NULL)
              cout << jFrag->name << jMatch.m_peakAnnotations[k].second
              << endl;

              if (kFrag != NULL)
              cout << kFrag->name << kMatch.m_peakAnnotations[k].second
              << endl;

              cout << "modified spectrum " << modifiedSpectrum[k][1] << endl;
#endif
              distinguishingIntensity += modifiedSpectrum[k][1];
            }
          }
          else if (jFrag != NULL || kFrag != NULL)
          {
#ifdef DEBUG
            if (jFrag != NULL)
            cout << jFrag->name << jMatch.m_peakAnnotations[k].second << endl;

            if (kFrag != NULL)
            cout << kFrag->name << kMatch.m_peakAnnotations[k].second << endl;

            cout << "modified spectrum " << modifiedSpectrum[k][1] << endl;
#endif
            distinguishingIntensity += modifiedSpectrum[k][1];
          }
        }
#ifdef DEBUG
        cout << "distinguishing intensity " << distinguishingIntensity << endl;
        cout << " totalIonCurrent " << totalIonCurrent << endl;
        cout << distinguishingIntensity / totalIonCurrent << endl;
#endif
        std::tr1::unordered_map<string, unsigned int>::iterator i_map =
            variantMap.find(kAnnotation);
        std::tr1::unordered_map<string, unsigned int>::iterator j_map =
            variantMap.find(jAnnotation);

        float currDistingishingIntensity = distinguishingIntensity
            / totalIonCurrent;

        if (currDistingishingIntensity < minimumGroupingPercentIntensity)
        {
          if (i_map != variantMap.end() && j_map != variantMap.end()) // both groups already in a group
          {
            if (i_map->second == j_map->second)
            {
              // do nothing
            }
            else
            {
              //merge groups
              vector < string > newGroup;
              for (int jIdx = 0; jIdx < variantGroups[j_map->second].size(); jIdx++)
              {
                newGroup.push_back(variantGroups[j_map->second][jIdx]);
              }
              for (int iIdx = 0; iIdx < variantGroups[i_map->second].size(); iIdx++)
              {
                newGroup.push_back(variantGroups[i_map->second][iIdx]);
              }
              //clear old variant vectors
              variantGroups[j_map->second].resize(0);
              variantGroups[i_map->second].resize(0);
              variantGroups.push_back(newGroup);

              //set new quantity
              float quantity = variantQuantities[i_map->second]
                  + variantQuantities[j_map->second];
              variantQuantities[i_map->second] = -1;
              variantQuantities[j_map->second] = -1;
              variantQuantities.push_back(quantity);

              //set new distinguishing intensity
              minDistinguishingIntensity[i_map->second] = 1.0;
              minDistinguishingIntensity[j_map->second] = 1.0;
              minDistinguishingIntensity.push_back(1.0);

              //map new location
              for (int gIdx = 0; gIdx < newGroup.size(); gIdx++)
              {
                variantMap[newGroup[gIdx]] = variantGroups.size() - 1;
              }
            }
          }
          else if (i_map != variantMap.end())
          {
            variantGroups[i_map->second].push_back(jAnnotation);
            variantQuantities[i_map->second] += quantities[j];

            if (minDistinguishingIntensity[i_map->second]
                > currDistingishingIntensity)
            {
              minDistinguishingIntensity[i_map->second]
                  = currDistingishingIntensity;
            }

            variantMap[jAnnotation] = variantMap[kAnnotation];
          }
          else if (j_map != variantMap.end())
          {
            variantGroups[j_map->second].push_back(kAnnotation);
            variantQuantities[j_map->second] += quantities[i];
            if (minDistinguishingIntensity[j_map->second]
                > currDistingishingIntensity)
            {
              minDistinguishingIntensity[j_map->second]
                  = currDistingishingIntensity;
            }

            variantMap[kAnnotation] = variantMap[jAnnotation];
          }
          else
          {
            vector < string > group;
            group.push_back(kAnnotation);
            group.push_back(jAnnotation);
            float quantity = quantities[i] + quantities[j];
            variantGroups.push_back(group);
            variantQuantities.push_back(quantity);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);

            variantMap[kAnnotation] = variantGroups.size() - 1;
            variantMap[jAnnotation] = variantGroups.size() - 1;
          }
        }
        else
        {
          if (i_map == variantMap.end())
          {
            vector < string > group1;
            group1.push_back(kAnnotation);
            variantGroups.push_back(group1);
            variantQuantities.push_back(quantities[i]);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);
            variantMap[kAnnotation] = variantGroups.size() - 1;
          }

          if (j_map == variantMap.end())
          {
            vector < string > group2;
            group2.push_back(jAnnotation);
            variantGroups.push_back(group2);
            variantQuantities.push_back(quantities[j]);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);
            variantMap[jAnnotation] = variantGroups.size() - 1;
          }
        }
      }
      tm.stop();
      //  DEBUG_MSG("Iteration " << i << " duration: " << tm.duration());
    }
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::groupVariants(Spectrum &modifiedSpectrum,
                                             Spectrum &allMassesSpectrum,
                                             vector<string> &variantSeq,
                                             const vector<float> &quantities,
                                             const vector<float> &massShifts,
                                             vector<vector<string> > &variantGroups,
                                             vector<float> &variantQuantities,
                                             vector<float> &minDistinguishingIntensity,
                                             float minimumGroupingPercentIntensity)
  {
    Timer tm;
    std::tr1::unordered_map<string, unsigned int> variantMap; //maps between variants and group indices

    string includeIons = "all";
    float totalIonCurrent = getTotalIonCurrent(modifiedSpectrum,
                                               allMassesSpectrum);

    for (int i = 0; i < variantSeq.size(); i++)
    {
      PeptideSpectrumMatch kMatch;
      kMatch.m_spectrum = &modifiedSpectrum;
      kMatch.annotate(variantSeq[i],
                      includeIons,
                      m_model,
                      0,
                      0,
                      m_peakTol,
                      false);

      tm.start();
      for (int j = i + 1; j < variantSeq.size(); j++)
      {
        float distinguishingIntensity = 0;
        PeptideSpectrumMatch jMatch;
        jMatch.m_spectrum = &modifiedSpectrum;
        jMatch.annotate(variantSeq[j],
                        includeIons,
                        m_model,
                        0,
                        0,
                        m_peakTol,
                        false);

#ifdef DEBUG
        cout << variantSeq[i] << endl;
        cout << variantSeq[j] << endl;
#endif

        for (int k = 0; k < modifiedSpectrum.size(); k++)
        {
          const ftIonFragment * jFrag = jMatch.m_peakAnnotations[k].first;
          const ftIonFragment * kFrag = kMatch.m_peakAnnotations[k].first;

          if (jFrag != NULL && kFrag != NULL)
          {
            //allow for peaks which have two possible ion interpretations to be considered in total intensity
            if (jFrag->name.compare(kFrag->name) != 0
                || jMatch.m_peakAnnotations[k].second
                    != kMatch.m_peakAnnotations[k].second)
            {
#ifdef DEBUG
              if (jFrag != NULL)
              cout << jFrag->name << jMatch.m_peakAnnotations[k].second
              << endl;

              if (kFrag != NULL)
              cout << kFrag->name << kMatch.m_peakAnnotations[k].second
              << endl;

              cout << "modified spectrum " << modifiedSpectrum[k][1] << endl;
#endif
              distinguishingIntensity += modifiedSpectrum[k][1];
            }
          }
          else if (jFrag != NULL || kFrag != NULL)
          {
#ifdef DEBUG
            if (jFrag != NULL)
            cout << jFrag->name << jMatch.m_peakAnnotations[k].second << endl;

            if (kFrag != NULL)
            cout << kFrag->name << kMatch.m_peakAnnotations[k].second << endl;

            cout << "modified spectrum " << modifiedSpectrum[k][1] << endl;
#endif
            distinguishingIntensity += modifiedSpectrum[k][1];
          }
        }
#ifdef DEBUG
        cout << "distinguishing intensity " << distinguishingIntensity << endl;
        cout << " totalIonCurrent " << totalIonCurrent << endl;
        cout << distinguishingIntensity / totalIonCurrent << endl;
#endif
        std::tr1::unordered_map<string, unsigned int>::iterator i_map =
            variantMap.find(variantSeq[i]);
        std::tr1::unordered_map<string, unsigned int>::iterator j_map =
            variantMap.find(variantSeq[j]);

        float currDistingishingIntensity = distinguishingIntensity
            / totalIonCurrent;

        if (currDistingishingIntensity < minimumGroupingPercentIntensity)
        {
          DEBUG_TRACE;
          if (i_map != variantMap.end() && j_map != variantMap.end()) // both groups already in a group
          {
            if (i_map->second == j_map->second)
            {
              // do nothing
            }
            else
            {
              //merge groups
              vector < string > newGroup;
              for (int jIdx = 0; jIdx < variantGroups[j_map->second].size(); jIdx++)
              {
                newGroup.push_back(variantGroups[j_map->second][jIdx]);
              }
              for (int iIdx = 0; iIdx < variantGroups[i_map->second].size(); iIdx++)
              {
                newGroup.push_back(variantGroups[i_map->second][iIdx]);
              }
              //clear old variant vectors
              variantGroups[j_map->second].resize(0);
              variantGroups[i_map->second].resize(0);
              variantGroups.push_back(newGroup);

              //set new quantity
              float quantity = variantQuantities[i_map->second]
                  + variantQuantities[j_map->second];
              variantQuantities[i_map->second] = -1;
              variantQuantities[j_map->second] = -1;
              variantQuantities.push_back(quantity);

              //set new distinguishing intensity
              minDistinguishingIntensity[i_map->second] = 1.0;
              minDistinguishingIntensity[j_map->second] = 1.0;
              minDistinguishingIntensity.push_back(1.0);

              //map new location
              for (int gIdx = 0; gIdx < newGroup.size(); gIdx++)
              {
                variantMap[newGroup[gIdx]] = variantGroups.size() - 1;
              }
            }
          }
          else if (i_map != variantMap.end())
          {
            variantGroups[i_map->second].push_back(variantSeq[j]);
            variantQuantities[i_map->second] += quantities[j];

            if (minDistinguishingIntensity[i_map->second]
                > currDistingishingIntensity)
            {
              minDistinguishingIntensity[i_map->second]
                  = currDistingishingIntensity;
            }

            variantMap[variantSeq[j]] = variantMap[variantSeq[i]];
          }
          else if (j_map != variantMap.end())
          {
            variantGroups[j_map->second].push_back(variantSeq[i]);
            variantQuantities[j_map->second] += quantities[i];
            if (minDistinguishingIntensity[j_map->second]
                > currDistingishingIntensity)
            {
              minDistinguishingIntensity[j_map->second]
                  = currDistingishingIntensity;
            }

            variantMap[variantSeq[i]] = variantMap[variantSeq[j]];
          }
          else
          {
            vector < string > group;
            group.push_back(variantSeq[i]);
            group.push_back(variantSeq[j]);
            float quantity = quantities[i] + quantities[j];
            variantGroups.push_back(group);
            variantQuantities.push_back(quantity);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);

            variantMap[variantSeq[i]] = variantGroups.size() - 1;
            variantMap[variantSeq[j]] = variantGroups.size() - 1;
          }
        }
        else
        {
          if (i_map == variantMap.end())
          {
            vector < string > group1;
            group1.push_back(variantSeq[i]);
            variantGroups.push_back(group1);
            variantQuantities.push_back(quantities[i]);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);
            variantMap[variantSeq[i]] = variantGroups.size() - 1;
          }

          if (j_map == variantMap.end())
          {
            vector < string > group2;
            group2.push_back(variantSeq[j]);
            variantGroups.push_back(group2);
            variantQuantities.push_back(quantities[j]);
            minDistinguishingIntensity.push_back(currDistingishingIntensity);
            variantMap[variantSeq[j]] = variantGroups.size() - 1;
          }
        }
        cout << distinguishingIntensity << "\t";
      }
      cout << endl;
      tm.stop();
      DEBUG_MSG("Iteration " << i << " duration:" << tm.duration());
    }
    for (int i = 0; i < variantSeq.size(); i++)
    {
      for (int j = i + 1; j < variantSeq.size(); j++)
      {
        cout << variantSeq[i] << "\t" << i;
        cout << variantSeq[j] << "\t" << j << endl;
      }
    }
  }

  // -------------------------------------------------------------------------

  bool compare(const ftIonFragment &i, const ftIonFragment &j)
  {
    return (i.prob < j.prob);
  }

  // -------------------------------------------------------------------------
  float totalModMass(vector<float> &modMass, vector<int> &modIndices)
  {
    float totalMass = 0;
    for (int i = 0; i < modIndices.size(); i++)
    {
      totalMass += modMass[modIndices[i]];
    }
    return totalMass;
  }
  // -------------------------------------------------------------------------
  // Helper function for generateFlrLP
  void FalseLocalizationRates::generateAllIonMasses(const string &peptide,
                                                    const vector<float> &modifications,
                                                    vector<float> &outputMasses)
  {
    vector<float> newMods = modifications;
    outputMasses.resize(0);

    //generate mod combos.
    vector < vector<int> > modCombos;

    for (int k = 1; k <= newMods.size(); k++)
    {
      MathUtils::combinations(newMods.size(), k, modCombos);
    }

    //generate srm and prm masses, no offsets.
    vector<float> prmMasses;
    vector<float> srmMasses;

    //generate total peptide mass (summed value of aa masses)
    float peptideMass;

    string unmodPeptide;
    PeptideSpectrumMatch::getUnmodifiedPeptide(peptide, unmodPeptide);
    DEBUG_VAR(unmodPeptide);

    m_jumps.getPRMandSRMMasses(unmodPeptide, prmMasses, srmMasses, peptideMass);

    short peptideLength = prmMasses.size(); //length of peptide

    vector<ftIonFragment> ionTypes;
    //copy fragments from m_model.
    ionTypes.resize(m_model.probs.size());
    copy(m_model.probs.begin(), m_model.probs.end(), ionTypes.begin());

    sort(ionTypes.begin(), ionTypes.end(), compare); //sort fragment types in ascending order so that low probability
    //annotations will be overwritten.

    short ionIdx;

    for (ionIdx = 0; ionIdx < peptideLength - 1; ionIdx++)
    {
      vector<ftIonFragment>::const_iterator currIonFrag;

      for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
      {
        if (currIonFrag->isIF)
        {
          for (int i = 0; i < modCombos.size(); i++)
          {
            //check to make sure ion length is enough
            if (ionIdx + 1 >= modCombos[i].size())
            {
              float currMass;

              if (currIonFrag->isNTerm)
              {
                float totalMass = totalModMass(newMods, modCombos[i]);

                currMass = (prmMasses[ionIdx] + currIonFrag->massOffset
                    + totalMass) / currIonFrag->charge;
              }
              else
              {
                float totalMass = totalModMass(newMods, modCombos[i]);

                currMass = (srmMasses[ionIdx] + currIonFrag->massOffset
                    + totalMass) / currIonFrag->charge;
              }
              outputMasses.push_back(currMass);
            }
          }
          //include unmodified ion
          float currMass;
          if (currIonFrag->isNTerm)
          {
            currMass = (prmMasses[ionIdx] + currIonFrag->massOffset)
                / currIonFrag->charge;
          }
          else
          {
            currMass = (srmMasses[ionIdx] + currIonFrag->massOffset)
                / currIonFrag->charge;
          }
          outputMasses.push_back(currMass);
        }
      }
    }
    sort(outputMasses.begin(), outputMasses.end());
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateAllIonMasses(const string &peptide,
                                                    const vector<vector<float> > &massShifts,
                                                    vector<float> &outputMasses)
  {
    vector<float> newMods;
    for (int i = 0; i < massShifts.size(); i++)
    {
      for (int j = 0; j < massShifts[i].size(); j++)
      {
        newMods.push_back(massShifts[i][j]);
      }
    }

    outputMasses.resize(0);

    //generate mod combos.
    vector < vector<int> > modCombos;

    for (int k = 1; k <= newMods.size(); k++)
    {
      MathUtils::combinations(newMods.size(), k, modCombos);
    }

    //generate srm and prm masses, no offsets.
    vector<float> prmMasses;
    vector<float> srmMasses;

    //generate total peptide mass (summed value of aa masses)
    float peptideMass;

    string unmodPeptide;
    PeptideSpectrumMatch::getUnmodifiedPeptide(peptide, unmodPeptide);
    DEBUG_VAR(unmodPeptide);

    m_jumps.getPRMandSRMMasses(unmodPeptide, prmMasses, srmMasses, peptideMass);

    short peptideLength = prmMasses.size(); //length of peptide

    vector<ftIonFragment> ionTypes;
    //copy fragments from m_model.
    ionTypes.resize(m_model.probs.size());
    copy(m_model.probs.begin(), m_model.probs.end(), ionTypes.begin());

    sort(ionTypes.begin(), ionTypes.end(), compare); //sort fragment types in ascending order so that low probability
    //annotations will be overwritten.

    short ionIdx;

    for (ionIdx = 0; ionIdx < peptideLength - 1; ionIdx++)
    {
      vector<ftIonFragment>::const_iterator currIonFrag;

      for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
      {
        if (currIonFrag->isIF)
        {
          for (int i = 0; i < modCombos.size(); i++)
          {
            //check to make sure ion length is enough
            if (ionIdx + 1 >= modCombos[i].size())
            {
              float currMass;

              if (currIonFrag->isNTerm)
              {
                float totalMass = totalModMass(newMods, modCombos[i]);

                currMass = (prmMasses[ionIdx] + currIonFrag->massOffset
                    + totalMass) / currIonFrag->charge;
              }
              else
              {
                float totalMass = totalModMass(newMods, modCombos[i]);

                currMass = (srmMasses[ionIdx] + currIonFrag->massOffset
                    + totalMass) / currIonFrag->charge;
              }
              outputMasses.push_back(currMass);
            }
          }
          //include unmodified ion
          float currMass;
          if (currIonFrag->isNTerm)
          {
            currMass = (prmMasses[ionIdx] + currIonFrag->massOffset)
                / currIonFrag->charge;
          }
          else
          {
            currMass = (srmMasses[ionIdx] + currIonFrag->massOffset)
                / currIonFrag->charge;
          }
          outputMasses.push_back(currMass);
        }
      }
    }
    sort(outputMasses.begin(), outputMasses.end());
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateFlrLP(Spectrum &unmodifiedSpectrum,
                                             Spectrum &modifiedSpectrum,
                                             const string &unmodifiedPeptide,
                                             const vector<float> &massShifts,
                                             const char * outputCplexFile,
                                             const char * outputAllPeaksSpectrum,
                                             const char * outputVariantMaps)
  {

    //Generate expected Intensities annotations;
    PeptideSpectrumMatch expectedIntensitiesPSM;
    expectedIntensitiesPSM.m_spectrum = &unmodifiedSpectrum;
    string includeIons = "all";
    expectedIntensitiesPSM.annotate(unmodifiedPeptide,
                                    includeIons,
                                    m_model,
                                    0,
                                    0,
                                    m_peakTol,
                                    false);

    std::tr1::unordered_map<string, int> unmodifiedSpectrumMap;
    vector < vector<int> > unmodifiedSpectrumIons;
    expectedIntensitiesPSM.mapIons(unmodifiedSpectrumIons,
                                   unmodifiedSpectrumMap);

    //generate variant sequences;
    vector < string > variants;
    DEBUG_VAR(massShifts.size());

    //make sure variants are from unmodified sequence (ignore fixed mods or SILAC mods)
    string unmodifiedSequence;
    PeptideSpectrumMatch::getUnmodifiedPeptide(unmodifiedPeptide,
                                               unmodifiedSequence);

    generateVariantSequences(unmodifiedSequence, massShifts, variants);
    DEBUG_VAR(variants.size());

    vector<float> peakMasses;

    generateAllIonMasses(unmodifiedPeptide, massShifts, peakMasses);

    Spectrum allPossiblePeaks;

    for (int i = 0; i < peakMasses.size(); i++)
    {
      TwoValues<float> currPeak;
      currPeak[0] = peakMasses[i];
      currPeak[1] = 0;
      allPossiblePeaks.push_back(currPeak);
    }

    vector < vector<float> > variantLHSConstraintMatrix;
    vector<float> variantRHSConstraints;

    variantLHSConstraintMatrix.resize(peakMasses.size());

    //initialize the matrix
    for (int i = 0; i < variantLHSConstraintMatrix.size(); i++)
    {
      variantLHSConstraintMatrix[i].resize(variants.size());
    }

    variantRHSConstraints.resize(peakMasses.size(), 0);

    int maxIndex = 0;
    int ionIndex = -1;
    float prevPeak = 0;
    vector<int> ionIndices(peakMasses.size());

    for (int i = 0; i < peakMasses.size(); i++)
    {
      if (peakMasses[i] - prevPeak >= m_peakTol)
      {
        ionIndex++;
      }
      ionIndices[i] = ionIndex;
      prevPeak = peakMasses[i];
    }
    maxIndex = ionIndex;

    for (int i = 0; i < variants.size(); i++)
    {
#ifdef DEBUG
      cout << "variant " << variants[i] << endl;
#endif
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allPossiblePeaks;

      string includeIons = "all";
      psm.annotate(variants[i], includeIons, m_model, 0, 0, m_peakTol, false);

      for (int j = 0; j < psm.m_peakAnnotations.size(); j++)
      {
        if (psm.m_peakAnnotations[j].first != NULL)
        {
          //initialize output variables.
          float expectedIntensity = 0;
          int startIndex = 0;

          //pull expected intensity
          const ftIonFragment * currIonFrag = psm.m_peakAnnotations[j].first;

          stringstream ss;
          ss << currIonFrag->name << psm.m_peakAnnotations[j].second;
          string ionName = ss.str();
#ifdef DEBUG
          cout << "ionName " << ionName << endl;
          cout << "mass " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (unmodifiedSpectrumMap.find(ionName)
              != unmodifiedSpectrumMap.end())
          {
            for (int k = 0; k
                < unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]].size(); k++)
            {

              int currPeakIndex =
                  unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]][k];
              expectedIntensity += unmodifiedSpectrum[currPeakIndex][1];
            }
          }
          else
          {
            expectedIntensity = 1;
          }
#ifdef DEBUG
          cout << "expectedIntensity " << expectedIntensity << endl;
#endif

          variantLHSConstraintMatrix[ionIndices[j]][i] = expectedIntensity;

          //get right hand side constraint
          vector<int> outputIndices;
#ifdef DEBUG
          cout << "allMasses " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (modifiedSpectrum.findMatches((*psm.m_spectrum)[j][0],
                                           m_peakTol,
                                           outputIndices,
                                           startIndex))
          {
            if (variantRHSConstraints[ionIndices[j]] == 0)
            {
              for (int matchIndex = 0; matchIndex < outputIndices.size(); matchIndex++)
              {
#ifdef DEBUG
                cout << "mod peak "
                << modifiedSpectrum[outputIndices[matchIndex]][0] << endl;
#endif

                variantRHSConstraints[ionIndices[j]]
                    += modifiedSpectrum[outputIndices[matchIndex]][1];
                startIndex = outputIndices[matchIndex];
              }
            }
#ifdef DEBUG
            cout << "variantRHSConstraints[j]"
            << variantRHSConstraints[ionIndices[j]] << endl;
#endif
          }
        }
      }
    }

    variantRHSConstraints.resize(maxIndex);
    variantLHSConstraintMatrix.resize(maxIndex);

    ofstream outputVariantHandle(outputVariantMaps, ios::binary);
    if (!outputVariantHandle.is_open() || !outputVariantHandle.good())
    {
      ERROR_MSG("Unable to open file " << outputVariantMaps);
      return;
    }

    DEBUG_TRACE;

    outputVariantHandle << "variantName\tlpVariantName" << endl;

    for (int i = 0; i < variants.size(); i++)
    {
      outputVariantHandle << variants[i] << "\t";
      stringstream ss;
      ss << setw(3) << setfill('0');
      ss << i + 1;

      string lpVariant = "Qv";
      lpVariant.append(ss.str());
      outputVariantHandle << lpVariant << endl;
    }

    LinearEquation lp;
    lp.m_lhsConstraints = variantLHSConstraintMatrix;
    lp.m_rhsConstraints = variantRHSConstraints;
    lp.saveCPLEXLP(outputCplexFile);
    SpecSet specs;
    specs.resize(1);
    specs[0] = allPossiblePeaks;
    specs.SaveSpecSet_mgf(outputAllPeaksSpectrum);
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateFlrLP(Spectrum &unmodifiedSpectrum,
                                             Spectrum &modifiedSpectrum,
                                             const string &unmodifiedPeptide,
                                             const vector<float> &modifications,
                                             const vector<vector<float> > &massShifts,
                                             const char * outputCplexFile,
                                             const char * outputAllPeaksSpectrum,
                                             const char * outputVariantMaps)
  {

    //Generate expected Intensities annotations;
    PeptideSpectrumMatch expectedIntensitiesPSM;
    expectedIntensitiesPSM.m_spectrum = &unmodifiedSpectrum;
    string includeIons = "all";
    expectedIntensitiesPSM.annotate(unmodifiedPeptide,
                                    includeIons,
                                    m_model,
                                    0,
                                    0,
                                    m_peakTol,
                                    false);

    std::tr1::unordered_map<string, int> unmodifiedSpectrumMap;
    vector < vector<int> > unmodifiedSpectrumIons;
    expectedIntensitiesPSM.mapIons(unmodifiedSpectrumIons,
                                   unmodifiedSpectrumMap);

    //generate variant sequences;
    vector < vector<string> > variants;
    DEBUG_VAR(massShifts.size());

    //make sure variants are from unmodified sequence (ignore fixed mods or SILAC mods)
    string unmodifiedSequence;
    PeptideSpectrumMatch::getUnmodifiedPeptide(unmodifiedPeptide,
                                               unmodifiedSequence);

    generateVariantSequences(unmodifiedSequence, massShifts, variants);
    DEBUG_VAR(variants.size());

    vector<float> peakMasses;

    generateAllIonMasses(unmodifiedPeptide, massShifts, peakMasses);

    Spectrum allPossiblePeaks;

    for (int i = 0; i < peakMasses.size(); i++)
    {
      TwoValues<float> currPeak;
      currPeak[0] = peakMasses[i];
      currPeak[1] = 0;
      allPossiblePeaks.push_back(currPeak);
    }

    vector < vector<float> > variantLHSConstraintMatrix;
    vector<float> variantRHSConstraints;

    variantLHSConstraintMatrix.resize(peakMasses.size());

    //initialize the matrix
    for (int i = 0; i < variantLHSConstraintMatrix.size(); i++)
    {
      variantLHSConstraintMatrix[i].resize(variants.size());
    }

    variantRHSConstraints.resize(peakMasses.size(), 0);

    int maxIndex = 0;
    int ionIndex = -1;
    float prevPeak = 0;
    vector<int> ionIndices(peakMasses.size());

    for (int i = 0; i < peakMasses.size(); i++)
    {
      if (peakMasses[i] - prevPeak >= m_peakTol)
      {
        ionIndex++;
      }
      ionIndices[i] = ionIndex;
      prevPeak = peakMasses[i];
    }
    maxIndex = ionIndex;

    for (int i = 0; i < variants.size(); i++)
    {
#ifdef DEBUG
      cout << "variant " << variants[i][0] << endl;
#endif
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allPossiblePeaks;

      string includeIons = "all";
      for (int j = 0; j < variants[i].size(); j++)
      {
        psm.annotate(variants[i][j],
                     includeIons,
                     m_model,
                     0,
                     0,
                     m_peakTol,
                     false,
                     false,
                     true);
      }

      for (int j = 0; j < psm.m_peakAnnotations.size(); j++)
      {
        if (psm.m_peakAnnotations[j].first != NULL)
        {
          //initialize output variables.
          float expectedIntensity = 0;
          int startIndex = 0;

          //pull expected intensity
          const ftIonFragment * currIonFrag = psm.m_peakAnnotations[j].first;

          stringstream ss;
          ss << currIonFrag->name << psm.m_peakAnnotations[j].second;
          string ionName = ss.str();
#ifdef DEBUG
          cout << "ionName " << ionName << endl;
          cout << "mass " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (unmodifiedSpectrumMap.find(ionName)
              != unmodifiedSpectrumMap.end())
          {
            for (int k = 0; k
                < unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]].size(); k++)
            {

              int currPeakIndex =
                  unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]][k];
              expectedIntensity += unmodifiedSpectrum[currPeakIndex][1];
            }
          }
          else
          {
            expectedIntensity = 1;
          }
#ifdef DEBUG
          cout << "expectedIntensity " << expectedIntensity << endl;
#endif

          variantLHSConstraintMatrix[ionIndices[j]][i] = expectedIntensity;

          //get right hand side constraint
          vector<int> outputIndices;
#ifdef DEBUG
          cout << "allMasses " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (modifiedSpectrum.findMatches((*psm.m_spectrum)[j][0],
                                           m_peakTol,
                                           outputIndices,
                                           startIndex))
          {
            if (variantRHSConstraints[ionIndices[j]] == 0)
            {
              for (int matchIndex = 0; matchIndex < outputIndices.size(); matchIndex++)
              {
#ifdef DEBUG
                cout << "mod peak "
                << modifiedSpectrum[outputIndices[matchIndex]][0] << endl;
#endif

                variantRHSConstraints[ionIndices[j]]
                    += modifiedSpectrum[outputIndices[matchIndex]][1];
                startIndex = outputIndices[matchIndex];
              }
            }
#ifdef DEBUG
            cout << "variantRHSConstraints[j]"
            << variantRHSConstraints[ionIndices[j]] << endl;
#endif
          }
        }
      }
    }

    variantRHSConstraints.resize(maxIndex);
    variantLHSConstraintMatrix.resize(maxIndex);

    ofstream outputVariantHandle(outputVariantMaps, ios::binary);
    if (!outputVariantHandle.is_open() || !outputVariantHandle.good())
    {
      ERROR_MSG("Unable to open file " << outputVariantMaps);
      return;
    }

    DEBUG_TRACE;

    outputVariantHandle << "variantName\tlpVariantName" << endl;

    for (int i = 0; i < variants.size(); i++)
    {
      string currVariant;
      stringJoin(currVariant, variants[i], ":");
      outputVariantHandle << currVariant << "\t";
      stringstream ss;
      ss << setw(3) << setfill('0');
      ss << i + 1;

      string lpVariant = "Qv";
      lpVariant.append(ss.str());
      outputVariantHandle << lpVariant << endl;
    }

    LinearEquation lp;
    lp.m_lhsConstraints = variantLHSConstraintMatrix;
    lp.m_rhsConstraints = variantRHSConstraints;
    lp.saveCPLEXLP(outputCplexFile);
    SpecSet specs;
    specs.resize(1);
    specs[0] = allPossiblePeaks;
    specs.SaveSpecSet_mgf(outputAllPeaksSpectrum);
  }
  // -------------------------------------------------------------------------
  void FalseLocalizationRates::generateFlrLPPhospho(Spectrum &unmodifiedSpectrum,
                                             Spectrum &modifiedSpectrum,
                                             const string &unmodifiedPeptide,
                                             const vector<float> &modifications,
                                             const vector<vector<float> > &massShifts,
                                             const char * outputCplexFile,
                                             const char * outputAllPeaksSpectrum,
                                             const char * outputVariantMaps)
  {

    //Generate expected Intensities annotations;
    PeptideSpectrumMatch expectedIntensitiesPSM;
    expectedIntensitiesPSM.m_spectrum = &unmodifiedSpectrum;
    string includeIons = "all";
    expectedIntensitiesPSM.annotate(unmodifiedPeptide,
                                    includeIons,
                                    m_model,
                                    0,
                                    0,
                                    m_peakTol,
                                    false);

    std::tr1::unordered_map<string, int> unmodifiedSpectrumMap;
    vector < vector<int> > unmodifiedSpectrumIons;
    expectedIntensitiesPSM.mapIons(unmodifiedSpectrumIons,
                                   unmodifiedSpectrumMap);

    //generate variant sequences;
    vector < vector<string> > variants;
    DEBUG_VAR(massShifts.size());

    //make sure variants are from unmodified sequence (ignore fixed mods or SILAC mods)
    string unmodifiedSequence;
    PeptideSpectrumMatch::getUnmodifiedPeptide(unmodifiedPeptide,
                                               unmodifiedSequence);

    generateVariantSequences(unmodifiedSequence, massShifts, variants);
    DEBUG_VAR(variants.size());

    vector<float> peakMasses;

    generateAllIonMasses(unmodifiedPeptide, massShifts, peakMasses);

    Spectrum allPossiblePeaks;

    for (int i = 0; i < peakMasses.size(); i++)
    {
      TwoValues<float> currPeak;
      currPeak[0] = peakMasses[i];
      currPeak[1] = 0;
      allPossiblePeaks.push_back(currPeak);
    }

    vector < vector<float> > variantLHSConstraintMatrix;
    vector<float> variantRHSConstraints;

    variantLHSConstraintMatrix.resize(peakMasses.size());

    //initialize the matrix
    for (int i = 0; i < variantLHSConstraintMatrix.size(); i++)
    {
      variantLHSConstraintMatrix[i].resize(variants.size());
    }

    variantRHSConstraints.resize(peakMasses.size(), 0);

    int maxIndex = 0;
    int ionIndex = -1;
    float prevPeak = 0;
    vector<int> ionIndices(peakMasses.size());

    for (int i = 0; i < peakMasses.size(); i++)
    {
      if (peakMasses[i] - prevPeak >= m_peakTol)
      {
        ionIndex++;
      }
      ionIndices[i] = ionIndex;
      prevPeak = peakMasses[i];
    }
    maxIndex = ionIndex;

    for (int i = 0; i < variants.size(); i++)
    {
#ifdef DEBUG
      cout << "variant " << variants[i][0] << endl;
#endif
      PeptideSpectrumMatch psm;
      psm.m_spectrum = &allPossiblePeaks;

      string includeIons = "all";

      vector<pair<const ftIonFragment*, short> > nonFlattenedPeakAnnotations;

      for (int j = 0; j < variants[i].size(); j++)
      {
        psm.annotate(variants[i][j],
                     includeIons,
                     m_model,
                     0,
                     0,
                     m_peakTol,
                     false,
                     false,
                     true);
        if (j ==0)
          nonFlattenedPeakAnnotations = psm.m_peakAnnotations;
      }

      for (int j = 0; j < psm.m_peakAnnotations.size(); j++)
      {
        if (psm.m_peakAnnotations[j].first != NULL)
        {
          //initialize output variables.
          float expectedIntensity = 0;
          int startIndex = 0;

          //pull expected intensity
          const ftIonFragment * currIonFrag = psm.m_peakAnnotations[j].first;

          stringstream ss;
          ss << currIonFrag->name << psm.m_peakAnnotations[j].second;
          string ionName = ss.str();
#ifdef DEBUG
          cout << "ionName " << ionName << endl;
          cout << "mass " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (unmodifiedSpectrumMap.find(ionName)
              != unmodifiedSpectrumMap.end() && nonFlattenedPeakAnnotations[j].first != NULL)
          {
            for (int k = 0; k
                < unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]].size(); k++)
            {

              int currPeakIndex =
                  unmodifiedSpectrumIons[unmodifiedSpectrumMap[ionName]][k];
              expectedIntensity += unmodifiedSpectrum[currPeakIndex][1];
            }
          }
          else
          {
            expectedIntensity = 1;
          }
#ifdef DEBUG
          cout << "expectedIntensity " << expectedIntensity << endl;
#endif

          variantLHSConstraintMatrix[ionIndices[j]][i] = expectedIntensity;

          //get right hand side constraint
          vector<int> outputIndices;
#ifdef DEBUG
          cout << "allMasses " << (*psm.m_spectrum)[j][0] << endl;
#endif

          if (modifiedSpectrum.findMatches((*psm.m_spectrum)[j][0],
                                           m_peakTol,
                                           outputIndices,
                                           startIndex))
          {
            if (variantRHSConstraints[ionIndices[j]] == 0)
            {
              for (int matchIndex = 0; matchIndex < outputIndices.size(); matchIndex++)
              {
#ifdef DEBUG
                cout << "mod peak "
                << modifiedSpectrum[outputIndices[matchIndex]][0] << endl;
#endif

                variantRHSConstraints[ionIndices[j]]
                    += modifiedSpectrum[outputIndices[matchIndex]][1];
                startIndex = outputIndices[matchIndex];
              }
            }
#ifdef DEBUG
            cout << "variantRHSConstraints[j]"
            << variantRHSConstraints[ionIndices[j]] << endl;
#endif
          }
        }
      }
    }

    variantRHSConstraints.resize(maxIndex);
    variantLHSConstraintMatrix.resize(maxIndex);

    ofstream outputVariantHandle(outputVariantMaps, ios::binary);
    if (!outputVariantHandle.is_open() || !outputVariantHandle.good())
    {
      ERROR_MSG("Unable to open file " << outputVariantMaps);
      return;
    }

    DEBUG_TRACE;

    outputVariantHandle << "variantName\tlpVariantName" << endl;

    for (int i = 0; i < variants.size(); i++)
    {
      string currVariant;
      stringJoin(currVariant, variants[i], ":");
      outputVariantHandle << currVariant << "\t";
      stringstream ss;
      ss << setw(3) << setfill('0');
      ss << i + 1;

      string lpVariant = "Qv";
      lpVariant.append(ss.str());
      outputVariantHandle << lpVariant << endl;
    }

    LinearEquation lp;
    lp.m_lhsConstraints = variantLHSConstraintMatrix;
    lp.m_rhsConstraints = variantRHSConstraints;
    lp.saveCPLEXLP(outputCplexFile);
    SpecSet specs;
    specs.resize(1);
    specs[0] = allPossiblePeaks;
    specs.SaveSpecSet_mgf(outputAllPeaksSpectrum);
  }
}
