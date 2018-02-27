#include "FdrPeptide.h"

namespace specnets
{
  void debugType(int fdrType)
  {
    string typeString = "";
    DEBUG_MSG("fdrType = " << std::hex << fdrType);
    switch (fdrType & FDR_TYPE_MASK) {
    case 0x0:
      typeString += "FDR_TYPE_NONE ";
      break;
    case 0x1:
      typeString += "FDR_TYPE_SPECTRUM ";
      break;
    case 0x2:
      typeString += "FDR_TYPE_PEPTIDE ";
      break;
    case 0x3:
      typeString += "FDR_TYPE_VARIANT ";
      break;
    case 0x4:
      typeString += "FDR_TYPE_PROTEIN ";
      break;
    default:
      typeString += "FDR_TYPE_UNK ";
      break;
    }
    switch (fdrType & FDR_REPLACE_MASK) {
    case 0x10:
      typeString += "FDR_REPLACE_SCORE ";
      break;
    case 0x20:
      typeString += "FDR_REPLACE_PVALUE ";
      break;
    default:
      typeString += "FDR_REPLACE_UNK ";
      break;
    }
    switch (fdrType & FDR_SORT_MASK) {
    case 0x100:
      typeString += "FDR_SORT_SCORE";
      break;
    case 0x200:
      typeString += "FDR_SORT_PVALUE";
      break;
    default:
      typeString += "FDR_SORT_UNK";
      break;
    }
    DEBUG_MSG("fdrType = " << std::hex << fdrType << "  " << typeString);
  }

  //---------------------------------------------------------------------------
  // helper function for calculatePValues, adds psm with maximum m_score to map.
  // -------------------------------------------------------------------------
  bool addToMaxMap(PeptideSpectrumMatchSet &inputSet,
                   std::tr1::unordered_map<string, psmPtr> &maxPsmMap,
                   int fdrType)
  {
    std::tr1::unordered_map<string, psmPtr>::const_iterator it;

    debugType(fdrType);
    int fdrTypeMasked = fdrType & FDR_TYPE_MASK;
    int fdrReplaceMasked = fdrType & FDR_REPLACE_MASK;
    DEBUG_MSG("fdrTypeMasked = " << std::hex << fdrTypeMasked);
    DEBUG_MSG("fdrReplaceMasked = " << std::hex << fdrReplaceMasked);

    for (unsigned int i = 0; i < inputSet.size(); i++)
    {
      psmPtr currMatch = inputSet[i];
      PeptideSpectrumMatch psm;

      stringstream ss;
      string key;
      if (fdrTypeMasked == FDR_TYPE_SPECTRUM) {
        ss << currMatch->m_scanNum << '_' << currMatch->m_spectrumFile;
        key = ss.str();
      } else if (fdrTypeMasked == FDR_TYPE_PEPTIDE) {
        string cleanAnnotation;
        currMatch->getUnmodifiedPeptide(cleanAnnotation);
        ss << cleanAnnotation;
        key = ss.str();
      } else if (fdrTypeMasked == FDR_TYPE_VARIANT) {
        ss << currMatch->m_variantGroup;
        key = ss.str();
      } else if (fdrTypeMasked == FDR_TYPE_PROTEIN) {
        ss << currMatch->m_protein;
        key = ss.str();
      } else if (fdrTypeMasked == FDR_TYPE_NONE) {
        ss << i;  // unique identifier
        key = ss.str();
      } else {
        ERROR_MSG("Unknown FDR type [" << hex << fdrTypeMasked << "]");
        return false;
      }

      it = maxPsmMap.find(key);

      if (it != maxPsmMap.end())
      {
        if (fdrReplaceMasked == FDR_REPLACE_PVALUE) {
          if (maxPsmMap[key]->m_pValue > currMatch->m_pValue)
          {
            maxPsmMap[key] = currMatch;
          }
        } else if (fdrReplaceMasked == FDR_REPLACE_SCORE) {
          if (maxPsmMap[key]->m_score < currMatch->m_score)
          {
            maxPsmMap[key] = currMatch;
          }
        } else {
          ERROR_MSG("Unknown replace type [" << hex << fdrReplaceMasked << "]");
          return false;
        }
      }
      else
      {
        maxPsmMap[key] = currMatch;
      }
    }
    return true;
  }

  //---------------------------------------------------------------------------
  bool FdrPeptide::compareFDR(psmPtr i, psmPtr j)
  {
    return i->m_score > j->m_score;
  }

  //---------------------------------------------------------------------------
  bool FdrPeptide::compareFDRPVal(psmPtr i, psmPtr j)
  {
    return i->m_pValue < j->m_pValue;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::concatenatedTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                           PeptideSpectrumMatchSet &outputPeptides,
                                           double scalingFactor,
                                           int fdrType /* = FDR_TYPE_SPECTRUM */)
  {
    debugType(fdrType);
    DEBUG_VAR(inputPeptides.size());
    //filter to top hit per scan
    PeptideSpectrumMatchSet temp;
    if (!filterToTopHit(inputPeptides, temp, fdrType))
    {
      ERROR_MSG("filterToTopHit failed!");
      return false;
    }
    DEBUG_VAR(temp.size());
    if (!calculatePValues(temp,
                          outputPeptides,
                          scalingFactor,
                          fdrType))
    {
      ERROR_MSG("calculatePValues failed!");
      return false;
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::separateTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                       PeptideSpectrumMatchSet &outputPeptides,
                                       double scalingFactor,
                                       int fdrType /* = FDR_TYPE_SPECTRUM */)
  {
    debugType(fdrType);
    //separate target and decoy hits
    PeptideSpectrumMatchSet target;
    PeptideSpectrumMatchSet decoy;
    for (unsigned int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr psm = inputPeptides[i];

      if (psm->m_isDecoy)
      {
        decoy.push_back(psm);
      }
      else
      {
        target.push_back(psm);
      }
    }
    PeptideSpectrumMatchSet targetFiltered;
    PeptideSpectrumMatchSet decoyFiltered;

    if (!filterToTopHit(target, targetFiltered, fdrType))
    {
      return false;
    }
    DEBUG_VAR(targetFiltered.size());

    if (!filterToTopHit(decoy, decoyFiltered, fdrType))
    {
      return false;
    }
    DEBUG_VAR(decoyFiltered.size());

    PeptideSpectrumMatchSet concatenatedFiltered;
    concatenateTargetDecoy(targetFiltered, decoyFiltered, concatenatedFiltered);

    if (!calculatePValues(concatenatedFiltered,
                          outputPeptides,
                          scalingFactor,
                          fdrType))
    {
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::calculatePValues(PeptideSpectrumMatchSet &inputPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides,
                                    double scalingFactor,
                                    int fdrType,
                                    bool calculatePepLevel)
  {
    DEBUG_VAR(inputPeptides.size());
    outputPeptides = inputPeptides;

    debugType(fdrType);
    int fdrSortMasked = fdrType & FDR_SORT_MASK;
    DEBUG_MSG("fdrSortMasked = " << std::hex << fdrSortMasked);
    if (fdrSortMasked == FDR_SORT_SCORE) {
      sort(outputPeptides.m_psmSet.begin(),
           outputPeptides.m_psmSet.end(),
           compareFDR);
    } else if (fdrSortMasked == FDR_SORT_PVALUE) {
      sort(outputPeptides.m_psmSet.begin(),
           outputPeptides.m_psmSet.end(),
           compareFDRPVal);
    } else {
      ERROR_MSG("Unknown sort type [" << hex << (fdrType & FDR_SORT_MASK) << "]");
      return false;
    }
    unsigned int correctHits = 0;
    unsigned int incorrectHits = 0;

    tr1::unordered_set<string> correctPeptides;
    tr1::unordered_set<string> incorrectPeptides;

    for (unsigned int i = 0; i < outputPeptides.size(); i++)
    {
      psmPtr psm = outputPeptides[i];
      if (psm->m_isDecoy)
      {
        incorrectHits++;
        if (calculatePepLevel)
          incorrectPeptides.insert(psm->m_annotation);
      }
      else
      {
        correctHits++;
        if (calculatePepLevel)
          correctPeptides.insert(psm->m_annotation);
      }
      if (correctHits > 0)
      {
        double FDR = ((double)incorrectHits * (double)scalingFactor)
            / (double)correctHits;
        if (FDR > 1)
        {
          outputPeptides[i]->m_fdr = 1;
        }
        else
        {
          outputPeptides[i]->m_fdr = FDR;
        }
      }
      else
      {
        outputPeptides[i]->m_fdr = 1;
      }

      if (calculatePepLevel)
      {
        if (correctPeptides.size() > 0)
        {
          double pepLevelFDR = ((double)incorrectPeptides.size()
              * (double)scalingFactor) / ((double)correctPeptides.size());
          if (pepLevelFDR > 1)
          {
            outputPeptides[i]->m_pepFdr = 1;
          }
          else
          {
            outputPeptides[i]->m_pepFdr = pepLevelFDR;
          }
        }
        else
        {
          outputPeptides[i]->m_pepFdr = 1;
        }
      }
    }

    float min_FDR = 1.0;
    for (int i = outputPeptides.size() - 1; i >= 0; i--)
    {
      //DEBUG_MSG(i<<"\t"<<outputPeptides.size());
      min_FDR = min(outputPeptides[i]->m_fdr, min_FDR);
      outputPeptides[i]->m_fdr = min_FDR;
    }

/*
    if (incorrectHits == 0)
    {
      ERROR_MSG("No decoys set on inputPeptides!");
      return false;
    }
*/
    //D
    cout << "outputPeptides.size() = " << outputPeptides.size() << endl;
    DEBUG_MSG("outputPeptides.size() = " << outputPeptides.size() << endl);

    return true;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::mergeTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                    PeptideSpectrumMatchSet &decoyPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides,
                                    int fdrType /* = FDR_TYPE_SPECTRUM */)
  {
    debugType(fdrType);
    std::tr1::unordered_map<string, psmPtr> maxPsmMap; //map between spectrum pointer and current max psm

    //cycle through target peptides first

    if (!addToMaxMap(targetPeptides, maxPsmMap, fdrType))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    if (!addToMaxMap(decoyPeptides, maxPsmMap, fdrType))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    std::tr1::unordered_map<string, psmPtr>::const_iterator it;

    for (it = maxPsmMap.begin(); it != maxPsmMap.end(); it++)
    {
      psmPtr currPtr = it->second;
      outputPeptides.push_back(currPtr);
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::filterToTopHit(PeptideSpectrumMatchSet &inputPeptides,
                                  PeptideSpectrumMatchSet &outputPeptides,
                                  int fdrType /* = FDR_TYPE_SPECTRUM */)
  {
    debugType(fdrType);
    std::tr1::unordered_map<string, psmPtr> maxPsmMap; //map between spectrum pointer and current max psm

    if (!addToMaxMap(inputPeptides, maxPsmMap, fdrType))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    std::tr1::unordered_map<string, psmPtr>::const_iterator it;

    for (it = maxPsmMap.begin(); it != maxPsmMap.end(); it++)
    {
      psmPtr currPtr = it->second;
      outputPeptides.push_back(currPtr);
    }
    DEBUG_VAR(outputPeptides.size());
    return true;
  }

  // -------------------------------------------------------------------------
  void FdrPeptide::concatenateTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                          PeptideSpectrumMatchSet &decoyPeptides,
                                          PeptideSpectrumMatchSet &outputPeptides)
  {
    for (unsigned int i = 0; i < targetPeptides.size(); i++)
    {
      psmPtr psm = targetPeptides[i];
      outputPeptides.push_back(psm);
    }

    for (unsigned int i = 0; i < decoyPeptides.size(); i++)
    {
      psmPtr psm = decoyPeptides[i];
      outputPeptides.push_back(psm);
    }
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::filterByPValue(PeptideSpectrumMatchSet &inputPeptides,
                                  double cutoff)
  {
    // Make sure we have something to filter
    if (inputPeptides.size() == 0)
    {
      return true;
    }

    int count = 0;
    PeptideSpectrumMatchSet filteredSet;
    int max_index = 0;
    for (int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr currMatch = inputPeptides[i];

      if (currMatch->m_fdr == -1)
      {
        ERROR_MSG("PValue not set on input!");
        return false;
      }

      if (currMatch->m_fdr <= cutoff)
      {
        max_index = i;
        //filteredSet.m_psmSet.push_back(currMatch);
      }
    }

    for (int i = 0; i <= max_index; i++)
    {
      psmPtr currMatch = inputPeptides[i];
      filteredSet.m_psmSet.push_back(currMatch);
    }

    inputPeptides = filteredSet;
    return true;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::removeDecoys(PeptideSpectrumMatchSet & inputPeptides)
  {
    // Make sure we have something to filter
    if (inputPeptides.size() == 0)
    {
      return true;
    }

    int count = 0;
    PeptideSpectrumMatchSet filteredSet;
    int max_index = 0;
    for (int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr currMatch = inputPeptides[i];
      if (!currMatch->m_isDecoy) {
        filteredSet.m_psmSet.push_back(currMatch);
      }
    }

    inputPeptides = filteredSet;
    return true;
  }

}
