/*
 * SpectrumLibrary.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: jsnedecor
 */

#include "SpectrumLibrary.h"

namespace specnets
{

  /*
   * Generates key for an associated annotation;
   */
  string getKey(string &annotation, int charge)
  {
    stringstream key;
    key << annotation << '_' << charge;
    return key.str();
  }

  // -------------------------------------------------------------------------
  void SpectrumLibrary::buildKeyMap(void)
  {
    m_map.clear();

    for (unsigned int i = 0; i < m_library.size(); i++)
    {
      string key = getKey(m_library[i]->m_annotation, m_library[i]->m_charge);
      pair<std::string, unsigned int> val;
      val.first = key;
      val.second = i;
      m_map.insert(val);
    }
  }

  // -------------------------------------------------------------------------
  psmPtr & SpectrumLibrary::operator[](unsigned int i)
  {
    return m_library[i];
  }

  // -------------------------------------------------------------------------
  const psmPtr & SpectrumLibrary::operator[](unsigned int i) const
  {
    return m_library[i];
  }

  // -------------------------------------------------------------------------
  SpectrumLibrary & SpectrumLibrary::operator=(SpectrumLibrary &other)
  {
    m_library.resize(other.m_library.size());
    m_library = other.m_library;
    buildKeyMap();
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumLibrary::size()
  {
    return (unsigned int)m_library.size();
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumLibrary::resize(unsigned int newSize)
  {
    m_library.resize(newSize);
    m_libraryScores.resize(newSize);
    buildKeyMap();
    return m_library.size();
  }

  // -------------------------------------------------------------------------
  // Helper function for addToLibrary
  void SpectrumLibrary::addSpectrumMatch(Spectrum &spectrum,
                                         PeptideSpectrumMatch &psm,
                                         string &key,
                                         float score)
  {
    psmPtr currPtr(new PeptideSpectrumMatch);
    *currPtr = psm;
    m_library.m_psmSet.push_back(currPtr);
    m_libraryScores.push_back(score);
    pair<std::string, unsigned int> val;
    val.first = key;
    val.second = m_library.size() - 1;
    m_map.insert(val);
  }

  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm,
                                     float score,
                                     bool scoreAscending)
  {
    string key = getKey(psm.m_annotation, psm.m_charge);
    std::tr1::unordered_multimap<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      psmPtr libMatch = m_library[it->second];

      if (scoreAscending)
      {
        // if new score is better than old score
        if (m_libraryScores[it->second] > score)
        {
          addSpectrumMatch(spectrum, psm, key, score);
          return true;
        }
      }
      else if (m_libraryScores[it->second] < score)
      {
        // if new score is better than old score
        addSpectrumMatch(spectrum, psm, key, score);
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      addSpectrumMatch(spectrum, psm, key, score);
      return true;
    }
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm)
  {
    string key = getKey(psm.m_annotation, psm.m_charge);
    addSpectrumMatch(spectrum, psm, key, 0);
    return true;
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm,
                                     string &key)
  {
    addSpectrumMatch(spectrum, psm, key, 0);
    return true;
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm,
                                     float score,
                                     string &key,
                                     bool scoreAscending)
  {
    std::tr1::unordered_multimap<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      psmPtr libMatch = m_library[it->second];

      if (scoreAscending)
      {
        // if new score is better than old score
        if (m_libraryScores[it->second] > score)
        {
          addSpectrumMatch(spectrum, psm, key, score);
          return true;
        }
      }
      else if (m_libraryScores[it->second] < score)
      {
        // if new score is better than old score
        addSpectrumMatch(spectrum, psm, key, score);
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      addSpectrumMatch(spectrum, psm, key, score);
      return true;
    }
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::getSpectrumLibraryMatch(string &annotation,
                                                int charge,
                                                psmPtr &outputPsm)
  {
    string key = getKey(annotation, charge);
    std::tr1::unordered_multimap<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      outputPsm = m_library[it->second];
      return true;
    }
    else
    {
      return false;
    }
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::pickBestHit(psmPtr &inputPsm,
                                    vector<psmPtr> &psmsToMatch,
                                    vector<string> &ionsToExtract,
                                    MS2ScoringModel &model,
                                    psmPtr &outputPsm)
  {
    if (psmsToMatch.size() == 1)
    {
      outputPsm = psmsToMatch[0];
      return true;
    }

    if (psmsToMatch.size() == 0)
    {
      WARN_MSG("No psms to match!");
      return false;
    }

    Spectrum * inputSpectrum = inputPsm->m_spectrum;

    if (!inputSpectrum)
    {
      ERROR_MSG("No spectrum defined for input Psm!");
      return false;
    }

    int bestMatch = 0; //index of best match
    float bestCosine = 0;

    for (int i = 1; i < psmsToMatch.size(); i++)
    {
      Spectrum * spectrumToMatch = psmsToMatch[i]->m_spectrum;

      if (!spectrumToMatch)
      {
        DEBUG_MSG("Spectrum not defined for psm!");
        continue;
      }

      CosineSimilarity c_similarity(ionsToExtract);
      float currSimilarity = c_similarity.similarity(*inputSpectrum,
                                                     *spectrumToMatch,
                                                     inputPsm->m_annotation,
                                                     psmsToMatch[i]->m_annotation,
                                                     model,
                                                     false);
      if (currSimilarity > bestCosine)
      {
        bestMatch = i;
        bestCosine = currSimilarity;
      }
    }
    outputPsm = psmsToMatch[bestMatch];
    return true;
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::getSpectrumLibraryMatches(string &annotation,
                                                  int charge,
                                                  vector<psmPtr> &outputPsms)
  {
    string key = getKey(annotation, charge);
    std::tr1::unordered_multimap<string, unsigned int>::const_iterator it;
    std::pair
        < std::tr1::unordered_multimap<string, unsigned int>::const_iterator, std::tr1::unordered_multimap<
        string, unsigned int>::const_iterator > ret;

    ret = m_map.equal_range(key);

    for (it = ret.first; it != ret.second; ++it)
    {
      outputPsms.push_back(m_library[it->second]);
    }

    if (outputPsms.size() == 0)
    {
      return false;
    }
    else
    {
      return true;
    }
  }
}

