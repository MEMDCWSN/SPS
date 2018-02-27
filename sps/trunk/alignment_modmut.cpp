#include "alignment_modmut.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <cstring>
#include <limits.h>

#include "Logger.h"

namespace specnets
{
  const static int MAX_AA_JUMP_REPORT = 2048;

  // ****************************************************************************************************
  //   AMM_match
  // ****************************************************************************************************

  AMM_match::AMM_match(const AMM_match &other)
  {
    proteinIdx = other.proteinIdx;
    matchScore = other.matchScore;
    specParentMass = other.specParentMass;
    orientationPRMs = other.orientationPRMs;
    aaStart = other.aaStart;
    aaEnd = other.aaEnd;
    modSize = other.modSize;
    exonAllele = other.exonAllele;
    predExon = other.predExon;

    deque<TwoValues<float> >::const_iterator p = other.matchedMasses.begin();
    while (p != other.matchedMasses.end())
      matchedMasses.push_back(*p++);

    deque<TwoValues<int> >::const_iterator q = other.matchedIndices.begin();
    while (q != other.matchedIndices.end())
      matchedIndices.push_back(*q++);
  }

  // -------------------------------------------------------------------------
  void AMM_match::setMatch(Spectrum spec,
                           Spectrum dbSpec,
                           AMM_peak_match start,
                           vector<vector<vector<AMM_peak_match> > > &matchMatrix,
                           float curModSize)
  {
    TwoValues<float> tmpMasses;
    TwoValues<int> tmpIndices;

    matchedMasses.clear();
    matchedIndices.clear();
    matchScore = start.score;
    aaEnd = start.predDbSpecIdx;
    modSize = curModSize;
    //cerr << "Best match with " << modIdx << " modifications:\n";
    //cerr << "start.predModsIdx = "<<start.predMods<<", predDbSpecIdx = "<<start.predDbSpecIdx<<", predSpecIdx = "<<start.predSpecIdx<<endl;

    while (start.predSpecIdx >= 0)
    {
      tmpMasses[0] = spec[start.predSpecIdx][0];
      tmpMasses[1] = dbSpec[start.predDbSpecIdx][0];
      tmpIndices[0] = start.predSpecIdx;
      tmpIndices[1] = start.predDbSpecIdx;
      matchedMasses.push_back(tmpMasses);
      matchedIndices.push_back(tmpIndices);
      aaStart = start.predDbSpecIdx;
      //cerr << "["<<dbSpec[start.predDbSpecIdx][0]<<"],["<<spec[start.predSpecIdx][0]<<"]\n";
      start
          = matchMatrix[start.predMods][start.predDbSpecIdx][start.predSpecIdx];
      //cerr << "start.predModsIdx = "<<start.predMods<<", predDbSpecIdx = "<<start.predDbSpecIdx<<", predSpecIdx = "<<start.predSpecIdx<<endl;
    }
    specParentMass = spec.parentMass;
    predExon = start.predExonMatch;
  }

  // -------------------------------------------------------------------------

  void AMM_match::setMatch(Spectrum spec,
                           Spectrum dbSpec,
                           AMM_peak_match start,
                           vector<vector<vector<vector<AMM_peak_match> > > > &matchMatrix,
                           float curModSize)
  {
    TwoValues<float> tmpMasses;
    TwoValues<int> tmpIndices;

    matchedMasses.clear();
    matchedIndices.clear();
    matchScore = start.score;
    aaEnd = start.predDbSpecIdx;
    modSize = curModSize;
    //cerr << "Best match with " << modIdx << " modifications:\n";
    //cerr << "start.predModsIdx = "<<start.predMods<<", predDbSpecIdx = "<<start.predDbSpecIdx<<", predSpecIdx = "<<start.predSpecIdx<<endl;

    while (start.predSpecIdx >= 0)
    {
      tmpMasses[0] = spec[start.predSpecIdx][0];
      tmpMasses[1] = dbSpec[start.predDbSpecIdx][0];
      tmpIndices[0] = start.predSpecIdx;
      tmpIndices[1] = start.predDbSpecIdx;
      matchedMasses.push_back(tmpMasses);
      matchedIndices.push_back(tmpIndices);
      aaStart = start.predDbSpecIdx;
      //cerr << "["<<dbSpec[start.predDbSpecIdx][0]<<"],["<<spec[start.predSpecIdx][0]<<"]\n";
      start
          = matchMatrix[start.predMods][start.predDbSpecIdx][start.predSpecIdx][start.predMPeaks];
      //cerr << "start.predModsIdx = "<<start.predMods<<", predDbSpecIdx = "<<start.predDbSpecIdx<<", predSpecIdx = "<<start.predSpecIdx<<endl;
    }
    specParentMass = spec.parentMass;
    predExon = start.predExonMatch;
  }

  // -------------------------------------------------------------------------

  AMM_match *AMM_match::output_csv(ostream &output,
                                   DB_fasta &db,
                                   float tolerance,
                                   char sep,
                                   short outputPRMs,
                                   int specIdx,
                                   int numMods)
  {
    char *buffer = (char *)malloc(2048), *buffer2 = (char *)malloc(2048);
    if (matchedMasses.size() == 0 or matchedIndices.size() == 0)
      return predExon;

    if (exonAllele[0] == -1)
      sprintf(buffer2, "%d", proteinIdx);
    else
      sprintf(buffer2, "<%d,%d>", exonAllele[0], exonAllele[1]);

    if (specIdx >= 0 and numMods >= 0)
      output << specIdx << sep << numMods << sep;
    if (orientationPRMs == 0)
      sprintf(buffer, "PRMs interpreted as b");
    else
      sprintf(buffer, "PRMs interpreted as y");
    if (outputPRMs)
      output << buffer << sep << matchScore << sep << buffer2 << sep << aaStart
          << sep << aaEnd << sep << aaEnd - aaStart << sep
          << db.getID(proteinIdx) << sep << db.getDesc(proteinIdx) << endl;
    else
      output << buffer << sep << matchScore << sep << buffer2 << sep << aaStart
          << sep << aaEnd << sep << aaEnd - aaStart << sep;
    char *protSeq = db[proteinIdx];
    Spectrum &protSpec = db.getMassesSpec(proteinIdx); // Must be indexed with idx-1 because the first element is not zero!

    deque<TwoValues<float> >::reverse_iterator p = matchedMasses.rbegin();
    deque<TwoValues<int> >::reverse_iterator q = matchedIndices.rbegin();
    string originalStr, matchStr;

    // Include previous 4 amino acids in originalStr
    if ((*q)[1] < 0 or (*q)[1] >= strlen(protSeq))
    {
      ERROR_MSG("ERROR in AMM_match::output_csv: Invalid amino acid index "
          << (*q)[1] << " (protein length is " << strlen(protSeq) << ")");
      return (AMM_match *)0;
    }
    int pos = max(0, (*q)[1] - 4);
    strncpy(buffer, &protSeq[pos], (*q)[1] - pos);
    pos = (*q)[1] - pos;
    buffer[pos] = '.';
    buffer[pos + 1] = '\0';
    originalStr += string(buffer);

    int prevIdx = -1;
    float prevSpecMass = -1;
    while (p != matchedMasses.rend())
    {
      //cerr<<" ---> q = ("<<(*q)[0]<<","<<(*q)[1]<<"), p = ("<<(*p)[0]<<","<<(*p)[1]<<")\n"; cerr.flush();
      if (outputPRMs)
        output << "[" << (*q)[0] << ", " << (*p)[0] << "],[" << (*q)[1] << ", "
            << (*p)[1] << "]\n";

      // Format of the output sequences:
      // (PE)  - Spectrum contains a jump of mass(P)+mass(E)
      // {DGH} - Spectrum contains an unmatched jump of mass(DGH) at the start/end of the spectrum
      // [n1,n2,n3] - n1 is size of the jump in the spectrum, n2 is the size of the jump in the sequence, n3 is |n1-n2|

      if (prevIdx >= 0)
      {
        if ((*q)[1] - prevIdx <= 0)
          WARN_MSG("WARNING: matched protein indices should be strictly increasing, found "
              << prevIdx << "->" << (*q)[1] << "!");
        if ((*q)[1] > prevIdx and (*q)[1] - prevIdx <= MAX_AA_JUMP_REPORT)
        {
          strncpy(buffer, &protSeq[prevIdx], (*q)[1] - prevIdx);
          buffer[(*q)[1] - prevIdx] = '\0';
        }
        else
          sprintf(buffer, "[%.1f]", protSpec[(*q)[1]][0] - protSpec[prevIdx][0]);
        originalStr += string(buffer);
        //            originalStr+=string(buffer)+"\n";
        if (fabs((*p)[1] - protSpec[prevIdx][0] - ((*p)[0] - prevSpecMass))
            <= tolerance)
        {
          if ((*q)[1] - prevIdx > 1 and (*q)[1] - prevIdx <= MAX_AA_JUMP_REPORT)
          {
            strcpy(buffer2, buffer);
            sprintf(buffer, "(%s)", buffer2);
          }
          //            } else sprintf(buffer,"[%.1f,%.1f,%.1f]",(*p)[0]-prevSpecMass,(*p)[1]-protSpec[prevIdx][0],((*p)[0]-prevSpecMass)-((*p)[1]-protSpec[prevIdx][0]));
        }
        else
        {
          strcpy(buffer2, buffer);
          sprintf(buffer,
                  "[%.1f,%s,%.1f]",
                  (*p)[0] - prevSpecMass,
                  buffer2,
                  ((*p)[0] - prevSpecMass) - ((*p)[1] - protSpec[prevIdx][0]));
        }
        matchStr += string(buffer);
        //            matchStr+=string(buffer)+"\n";
      }
      else if ((*p)[0] > 0)
      {
        sprintf(buffer, "{%.1f}", (*p)[0]);
        matchStr += string(buffer);
      }
      prevIdx = (*q)[1];
      prevSpecMass = (*p)[0];
      p++;
      q++;
    }

    if (specParentMass - 19 - prevSpecMass > 2 * tolerance) // Trailing mass not matched to the DB sequence?
    {
      sprintf(buffer, "{%.1f}", specParentMass - 19 - prevSpecMass);
      matchStr += string(buffer);
    }

    // Include the next 4 amino acids in originalStr
    if (prevIdx < strlen(protSeq))
    {
      pos = min((int)strlen(protSeq) - 1, prevIdx + 4);
      strncpy(buffer, &protSeq[prevIdx], pos - prevIdx);
      buffer[pos - prevIdx] = '\0';
      originalStr += string(".");
      originalStr += string(buffer);
    }
    else
      originalStr += string(".");

    if (outputPRMs)
      output << "Sequence: " << originalStr << sep << matchStr << endl;
    else
      output << originalStr << sep << matchStr << sep << db.getID(proteinIdx)
          << sep << db.getDesc(proteinIdx) << endl;
    free(buffer);
    free(buffer2);
    return predExon;
  }

  // ****************************************************************************************************
  //   AMM_match_spec
  // ****************************************************************************************************

  void AMM_match_spec::addMatch(short numMods, AMM_match match)
  {
    bool sameRegion;

    if (match.matchScore <= lowestScore[numMods])
      return;
    deque<AMM_match>::iterator p = matches[numMods].begin();
    // Check if match will be added or remove any overlapping matches of lower score
    while (p != matches[numMods].end())
    {
      sameRegion = (match.proteinIdx == (*p).proteinIdx) and ((match.aaStart
          >= (*p).aaStart and match.aaStart <= (*p).aaEnd) || (match.aaEnd
          >= (*p).aaStart and match.aaEnd <= (*p).aaEnd));
      if (sameRegion)
      {
        if ((*p).matchScore - match.matchScore > 0.0001)
          return;
        else
        {
          if (match.matchScore > (*p).matchScore or fabs(match.modSize)
              <= fabs((*p).modSize))
            p = matches[numMods].erase(p);
          else
            return;
        }
      }
      else
        p++;
    }
    for (p = matches[numMods].begin(); p != matches[numMods].end()
        and match.matchScore > (*p).matchScore; p++)
      ;

    /* Old 05.05.18
     while (p!=matches[numMods].end() and match.matchScore<(*p).matchScore) p++;
     if(matches[numMods].size()>=numMatchesToKeep) matches[numMods].pop_back();
     deque<AMM_match>::iterator p = matches[numMods].begin();
     while (p!=matches[numMods].end() and match.matchScore<(*p).matchScore) p++;
     */
    matches[numMods].insert(p, match);
    if (matches[numMods].size() > numMatchesToKeep)
      matches[numMods].pop_front();
    lowestScore[numMods] = matches[numMods].front().matchScore;
  }

  // -------------------------------------------------------------------------
  void AMM_match_spec::output_csv(ostream &output,
                                  DB_fasta &db,
                                  float tolerance,
                                  char separator,
                                  short outputPRMs,
                                  int specIdx)
  {
    for (unsigned int i = 0; i < matches.size(); i++)
    {
      if (specIdx < 0)
        output << "Num mod/muts = " << i << endl;
      deque<AMM_match>::reverse_iterator p = matches[i].rbegin();
      while (p != matches[i].rend())
      {
        if (specIdx >= 0)
          p->output_csv(output,
                        db,
                        tolerance,
                        separator,
                        outputPRMs,
                        specIdx,
                        i);
        else
          p->output_csv(output, db, tolerance, separator, outputPRMs);
        p++;
      }
    }
  }

  // ****************************************************************************************************
  //   AMM_match_set
  // ****************************************************************************************************

  AMM_match_set::AMM_match_set(DB_fasta &in_db, unsigned int szSpecSet)
  {
    db = &in_db;
    matches.resize(db->size());
    resize(szSpecSet);
  }

  // -------------------------------------------------------------------------
  unsigned int AMM_match_set::resize(unsigned int szSpecSet)
  {
    unsigned int oldSize = matchedPeaks.size();

    matchedPeaks.resize(szSpecSet);
    specProtMatches.resize(szSpecSet);
    for (unsigned int i = oldSize; i < szSpecSet; i++)
    {
      specProtMatches[i].resize(3);
      specProtMatches[i][0] = -1;
      specProtMatches[i][1] = 0;
      specProtMatches[i][2] = 0;
    }
    return szSpecSet;
  }

  // -------------------------------------------------------------------------
  bool AMM_match_set::set(unsigned int specIdx, AMM_match &match)
  {

    //DEBUG_MSG(">>> specIdx = " << specIdx << ", matchedPeaks.size() = "
    //    << matchedPeaks.size() << ", matches.size() = " << matches.size()
    //    << ", specProtMatches.size() = " << specProtMatches.size());
    if (specIdx >= matchedPeaks.size() or specIdx >= matches.size() or specIdx
        >= specProtMatches.size())
      return false;

    specProtMatches[specIdx][0] = match.proteinIdx;
    specProtMatches[specIdx][1] = 0; // Can't tell from AMM_match
    specProtMatches[specIdx][2] = (int)match.orientationPRMs;
    //DEBUG_MSG(">>> specProtMatches[" << specIdx << "] = "
    //    << specProtMatches[specIdx][0] << ", " << specProtMatches[specIdx][1]
    //    << ", " << specProtMatches[specIdx][2]);

    unsigned int pairIdx = 0;
    matchedPeaks[specIdx].resize(match.matchedMasses.size());
    //DEBUG_MSG(">>> matchedPeaks[" << specIdx << "].size() = "
    //    << matchedPeaks[specIdx].size());
    for (deque<TwoValues<int> >::reverse_iterator iter =
        match.matchedIndices.rbegin(); iter != match.matchedIndices.rend(); iter++, pairIdx++)
      matchedPeaks[specIdx][pairIdx].set((float)(*iter)[0], (float)(*iter)[1]);

    TwoValues<int> curMatch;
    if (match.proteinIdx >= 0 and match.proteinIdx < matches.size())
      if (matchedPeaks[specIdx].size() > 0)
      {
        curMatch.set((int)round(matchedPeaks[specIdx][0][1]), (int)specIdx);
        matches[match.proteinIdx].push_back(curMatch);
        matches[match.proteinIdx].sort();
      }

    return true;
  }

  // -------------------------------------------------------------------------
  bool AMM_match_set::SetMatches(vector<vector<int> > &in_matchedProts,
                                 SpecSet &in_matchedPeaks)
  {
    // Set member variables if needed
    if (&in_matchedProts != &specProtMatches)
    {
      specProtMatches.resize(in_matchedProts.size());
      for (unsigned int i = 0; i < in_matchedProts.size(); i++)
      {
        specProtMatches[i].resize(in_matchedProts[i].size());
        for (unsigned int j = 0; j < in_matchedProts[i].size(); j++)
          specProtMatches[i][j] = in_matchedProts[i][j];
      }
    }

    if (&in_matchedPeaks != &matchedPeaks)
    {
      matchedPeaks.resize(in_matchedPeaks.size());
      for (unsigned int i = 0; i < in_matchedProts.size(); i++)
        matchedPeaks[i] = in_matchedPeaks[i];
    }

    // Read all matches into matches
    TwoValues<int> curMatch;
    for (unsigned int matchIdx = 0; matchIdx < matchedPeaks.size(); matchIdx++)
    {
      if (specProtMatches[matchIdx][0] >= 0 and specProtMatches[matchIdx][0]
          < matches.size())
      {
        if (matchedPeaks[matchIdx].size() > 0)
        {
          curMatch.set((int)round(matchedPeaks[matchIdx][0][1]), (int)matchIdx);
          matches[specProtMatches[matchIdx][0]].push_back(curMatch);
        }
      }
    }

    // Sort every matches[i] by increasing start position on the protein
    for (unsigned int pIdx = 0; pIdx < matches.size(); pIdx++)
      matches[pIdx].sort();

    return true;
  }

  bool AMM_match_set::SetMatches(SpecSet * specset)
  {
    specProtMatches.resize(specset->size());
    for (unsigned int i = 0; i < specset->size(); i++)
    {
      specProtMatches[i].resize(3,0);
      if ((*specset)[i].psmList.size() > 0)
      {
        specProtMatches[i][0] = (*specset)[i].psmList.front()->m_dbIndex;
        specProtMatches[i][1] = (*specset)[i].psmList.front()->m_numMods;
        specProtMatches[i][2] = (*specset)[i].psmList.front()->m_matchOrientation;
      }
    }

    matchedPeaks.resize(specset->size());
    for (unsigned int i = 0; i < specset->size(); i++)
    {
      if ((*specset)[i].psmList.size() == 0)
      {
        continue;
      }
      int peakListSize = (*specset)[i].psmList.front()->m_matchedPeaks.size();
      matchedPeaks[i].resize(peakListSize);
      for (unsigned int j = 0; j < peakListSize; j++)
      {
        matchedPeaks[i][j][0]
            = (*specset)[i].psmList.front()->m_matchedPeaks[j][0];
        matchedPeaks[i][j][1]
            = (*specset)[i].psmList.front()->m_matchedPeaks[j][1];
      }
    }

    // Read all matches into matches
    TwoValues<int> curMatch;
    for (unsigned int matchIdx = 0; matchIdx < matchedPeaks.size(); matchIdx++)
    {
      if (specProtMatches[matchIdx][0] >= 0 and specProtMatches[matchIdx][0]
          < matches.size())
      {
        if (matchedPeaks[matchIdx].size() > 0)
        {
          curMatch.set((int)round(matchedPeaks[matchIdx][0][1]), (int)matchIdx);
          matches[specProtMatches[matchIdx][0]].push_back(curMatch);
        }
      }
    }

    // Sort every matches[i] by increasing start position on the protein
    for (unsigned int pIdx = 0; pIdx < matches.size(); pIdx++)
    {
      matches[pIdx].sort();
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool AMM_match_set::LoadMatches(const char *matchedProts,
                                  const char *matchedPeaksIdx)
  {

    Load_binArray<int> (matchedProts, specProtMatches);
    if (specProtMatches.size() == 0 or specProtMatches[0].size() != 3)
    {
      ERROR_MSG("ERROR loading " << matchedProts << "!");
      return false;
    }
    if (matchedPeaks.loadPklBin(matchedPeaksIdx) <= 0
        or specProtMatches.size() != matchedPeaks.size())
    {
      ERROR_MSG("ERROR loading " << matchedPeaksIdx << "!");
      return false;
    }
    return SetMatches(specProtMatches, matchedPeaks);
  }

  // -------------------------------------------------------------------------
  void AMM_match_set::GenerateGlues(SpectrumPairSet &pairsPA, vector<vector<
      TwoValues<int> > > &vPairMatchedPeaks)
  {

    list < vector<TwoValues<int> > > pairMatchedPeaks;
    list<TwoValues<int> > specPairs;
    TwoValues<int> curPair;
    unsigned int lastProtMatch, idxSpec1, idxSpec2, peakIdx1, peakIdx2,
        numMatches;
    list<TwoValues<int> >::iterator specsIter, otherIter;
    vector<TwoValues<int> > curMatches;
    curMatches.reserve(1000);

    pairsPA.resize(0);
    pairMatchedPeaks.resize(0);
    for (unsigned int pIdx = 0; pIdx < matches.size(); pIdx++)
      for (specsIter = matches[pIdx].begin(); specsIter != matches[pIdx].end(); specsIter++)
      {
        idxSpec1 = (*specsIter)[1];
        lastProtMatch
            = (int)round(matchedPeaks[idxSpec1][matchedPeaks[idxSpec1].size()
                - 1][1]);
        otherIter = specsIter;
        otherIter++;

        for (; otherIter != matches[pIdx].end() and (*otherIter)[0]
            <= lastProtMatch; otherIter++)
        { // Spectra overlap, look for matches
          idxSpec2 = (*otherIter)[1];
          peakIdx1 = 0;
          peakIdx2 = 0;
          numMatches = 0;
          curMatches.resize(min(matchedPeaks[idxSpec1].size(),
                                matchedPeaks[idxSpec2].size()));
          while (peakIdx1 < matchedPeaks[idxSpec1].size() and peakIdx2
              < matchedPeaks[idxSpec2].size())
          {
            if (matchedPeaks[idxSpec1][peakIdx1][1]
                == matchedPeaks[idxSpec2][peakIdx2][1])
            {
              curMatches[numMatches++].set((int)round(matchedPeaks[idxSpec1][peakIdx1++][0]),
                                           (int)round(matchedPeaks[idxSpec2][peakIdx2++][0]));
            }
            else
            {
              if (matchedPeaks[idxSpec1][peakIdx1][1]
                  < matchedPeaks[idxSpec2][peakIdx2][1])
                peakIdx1++;
              else
                peakIdx2++;
            }
          }

          if (numMatches > 0)
          {
            curPair.set(idxSpec1, idxSpec2);
            specPairs.push_back(curPair);
            curMatches.resize(numMatches);
            pairMatchedPeaks.push_back(curMatches);
          }
        }
      }

    unsigned int pairIdx;
    pairsPA.resize(specPairs.size());
    otherIter = specPairs.begin();
    for (pairIdx = 0; pairIdx < pairsPA.size(); pairIdx++, otherIter
        = specPairs.erase(otherIter))
    {
      pairsPA[pairIdx].spec1 = (*otherIter)[0];
      pairsPA[pairIdx].spec2 = (*otherIter)[1];
    }
    vPairMatchedPeaks.resize(pairMatchedPeaks.size());
    list<vector<TwoValues<int> > >::iterator mpIter = pairMatchedPeaks.begin();
    for (pairIdx = 0; pairIdx < vPairMatchedPeaks.size(); pairIdx++, mpIter
        = pairMatchedPeaks.erase(mpIter))
      vPairMatchedPeaks[pairIdx].swap(*mpIter);
  }

  // -------------------------------------------------------------------------
  //  MergeIntoReference - Merges the spectrum matches to protein otherIdx into protein refIdx using
  //   the protein-protein correspondence specified in matchedIndices. Matches to otherIdx with no
  //   correspondence in refIdx (via matchedIndices) remain associated with protein otherIdx.
  //
  //   refIdx / otherIdx - Index of the protein where the spectra matches are imported to/from, respectively
  //   matchedIndices - matched mass indices between proteins refIdx and otherIdx (cols 0/1)
  // -------------------------------------------------------------------------
  void AMM_match_set::MergeIntoReference(int refIdx, int otherIdx, vector<
      TwoValues<int> > &matchedIndices)
  {
    list<TwoValues<int> >::iterator otherIter = matches[otherIdx].begin();
    unsigned int specIdx, peakIdx, mPeaksIdx, // Iterator over contig/other-protein matches (AMM_match_set::matchedPeaks)
        mProtsIdx; // Iterator over ref-protein/other-protein matches (e.g. clustalw)
    TwoValues<int> curMatch;

    mProtsIdx = 0;
    while (otherIter != matches[otherIdx].end())
    {
      peakIdx = (*otherIter)[0]; // First matched protein mass
      specIdx = (*otherIter)[1]; // Index of matched contig

      // Find aligned protein-protein mass also matched to the contig
      if (mProtsIdx >= matchedIndices.size())
        mProtsIdx = matchedIndices.size() - 1;
      while (mProtsIdx > 1 and matchedIndices[mProtsIdx][1] > peakIdx)
        mProtsIdx--;
      mPeaksIdx = 0;
      while (mPeaksIdx < matchedPeaks[specIdx].size() and mProtsIdx
          < matchedIndices.size() and matchedPeaks[specIdx][mPeaksIdx][1]
          != matchedIndices[mProtsIdx][1])
        if (matchedPeaks[specIdx][mPeaksIdx][1] < matchedIndices[mProtsIdx][1])
          mPeaksIdx++;
        else
          mProtsIdx++;

      // Transfer contig matches to reference protein
      if (mPeaksIdx < matchedPeaks[specIdx].size() and mProtsIdx
          < matchedIndices.size() and matchedPeaks[specIdx][mPeaksIdx][1]
          == matchedIndices[mProtsIdx][1])
      {
        curMatch.set(matchedIndices[mProtsIdx][0], specIdx);
        matches[refIdx].push_back(curMatch);
        peakIdx = 0; // Iterator over matched contig/ref-protein matches
        mPeaksIdx = 0; // Iterator over matched contig/other-protein matches
        while (mPeaksIdx < matchedPeaks[specIdx].size() and mProtsIdx
            < matchedIndices.size())
          if (matchedPeaks[specIdx][mPeaksIdx][1]
              == matchedIndices[mProtsIdx][1])
          {
            matchedPeaks[specIdx][peakIdx++].set(matchedPeaks[specIdx][mPeaksIdx++][0],
                                                 matchedIndices[mProtsIdx++][0]);
          }
          else
          {
            if (matchedPeaks[specIdx][mPeaksIdx][1]
                < matchedIndices[mProtsIdx][1])
              mPeaksIdx++;
            else
              mProtsIdx++;
          }
        matchedPeaks[specIdx].resize(peakIdx);
        specProtMatches[specIdx][0] = refIdx;
        otherIter = matches[otherIdx].erase(otherIter);
      }
      else
        otherIter++;
    }

    matches[refIdx].sort();
  }

  // ****************************************************************************************************
  //   Functions - scoreOverlapAMM
  //
  //   curMatch - information on current match: orientation of PRMs/endpoints, proteinIdx
  //   minSpecDist - minimum distance between matched peaks in the PRM spectrum
  //   maxDbSpecMod - maximum mass jump of amy single modification on the database masses
  //   enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks
  //                       (consuming available modifications as necessary)
  // ****************************************************************************************************

  //vector<vector<AMM_match> > scoreOverlapAMM(SpecSet matchSpecs, DB_fasta db, Clusters clst, SpecSet allSpecs, short maxNumMods, short keepTopK, float pmTolerance, float tolerance) {
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       short maxNumMods,
                       AMM_match_spec &topKmatches,
                       AMM_match &curMatch,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks)
  {
    int specIdx, dbSpecIdx, modIdx, predSpecIdx, predDbSpecIdx;

    vector<vector<vector<AMM_peak_match> > > matchMatrix, bestSoFar;
    AMM_peak_match start; // Used at the end to submit entries to topKmatches

    // Compute the match between the two spectra
    float usedModSize; // Keeps track of the size of the mods used - lower mods are considered better
    matchMatrix.resize(maxNumMods + 1);
    bestSoFar.resize(maxNumMods + 1);
    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      matchMatrix[modIdx].resize(dbSpec.size());
      bestSoFar[modIdx].resize(dbSpec.size());
      for (dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++)
      {
        matchMatrix[modIdx][dbSpecIdx].resize(spec.size());
        bestSoFar[modIdx][dbSpecIdx].resize(spec.size());

        for (specIdx = 0; specIdx < spec.size(); specIdx++)
        {
          // Line below: Matching can start on any mass, assume no predecessor (if enforceEnpeaks is false)
          if (!enforceEndpeaks or specIdx == 0)
            matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1],
                                                                   modIdx,
                                                                   -1,
                                                                   dbSpecIdx
                                                                       - 1);
          //                if(specIdx==0) matchMatrix[modIdx][dbSpecIdx][0].replaceWithMax(spec[0][1],modIdx,-1,max(0,dbSpecIdx-1));  // max(0,dbSpecIdx-1) because position 0 in dbSpec is already the mass of the first AA
          //cerr << "dbSpecIdx = " << dbSpecIdx << ", specIdx = " << specIdx << endl;

          // 0 Mods - find other peaks at the same distance in both spectra
          predSpecIdx = specIdx - 1;
          predDbSpecIdx = dbSpecIdx - 1;
          usedModSize = 0;
          // Zero mods case also needs to search for best predecessor because the previous peaks may have been
          //   missed due to peak tolerance problems. Search for allowable jump with maximum total score
          while (predSpecIdx >= 0 && predDbSpecIdx >= 0)
          {
            // Find an eligible predecessor
            while (predSpecIdx >= 0 && predDbSpecIdx >= 0
                && (fabs((spec[specIdx][0] - spec[predSpecIdx][0])
                    - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                    > tolerance))
            {
              if ((spec[specIdx][0] - spec[predSpecIdx][0])
                  > (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                predDbSpecIdx--;
              else
                predSpecIdx--;
            }
            // Use predecessor if score increases
            if (predSpecIdx >= 0 && predDbSpecIdx >= 0
                && matchMatrix[modIdx][dbSpecIdx][specIdx].score
                    < spec[specIdx][1]
                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx].score)
              matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1]
                                                                         + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx].score,
                                                                     modIdx,
                                                                     predSpecIdx,
                                                                     predDbSpecIdx);
            predSpecIdx--; //predDbSpecIdx--;
          }

          // 1 Mods - any pair of peaks with one less mutation is fair game
          if (modIdx > 0 && dbSpecIdx > 0 && specIdx > 0 && (bestSoFar[modIdx
              - 1][dbSpecIdx - 1][specIdx - 1].score + spec[specIdx][1]
              > matchMatrix[modIdx][dbSpecIdx][specIdx].score))
          {
            // Find first eligible predecessor
            predSpecIdx = specIdx - 1;
            predDbSpecIdx = dbSpecIdx - 1;
            // Find first eligible predSpecIdx
            while (predSpecIdx >= 0 && spec[specIdx][0] - spec[predSpecIdx][0]
                < minSpecDist - 2 * tolerance)
              predSpecIdx--;
            // Find first eligible predSpecIdx/predDbSpecIdx pair
            //                  while(predSpecIdx>=0 && predDbSpecIdx>=0 && fabs(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]-(spec[specIdx][0]-spec[predSpecIdx][0]))>=maxDbSpecMod+2*tolerance)
            while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                and ((spec[specIdx][0] - spec[predSpecIdx][0]
                    - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                    > maxDbSpecMod + tolerance or (spec[specIdx][0]
                    - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                    - dbSpec[predDbSpecIdx][0])) < minDbSpecMod - tolerance))
            {
              if (predDbSpecIdx > 0 and spec[specIdx][0] - spec[predSpecIdx][0]
                  - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0])
                  > minDbSpecMod - tolerance)
                predDbSpecIdx--;
              else
              {
                predSpecIdx--;
                predDbSpecIdx = dbSpecIdx - 1;
              }
            }
            // Initialize with first eligible predecessor pair
            if (predSpecIdx >= 0 && predDbSpecIdx >= 0)
            {
              matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1]
                                                                         + matchMatrix[modIdx
                                                                             - 1][predDbSpecIdx][predSpecIdx].score,
                                                                     modIdx - 1,
                                                                     predSpecIdx,
                                                                     predDbSpecIdx);
              usedModSize = fabs(dbSpec[dbSpecIdx][0]
                  - dbSpec[predDbSpecIdx][0] - (spec[specIdx][0]
                  - spec[predSpecIdx][0]));
            }

            // Try all eligible predecessors
            //                  while(predSpecIdx>=0 && predDbSpecIdx>=0 && dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]-(spec[specIdx][0]-spec[predSpecIdx][0])<=maxDbSpecMod+2*tolerance) {
            while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
            {
              float tmpScore =
                  matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx].score
                      + spec[specIdx][1], tmpModSize =
                  fabs(dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]
                      - (spec[specIdx][0] - spec[predSpecIdx][0]));

              if (tmpScore > matchMatrix[modIdx][dbSpecIdx][specIdx].score
                  || (matchMatrix[modIdx][dbSpecIdx][specIdx].score - tmpScore
                      < 0.0001 && tmpModSize < usedModSize))
              {
                matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(tmpScore,
                                                                       modIdx
                                                                           - 1,
                                                                       predSpecIdx,
                                                                       predDbSpecIdx);
                usedModSize = tmpModSize;
              }
              // Find next eligible predecessor pair
              predDbSpecIdx--;
              while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                  and ((spec[specIdx][0] - spec[predSpecIdx][0]
                      - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                      > maxDbSpecMod + tolerance or (spec[specIdx][0]
                      - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                      - dbSpec[predDbSpecIdx][0])) < minDbSpecMod - tolerance))
                if (predDbSpecIdx > 0 and spec[specIdx][0]
                    - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                    - dbSpec[predDbSpecIdx][0]) > minDbSpecMod - tolerance)
                  predDbSpecIdx--;
                else
                {
                  predSpecIdx--;
                  predDbSpecIdx = dbSpecIdx - 1;
                }
            }
          }

          // Update bestSoFar; give precedence to options with same score and less mutations
          if (modIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx].replaceWithMax(bestSoFar[modIdx
                - 1][dbSpecIdx][specIdx]);
          if (specIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx][specIdx
                - 1]);
          if (dbSpecIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                - 1][specIdx]);
          bestSoFar[modIdx][dbSpecIdx][specIdx].replaceWithMax(matchMatrix[modIdx][dbSpecIdx][specIdx].score,
                                                               modIdx,
                                                               specIdx,
                                                               dbSpecIdx);

        } // for specIdx

        if (enforceEndpeaks)
        { // Enforce that matches in bestSoFar[][][spec.size()-1] must end at spec.size()-1
          bestSoFar[modIdx][dbSpecIdx][specIdx - 1]
              = matchMatrix[modIdx][dbSpecIdx][specIdx - 1];
          if (modIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1].replaceWithMax(bestSoFar[modIdx
                - 1][dbSpecIdx][specIdx - 1]);
          if (dbSpecIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                - 1][specIdx - 1]);
        }
      } // for dbSpecIdx
    } // for modIdx

    // Update topKmatches
    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      // Keep track of topKmatches
      if (bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1].score
          > topKmatches.lowestScore[modIdx])
      {
        curMatch.setMatch(spec,
                          dbSpec,
                          bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1],
                          matchMatrix,
                          usedModSize);
        topKmatches.addMatch(modIdx, curMatch);
      }
    }
  }

  // -------------------------------------------------------------------------
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       short maxNumMods,
                       int minMatchedPeaks,
                       AMM_match_spec &topKmatches,
                       AMM_match &curMatch,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks)
  {
    int specIdx, dbSpecIdx, modIdx, mPeaksIdx, predSpecIdx, predDbSpecIdx,
        predMPeaksIdx;
    if (spec.size() == 0 or spec.size() < minMatchedPeaks)
      return;

    // matchMatrix[i,j,k,p] where
    //  i - number of used mod/mut mass offsets
    //  j - index in dbSpec
    //  k - index in spec
    //  p - at least p+1 matched peaks in spec
    vector<vector<vector<vector<AMM_peak_match> > > > matchMatrix, // Each [i,j,k,p] position contains the score and predecessors for the highest-scoring
        //   match up to and including masses j,k in dbSpec and spec, respectively, with i mod/mut
        //   mass offsets and at least p matched peaks
        bestSoFar; // Each [i,j,k,p] position contains the coordinates of the highest scoring match
    //   in matchMatrix with coordinates up to and possibly including [i,j,k,p]

    // Compute the match between the two spectra
    AMM_peak_match tmp;
    float usedModSize; // Keeps track of the size of the mods used - lower mods are considered better
    matchMatrix.resize(maxNumMods + 1);
    bestSoFar.resize(maxNumMods + 1);
    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      matchMatrix[modIdx].resize(dbSpec.size());
      bestSoFar[modIdx].resize(dbSpec.size());
      for (dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++)
      {
        matchMatrix[modIdx][dbSpecIdx].resize(spec.size());
        bestSoFar[modIdx][dbSpecIdx].resize(spec.size());

        for (specIdx = 0; specIdx < spec.size(); specIdx++)
        {
          matchMatrix[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);
          bestSoFar[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);
          // Line below: Matching can start on any mass, assume no predecessor (if enforceEnpeaks is false)
          if (!enforceEndpeaks or specIdx == 0)
            matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],
                                                                      modIdx,
                                                                      -1,
                                                                      -1,
                                                                      -1);
          //          matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],modIdx,-1,dbSpecIdx-1,-1);

          unsigned int predMPeaksIdx;
          for (mPeaksIdx = 0; mPeaksIdx < min(minMatchedPeaks, specIdx + 1); mPeaksIdx++)
          {
            //cerr<<" >> start ("<<modIdx<<","<<dbSpecIdx<<","<<specIdx<<","<<mPeaksIdx<<"), score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<", peak score = "<<spec[specIdx][1]<<"\n";
            // 0 Mods - find other peaks at the same distance in both spectra
            predSpecIdx = specIdx - 1;
            while (predSpecIdx >= 0 and fabs(spec[specIdx][0]
                - spec[predSpecIdx][0]) < minSpecDist - tolerance)
              predSpecIdx--;
            predDbSpecIdx = dbSpecIdx - 1;
            usedModSize = 0;
            // Zero mods case also needs to search for best predecessor because the previous peaks may have been
            //   missed due to peak tolerance problems. Search for allowable jump with maximum total score
            if (mPeaksIdx > 0 or mPeaksIdx == minMatchedPeaks - 1)
            {
              while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
              {
                // Find an eligible predecessor
                while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and (fabs((spec[specIdx][0] - spec[predSpecIdx][0])
                        - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                        > tolerance))
                {
                  //if(modIdx==0 and dbSpecIdx==11 and specIdx==9)
                  //  cerr<<" -0-->> skipped pred ("<<predSpecIdx<<","<<predDbSpecIdx<<"): fabs("<<(spec[specIdx][0]-spec[predSpecIdx][0])<<","<<(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0])<<")="<<fabs((spec[specIdx][0]-spec[predSpecIdx][0])-(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]))<<" - ";
                  if ((spec[specIdx][0] - spec[predSpecIdx][0])
                      > (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                    predDbSpecIdx--;
                  else
                    predSpecIdx--;
                  //if(modIdx==0 and dbSpecIdx==11 and specIdx==9)
                  //  cerr<<"next putative pred is ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";
                }
                //cerr<<" -0- pred zero-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), cur score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score;
                //if (mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                //else cerr<<" (no pred found)\n";
                // Use predecessor if score increases
                if (mPeaksIdx > 0 and predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                - 1].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                                                                            - 1].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx
                                                                                        - 1);
                }
                //cerr<<" -0- pred zero-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), mPeaksIdx = "<<mPeaksIdx<<", minMatchedPeaks = "<<minMatchedPeaks;
                //if (mPeaksIdx==minMatchedPeaks-1 and mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                //else cerr<<" (no pred found)\n";
                if (mPeaksIdx == minMatchedPeaks - 1 and predSpecIdx >= 0
                    and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx);
                  //cerr<<" -0--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                }
                predSpecIdx--; //predDbSpecIdx--;
              }

              // 1+ Mods - any pair of peaks with one less mutation and within mod mass limits is fair game
              float maxPredScore = -3.4e37, tmpScore, tmpModSize;
              if (modIdx > 0 and dbSpecIdx > 0 and specIdx > 0)
              {
                if (mPeaksIdx > 0)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx
                                - 1].score);
                if (mPeaksIdx == minMatchedPeaks - 1)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx].score);
                if (maxPredScore + spec[specIdx][1]
                    > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score)
                {

                  // Find first eligible predecessor
                  predSpecIdx = specIdx - 1;
                  predDbSpecIdx = dbSpecIdx - 1;
                  // Find first eligible predSpecIdx
                  while (predSpecIdx >= 0 && spec[specIdx][0]
                      - spec[predSpecIdx][0] < minSpecDist - 2 * tolerance)
                    predSpecIdx--;
                  // Find first eligible predSpecIdx/predDbSpecIdx pair
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                      and ((spec[specIdx][0] - spec[predSpecIdx][0]
                          - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                          > maxDbSpecMod + tolerance or (spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                          - tolerance))
                  {
                    if (predDbSpecIdx > 0 and spec[specIdx][0]
                        - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0]) > minDbSpecMod - tolerance)
                      predDbSpecIdx--;
                    else
                    {
                      predSpecIdx--;
                      predDbSpecIdx = dbSpecIdx - 1;
                    }
                  }
                  // Initialize with first eligible predecessor pair
                  //                if(predSpecIdx>=0 and predDbSpecIdx>=0) {
                  //                  matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1]+matchMatrix[modIdx-1][predDbSpecIdx][predSpecIdx].score,modIdx-1,predSpecIdx,predDbSpecIdx);
                  //                  usedModSize = fabs(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]-(spec[specIdx][0]-spec[predSpecIdx][0]));
                  //                }

                  // Try all eligible predecessors
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
                  {
                    if (mPeaksIdx > 0)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                  - 1].score;
                      predMPeaksIdx = mPeaksIdx - 1;
                    }
                    if (mPeaksIdx == minMatchedPeaks - 1
                        and tmpScore < spec[specIdx][1] + matchMatrix[modIdx
                            - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score;
                      predMPeaksIdx = mPeaksIdx;
                    }
                    tmpModSize = fabs(dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0] - (spec[specIdx][0]
                        - spec[predSpecIdx][0]));

                    //cerr<<" -1- pred ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";
                    //cerr<<" -1- pred one-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), cur score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score;
                    //if (mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                    //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                    //else cerr<<" (no pred found)\n";

                    if (tmpScore
                        > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        or (matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                            - tmpScore < 0.0001 and tmpModSize < usedModSize))
                    {
                      matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(tmpScore,
                                                                                        modIdx
                                                                                            - 1,
                                                                                        predSpecIdx,
                                                                                        predDbSpecIdx,
                                                                                        predMPeaksIdx);
                      usedModSize = tmpModSize;
                      //cerr<<" -1--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                    }
                    // Find next eligible predecessor pair
                    predDbSpecIdx--;
                    while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                        and ((spec[specIdx][0] - spec[predSpecIdx][0]
                            - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                            > maxDbSpecMod + tolerance or (spec[specIdx][0]
                            - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                            - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                            - tolerance))
                      if (predDbSpecIdx > 0 and spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0]) > minDbSpecMod
                          - tolerance)
                        predDbSpecIdx--;
                      else
                      {
                        predSpecIdx--;
                        predDbSpecIdx = dbSpecIdx - 1;
                      }
                  }
                } // 1+ mods, eligible predecessor
              } // 1+ mods
            } // mPeaksIdx>0 or mPeaksIdx==minMatchedPeaks-1

            // Update bestSoFar; give precedence to options with same score and less mutations or more matched peaks
            if (modIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx
                  - 1][dbSpecIdx][specIdx][mPeaksIdx]);
            if (specIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx][specIdx
                  - 1][mPeaksIdx]);
            if (dbSpecIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                  - 1][specIdx][mPeaksIdx]);
            bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score,
                                                                            modIdx,
                                                                            specIdx,
                                                                            dbSpecIdx,
                                                                            mPeaksIdx);

            //tmp = matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx];
            //cerr<<" >> end: (mod="<<modIdx<<",db="<<dbSpecIdx<<",spec="<<specIdx<<",mp="<<mPeaksIdx<<"): score = "<<tmp.score<<", pred=("<<tmp.predMods<<","<<tmp.predDbSpecIdx<<","<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<")\n"; cerr.flush();
          } // for mPeaksIdx
        } // for specIdx

        //tmp = bestSoFar[modIdx][dbSpecIdx][specIdx-1][minMatchedPeaks-1];
        //cerr<<"(mod="<<modIdx<<",db="<<dbSpecIdx<<"): score = "<<tmp.score<<", mod = "<<tmp.predMods<<", db = "<<tmp.predDbSpecIdx<<", spec = "<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<"\n"; cerr.flush();

        if (enforceEndpeaks)
        { // Enforce that matches in bestSoFar[][][spec.size()-1] must end at spec.size()-1
          bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]
              = matchMatrix[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1];
          if (modIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx
                - 1][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]);
          if (dbSpecIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                - 1][specIdx - 1][minMatchedPeaks - 1]);
        }
      } // for dbSpecIdx
    } // for modIdx

    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      // Keep track of topKmatches
      if (bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1][minMatchedPeaks - 1].score
          > topKmatches.lowestScore[modIdx])
      {
        curMatch.setMatch(spec,
                          dbSpec,
                          bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1][minMatchedPeaks
                              - 1],
                          matchMatrix,
                          usedModSize);
        topKmatches.addMatch(modIdx, curMatch);
      }
    }
  }

  // ****************************************************************************************************
  //   Functions - scoreOverlapAMM
  //
  //   minSpecDist - minimum distance between matched peaks in the PRM spectrum
  //   maxDbSpecMod - maximum mass jump of amy single modification on the database masses
  //   enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks
  //                       (consuming available modifications as necessary)
  // ****************************************************************************************************
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       int dbIndex,
                       int matchOrientation,
                       set<float> & startRange,
                       short maxNumMods,
                       int minMatchedPeaks,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks)
  {
    int specIdx, dbSpecIdx, modIdx, mPeaksIdx, predSpecIdx, predDbSpecIdx,
        predMPeaksIdx;
    if (spec.size() == 0 or spec.size() < minMatchedPeaks)
    {
      WARN_MSG("Spectrum size [" << spec.size() << "] is less than minimum matched peaks [" << minMatchedPeaks << "]");
      return;
    }

    // matchMatrix[i,j,k,p] where
    //  i - number of used mod/mut mass offsets
    //  j - index in dbSpec
    //  k - index in spec
    //  p - at least p+1 matched peaks in spec
    // Each [i,j,k,p] position contains the score and predecessors for the highest-scoring
    vector<vector<vector<vector<AMM_peak_match> > > > matchMatrix,
    //   match up to and including masses j,k in dbSpec and spec, respectively, with i mod/mut
        //   mass offsets and at least p matched peaks
        bestSoFar; // Each [i,j,k,p] position contains the coordinates of the highest scoring match
    //   in matchMatrix with coordinates up to and possibly including [i,j,k,p]

    // Initialize starting flag array
    vector < vector<char> > startFlags(spec.size());
    for (int i = 0; i < spec.size(); i++)
    {
      vector<char> newArray(dbSpec.size());
      startFlags[i] = newArray;
      for (int j = 0; j < dbSpec.size(); j++)
      {
        // If the startRange set is empty then allow any start position
        if (startRange.size() == 0)
        {
          startFlags[i][j] = 1;
        }
        else
        {
          startFlags[i][j] = 0;
        }
      }
    }

    // Find all the valid starting points in the matrix    
    set<float>::iterator itr = startRange.begin();
    set<float>::iterator itrEnd = startRange.end();
    for (; itr != itrEnd; itr++)
    {
      int minIdx1, maxIdx1, minIdx2, maxIdx2;
      float minStartMass = *itr + minDbSpecMod * maxNumMods;
      float maxStartMass = *itr + maxDbSpecMod * maxNumMods;
      float lastMass = 0;
      for (int specIdx = 0; specIdx < spec.size(); specIdx++)
      {
        float massDiff = spec[specIdx][0] - lastMass;
        minStartMass += massDiff;
        maxStartMass += massDiff;
        lastMass = spec[specIdx][0];
        // Find peaks close to desired bounds
        float T = 200.0;

        list<int> matches1;
        list<int> matches2;
        dbSpec.setPeakTolerance(0);
        dbSpec.findPeaks(minStartMass, T, &matches1);
        dbSpec.findPeaks(maxStartMass, T, &matches2);
        minIdx1 = (matches1.size() == 0) ? -1 : matches1.front();
        maxIdx1 = (matches1.size() == 0) ? -1 : matches1.back();
        minIdx2 = (matches2.size() == 0) ? -1 : matches2.front();
        maxIdx2 = (matches2.size() == 0) ? -1 : matches2.back();

        // If peak is not found, index could be -1.. so make sure at least 0
        minIdx1 = max<int> (minIdx1, 0);
        minIdx2 = max<int> (minIdx2, 0);
        // Make sure min peak is INSIDE the desired range
        while (minIdx1 < dbSpec.size() && dbSpec[minIdx1][0] < minStartMass
            - tolerance - AAJumps::massH2O)
          minIdx1++;

        // Make sure max peak is OUTSIDE (above) the desired range
        while (minIdx2 < dbSpec.size() && dbSpec[minIdx2][0] < maxStartMass
            + tolerance + AAJumps::massH2O)
          minIdx2++;

        // Mark all peaks in range as valid starting points in matrix        
        for (int i = minIdx1; i < minIdx2 && i < dbSpec.size(); i++)
        {
          startFlags[specIdx][i] = 1;
        }
      } // for (int specIdx = 0; specIdx < spec.size(); specIdx++)
    } // for (; itr != itrEnd; itr++) {

    // Compute the match between the two spectra
    AMM_peak_match tmp;
    float usedModSize; // Keeps track of the size of the mods used - lower mods are considered better
    matchMatrix.resize(maxNumMods + 1);
    bestSoFar.resize(maxNumMods + 1);

    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {

      matchMatrix[modIdx].resize(dbSpec.size());
      bestSoFar[modIdx].resize(dbSpec.size());

      for (dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++)
      {

        matchMatrix[modIdx][dbSpecIdx].resize(spec.size());
        bestSoFar[modIdx][dbSpecIdx].resize(spec.size());

        for (specIdx = 0; specIdx < spec.size(); specIdx++)
        {

          matchMatrix[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);
          bestSoFar[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);

          //if location is not a valid start, make score infinitely small and skip all other scoring
          if (!startFlags[specIdx][dbSpecIdx])
          {
            for (mPeaksIdx = 0; mPeaksIdx < min(minMatchedPeaks, specIdx + 1); mPeaksIdx++)
            {

              matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(-(float)INT_MAX,
                                                                                modIdx,
                                                                                -1,
                                                                                -1,
                                                                                -1);
              // Update bestSoFar; give precedence to options with same score and less mutations or more matched peaks
              if (modIdx > 0)
                bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx
                    - 1][dbSpecIdx][specIdx][mPeaksIdx]);
              if (specIdx > 0)
                bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx][specIdx
                    - 1][mPeaksIdx]);
              if (dbSpecIdx > 0)
                bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                    - 1][specIdx][mPeaksIdx]);
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score,
                                                                              modIdx,
                                                                              specIdx,
                                                                              dbSpecIdx,
                                                                              0);
            }

            continue;
          }

          // Matching can start on any mass, assume no predecessor (if enforceEnpeaks is false)
          if (!enforceEndpeaks or specIdx == 0)
          {
            matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],
                                                                      modIdx,
                                                                      -1,
                                                                      -1,
                                                                      -1);
          }
          // matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],modIdx,-1,dbSpecIdx-1,-1);

          unsigned int predMPeaksIdx;
          for (mPeaksIdx = 0; mPeaksIdx < min(minMatchedPeaks, specIdx + 1); mPeaksIdx++)
          {
            //cerr<<" >> start ("<<modIdx<<","<<dbSpecIdx<<","<<specIdx<<","<<mPeaksIdx<<"), score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<", peak score = "<<spec[specIdx][1]<<"\n";
            // 0 Mods - find other peaks at the same distance in both spectra
            predSpecIdx = specIdx - 1;
            predDbSpecIdx = dbSpecIdx - 1;
            usedModSize = 0;
            // Zero mods case also needs to search for best predecessor because the previous peaks may have been
            //   missed due to peak tolerance problems. Search for allowable jump with maximum total score
            if (mPeaksIdx > 0 or mPeaksIdx == minMatchedPeaks - 1)
            {
              while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
              {
                // Find an eligible predecessor
                while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and (fabs((spec[specIdx][0] - spec[predSpecIdx][0])
                        - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                        > tolerance))
                {
                  //if(modIdx==0 and dbSpecIdx==11 and specIdx==9)
                  //  cerr<<" -0-->> skipped pred ("<<predSpecIdx<<","<<predDbSpecIdx<<"): fabs("<<(spec[specIdx][0]-spec[predSpecIdx][0])<<","<<(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0])<<")="<<fabs((spec[specIdx][0]-spec[predSpecIdx][0])-(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]))<<" - ";
                  if ((spec[specIdx][0] - spec[predSpecIdx][0])
                      > (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                    predDbSpecIdx--;
                  else
                    predSpecIdx--;
                  //if(modIdx==0 and dbSpecIdx==11 and specIdx==9)
                  //cerr<<"next putative pred is ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";
                }
                //cerr<<" -0- pred zero-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), cur score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score;
                //if (mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                //else 
                //  cerr<<" (no pred found)\n";
                // Use predecessor if score increases
                if (mPeaksIdx > 0 and predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                - 1].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                                                                            - 1].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx
                                                                                        - 1);
                  //DEBUG_MSG(modIdx << "," << dbSpecIdx << "," << specIdx << "," << mPeaksIdx << " = " << spec[specIdx][1]
                  //                                                                      + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                  //                                                                          - 1].score);
                }

                //cerr<<" -0- pred zero-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), mPeaksIdx = "<<mPeaksIdx<<", minMatchedPeaks = "<<minMatchedPeaks;
                //if (mPeaksIdx==minMatchedPeaks-1 and mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                //else cerr<<" (no pred found)\n";

                if (mPeaksIdx == minMatchedPeaks - 1 and predSpecIdx >= 0
                    and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx);
                  //DEBUG_MSG(modIdx << "," << dbSpecIdx << "," << specIdx << "," << mPeaksIdx << " = " << spec[specIdx][1]
                  //                                                                      + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score);
                  //cerr<<" -0--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                }
                predSpecIdx--; //predDbSpecIdx--;
              }

#if 1
              // 1+ Mods - any pair of peaks with one less mutation and within mod mass limits is fair game
              float maxPredScore = -3.4e37, tmpScore, tmpModSize;
              if (modIdx > 0 and dbSpecIdx > 0 and specIdx > 0)
              {
                if (mPeaksIdx > 0)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx
                                - 1].score);
                if (mPeaksIdx == minMatchedPeaks - 1)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx].score);
                if (maxPredScore + spec[specIdx][1]
                    > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score)
                {

                  // Find first eligible predecessor
                  predSpecIdx = specIdx - 1;
                  predDbSpecIdx = dbSpecIdx - 1;
                  // Find first eligible predSpecIdx
                  while (predSpecIdx >= 0 && spec[specIdx][0]
                      - spec[predSpecIdx][0] < minSpecDist - 2 * tolerance)
                    predSpecIdx--;
                  // Find first eligible predSpecIdx/predDbSpecIdx pair
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                      and ((spec[specIdx][0] - spec[predSpecIdx][0]
                          - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                          > maxDbSpecMod + tolerance or (spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                          - tolerance))
                  {
                    if (predDbSpecIdx > 0 and spec[specIdx][0]
                        - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0]) > minDbSpecMod - tolerance)
                      predDbSpecIdx--;
                    else
                    {
                      predSpecIdx--;
                      predDbSpecIdx = dbSpecIdx - 1;
                    }
                  }
                  // Initialize with first eligible predecessor pair
                  //                if(predSpecIdx>=0 and predDbSpecIdx>=0) {
                  //                  matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1]+matchMatrix[modIdx-1][predDbSpecIdx][predSpecIdx].score,modIdx-1,predSpecIdx,predDbSpecIdx);
                  //                  usedModSize = fabs(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]-(spec[specIdx][0]-spec[predSpecIdx][0]));
                  //                }

                  // Try all eligible predecessors
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
                  {
                    if (mPeaksIdx > 0)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                  - 1].score;
                      predMPeaksIdx = mPeaksIdx - 1;
                    }
                    if (mPeaksIdx == minMatchedPeaks - 1
                        and tmpScore < spec[specIdx][1] + matchMatrix[modIdx
                            - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score;
                      predMPeaksIdx = mPeaksIdx;
                    }
                    tmpModSize = fabs(dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0] - (spec[specIdx][0]
                        - spec[predSpecIdx][0]));

                    //cerr<<" -1- pred ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";
                    //cerr<<" -1- pred one-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<"), cur score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score;
                    //if (mPeaksIdx>0 and predSpecIdx>=0 and predDbSpecIdx>=0)
                    //  cerr<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<", score-with-pred = "<<spec[specIdx][1]+matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx-1].score<<"\n";
                    //else cerr<<" (no pred found)\n";

                    if (tmpScore
                        > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        or (matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                            - tmpScore < 0.0001 and tmpModSize < usedModSize))
                    {
                      matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(tmpScore,
                                                                                        modIdx
                                                                                            - 1,
                                                                                        predSpecIdx,
                                                                                        predDbSpecIdx,
                                                                                        predMPeaksIdx);
                      usedModSize = tmpModSize;
                      //cerr<<" -1--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                    }
                    // Find next eligible predecessor pair
                    predDbSpecIdx--;
                    while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                        and ((spec[specIdx][0] - spec[predSpecIdx][0]
                            - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                            > maxDbSpecMod + tolerance or (spec[specIdx][0]
                            - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                            - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                            - tolerance))
                      if (predDbSpecIdx > 0 and spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0]) > minDbSpecMod
                          - tolerance)
                        predDbSpecIdx--;
                      else
                      {
                        predSpecIdx--;
                        predDbSpecIdx = dbSpecIdx - 1;
                      }
                  }
                } // 1+ mods, eligible predecessor
              } // 1+ mods
#endif

            } // mPeaksIdx>0 or mPeaksIdx==minMatchedPeaks-1

            // Update bestSoFar; give precedence to options with same score and less mutations or more matched peaks
            if (modIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx
                  - 1][dbSpecIdx][specIdx][mPeaksIdx]);
            if (specIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx][specIdx
                  - 1][mPeaksIdx]);
            if (dbSpecIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                  - 1][specIdx][mPeaksIdx]);
            bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score,
                                                                            modIdx,
                                                                            specIdx,
                                                                            dbSpecIdx,
                                                                            mPeaksIdx);
            //DEBUG_MSG("bestSoFar: " << modIdx << "," << dbSpecIdx << "," << specIdx << "," << mPeaksIdx << " = " << matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score);

            //tmp = matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx];
            //cerr<<" >> end: (mod="<<modIdx<<",db="<<dbSpecIdx<<",spec="<<specIdx<<",mp="<<mPeaksIdx<<"): score = "<<tmp.score<<", pred=("<<tmp.predMods<<","<<tmp.predDbSpecIdx<<","<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<")\n"; cerr.flush();

          } // for mPeaksIdx

        } // for specIdx

        //tmp = bestSoFar[modIdx][dbSpecIdx][specIdx-1][minMatchedPeaks-1];
        //cerr<<"(mod="<<modIdx<<",db="<<dbSpecIdx<<"): score = "<<tmp.score<<", mod = "<<tmp.predMods<<", db = "<<tmp.predDbSpecIdx<<", spec = "<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<"\n"; cerr.flush();

        if (enforceEndpeaks)
        { // Enforce that matches in bestSoFar[][][spec.size()-1] must end at spec.size()-1
          bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]
              = matchMatrix[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1];
          if (modIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx
                - 1][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]);
          if (dbSpecIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                - 1][specIdx - 1][minMatchedPeaks - 1]);
        }
      } // for dbSpecIdx
    } // for modIdx

    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      /*
       // Keep track of topKmatches
       if (bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1][minMatchedPeaks - 1].score > topKmatches.lowestScore[modIdx]) {
       curMatch.setMatch(spec,
       dbSpec,
       bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1][minMatchedPeaks - 1],
       matchMatrix,
       usedModSize);
       topKmatches.addMatch(modIdx, curMatch);
       }
       */
      float
          score = bestSoFar[modIdx][dbSpecIdx - 1][specIdx - 1][minMatchedPeaks
              - 1].score;
      //DEBUG_VAR(score);
      if (score != -(float)INT_MAX && score != -3.4e37)
      {
        // Create a PSM for this match      
        psmPtr p(new PeptideSpectrumMatch);
        p->m_score = score;
        p->m_dbIndex = dbIndex;
        p->m_matchOrientation = matchOrientation;
        p->m_spectrum = &spec;

        // Copy the peak matches to the PSM
        AMM_peak_match & matchPtr = bestSoFar[modIdx][dbSpecIdx - 1][specIdx
            - 1][minMatchedPeaks - 1];
        while (matchPtr.predSpecIdx >= 0)
        {
          TwoValues<int> tmpIndices;
          tmpIndices[0] = matchPtr.predSpecIdx;
          tmpIndices[1] = matchPtr.predDbSpecIdx;
          p->m_matchedPeaks.push_back(tmpIndices);
          matchPtr
              = matchMatrix[matchPtr.predMods][matchPtr.predDbSpecIdx][matchPtr.predSpecIdx][matchPtr.predMPeaks];
        }
        //DEBUG_VAR(p->m_matchedPeaks.size());
        reverse(p->m_matchedPeaks.begin(), p->m_matchedPeaks.end());
        if (p->m_matchedPeaks.size() != 0) {
          p->m_startMass = dbSpec[p->m_matchedPeaks[0][1]][0];
          spec.psmList.push_back(p);
        }
      }

    } // for (modIdx = 0; modIdx <= maxNumMods; modIdx++)

    return;
  }

  // -------------------------------------------------------------------------
  void scoreOverlapAMMme(Spectrum &spec,
                         Spectrum &dbSpec,
                         short maxNumMods,
                         int minMatchedPeaks,
                         AMM_match_spec &topKmatches,
                         AMM_match &curMatch,
                         vector<vector<AMM_match> > *prevMaxMatch,
                         vector<vector<AMM_match> > *curMaxMatch,
                         float pmTolerance,
                         float tolerance,
                         float minSpecDist,
                         float maxDbSpecMod,
                         float minDbSpecMod,
                         bool enforceEndpeaks)
  {
    int specIdx, dbSpecIdx, modIdx, mPeaksIdx, predSpecIdx, predDbSpecIdx,
        predMPeaksIdx;
    if (spec.size() == 0 or spec.size() < minMatchedPeaks)
      return;

    // matchMatrix[i,j,k,p] where
    //  i - number of used mod/mut mass offsets
    //  j - index in dbSpec
    //  k - index in spec
    //  p - at least p+1 matched peaks in spec
    vector<vector<vector<vector<AMM_peak_match> > > > matchMatrix, // Each [i,j,k,p] position contains the score and predecessors for the highest-scoring
        //   match up to and including masses j,k in dbSpec and spec, respectively, with i mod/mut
        //   mass offsets and at least p matched peaks
        bestSoFar; // Each [i,j,k,p] position contains the coordinates of the highest scoring match
    //   in matchMatrix with coordinates up to and possibly including [i,j,k,p]

    // Compute the match between the two spectra
    AMM_peak_match tmp;
    float usedModSize; // Keeps track of the size of the mods used - lower mods are considered better
    matchMatrix.resize(maxNumMods + 1);
    bestSoFar.resize(maxNumMods + 1);
    for (modIdx = 0; modIdx <= maxNumMods; modIdx++)
    {
      matchMatrix[modIdx].resize(dbSpec.size());
      bestSoFar[modIdx].resize(dbSpec.size());
      for (dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++)
      {
        matchMatrix[modIdx][dbSpecIdx].resize(spec.size());
        bestSoFar[modIdx][dbSpecIdx].resize(spec.size());

        for (specIdx = 0; specIdx < spec.size(); specIdx++)
        {
          matchMatrix[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);
          bestSoFar[modIdx][dbSpecIdx][specIdx].resize(minMatchedPeaks);
          // Line below: Matching can start on any mass, assume no predecessor (if enforceEnpeaks is false)
          if (!enforceEndpeaks or specIdx == 0)
          {
            if (prevMaxMatch and (*prevMaxMatch)[modIdx][specIdx].matchScore
                > spec[specIdx][1] + 0.0001)
            {
              matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax((*prevMaxMatch)[modIdx][specIdx].matchScore,
                                                                        modIdx,
                                                                        -1,
                                                                        -1,
                                                                        -1,
                                                                        &(*prevMaxMatch)[modIdx][specIdx]);
              cerr << "[mod=" << modIdx << ",spec=" << specIdx
                  << "] selected pred <"
                  << (*prevMaxMatch)[modIdx][specIdx].exonAllele[0] << ","
                  << (*prevMaxMatch)[modIdx][specIdx].exonAllele[1]
                  << "> with score "
                  << (*prevMaxMatch)[modIdx][specIdx].matchScore
                  << ", predExon = "
                  << (*prevMaxMatch)[modIdx][specIdx].predExon << endl;
            }
            else
              matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],
                                                                        modIdx,
                                                                        -1,
                                                                        -1,
                                                                        -1);
            //          matchMatrix[modIdx][dbSpecIdx][specIdx][0].replaceWithMax(spec[specIdx][1],modIdx,-1,dbSpecIdx-1,-1);
          }

          unsigned int predMPeaksIdx;
          for (mPeaksIdx = 0; mPeaksIdx < min(minMatchedPeaks, specIdx + 1); mPeaksIdx++)
          {
            //cerr<<" >> start ("<<modIdx<<","<<dbSpecIdx<<","<<specIdx<<","<<mPeaksIdx<<"), score = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
            // 0 Mods - find other peaks at the same distance in both spectra
            predSpecIdx = specIdx - 1;
            predDbSpecIdx = dbSpecIdx - 1;
            usedModSize = 0;
            // Zero mods case also needs to search for best predecessor because the previous peaks may have been
            //   missed due to peak tolerance problems. Search for allowable jump with maximum total score
            if (mPeaksIdx > 0 or mPeaksIdx == minMatchedPeaks - 1)
            {
              while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
              {
                // Find an eligible predecessor
                while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and (fabs((spec[specIdx][0] - spec[predSpecIdx][0])
                        - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                        > tolerance))
                {
                  if ((spec[specIdx][0] - spec[predSpecIdx][0])
                      > (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                    predDbSpecIdx--;
                  else
                    predSpecIdx--;
                }
                //cerr<<" -0- pred zero-mod ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";
                // Use predecessor if score increases
                if (mPeaksIdx > 0 and predSpecIdx >= 0 and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                - 1].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                                                                            - 1].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx
                                                                                        - 1);
                }
                if (mPeaksIdx == minMatchedPeaks - 1 and predSpecIdx >= 0
                    and predDbSpecIdx >= 0
                    and matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        < spec[specIdx][1]
                            + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                {
                  matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(spec[specIdx][1]
                                                                                        + matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score,
                                                                                    modIdx,
                                                                                    predSpecIdx,
                                                                                    predDbSpecIdx,
                                                                                    mPeaksIdx);
                  //cerr<<" -0--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                }
                predSpecIdx--; //predDbSpecIdx--;
              }

              // 1+ Mods - any pair of peaks with one less mutation and within mod mass limits is fair game
              float maxPredScore = -3.4e37, tmpScore, tmpModSize;
              if (modIdx > 0 and dbSpecIdx > 0 and specIdx > 0)
              {
                if (mPeaksIdx > 0)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx
                                - 1].score);
                if (mPeaksIdx == minMatchedPeaks - 1)
                  maxPredScore
                      = max(maxPredScore,
                            bestSoFar[modIdx - 1][dbSpecIdx - 1][specIdx - 1][mPeaksIdx].score);
                if (maxPredScore + spec[specIdx][1]
                    > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score)
                {

                  // Find first eligible predecessor
                  predSpecIdx = specIdx - 1;
                  predDbSpecIdx = dbSpecIdx - 1;
                  // Find first eligible predSpecIdx
                  while (predSpecIdx >= 0 && spec[specIdx][0]
                      - spec[predSpecIdx][0] < minSpecDist - 2 * tolerance)
                    predSpecIdx--;
                  // Find first eligible predSpecIdx/predDbSpecIdx pair
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                      and ((spec[specIdx][0] - spec[predSpecIdx][0]
                          - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                          > maxDbSpecMod + tolerance or (spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                          - tolerance))
                  {
                    if (predDbSpecIdx > 0 and spec[specIdx][0]
                        - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0]) > minDbSpecMod - tolerance)
                      predDbSpecIdx--;
                    else
                    {
                      predSpecIdx--;
                      predDbSpecIdx = dbSpecIdx - 1;
                    }
                  }
                  // Initialize with first eligible predecessor pair
                  //                if(predSpecIdx>=0 and predDbSpecIdx>=0) {
                  //                  matchMatrix[modIdx][dbSpecIdx][specIdx].replaceWithMax(spec[specIdx][1]+matchMatrix[modIdx-1][predDbSpecIdx][predSpecIdx].score,modIdx-1,predSpecIdx,predDbSpecIdx);
                  //                  usedModSize = fabs(dbSpec[dbSpecIdx][0]-dbSpec[predDbSpecIdx][0]-(spec[specIdx][0]-spec[predSpecIdx][0]));
                  //                }

                  // Try all eligible predecessors
                  while (predSpecIdx >= 0 and predDbSpecIdx >= 0)
                  {
                    //cerr<<" -1- pred ("<<predSpecIdx<<","<<predDbSpecIdx<<")\n";

                    if (mPeaksIdx > 0)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx
                                  - 1].score;
                      predMPeaksIdx = mPeaksIdx - 1;
                    }
                    if (mPeaksIdx == minMatchedPeaks - 1
                        and tmpScore < spec[specIdx][1] + matchMatrix[modIdx
                            - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score)
                    {
                      tmpScore
                          = spec[specIdx][1]
                              + matchMatrix[modIdx - 1][predDbSpecIdx][predSpecIdx][mPeaksIdx].score;
                      predMPeaksIdx = mPeaksIdx;
                    }
                    tmpModSize = fabs(dbSpec[dbSpecIdx][0]
                        - dbSpec[predDbSpecIdx][0] - (spec[specIdx][0]
                        - spec[predSpecIdx][0]));

                    if (tmpScore
                        > matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                        or (matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score
                            - tmpScore < 0.0001 and tmpModSize < usedModSize))
                    {
                      matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(tmpScore,
                                                                                        modIdx
                                                                                            - 1,
                                                                                        predSpecIdx,
                                                                                        predDbSpecIdx,
                                                                                        predMPeaksIdx);
                      usedModSize = tmpModSize;
                      //cerr<<" -1--- mPeaksIdx = "<<mPeaksIdx<<", pred score = "<<matchMatrix[modIdx][predDbSpecIdx][predSpecIdx][mPeaksIdx].score<<", curScore = "<<matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score<<"\n";
                    }
                    // Find next eligible predecessor pair
                    predDbSpecIdx--;
                    while (predSpecIdx >= 0 and predDbSpecIdx >= 0
                        and ((spec[specIdx][0] - spec[predSpecIdx][0]
                            - (dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0]))
                            > maxDbSpecMod + tolerance or (spec[specIdx][0]
                            - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                            - dbSpec[predDbSpecIdx][0])) < minDbSpecMod
                            - tolerance))
                      if (predDbSpecIdx > 0 and spec[specIdx][0]
                          - spec[predSpecIdx][0] - (dbSpec[dbSpecIdx][0]
                          - dbSpec[predDbSpecIdx][0]) > minDbSpecMod
                          - tolerance)
                        predDbSpecIdx--;
                      else
                      {
                        predSpecIdx--;
                        predDbSpecIdx = dbSpecIdx - 1;
                      }
                  }
                } // 1+ mods, eligible predecessor
              } // 1+ mods
            } // mPeaksIdx>0 or mPeaksIdx==minMatchedPeaks-1

            // Update bestSoFar; give precedence to options with same score and less mutations or more matched peaks
            if (modIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx
                  - 1][dbSpecIdx][specIdx][mPeaksIdx]);
            if (specIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx][specIdx
                  - 1][mPeaksIdx]);
            if (dbSpecIdx > 0)
              bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                  - 1][specIdx][mPeaksIdx]);
            //cerr<<"modIdx = "<<modIdx<<", dbSpecIdx = "<<dbSpecIdx<<", specIdx = "<<specIdx<<", mPeaksIdx = "<<mPeaksIdx<<"\n"; cerr.flush();
            bestSoFar[modIdx][dbSpecIdx][specIdx][mPeaksIdx].replaceWithMax(matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx].score,
                                                                            modIdx,
                                                                            specIdx,
                                                                            dbSpecIdx,
                                                                            mPeaksIdx);

            //tmp = matchMatrix[modIdx][dbSpecIdx][specIdx][mPeaksIdx];
            //cerr<<" >> end: (mod="<<modIdx<<",db="<<dbSpecIdx<<",spec="<<specIdx<<",mp="<<mPeaksIdx<<"): score = "<<tmp.score<<", pred=("<<tmp.predMods<<","<<tmp.predDbSpecIdx<<","<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<")\n"; cerr.flush();
          } // for mPeaksIdx
        } // for specIdx

        //tmp = bestSoFar[modIdx][dbSpecIdx][specIdx-1][minMatchedPeaks-1];
        //cerr<<"(mod="<<modIdx<<",db="<<dbSpecIdx<<"): score = "<<tmp.score<<", mod = "<<tmp.predMods<<", db = "<<tmp.predDbSpecIdx<<", spec = "<<tmp.predSpecIdx<<", mPeaks = "<<tmp.predMPeaks<<"\n"; cerr.flush();

        if (enforceEndpeaks)
        { // Enforce that matches in bestSoFar[][][spec.size()-1] must end at spec.size()-1
          bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]
              = matchMatrix[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1];
          if (modIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx
                - 1][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1]);
          if (dbSpecIdx > 0)
            bestSoFar[modIdx][dbSpecIdx][specIdx - 1][minMatchedPeaks - 1].replaceWithMax(bestSoFar[modIdx][dbSpecIdx
                - 1][specIdx - 1][minMatchedPeaks - 1]);
        }
      } // for dbSpecIdx

      // Keep track of top matches per exon type
      AMM_peak_match none;
      if (curMaxMatch)
      {
        for (specIdx = 0; specIdx < spec.size(); specIdx++)
          if (specIdx < minMatchedPeaks - 1)
            (*curMaxMatch)[modIdx][specIdx].setMatch(spec,
                                                     dbSpec,
                                                     none,
                                                     matchMatrix,
                                                     0);
          else if ((*curMaxMatch)[modIdx][specIdx].matchScore
              < bestSoFar[modIdx][dbSpecIdx - 1][specIdx][minMatchedPeaks - 1].score)
          {
            (*curMaxMatch)[modIdx][specIdx].proteinIdx = curMatch.proteinIdx;
            (*curMaxMatch)[modIdx][specIdx].exonAllele = curMatch.exonAllele;
            (*curMaxMatch)[modIdx][specIdx].orientationPRMs
                = curMatch.orientationPRMs;
            (*curMaxMatch)[modIdx][specIdx].setMatch(spec,
                                                     dbSpec,
                                                     bestSoFar[modIdx][dbSpecIdx
                                                         - 1][specIdx][minMatchedPeaks
                                                         - 1],
                                                     matchMatrix,
                                                     usedModSize);
          }
      }

      // Keep track of overall topKmatches
      if (bestSoFar[modIdx][dbSpecIdx - 1][spec.size() - 1][minMatchedPeaks - 1].score
          > topKmatches.lowestScore[modIdx])
      {
        curMatch.setMatch(spec,
                          dbSpec,
                          bestSoFar[modIdx][dbSpecIdx - 1][spec.size() - 1][minMatchedPeaks
                              - 1],
                          matchMatrix,
                          usedModSize);
        cerr << " -- Added curMatch with matchScore = " << curMatch.matchScore
            << ", exonAllele = (" << curMatch.exonAllele[0] << ","
            << curMatch.exonAllele[1] << ") and predExon = "
            << curMatch.predExon << "\n";
        cerr.flush();
        topKmatches.addMatch(modIdx, curMatch);
      }

    } // for modIdx
  }

}

