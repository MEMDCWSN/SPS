/*
 * ExecIntersectingSpecProb.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: aguthals
 */

#include "ExecIntersectingSpecProb.h"
#include "GFTableIntersection.h"
#include "alignment_scoring.h"
#include "AlignmentUtils.h"
#include "FileUtils.h"
#include <time.h>

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

using namespace std;

namespace specnets
{

  typedef tr1::unordered_map<string, set<int> > peptide_index_map;

  void ExecIntersectingSpecProb::generatePSMAligns(const ParameterList &inputParams,
                                                   const AAJumps &inputJumps,
                                                   const SpecSet &inputSpectra,
                                                   const ClusterSet &inputClusters,
                                                   SpectrumPairSet &outputPairs)
  {
    const int minOverlapPeaks = inputParams.getValueInt("MIN_OVERLAP_PEAKS", 3);

    const int minMatchedPeaks = inputParams.getValueInt("MIN_MATCHED_PEAKS", 3);

    const float peakTol = inputParams.getValueFloat("TOLERANCE_PEAK", 0.5);

    const float pmTol = inputParams.getValueFloat("TOLERANCE_PM", 0.1);

    const int minAAOverlap = inputParams.getValueInt("MIN_MATCHED_AA", 7);

    const int topAlignsPerSpec = inputParams.getValueInt("TOP_ALIGNS_PER_SPEC",
                                                         100);

    const double maxSpecProb = inputParams.getValueDouble("MAX_PVALUE", 1e-05);

    const float alignFdr = inputParams.getValueFloat("ALIGN_FDR", 0.01000);

    const string decoyID = inputParams.getValue("DECOY_PROTEIN_ID", "DECOY");

    const bool useAlignGF = inputParams.getValueBool("USE_ALIGNGF", false);

    const double maxAlignGFPval = inputParams.getValueDouble("MAX_ALIGNGF_PVAL",
                                                             1);

    const bool offRankPairs = inputParams.getValueBool("OFFRANK_PAIRS", false);

    DEBUG_VAR(alignFdr);
    DEBUG_VAR(useAlignGF);
    DEBUG_VAR(maxAlignGFPval);

    int debugIdx1 = inputParams.getValueInt("DEBUG_IDX1", -1);
    int debugIdx2 = inputParams.getValueInt("DEBUG_IDX2", -1);

    float normalizedSpecScore =
        inputParams.getValueFloat("NORMALIZED_SPECTRUM_SCORE", -1.0);

    int indexLen = min(minAAOverlap, 7);

    vector<string> peptidesUse(inputSpectra.size(), "");
    vector<int> psmRanks(inputSpectra.size(), -1);
    vector<int> specToClust(inputSpectra.size(), -1);

    for (int c = 0; c < inputClusters.size(); c++)
    {
      for (int r = 0; r < inputClusters[c].size(); r++)
      {
        psmRanks[inputClusters[c][r].m_index] = r + 1;
        specToClust[inputClusters[c][r].m_index] = inputClusters[c].m_index;
      }
    }

    for (int i = 0; i < inputSpectra.size(); i++)
    {
      const Spectrum &spec = inputSpectra[i];
      if (spec.psmList.size() == 0)
      {
        continue;
      }

      const string &peptide = spec.psmList.front()->m_annotation;

      if (specToClust[i] == debugIdx1 || specToClust[i] == debugIdx2)
      {
        DEBUG_VAR(i);
        DEBUG_VAR(specToClust[i]);
        DEBUG_VAR(inputSpectra[i].scan);
        DEBUG_VAR(peptide);
        DEBUG_VAR(spec.psmList.front()->m_pValue);
      }

      if (spec.psmList.front()->m_pValue > maxSpecProb)
      {
        continue;
      }

      peptidesUse[i] = AAJumps::stripMods(peptide);

      //replaceAll(peptidesUse[i], "I", "L");
      //replaceAll(peptidesUse[i], "K", "Q");
    }

    vector<vector<float> > matchedPRMs(inputSpectra.size());

    peptide_index_map tagIndex(inputSpectra.size());
    outputPairs.resize(0);

    SpecSet specsUse(inputSpectra.size());

    float stepProg = 0.1;
    float curProg = 0;

    for (int i = 0; i < inputSpectra.size(); i++)
    {
      /*if (DEBUG_IDX1 >= 0 && DEBUG_IDX2 >= 0 && i != DEBUG_IDX1
       && i != DEBUG_IDX2)
       {
       continue;
       }*/

      float prog = ((float)i) / ((float)inputSpectra.size());
      if (prog >= curProg)
      {
        DEBUG_MSG("Indexing peptide " << indexLen << "-mers ... " << parseFloat(100.0 * prog, 1) << "%");
        curProg += stepProg;
      }

      if (debugIdx1 >= 0 && debugIdx2 >= 0 && specToClust[i] != debugIdx1
          && specToClust[i] != debugIdx2)
      {
        continue;
      }

      const Spectrum &spec = inputSpectra[i];

      const string &peptide = peptidesUse[i];

      if (peptide.length() == 0)
      {
        continue;
      }

      bool debug = false;
      if (specToClust[i] == debugIdx1 || specToClust[i] == debugIdx2)
      {
        DEBUG_VAR(i);
        DEBUG_VAR(inputSpectra[i].scan);
        DEBUG_VAR(peptide);
        //DEBUG_TRACE;
        //spec.output(cerr);
        DEBUG_TRACE;
        debug = true;
      }

      specsUse[i] = spec;
      specsUse[i].setPeakTolerance(0);

      if (normalizedSpecScore > 0
          && spec.getTotalIonCurrent() > normalizedSpecScore)
      {
        specsUse[i].normalize(normalizedSpecScore, 0);
      }

      Spectrum prmMasses;
      inputJumps.getRoundedPRMMasses(peptide, prmMasses);

      if (debug)
      {
        DEBUG_TRACE;
        //prmMasses.output(cerr);
        DEBUG_TRACE;
      }

      matchedPRMs[i].assign(prmMasses.size(), 0);

      vector<int> idxMatched1;
      vector<int> idxMatched2;
      FindMatchPeaksAll2(specsUse[i],
                         prmMasses,
                         0,
                         peakTol,
                         idxMatched1,
                         idxMatched2);

      int numMatchedPks = 0;
      for (int j = 0; j < idxMatched2.size(); j++)
      {
        if (specsUse[i][idxMatched1[j]][0] < 10.0
            || specsUse[i][idxMatched1[j]][0] > specsUse[i].parentMass - 37.0)
        {
          continue;
        }
        numMatchedPks++;
        matchedPRMs[i][idxMatched2[j]] += specsUse[i][idxMatched1[j]][1];
        if (debug)
        {
          DEBUG_MSG("Matched at index = " << idxMatched2[j] << ", mass = " << specsUse[i][idxMatched1[j]][0] << ", score = " << specsUse[i][idxMatched1[j]][1]);
        }
      }
      //matchedPRMs[i][matchedPRMs[i].size() - 1] = 1.0;
      //matchedPRMs[i][0] = 1.0;

      if (debug)
      {
        DEBUG_VAR(minMatchedPeaks);
        DEBUG_VAR(numMatchedPks);
        DEBUG_VAR(peptide.length());
      }

      if (peptide.length() < minAAOverlap || numMatchedPks < minMatchedPeaks)
      {
        continue;
      }

      for (int pos = 0; pos <= peptide.length() - indexLen; pos++)
      {
        string tag = peptide.substr(pos, indexLen);
        if (tagIndex.count(tag) > 0)
        {
          tagIndex[tag].insert(i);
        }
        else
        {
          set<int> newSet;
          newSet.insert(i);
          tagIndex[tag] = newSet;
        }
      }
    }

    DEBUG_MSG("Successfully indexed " << tagIndex.size() << " tags");
    curProg = 0;
    float curTagNum = 0;

    vector<set<int> > seenPairs(inputSpectra.size());

    vector<vector<TwoValues<int> > > scoredPairs(inputSpectra.size());

    vector<TwoValues<float> > rankedPairs(inputSpectra.size());
    int rankPairIdx = 0;

    SpectrumPair nextPair;
    SpectrumPairSet tempPairs;
    for (peptide_index_map::const_iterator tagIt = tagIndex.begin();
        tagIt != tagIndex.end(); tagIt++)
    {
      float prog = curTagNum / ((float)tagIndex.size());
      if (prog >= curProg)
      {
        DEBUG_MSG("Extracting PSM pairs ... " << parseFloat(100.0 * prog, 1) << "%");
        curProg += stepProg;
      }
      curTagNum += 1.0;

      if (tagIt->second.size() < 2)
      {
        continue;
      }
      for (set<int>::const_reverse_iterator idxIt1 = tagIt->second.rbegin();
          idxIt1 != tagIt->second.rend(); idxIt1++)
      {
        int idx1 = *idxIt1;

        /*if (DEBUG_IDX2 >= 0 && idx1 != DEBUG_IDX2)
         {
         continue;
         }*/
        const string &pep1 = peptidesUse[idx1];

        bool pep1IsDecoy =
            inputSpectra[idx1].psmList.front()->m_protein.find(decoyID)
                != string::npos;

        set<int>::const_reverse_iterator idxIt2 = idxIt1;
        idxIt2++;
        for (; idxIt2 != tagIt->second.rend(); idxIt2++)
        {
          int idx2 = *idxIt2;

          /*if (DEBUG_IDX1 >= 0 && idx2 != DEBUG_IDX1)
           {
           continue;
           }*/

          if (specToClust[idx1] == specToClust[idx2]
              || seenPairs[idx1].count(idx2) > 0)
          {
            continue;
          }
          seenPairs[idx1].insert(idx2);

          if ((!offRankPairs) && psmRanks[idx1] != 1 && psmRanks[idx2] != 1)
          {
            continue;
          }

          const string &pep2 = peptidesUse[idx2];

          bool pep2IsDecoy =
              inputSpectra[idx2].psmList.front()->m_protein.find(decoyID)
                  != string::npos;

          bool debug = false;
          if (specToClust[idx1] == debugIdx1 || specToClust[idx1] == debugIdx2)
          {
            DEBUG_VAR(inputSpectra[idx1].scan);
            DEBUG_VAR(inputSpectra[idx2].scan);
            DEBUG_VAR(pep1);
            DEBUG_VAR(pep2);
            debug = true;
          }

          /*
           if (pep1 == pep2
           && inputSpectra[idx1].parentCharge
           == inputSpectra[idx2].parentCharge)
           {
           continue;
           }
           */

          int daShift;
          int resShift;
          float peakShift;

          if (inputJumps.getPeptideShift(pep1,
                                         pep2,
                                         minAAOverlap,
                                         &daShift,
                                         &peakShift,
                                         &resShift))
          {
            int numPks1 = 0, numPks2 = 0, numMatchedPks = 0, recScore1 = 0,
                recScore2 = 0, ovlpScore1 = 0, ovlpScore2 = 0;

            int start1 = (resShift >= 0) ? resShift : 0;
            int start2 = (resShift < 0) ? 0 - resShift : 0;
            int maxOverlap = min(pep1.length() - start1,
                                 pep2.length() - start2);

            if (debug)
              DEBUG_MSG("Have alignment start1=" << start1 << "; start2=" << start2 << "; maxOverlap=" << maxOverlap << "; daShift=" << daShift << "; resShift = " << resShift);

            for (int aaI = 0; aaI < maxOverlap; aaI++)
            {
              int rootIdx1 = aaI + start1;
              int rootIdx2 = aaI + start2;

              float matchScore1 = matchedPRMs[idx1][rootIdx1];
              float matchScore2 = matchedPRMs[idx2][rootIdx2];
              if (matchScore1 > 0)
              {
                if (debug)
                  DEBUG_MSG("1 matches at " << rootIdx1);
                numPks1++;
              }
              if (matchScore2 > 0)
              {
                if (debug)
                  DEBUG_MSG("2 matches at " << rootIdx2);
                numPks2++;
              }
              if (matchScore1 > 0 && matchScore2 > 0)
              {
                if (debug)
                  DEBUG_MSG("1 matches at " << rootIdx1 << " and 2 matches at " << rootIdx2);
                numMatchedPks++;
                ovlpScore1 += floatToInt(matchScore2);
                ovlpScore2 += floatToInt(matchScore1);
              }

              if (matchScore1 <= 0 && matchScore2 > 0)
              {
                recScore1 += floatToInt(matchScore2);
              }
              else if (matchScore1 > 0 && matchScore2 <= 0)
              {
                recScore2 += floatToInt(matchScore1);
              }
            }

            if (debug)
            {
              DEBUG_MSG("numPks1 = " << numPks1 << ", numPks2 = " << numPks2 << ", numMatchedPks = " << numMatchedPks);
            }

            if (numPks1 < minOverlapPeaks || numPks2 < minOverlapPeaks
                || numMatchedPks < minMatchedPeaks)
            {
              continue;
            }

            const Spectrum &spec1 = specsUse[idx1];
            const Spectrum &spec2 = specsUse[idx2];
            const float floatShift = daShift;

            float totalOverlap1 = 0;
            float totalOverlap2 = 0;
            float totalMatched1 = 0;
            float totalMatched2 = 0;

            for (int pIdx = 0; pIdx < spec2.size(); pIdx++)
            {
              const float offsetMass = spec2[pIdx][0] + floatShift;
              if (offsetMass > 0 && offsetMass < spec1[spec1.size() - 1][0])
              {
                totalOverlap2 += spec2[pIdx][1];

                if (spec1.findPeaks(offsetMass, peakTol * 2) >= 0)
                {
                  totalMatched2 += spec2[pIdx][1];
                }
              }
            }

            for (int pIdx = 0; pIdx < spec1.size(); pIdx++)
            {
              const float offsetMass = spec1[pIdx][0] - floatShift;
              if (offsetMass > 0 && offsetMass < spec2[spec2.size() - 1][0])
              {
                totalOverlap1 += spec1[pIdx][1];

                if (spec2.findPeaks(offsetMass, peakTol * 2) >= 0)
                {
                  totalMatched1 += spec1[pIdx][1];
                }
              }
            }

            /*
             if (idx1 == 385256 && idx2 == 215075)
             {
             DEBUG_MSG("(" << idx1 << "," << idx2 << ") shift=" << daShift << " peps=( " << pep1 << " , " << pep2 << " )");
             }
             */

            if (debug)
              DEBUG_MSG("Have pair " << idx1 << " , " << idx2);

            recScore1 = (recScore1 * 1000) + ovlpScore1;
            recScore2 = (recScore2 * 1000) + ovlpScore2;

            /*
             if (pep1 == pep2
             && inputSpectra[idx1].parentCharge
             == inputSpectra[idx2].parentCharge)
             {
             continue;
             }
             */

            nextPair.spec1 = idx1;
            nextPair.spec2 = idx2;
            nextPair.shift1 = daShift;
            nextPair.shift2 = peakShift;
            nextPair.score1 = totalMatched1 / totalOverlap1; //float(recScore1);
            nextPair.score2 = totalMatched2 / totalOverlap2; //float(recScore2);
            nextPair.specC = (pep1IsDecoy || pep2IsDecoy) ? 1 : 0;

            //nextPair.spec2rev = false;
            const int pairIdx = tempPairs.size();

            tempPairs.push_back(nextPair);

            if (rankPairIdx >= rankedPairs.size())
            {
              rankedPairs.resize(rankPairIdx + 10000);
            }

            const float alignmentScore = ((float)numMatchedPks)
                * min(totalMatched1 / totalOverlap1,
                      totalMatched2 / totalOverlap2);

            scoredPairs[idx1].push_back(TwoValues<int>(floatToInt(alignmentScore
                                                           * 1000),
                                                       pairIdx));
            scoredPairs[idx2].push_back(TwoValues<int>(floatToInt(alignmentScore
                                                           * 1000),
                                                       pairIdx));

            /*
             const float alignmentScore = 1.0
             - max(inputSpectra[idx1].psmList.front()->m_pValue,
             inputSpectra[idx2].psmList.front()->m_pValue);
             */

            TwoValues<float> rankPr(0 - alignmentScore, (float)pairIdx);
            rankedPairs[rankPairIdx++] = rankPr;
          }
          else
          {
            if (debug)
              DEBUG_MSG("Skipped alignment");
          }
        }
      }
    }
    rankedPairs.resize(rankPairIdx);

    /*
     if (useAlignGF)
     {
     DEBUG_TRACE;

     vector<TwoValues<double> > pvalues;
     SpecSet specsCopy;
     specsCopy = inputSpectra;

     DEBUG_VAR(tempPairs.size());

     getPairAlignGFPValues22(specsCopy,
     tempPairs,
     pvalues,
     pmTol,
     peakTol,
     inputJumps,
     false);

     DEBUG_VAR(pvalues.size());

     OutputTable pairTable;
     int row = 0, col = 0;

     const string outputPvalFile = inputParams.getValue("OUTPUT_ALIGNGF_PVALS", "");

     if (outputPvalFile.length() > 0)
     {
     pairTable.values.resize(pvalues.size()+1);

     pair<string, bool> defaultHeader("", true);
     pairTable.values[row].assign(14, defaultHeader);
     pairTable.values[row][col++].first = "Scan 1";
     pairTable.values[row][col++].first = "Scan 2";
     pairTable.values[row][col++].first = "Charge 1";
     pairTable.values[row][col++].first = "Charge 2";
     pairTable.values[row][col++].first = "PSM Rank 1";
     pairTable.values[row][col++].first = "PSM Rank 2";
     pairTable.values[row][col++].first = "Peptide 1";
     pairTable.values[row][col++].first = "Peptide 2";
     pairTable.values[row][col++].first = "SpecProb 1";
     pairTable.values[row][col++].first = "SpecProb 2";
     pairTable.values[row][col++].first = "AlignGF Pval 1";
     pairTable.values[row][col++].first = "AlignGF Pval 2";
     pairTable.values[row][col++].first = "Max Pval";
     pairTable.values[row][col++].first = "Align Score";
     col = 0;
     row++;
     }

     pair<string, bool> defaultStat("", false);
     for (int i = 0; i < rankedPairs.size(); i++)
     {
     int pairIdx = floatToInt(rankedPairs[i][1]);
     const double alignGFPval = max(pvalues[pairIdx][0], pvalues[pairIdx][1]);

     if (outputPvalFile.length() > 0)
     {
     pairTable.values[row].assign(14, defaultStat);
     pairTable.values[row][col++].first = parseInt(specsUse[tempPairs[pairIdx].spec1].scan);
     pairTable.values[row][col++].first = parseInt(specsUse[tempPairs[pairIdx].spec2].scan);
     pairTable.values[row][col++].first = parseInt(specsUse[tempPairs[pairIdx].spec1].parentCharge);
     pairTable.values[row][col++].first = parseInt(specsUse[tempPairs[pairIdx].spec2].parentCharge);
     pairTable.values[row][col++].first = parseInt(psmRanks[tempPairs[pairIdx].spec1]);
     pairTable.values[row][col++].first = parseInt(psmRanks[tempPairs[pairIdx].spec2]);
     pairTable.values[row][col++].first = peptidesUse[tempPairs[pairIdx].spec1];
     pairTable.values[row][col++].first = peptidesUse[tempPairs[pairIdx].spec2];
     pairTable.values[row][col++].first = parseDoubleSci(specsUse[tempPairs[pairIdx].spec1].psmList.front()->m_pValue, 4);
     pairTable.values[row][col++].first = parseDoubleSci(specsUse[tempPairs[pairIdx].spec2].psmList.front()->m_pValue, 4);
     pairTable.values[row][col++].first = parseDoubleSci(pvalues[pairIdx][0], 4);
     pairTable.values[row][col++].first = parseDoubleSci(pvalues[pairIdx][1], 4);
     pairTable.values[row][col++].first = parseDoubleSci(alignGFPval, 4);
     pairTable.values[row][col++].first = parseFloat(rankedPairs[i][0], 4);
     col = 0;
     row++;
     }

     if (alignGFPval > maxAlignGFPval || alignGFPval == 0)
     {
     rankedPairs[i][0] = 1001.0;
     continue;
     }

     rankedPairs[i][0] = alignGFPval;
     }

     if (outputPvalFile.length() > 0)
     {
     pairTable.values.resize(row);

     if (!pairTable.printToCSV(outputPvalFile.c_str(), "\t"))
     {
     ERROR_MSG("Failed to save AlignGF p-values to \'" << outputPvalFile << "\'");
     }
     }
     }
     */

    DEBUG_VAR(tempPairs.size());
    DEBUG_MSG("Filtering pairs by " << alignFdr << " FDR ...");
    sort(rankedPairs.begin(), rankedPairs.end());

    float numTP = 0.001, numFP = 0;
    int rankCutoff = rankedPairs.size() - 1;

    float numSamePep = 0, numSameCharge = 0;
    for (int i = 0; i < rankedPairs.size(); i++)
    {
      int pairIdx = floatToInt(rankedPairs[i][1]);

      if (rankedPairs[i][0] > 1000.0)
      {
        break;
      }

      if (tempPairs[pairIdx].specC > 0)
      {
        numFP += 1.0;

        int idx1 = tempPairs[pairIdx].spec1;
        int idx2 = tempPairs[pairIdx].spec2;

        if (inputSpectra[idx1].psmList.front()->m_annotation
            == inputSpectra[idx2].psmList.front()->m_annotation)
        {
          numSamePep += 1.0;

          if (inputSpectra[idx1].parentCharge
              == inputSpectra[idx2].parentCharge)
          {
            numSameCharge += 1.0;
          }
        }

        //DEBUG_VAR(tempPairs[pairIdx].specC);
        //DEBUG_MSG("False-positive pair ( pep1=" << inputSpectra[idx1].psmList.front()->m_annotation << ", pval=" << rankedPairs[i][0]);
        //DEBUG_MSG("False-positive pair ( pep2=" << inputSpectra[idx2].psmList.front()->m_annotation << ", pval=" << rankedPairs[i][0]);
      }
      else
      {
        numTP += 1.0;
      }

      rankCutoff = i;
      if ((numFP / numTP) > alignFdr)
      {
        break;
      }
    }

    DEBUG_MSG("Found " << rankCutoff+1 << " pairs at " << (numFP / numTP) << " FDR ...");

    DEBUG_VAR(rankedPairs.size());
    DEBUG_VAR(numFP);
    DEBUG_VAR(numSamePep);
    DEBUG_VAR(numSameCharge);

    //abort();

    for (int i = rankCutoff + 1; i < rankedPairs.size(); i++)
    {
      int pairIdx = floatToInt(rankedPairs[i][1]);
      tempPairs[pairIdx].specC = 0;
      tempPairs[pairIdx].score1 = rankedPairs[i][0];
      tempPairs[pairIdx].score2 = rankedPairs[i][0];
    }

    for (int i = 0; i <= rankCutoff; i++)
    {
      int pairIdx = floatToInt(rankedPairs[i][1]);
      tempPairs[pairIdx].specC = 1;
      tempPairs[pairIdx].score1 = rankedPairs[i][0];
      tempPairs[pairIdx].score2 = rankedPairs[i][0];
    }

    DEBUG_MSG("Extracting top " << topAlignsPerSpec << " pairs per spectrum ...");

    vector<vector<TwoValues<int> > > sortedPairs(inputSpectra.size());

    for (int i = 0; i < scoredPairs.size(); i++)
    {
      sort(scoredPairs[i].begin(), scoredPairs[i].end());

      int pairCount = 0;
      //DEBUG_VAR(i);
      for (int j = scoredPairs[i].size() - 1; j >= 0; j--)
      {
        if (pairCount >= topAlignsPerSpec)
        {
          //DEBUG_VAR(scoredPairs[i][j][0]);
          break;
        }
        if (tempPairs[scoredPairs[i][j][1]].specC > 0)
        {
          //DEBUG_VAR(scoredPairs[i][j][0]);
          int pairIdx = scoredPairs[i][j][1];
          int spec1 = tempPairs[pairIdx].spec1;
          int spec2 = tempPairs[pairIdx].spec2;

          if (spec1 < spec2)
          {
            tempPairs[pairIdx].spec1 = spec2;
            tempPairs[pairIdx].spec2 = spec1;
            tempPairs[pairIdx].shift1 = 0 - tempPairs[pairIdx].shift1;

            float temp = tempPairs[pairIdx].score1;
            tempPairs[pairIdx].score1 = tempPairs[pairIdx].score2;
            tempPairs[pairIdx].score2 = temp;
          }

          sortedPairs[tempPairs[pairIdx].spec1].push_back(TwoValues<int>(tempPairs[pairIdx].spec2,
                                                                         pairIdx));
          tempPairs[scoredPairs[i][j][1]].specC = 0;
        }
        pairCount++;
      }
    }

    DEBUG_MSG("Sorting pairs ...");

    for (int i = 0; i < sortedPairs.size(); i++)
    {
      if (sortedPairs[i].size() == 0)
      {
        continue;
      }

      sort(sortedPairs[i].begin(), sortedPairs[i].end());

      for (int j = 0; j < sortedPairs[i].size(); j++)
      {
        outputPairs.push_back(tempPairs[sortedPairs[i][j][1]]);
      }
    }
    DEBUG_VAR(outputPairs.size());

    /*string pairsFilename = "input_pairs_dancik.bin";
     if (!m_inputPairs->saveToBinaryFile(pairsFilename))
     {
     ERROR_MSG("Failed to save pairs to \'" << pairsFilename << "\'");
     }*/
  }
  //const int ExecIntersectingSpecProb::DEBUG_IDX1 = -1; //80537;
  //const int ExecIntersectingSpecProb::DEBUG_IDX2 = -1; //837986;

  ExecIntersectingSpecProb::ExecIntersectingSpecProb(void) :
      ExecBase(), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputPSMs(0x0), m_inputPairs(0x0), m_inputClusters(0x0), m_inputErrorModelTarget(0x0), m_inputErrorModelDecoy(0x0), m_inputMsSpectra(0x0), m_outputPairs(0x0), m_outputProbs(0x0), m_jumps(0x0)
  {
    m_name = "ExecIntersectingSpecProb";
    m_type = "ExecIntersectingSpecProb";
    m_inputSpectra = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputPairs = new SpectrumPairSet;
    m_inputClusters = new ClusterSet;
    m_inputMsSpectra = new SpecSet;
    m_inputErrorModelTarget = new MassErrorModel;
    m_inputErrorModelDecoy = new MassErrorModel;
    m_outputPairs = new SpectrumPairSet;
    m_outputProbs = new vector<vector<double> >;
    m_jumps = new AAJumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);
    m_jumps->multiplyMasses(AA_ROUNDING);
  }

  ExecIntersectingSpecProb::ExecIntersectingSpecProb(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputPairs(0x0), m_inputPSMs(0x0), m_inputClusters(0x0), m_inputErrorModelTarget(0x0), m_inputErrorModelDecoy(0x0), m_inputMsSpectra(0x0), m_outputPairs(0x0), m_outputProbs(0x0), m_jumps(0x0)
  {
    m_name = "ExecIntersectingSpecProb";
    m_type = "ExecIntersectingSpecProb";
    m_inputSpectra = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputPairs = new SpectrumPairSet;
    m_inputClusters = new ClusterSet;
    m_inputMsSpectra = new SpecSet;
    m_inputErrorModelTarget = new MassErrorModel;
    m_inputErrorModelDecoy = new MassErrorModel;
    m_outputPairs = new SpectrumPairSet;
    m_outputProbs = new vector<vector<double> >;
    m_jumps = new AAJumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);
    m_jumps->multiplyMasses(AA_ROUNDING);
  }

  ExecIntersectingSpecProb::ExecIntersectingSpecProb(const ParameterList &inputParams,
                                                     SpecSet *inputSpectra,
                                                     PeptideSpectrumMatchSet *inputPSMs,
                                                     SpectrumPairSet *inputPairs,
                                                     ClusterSet* inputClusters,
                                                     MassErrorModel *inputErrorModelTarget,
                                                     MassErrorModel *inputErrorModelDecoy,
                                                     SpecSet *inputMsSpectra,
                                                     AAJumps *jumps) :
      ExecBase(inputParams), ownInput(false), ownOutput(true), m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs), m_inputPairs(inputPairs), m_inputClusters(inputClusters), m_inputErrorModelTarget(inputErrorModelTarget), m_inputErrorModelDecoy(inputErrorModelDecoy), m_inputMsSpectra(inputMsSpectra), m_outputPairs(0x0), m_outputProbs(0x0), m_jumps(jumps)
  {
    m_name = "ExecIntersectingSpecProb";
    m_type = "ExecIntersectingSpecProb";
    m_outputPairs = new SpectrumPairSet;
    m_outputProbs = new vector<vector<double> >;
  }

  ExecIntersectingSpecProb::ExecIntersectingSpecProb(const ParameterList &inputParams,
                                                     SpecSet *inputSpectra,
                                                     PeptideSpectrumMatchSet *inputPSMs,
                                                     SpectrumPairSet *inputPairs,
                                                     ClusterSet* inputClusters,
                                                     MassErrorModel *inputErrorModelTarget,
                                                     MassErrorModel *inputErrorModelDecoy,
                                                     SpecSet *inputMsSpectra,
                                                     AAJumps *jumps,
                                                     SpectrumPairSet *outputPairs,
                                                     vector<vector<double> > *outputProbs) :
      ExecBase(inputParams), ownInput(false), ownOutput(false), m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs), m_inputPairs(inputPairs), m_inputClusters(inputClusters), m_inputErrorModelTarget(inputErrorModelTarget), m_inputErrorModelDecoy(inputErrorModelDecoy), m_inputMsSpectra(inputMsSpectra), m_outputPairs(outputPairs), m_outputProbs(outputProbs), m_jumps(jumps)
  {
    m_name = "ExecIntersectingSpecProb";
    m_type = "ExecIntersectingSpecProb";
  }

  ExecIntersectingSpecProb::~ExecIntersectingSpecProb(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
      delete m_inputPSMs;
      delete m_inputPairs;
      delete m_inputClusters;
      delete m_inputErrorModelTarget;
      delete m_inputErrorModelDecoy;
      delete m_inputMsSpectra;
      delete m_jumps;
    }
    if (ownOutput)
    {
      delete m_outputPairs;
      delete m_outputProbs;
    }
  }

  ExecBase * ExecIntersectingSpecProb::clone(const ParameterList & inputParams) const
  {
    return new ExecIntersectingSpecProb(inputParams);
  }

  bool ExecIntersectingSpecProb::invoke(void)
  {

    //double minPValueCompute = m_params.getValueDouble("MIN_PVALUE_COMPUTE", -1.0);

    const int numComputationsPerCompress = 60;

    const float peakTol = 0.5; //m_params.getValueFloat("TOLERANCE_PEAK");
    const float pmTol = m_params.getValueFloat("TOLERANCE_PM");

    int pairIdxEnd = m_params.getValueInt("PAIR_IDX_END",
                                          ((int)m_inputPairs->size()) - 1);
    int pairIdxStart = m_params.getValueInt("PAIR_IDX_START", 0);

    const int debugIdx1 = m_params.getValueInt("DEBUG_IDX1", -1);
    const int debugIdx2 = m_params.getValueInt("DEBUG_IDX2", -1);

    const bool verbose = m_params.getValueBool("VERBOSE", false);

    const bool mergePrms = m_params.getValueBool("MERGE_PRMS", false);

    const bool mergeAlignedPrms = m_params.getValueBool("MERGE_ALIGNED_PRMS",
                                                        false);

    const double minPvalCompute = m_params.getValueDouble("MIN_PVAL_COMPUTE",
                                                          -1.0);

    const float ppmEdgeError = m_params.getValueFloat("TOLERANCE_PEAK_PPM")
        * 2.0;

    const float normalizedSpectrumScore =
        m_params.getValueFloat("NORMALIZED_SPECTRUM_SCORE", -1.0);

    const bool ignoreMS2ErrorModel = (m_inputErrorModelTarget->size() == 0
        || m_inputErrorModelDecoy->size() == 0 || m_inputMsSpectra->size() == 0);

    Spectrum emptySpectrum;

    if (m_inputPairs->size() == 0)
    {
      generatePSMAligns(m_params,
                        *m_jumps,
                        *m_inputSpectra,
                        *m_inputClusters,
                        *m_inputPairs);
      pairIdxEnd = ((int)m_inputPairs->size()) - 1;
      pairIdxStart = 0;
    }

#pragma omp critical
    {
      DEBUG_VAR(m_inputPairs->size());
      DEBUG_VAR(normalizedSpectrumScore);
      DEBUG_VAR(mergePrms);
      DEBUG_VAR(mergeAlignedPrms);
      DEBUG_VAR(minPvalCompute);
      DEBUG_VAR(pairIdxStart);
      DEBUG_VAR(pairIdxEnd);
    }

    int lastIdx1 = -1;
    GFTable tab1;
    GFTable tab2;
    GFTableIntersection intersection;

    m_outputPairs->resize(0);
    m_outputProbs->resize(0);

    clock_t init, final;
    int numComputed = 0;
    init = clock();

    int lastProg = -1;

    int matchScore1, matchScore2, matchScoreOvlp;

    bool debug = debugIdx1 >= 0 && debugIdx2 >= 0;

    Spectrum specUse1;
    Spectrum specUse2;

    vector<int> specToClust(m_inputSpectra->size(), -1);
    for (int i = 0; i < m_inputClusters->size(); i++)
    {

      for (int j = 0; j < (*m_inputClusters)[i].size(); j++)
      {
        specToClust[(*m_inputClusters)[i][j].m_index] = i;
      }
    }

    vector<int> msScanToIdx(m_inputMsSpectra->size(), -1);
    for (int i = 0; i < m_inputMsSpectra->size(); i++)
    {
      const int msScan = (*m_inputMsSpectra)[i].scan;

      if (msScan >= msScanToIdx.size())
      {
        msScanToIdx.resize(msScan + 1, -1);
      }

      if (msScanToIdx[msScan] > -1)
      {
        ERROR_MSG("Duplicate scans " << msScanToIdx[msScan]);
        abort();
      }
      msScanToIdx[msScan] = i;
    }

    double starP1, starP2;

    vector<double> starProbs(m_inputSpectra->size(), 1.0);

    for (int i = pairIdxEnd; i >= pairIdxStart && i >= 0; i--)
    {
      int idx1 = (*m_inputPairs)[i].spec1;
      int idx2 = (*m_inputPairs)[i].spec2;

      starP1 = starProbs[idx1];
      starP2 = starProbs[idx2];

      if (minPvalCompute > 0 && starP1 < minPvalCompute
          && starP2 < minPvalCompute)
      {
        continue;
      }

      const Spectrum &spec1 = (*m_inputSpectra)[idx1];
      const Spectrum &spec2 = (*m_inputSpectra)[idx2];

      if (spec1.size() == 0 || spec2.size() == 0)
      {
        continue;
      }
      string pepUse1 = (*m_inputSpectra)[idx1].psmList.front()->m_annotation;
      string pepUse2 = (*m_inputSpectra)[idx2].psmList.front()->m_annotation;

      /*
       replaceAll(pepUse1, "I", "L");
       replaceAll(pepUse1, "K", "Q");

       replaceAll(pepUse2, "I", "L");
       replaceAll(pepUse2, "K", "Q");
       */

      float shiftUse = (*m_inputPairs)[i].shift1;
      int intShift = floatToInt((*m_inputPairs)[i].shift1);

      if (debug)
      {
        DEBUG_MSG("Pair (" << idx1 << " , " << idx2 << "), scans = ( " << (*m_inputSpectra)[idx1].scan << " , " << (*m_inputSpectra)[idx2].scan << " )");
        DEBUG_VAR(pepUse1);
        DEBUG_VAR(pepUse2);
        DEBUG_VAR(intShift);
        DEBUG_VAR(spec1.psmList.front()->m_score);
        DEBUG_VAR(spec1.psmList.front()->m_pValue);
        DEBUG_VAR(spec2.psmList.front()->m_score);
        DEBUG_VAR(spec2.psmList.front()->m_pValue);
      }

      specUse1 = spec1;

      if (normalizedSpectrumScore > 0
          && specUse1.getTotalIonCurrent() > normalizedSpectrumScore)
      {
        specUse1.normalize(normalizedSpectrumScore);
      }

      specUse2 = spec2;

      if (normalizedSpectrumScore > 0
          && specUse2.getTotalIonCurrent() > normalizedSpectrumScore)
      {
        specUse2.normalize(normalizedSpectrumScore);
      }

      /*
       tab1.initialize(specUse1, *m_jumps, specUse1);
       //tab1.addMatchingPRMScores(intShift, specUse2);
       tab2.initialize(specUse2, *m_jumps, specUse2);
       //tab2.addMatchingPRMScores(0 - intShift, specUse1);
       int indMatchScore1 = tab1.getMatchScore(pepUse1);
       int indMatchScore2 = tab2.getMatchScore(pepUse2);
       tab1.addPRMScores(intShift, specUse2);

       matchScore1 = tab1.getMatchScore(pepUse1);
       //GFTable::DEBUG = true;
       tab1.computeGeneratingFunction(matchScore1);
       //GFTable::DEBUG = false;

       DEBUG_VAR(indMatchScore1);
       DEBUG_VAR(indMatchScore2);
       DEBUG_VAR(matchScore1);

       DEBUG_MSG("Pair (" << idx1 << " , " << idx2 << "), scans = ( " << (*m_inputSpectra)[idx1].scan << " , " << (*m_inputSpectra)[idx2].scan << " )");
       DEBUG_VAR(pepUse1);
       DEBUG_VAR(pepUse2);
       DEBUG_VAR(intShift);
       DEBUG_VAR(tab1.getPValue(matchScore1));

       tab1.initialize(specUse1, *m_jumps, specUse1);
       DEBUG_TRACE;
       tab1.computeGeneratingFunction(indMatchScore1);
       DEBUG_TRACE;
       tab2.computeGeneratingFunction(indMatchScore2);
       DEBUG_TRACE;

       intersection.initialize(tab1, tab2, intShift);
       DEBUG_TRACE;
       intersection.computeGeneratingFunction(false);
       DEBUG_TRACE;
       pair<double, double> pvals = intersection.getIntersectingPValue(pepUse1,
       pepUse2);
       DEBUG_MSG("intersecting pval1 = " << pvals.first);
       DEBUG_MSG("intersecting pval2 = " << pvals.second);

       GFTable mergedTable;
       mergedTable.initialize(specUse1, *m_jumps, specUse1);
       mergedTable.addPRMScores(intShift, specUse2);
       int mergedScore1 = mergedTable.getMatchScore(pepUse1);
       mergedTable.computeGeneratingFunction(mergedScore1);
       double mergedPval1 = mergedTable.getPValue(mergedScore1);

       mergedTable.initialize(specUse2, *m_jumps, specUse2);
       mergedTable.addPRMScores(0 - intShift, specUse1);
       int mergedScore2 = mergedTable.getMatchScore(pepUse2);
       mergedTable.computeGeneratingFunction(mergedScore2);
       double mergedPval2 = mergedTable.getPValue(mergedScore2);
       DEBUG_VAR(mergedScore1);
       DEBUG_VAR(mergedScore2);

       DEBUG_MSG("merged pval1 = " << mergedPval1);
       DEBUG_MSG("merged pval2 = " << mergedPval2);
       abort();
       */

      if (debug)
      {
        DEBUG_VAR(msScanToIdx[(*m_inputClusters)[specToClust[idx1]].m_scan]);
        DEBUG_VAR(msScanToIdx[(*m_inputClusters)[specToClust[idx2]].m_scan]);
      }

      const Spectrum &msSpectrum1 =
          (ignoreMS2ErrorModel) ? emptySpectrum :
              (*m_inputMsSpectra)[msScanToIdx[(*m_inputClusters)[specToClust[idx1]].m_scan]];

      const Spectrum &msSpectrum2 =
          (ignoreMS2ErrorModel) ? emptySpectrum :
              (*m_inputMsSpectra)[msScanToIdx[(*m_inputClusters)[specToClust[idx2]].m_scan]];

      if (debug)
      {
        DEBUG_VAR((*m_inputClusters)[specToClust[idx1]].m_scan);
        DEBUG_VAR(msSpectrum1.scan);
        DEBUG_VAR((*m_inputClusters)[specToClust[idx2]].m_scan);
        DEBUG_VAR(msSpectrum2.scan);
      }

      if (mergePrms)
      {
        tab1.initialize(specUse1,
                        *m_jumps,
                        msSpectrum1,
                        *m_inputErrorModelTarget,
                        *m_inputErrorModelDecoy,
                        ppmEdgeError);

        tab1.addPRMScores(intShift, specUse2);

        matchScore1 = tab1.getMatchScore(pepUse1);
        tab1.computeGeneratingFunction(matchScore1);

        vector<double> probs(2, 0);
        probs[0] = tab1.getPValue(matchScore1);

        if (intShift == 0 && pepUse1 == pepUse2)
        {
          probs[1] = probs[0];
        }
        else
        {
          tab2.initialize(specUse2,
                          *m_jumps,
                          msSpectrum2,
                          *m_inputErrorModelTarget,
                          *m_inputErrorModelDecoy,
                          ppmEdgeError);

          tab2.addPRMScores(0 - intShift, specUse1);

          matchScore2 = tab2.getMatchScore(pepUse2);
          tab2.computeGeneratingFunction(matchScore2);

          probs[1] = tab2.getPValue(matchScore2);
        }

        (*m_outputPairs).push_back((*m_inputPairs)[i]);
        (*m_outputProbs).push_back(probs);

        starProbs[idx1] = min(starProbs[idx1], probs[0]);
        starProbs[idx2] = min(starProbs[idx2], probs[1]);

      }
      else
      {

        if (idx1 != lastIdx1 || mergeAlignedPrms)
        {
          if (debug)
            DEBUG_MSG("Initializing table for idx " << idx1);

          tab1.initialize(specUse1,
                          *m_jumps,
                          msSpectrum1,
                          *m_inputErrorModelTarget,
                          *m_inputErrorModelDecoy,
                          ppmEdgeError);

          if (mergeAlignedPrms)
          {
            tab1.addMatchingPRMScores(intShift, specUse2);
          }

          matchScore1 = tab1.getMatchScore(pepUse1);

          if (debug)
          {
            tab1.computeGeneratingFunction(matchScore1);
          }
          else
          {
            tab1.optimizeDimensions(matchScore1);
          }

          lastIdx1 = idx1;
        }

        if (intShift == 0)
        {
          matchScoreOvlp = matchScore1;
        }
        else if (intShift > 0)
        {
          string pepOvlp = m_jumps->getPeptideSuffix(pepUse1, intShift);
          matchScoreOvlp = tab1.lookupPeptide(intShift, 0, pepOvlp).second;
        }

        tab2.initialize(specUse2,
                        *m_jumps,
                        msSpectrum2,
                        *m_inputErrorModelTarget,
                        *m_inputErrorModelDecoy,
                        ppmEdgeError);

        if (mergeAlignedPrms)
        {
          tab2.addMatchingPRMScores(0 - intShift, specUse1);
        }

        if (debug)
          DEBUG_MSG("Initializing table for idx " << idx2);

        matchScore2 = tab2.getMatchScore(pepUse2);

        if (intShift < 0)
        {
          string pepOvlp = m_jumps->getPeptideSuffix(pepUse2, 0 - intShift);
          matchScoreOvlp = tab2.lookupPeptide(0 - intShift, 0, pepOvlp).second;
        }

        if (debug)
        {
          tab2.computeGeneratingFunction(matchScore2);
        }
        else
        {
          tab2.optimizeDimensions(matchScore2);
        }

        if (debug)
        {
          DEBUG_VAR(matchScore1);
          DEBUG_VAR(matchScore2);
          DEBUG_VAR(tab1.getPValue(matchScore1));
          DEBUG_VAR(tab2.getPValue(matchScore2));
          DEBUG_VAR(tab1.getPValue(spec1.psmList.front()->m_score));
          DEBUG_VAR(tab2.getPValue(spec2.psmList.front()->m_score));
        }

        vector<double> probs(2, 0);

        if (intShift < 0)
        {
          if (debug)
          {
            DEBUG_TRACE;
          }
          intersection.initialize(tab2, tab1, 0 - intShift, matchScoreOvlp);
          if (debug)
          {
            DEBUG_TRACE;
          }
          intersection.computeGeneratingFunction(false);
          if (debug)
          {
            DEBUG_TRACE;
          }
          pair<double, double> pvals =
          intersection.getIntersectingPValue(pepUse2, pepUse1);
          probs[0] = pvals.second;
          probs[1] = pvals.first;
        }
        else
        {
          if (debug)
          {
            DEBUG_TRACE;
          }
          intersection.initialize(tab1, tab2, intShift, matchScoreOvlp);
          if (debug)
          {
            DEBUG_TRACE;
            //GFTableIntersection::DEBUG = true;
          }
          intersection.computeGeneratingFunction(false);
          if (debug)
          {
            DEBUG_TRACE;
          }
          pair<double, double> pvals =
          intersection.getIntersectingPValue(pepUse1, pepUse2);

          probs[0] = pvals.first;
          probs[1] = pvals.second;
        }

        SpectrumPair newPair;
        newPair = (*m_inputPairs)[i];

        if (debug)
        {
          DEBUG_VAR(probs[0]);
          DEBUG_VAR(probs[1]);
        }

        if ((probs[0] / 10.0) > spec1.psmList.front()->m_pValue)
        {
          ERROR_MSG("Higher pvalues for scans (" << msSpectrum1.scan << "," << msSpectrum2.scan << ") = " << probs[0] << " > " << spec1.psmList.front()->m_pValue);
          abort();
        }

        if ((probs[1] / 10.0) > spec2.psmList.front()->m_pValue)
        {
          ERROR_MSG("Higher pvalues for scans (" << msSpectrum2.scan << "," << msSpectrum1.scan << ") = " << probs[1] << " > " << spec2.psmList.front()->m_pValue);
          abort();
        }

        if (probs[0] == 0 || probs[1] == 0)
        {
#pragma omp critical
          {
            WARN_MSG("Failed to compute intersection for (" << idx1 << "," << idx2 << "), shift = " << shiftUse << ", peptides = ( " << pepUse1 << " , " << pepUse2 << " )");
          }
          //abort();
        }
        /*else if (probs[0] > (*m_inputSpectra)[idx1].psmList.front()->m_pValue * 1.001
         || probs[1] > (*m_inputSpectra)[idx2].psmList.front()->m_pValue * 1.001)
         {
         #pragma omp critical
         {
         WARN_MSG("Invalid intersecting probabilities for (" << idx1 << "," << idx2 << "), shift = " << shiftUse << ", peptides = ( " << pepUse1 << " , " << pepUse2 << " )");
         WARN_MSG((*m_inputSpectra)[idx1].psmList.front()->m_pValue << " -> " << probs[0]);
         WARN_MSG((*m_inputSpectra)[idx2].psmList.front()->m_pValue << " -> " << probs[1]);
         }
         }
         */else
        {
          (*m_outputPairs).push_back((*m_inputPairs)[i]);
          (*m_outputProbs).push_back(probs);

          starProbs[idx1] = min(starProbs[idx1], probs[0]);
          starProbs[idx2] = min(starProbs[idx2], probs[1]);

          if (verbose)
          {
#pragma omp critical
            {
              DEBUG_MSG("Intersection (" << idx1 << "," << idx2 << ") has probabilities ( " << (*m_inputSpectra)[idx1].psmList.front()->m_pValue << " -> " << probs[0] << " , " << (*m_inputSpectra)[idx2].psmList.front()->m_pValue << " -> " << probs[1] << " )");
            }
          }
        }
      }

      numComputed++;

      if (numComputed % numComputationsPerCompress == 0)
      {
        intersection.compress();
      }

      double curProg = ((((double)pairIdxEnd) - ((double)i)) * 100.0)
          / (((double)pairIdxEnd) - ((double)pairIdxStart));

      if ((int)curProg > lastProg)
      {
        lastProg = (int)curProg;

        final = clock();
        double seconds = (((double)final) - ((double)init))
            / ((double)CLOCKS_PER_SEC);
        double secondsPerIntersection = seconds / ((double)numComputed);

#pragma omp critical
        {
          DEBUG_MSG(m_name << " ... " << lastProg << "%, sec/pair = " << secondsPerIntersection);
        }
      }

    }

    DEBUG_VAR(m_outputPairs->size());
    DEBUG_VAR(m_outputProbs->size());

    return true;
  }

  bool ExecIntersectingSpecProb::loadInputData(void)
  {

    int debugIdx1 = m_params.getValueInt("DEBUG_IDX1", -1);
    int debugIdx2 = m_params.getValueInt("DEBUG_IDX2", -1);

    if (m_params.exists("INPUT_SPECTRA"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA from \'" << m_params.getValue("INPUT_SPECTRA") << "\' ...");
      if (!m_inputSpectra->Load(m_params.getValue("INPUT_SPECTRA").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA") << "\'");
        return false;
      }
      for (int i = 0; i < m_inputSpectra->size(); i++)
      {
        (*m_inputSpectra)[i].psmList.clear();
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!!!");
      return false;
    }

    if (m_params.exists("INPUT_PAIRS"))
    {
      DEBUG_MSG("Loading INPUT_PAIRS from \'" << m_params.getValue("INPUT_PAIRS") << "\' ...");
      if (!m_inputPairs->loadFromBinaryFile(m_params.getValue("INPUT_PAIRS")))
      {
        ERROR_MSG("Failed to load spectrum pairs from \'" << m_params.getValue("INPUT_PAIRS") << "\'");
        return false;
      }
    }

    if (m_params.exists("INPUT_CLUSTERS"))
    {
      DEBUG_MSG("Loading INPUT_CLUSTERS from \'" << m_params.getValue("INPUT_CLUSTERS") << "\' ...");
      if (!m_inputClusters->loadBinaryFile(m_params.getValue("INPUT_CLUSTERS")))
      {
        ERROR_MSG("Failed to load clusters from \'" << m_params.getValue("INPUT_CLUSTERS") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputClusters->size());
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("No input spectra!!!");
      return false;
    }

    if (m_params.exists("INPUT_SPECTRA_MS"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA_MS from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\' ...");
      if (!m_inputMsSpectra->Load(m_params.getValue("INPUT_SPECTRA_MS").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputMsSpectra->size());

      DEBUG_MSG("Pre-processing MS/MS spectra ...");
      for (int i = 0; i < m_inputClusters->size(); i++)
      {
        int clustIdx = (*m_inputClusters)[i].m_index;
        if ((*m_inputClusters)[i].size() == 0
            || (debugIdx1 >= 0 && debugIdx2 >= 0 && clustIdx != debugIdx1
                && clustIdx != debugIdx2))
        {
          continue;
        }
        (*m_inputMsSpectra)[clustIdx].rankFilterPeaks(10);
        (*m_inputMsSpectra)[clustIdx].setPeakTolerance(0);
        (*m_inputMsSpectra)[clustIdx].roundPeaks(AA_ROUNDING, true, false);

      }
      DEBUG_MSG("finished");
    }

    if (m_params.exists("INPUT_ERROR_MODEL_TARGET"))
    {
      string param = m_params.getValue("INPUT_ERROR_MODEL_TARGET");
      DEBUG_MSG("Loading INPUT_ERROR_MODEL_TARGET from \'" << param << "\' ...");
      if (!m_inputErrorModelTarget->loadBinaryFile(param))
      {
        ERROR_MSG("Failed to load target error model from \'" << param << "\'");
        return false;
      }
      DEBUG_VAR(m_inputErrorModelTarget->size());
    }

    if (m_params.exists("INPUT_ERROR_MODEL_DECOY"))
    {
      string param = m_params.getValue("INPUT_ERROR_MODEL_DECOY");
      DEBUG_MSG("Loading INPUT_ERROR_MODEL_DECOY from \'" << param << "\' ...");
      if (!m_inputErrorModelDecoy->loadBinaryFile(param))
      {
        ERROR_MSG("Failed to load target error model from \'" << param << "\'");
        return false;
      }
      DEBUG_VAR(m_inputErrorModelDecoy->size());
    }

    if (m_params.exists("INPUT_PSMS"))
    {
      const string psmType = m_params.getValue("INPUT_PSMS_TYPE", "");
      const string psmFile = m_params.getValue("INPUT_PSMS");

      DEBUG_MSG("Loading INPUT_PSMS from \'" << m_params.getValue("INPUT_PSMS") << "\' ...");
      if (!m_inputPSMs->Load(psmFile, psmType))
      {
        ERROR_MSG("Failed to load PSMs from \'" << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }

      int numAdded = m_inputPSMs->addMostSpectra(m_inputSpectra);

      if (numAdded == 0)
      {
        ERROR_MSG("No PSMs could be matched to input spectra!!!");
        return false;
      }
    }

    if (m_params.exists("INPUT_MSGFDB_PSM_PAIRS"))
    {
      const string psmFile = m_params.getValue("INPUT_MSGFDB_PSM_PAIRS");

      OutputTable tempPsms;

      DEBUG_MSG("Loading INPUT_MSGFDB_PSM_PAIRS from \'" << m_params.getValue("INPUT_MSGFDB_PSM_PAIRS") << "\' ...");
      if (!tempPsms.loadFromCSV(psmFile.c_str(), "\t"))
      {
        ERROR_MSG("Failed to load PSMs from \'" << psmFile << "\'");
        return false;
      }

      if (tempPsms.values.size() <= 1)
      {
        ERROR_MSG("No rows found in PSM pair table " << psmFile << "!");
        return false;
      }

      DEBUG_MSG("Locating pairs ...");

      int specFileIdx = -1, specIdxIdx = -1, specScanIdx = -1;
      for (int i = 0; i < tempPsms.values[0].size(); i++)
      {
        if (tempPsms.values[0][i].first == "#SpecFile")
        {
          specFileIdx = i;
        }
        else if (tempPsms.values[0][i].first == "SpecIndex")
        {
          specIdxIdx = i;
        }
        else if (tempPsms.values[0][i].first == "Scan#")
        {
          specScanIdx = i;
        }
        else if (specFileIdx >= 0 && specIdxIdx >= 0 && specScanIdx >= 0)
        {
          break;
        }
        else
        {
          continue;
        }
      }

      const char* clustSep = "/";
      string scanStr("");
      vector<string> clusteredScans(0);
      string specIdxStr("");
      string specFilename("");
      vector<string> clusteredIdxs(0);

      m_inputSpectra->index();

      vector<set<int> > seenPairs(m_inputSpectra->size());

      vector<vector<TwoValues<int> > > sortedPairs(m_inputSpectra->size());
      SpectrumPairSet tempPairs;

      for (int i = 1; i < tempPsms.values.size(); i++)
      {

        scanStr = tempPsms.values[i][specScanIdx].first;
        specIdxStr = tempPsms.values[i][specIdxIdx].first;
        specFilename = tempPsms.values[i][specFileIdx].first;

        FilenameManager fm(specFilename);

        splitText(scanStr.c_str(), clusteredScans, clustSep);

        for (int j1 = 0; j1 < clusteredScans.size(); j1++)
        {
          string id1 = clusteredScans[j1];
          id1 += "$";
          id1 += fm.filename;
          int idx1 = m_inputSpectra->getIndexFromID(id1);

          if (idx1 < 0)
          {
            //WARN_MSG("Failed to locate spectrum ID " << id1);
            continue;
          }

          for (int j2 = j1 + 1; j2 < clusteredScans.size(); j2++)
          {
            string id2 = clusteredScans[j2];
            id2 += "$";
            id2 += fm.filename;
            int idx2 = m_inputSpectra->getIndexFromID(id2);

            if (idx2 < 0)
            {
              //WARN_MSG("Failed to locate spectrum ID " << id2);
              continue;
            }

            if (seenPairs[idx1].count(idx2) > 0)
            {
              continue;
            }
            seenPairs[idx1].insert(idx2);
            seenPairs[idx2].insert(idx1);

            SpectrumPair newPair;
            newPair.spec1 = idx2;
            newPair.spec2 = idx1;
            newPair.shift1 = 0;
            newPair.shift2 = 0;

            int pairIdx = tempPairs.size();
            tempPairs.push_back(newPair);

            sortedPairs[newPair.spec1].push_back(TwoValues<int>(newPair.spec2,
                                                                pairIdx));
          }
        }
      }

      DEBUG_MSG("Sorting pairs ...");

      for (int i = 0; i < sortedPairs.size(); i++)
      {
        if (sortedPairs[i].size() == 0)
        {
          continue;
        }

        sort(sortedPairs[i].begin(), sortedPairs[i].end());

        for (int j = 0; j < sortedPairs[i].size(); j++)
        {
          m_inputPairs->push_back(tempPairs[sortedPairs[i][j][1]]);
        }
      }

      DEBUG_MSG("Added " << tempPairs.size() << " pairs from MSGFDB PSMs");
    }

    if (m_inputPSMs->size() == 0)
    {
      ERROR_MSG("Input PSMs size is 0!!!");
      return false;
    }

    return true;
  }

  bool ExecIntersectingSpecProb::saveOutputData(void)
  {
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_PAIRS"))
    {
      if (!m_outputPairs->saveToBinaryFile(m_params.getValue("OUTPUT_PAIRS")))
      {
        ERROR_MSG("Failed to save spectrum pairs to \'" << m_params.getValue("OUTPUT_PAIRS") << "\'");
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PROBS"))
    {
      string outProbsFile = m_params.getValue("OUTPUT_PROBS");
      if (Save_binArray<double>(outProbsFile.c_str(), *m_outputProbs) != 1)
      {
        ERROR_MSG("Failed to save intersecting probabilities to \'" << outProbsFile << "\'");
        return false;
      }
    }
    return true;
  }

  bool ExecIntersectingSpecProb::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();

    string spectraFilename = baseDirectory + "input_spectra_intersect.pklbin";
    if (!fileExists(spectraFilename))
    {
      if (m_inputSpectra->savePklBin(spectraFilename.c_str()) != 1)
      {
        ERROR_MSG("Failed to save spectra to \'" << spectraFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << spectraFilename << " (already exists)");
    }
    m_params.setValue("INPUT_SPECTRA", spectraFilename);

    string clustersFilename = baseDirectory + "input_spectra_intersect.clust";
    if (!fileExists(clustersFilename))
    {
      if (!m_inputClusters->saveBinaryFile(clustersFilename.c_str()))
      {
        ERROR_MSG("Failed to save clusters to \'" << clustersFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << clustersFilename << " (already exists)");
    }
    m_params.setValue("INPUT_CLUSTERS", clustersFilename);

    string psmsFilename = baseDirectory + "input_psms_intersect.txt";
    if (!fileExists(psmsFilename))
    {
      if (!m_inputPSMs->saveToFile(psmsFilename.c_str(), true))
      {
        ERROR_MSG("Failed to save psms to \'" << psmsFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << psmsFilename << " (already exists)");
    }
    m_params.setValue("INPUT_PSMS", psmsFilename);

    string pairsFilename = baseDirectory + "input_pairs_intersect.bin";
    if (!fileExists(pairsFilename))
    {
      if (!m_inputPairs->saveToBinaryFile(pairsFilename))
      {
        ERROR_MSG("Failed to save pairs to \'" << pairsFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << pairsFilename << " (already exists)");
    }
    m_params.setValue("INPUT_PAIRS", pairsFilename);

    string massErrorTgtFilename = baseDirectory + "error_model_target.merr";
    if (!fileExists(massErrorTgtFilename))
    {
      if (!m_inputErrorModelTarget->saveBinaryFile(massErrorTgtFilename))
      {
        ERROR_MSG("Failed to save target error model to \'" << massErrorTgtFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << massErrorTgtFilename << " (already exists)");
    }
    m_params.setValue("INPUT_ERROR_MODEL_TARGET", massErrorTgtFilename);

    string massErrorDecFilename = baseDirectory + "error_model_decoy.merr";
    if (!fileExists(massErrorDecFilename))
    {
      if (!m_inputErrorModelDecoy->saveBinaryFile(massErrorDecFilename))
      {
        ERROR_MSG("Failed to save decoy error model to \'" << massErrorDecFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << massErrorDecFilename << " (already exists)");
    }
    m_params.setValue("INPUT_ERROR_MODEL_DECOY", massErrorDecFilename);

    string spectraMsFilename = baseDirectory + "input_spectra_ms.pklbin";
    if (!fileExists(spectraMsFilename))
    {
      if (m_inputMsSpectra->savePklBin(spectraMsFilename.c_str()) != 1)
      {
        ERROR_MSG("Failed to save MS/MS spectra to \'" << spectraMsFilename << "\'");
        return false;
      }
    }
    else
    {
      DEBUG_MSG("Not Saving " << spectraMsFilename << " (already exists)");
    }
    m_params.setValue("INPUT_SPECTRA_MS", spectraMsFilename);

    //SpecSet m_inputSpectra; // the input spectra
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(spectraFilename);
    filenames.push_back(clustersFilename);
    filenames.push_back(psmsFilename);
    filenames.push_back(pairsFilename);
    filenames.push_back(massErrorTgtFilename);
    filenames.push_back(massErrorDecFilename);
    filenames.push_back(spectraMsFilename);

    return true;
  }

  bool ExecIntersectingSpecProb::loadOutputData(void)
  {
    if (m_params.exists("OUTPUT_PAIRS"))
    {
      const string filename = m_params.getValue("OUTPUT_PAIRS");
      DEBUG_MSG("Loading output pairs from \'" << filename << "\'");
      if (!m_outputPairs->loadFromBinaryFile(filename))
      {
        ERROR_MSG("Failed to load spectrum pairs from \'" << filename << "\'");
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PROBS"))
    {
      string filename = m_params.getValue("OUTPUT_PROBS");
      DEBUG_MSG("Loading output intersecting probabilities from \'" << filename << "\'");
      if (Load_binArray<double>(filename.c_str(), *m_outputProbs) != 1)
      {
        ERROR_MSG("Failed to load intersecting probabilities from \'" << filename << "\'");
        return false;
      }
    }
    return true;
  }

// -------------------------------------------------------------------------
  vector<ExecBase*>
  const & ExecIntersectingSpecProb::split(int numSplit)
  {

    if (m_inputSpectra->size() == 0 || m_inputPSMs->size() == 0)
    {
      if (!loadInputData())
      {
        return m_subModules;
      }
    }

    int pairIdxEnd = m_params.getValueInt("PAIR_IDX_END",
                                          ((int)m_inputPairs->size()) - 1);
    int pairIdxStart = m_params.getValueInt("PAIR_IDX_START", 0);

    DEBUG_VAR(m_inputPairs->size());
    if (m_inputPairs->size() == 0)
    {
      generatePSMAligns(m_params,
                        *m_jumps,
                        *m_inputSpectra,
                        *m_inputClusters,
                        *m_inputPairs);
      pairIdxEnd = ((int)m_inputPairs->size()) - 1;
      pairIdxStart = 0;
    }
    DEBUG_VAR(m_inputPairs->size());

    DEBUG_VAR(pairIdxStart);
    DEBUG_VAR(pairIdxEnd);

    m_inputPairs->sort_pairs_by_index();

    DEBUG_VAR(numSplit);

    double totalComputations = 0;
    vector<double> computationsPerPair(m_inputPairs->size(), 0);
    for (int i = pairIdxStart; i <= pairIdxEnd; i++)
    {
      int idx1 = (*m_inputPairs)[i].spec1;
      int idx2 = (*m_inputPairs)[i].spec2;

      const Spectrum &spec1 = (*m_inputSpectra)[idx1];
      const Spectrum &spec2 = (*m_inputSpectra)[idx2];

      float shift = (*m_inputPairs)[i].shift1;

      double numMassBins = (double)max(spec1.parentMass - shift,
                                       spec2.parentMass + shift);
      double score1 = spec1.getTotalIonCurrent();
      double score2 = spec2.getTotalIonCurrent();

      computationsPerPair[i] = numMassBins * score1 * score2;
      totalComputations += computationsPerPair[i];
    }
    double computationsPerModule = totalComputations / ((double)numSplit);

    DEBUG_VAR(totalComputations);
    DEBUG_VAR(computationsPerModule);

    int lastIdx1 = -1;
    int curPairIdxEnd = pairIdxEnd;
    int curSplit = 0;
    m_subModules.resize(numSplit, (ExecBase*)0);
    double runningComputations = 0;
    double runningComputationsPerModule = 0;

    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }

    for (int i = pairIdxEnd; i >= pairIdxStart; i--)
    {
      int idx1 = (*m_inputPairs)[i].spec1;
      int idx2 = (*m_inputPairs)[i].spec2;

      if (lastIdx1 < 0)
      {
        lastIdx1 = idx1;
      }

      runningComputations += computationsPerPair[i];
      runningComputationsPerModule += computationsPerPair[i];

      if ((runningComputationsPerModule >= computationsPerModule
          && curSplit < numSplit + 1) || i == pairIdxStart)
      {

        ParameterList childParams(m_params);
        childParams.setValue("PAIR_IDX_END", parseInt(curPairIdxEnd));
        childParams.setValue("PAIR_IDX_START", parseInt(i));

        DEBUG_MSG("Job " << curSplit << " has range [" << i << "," << curPairIdxEnd << "]");

        ExecBase * theClone =
            new ExecIntersectingSpecProb(childParams,
                                         m_inputSpectra,
                                         m_inputPSMs,
                                         m_inputPairs,
                                         m_inputClusters,
                                         m_inputErrorModelTarget,
                                         m_inputErrorModelDecoy,
                                         m_inputMsSpectra,
                                         m_jumps);

        theClone->setName(makeName(m_name, curSplit));

        string baseName = dataDir + "/" + theClone->getName();
        theClone->m_params.setValue("OUTPUT_PAIRS",
                                    baseName + "_output_aligns.bin");
        theClone->m_params.setValue("OUTPUT_PROBS",
                                    baseName + "_output_probs.bin");

        m_subModules[curSplit++] = theClone;
        curPairIdxEnd = i - 1;

        runningComputationsPerModule = 0;

        if (curSplit < numSplit)
        {
          computationsPerModule = (totalComputations - runningComputations)
              / (((double)numSplit) - ((double)curSplit));
        }
      }
    }

    DEBUG_VAR(curSplit);

    return m_subModules;
  }

// -------------------------------------------------------------------------
  bool ExecIntersectingSpecProb::merge(void)
  {
    m_outputPairs->resize(0);
    m_outputProbs->resize(0);

    for (int nSplit = 0; nSplit < m_subModules.size(); nSplit++)
    {
      if (m_subModules[nSplit] == 0)
      {
        ERROR_MSG("ExecBase module at index " << nSplit << " is zero");
        return false;
      }
      ExecIntersectingSpecProb* module =
          dynamic_cast<ExecIntersectingSpecProb*>(m_subModules[nSplit]);

      if (module == 0)
      {
        ERROR_MSG("Failed to cast module \'" << m_subModules[nSplit]->getName() << "\' to ExecIntersectingSpecProb");
        return false;
      }
      for (int i = 0; i < module->m_outputPairs->size(); i++)
      {
        m_outputPairs->push_back((*module->m_outputPairs)[i]);
        m_outputProbs->push_back((*module->m_outputProbs)[i]);
      }
    }
    return true;
  }

  bool ExecIntersectingSpecProb::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
    VALIDATE_PARAM_EXIST("INPUT_PSMS");

    m_isValid = true;
    return true;
  }
}
