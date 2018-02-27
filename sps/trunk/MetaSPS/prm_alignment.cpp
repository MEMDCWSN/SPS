/*
 * prm_alignment.cpp
 *
 *  Created on: Feb 18, 2011
 *      Author: aguthals
 */

#include "prm_alignment.h"

namespace specnets
{
  /**
   * Computes the highest-scoring shift of spectrum2 in relation to spectrum1
   *   and reversed(spectrum2) in relation to spectrum1. This alignment is
   *   stored in the SpectrumPair class variables
   * @param minMatchedPeaks minimum number of matching peaks allowed
   * @param minScore minimum alignment score allowed
   * @param score_type see method gettShiftScore
   * @return true if alignment was found, false if none met criteria
   */
  bool PRMAlignment::computeShiftNoModFR(int minMatchedPeaks,
                                         float minScore,
                                         short score_type)
  {
    SpectrumPair alignF;
    SpectrumPair alignR;

    reverseSpec2(false);
    bool resF = getShiftNoMod(minMatchedPeaks, minScore, score_type, alignF);
    float scoreF = min(alignF.score1, alignF.score2);

    reverseSpec2(true);
    bool resR = getShiftNoMod(minMatchedPeaks, minScore, score_type, alignR);
    float scoreR = min(alignR.score1, alignR.score2);

    if ((resF && !resR) || (resF && resR && scoreF >= scoreR))
    {
      score1 = alignF.score1;
      score2 = alignF.score2;
      shift1 = alignF.shift1;
      shift2 = alignF.shift2;
      spec2rev = false;
      return true;
    }
    else if ((!resF && resR) || (resF && resR && scoreF < scoreR))
    {
      score1 = alignR.score1;
      score2 = alignR.score2;
      shift1 = alignR.shift1;
      shift2 = alignR.shift2;
      spec2rev = true;
      return true;
    }
    else
    {
      return false;
    }
  }

  /**
   * Computes the highest-scoring shift of spectrum2 in relation to spectrum1
   *   and stores the alignment result in "alignment"
   * @param minMatchedPeaks minimum number of matching peaks allowed
   * @param minScore minimum alignment score allowed
   * @param score_type see method computShiftScore
   * @return true if alignment was found, false if none met criteria
   */
  bool PRMAlignment::getShiftNoMod(int minMatchedPeaks,
                                   float minScore,
                                   short score_type,
                                   SpectrumPair& alignment)
  {
    computeShifts(minMatchedPeaks);

    float max_score = -2;
    float cur_score;

    pair<int, pair<float, float> > shift_scores;
    pair<int, pair<float, float> > best_shift_scores;
    MZRange bestShift;
    MZRange curShift;

    short mergeType = (score_type == 0) ? 1 : 0;

    for (set<MZRange>::iterator shiftIt = bestShifts->begin(); shiftIt
        != bestShifts->end(); shiftIt++)
    {

      curShift.MergeMZRanges(&(*shift_mp)[*shiftIt].first, mergeType);

      shift_scores = getShiftScore(curShift.getMass(),
                                   curShift.getTolerance(),
                                   score_type);

      cur_score = min(shift_scores.second.first, shift_scores.second.second);

      if (cur_score > max_score)
      {
        max_score = cur_score;
        best_shift_scores = shift_scores;
        bestShift = curShift;
      }
    }

    if (max_score < minScore)
    {
      return false;
    }

    alignment.shift1 = bestShift.getMass();
    alignment.score1 = best_shift_scores.second.first;
    alignment.score2 = best_shift_scores.second.second;

    return true;

  }

  /**
   * Computes the alignment scores of spectrum2 shifted wrt spectrum1 and returns number of matching peaks
   * @param shift
   * @param shift_tol Da tolerance of shift
   * @param alignedScores optional output data structure containing:
   *   first[0] = matched score of spectrum1
   *   first[1] = overlapping score of spectrum1
   *   first[2] = total score of spectrum1
   *   second[0] = matched score of spectrum2
   *   second[1] = overlapping score of spectrum2
   *   second[2] = total score of spectrum2
   * @param score_type which type of scoring metric to use
   *   0 - overlapScore
   *   1 - totalScore
   * @return first = matched peaks, second.first = score of spectrum1, second.second = score of spectrum2
   */
  pair<int, pair<float, float> > PRMAlignment::getShiftScore(float shift,
                                                             float shift_tol,
                                                             short score_type,
                                                             pair<
                                                                 vector<float> ,
                                                                 vector<float> >* alignedScores) const
  {

    if (spectrum1->size() == 0 || spectrum2->size() == 0)
    {
      return pair<int, pair<float, float> > (0, pair<float, float> (0, 0));
    }

    int numMP = 0;
    float matched1 = 0, matched2 = 0;
    float overlap1 = 0, overlap2 = 0;
    float total1 = 0, total2 = 0;

    MZRange shift_range(shift, 0, shift_tol);

    MZRange endPt1(spectrum1->parentMass - AAJumps::massHion,
                   0,
                   spectrum1->getTolerance(0));
    MZRange endPt2(spectrum2->parentMass - AAJumps::massHion,
                   0,
                   spectrum2->getTolerance(0));

    endPt2 = shift_range + endPt2;

    MZRange range1;
    /*
     cout << "spec1:\n";
     spectrum1->output(cout);

     cout << "\nspec2:\n";
     spectrum2->output(cout);
     */
    int j = 0;
    bool finished2 = false;
    MZRange range2((*spectrum2)[j][0] + shift_range.getMass(),
                   (*spectrum2)[j][1],
                   spectrum2->getTolerance(j));

    while (range2 < 0.0 && !finished2)
    {
      total2 += range2.getIntensity();

      j++;
      if (j >= szSpec2)
      {
        finished2 = true;
        break;
      }

      range2.set((*spectrum2)[j][0] + shift_range.getMass(),
                 (*spectrum2)[j][1],
                 spectrum2->getTolerance(j));
    }

    for (int i = 0; i < szSpec1; i++)
    {
      range1.set((*spectrum1)[i][0],
                 (*spectrum1)[i][1],
                 spectrum1->getTolerance(i));
      total1 += range1.getIntensity();

      if (range1 >= shift_range && range1 <= endPt2)
      {
        overlap1 += range1.getIntensity();
      }

      while (range2 < range1 && range2 <= endPt1 && !finished2)
      {
        total2 += range2.getIntensity();
        overlap2 += range2.getIntensity();

        j++;
        if (j >= szSpec2)
        {
          finished2 = true;
          break;
        }
        range2.set((*spectrum2)[j][0] + shift_range.getMass(),
                   (*spectrum2)[j][1],
                   spectrum2->getTolerance(j));
      }

      while (range2 == range1 && !finished2)
      {
        matched1 += range1.getIntensity();

        matched2 += range2.getIntensity();
        overlap2 += range2.getIntensity();
        total2 += range2.getIntensity();

        numMP++;

        j++;
        if (j >= szSpec2)
        {
          finished2 = true;
          break;
        }
        range2.set((*spectrum2)[j][0] + shift_range.getMass(),
                   (*spectrum2)[j][1],
                   spectrum2->getTolerance(j));
      }
    }

    while (!finished2)
    {
      total2 += range2.getIntensity();

      if (range2 >= 0.0 && range2 <= endPt1)
      {
        overlap2 += range2.getIntensity();
      }

      j++;
      if (j >= szSpec2)
      {
        finished2 = true;
        break;
      }
      range2.set((*spectrum2)[j][0] + shift_range.getMass(),
                 (*spectrum2)[j][1],
                 spectrum2->getTolerance(j));
    }

    float score1, score2;

    switch (score_type)
    {
    case 0:
      score1 = overlapScore(numMP, matched1, overlap1);
      score2 = overlapScore(numMP, matched2, overlap2);
      break;
    case 1:
      score1 = totalScore(matched1, total1);
      score2 = totalScore(matched2, total2);
      break;
    }

    if (alignedScores != 0)
    {
      alignedScores->first.resize(3);
      alignedScores->first[0] = matched1;
      alignedScores->first[1] = overlap1;
      alignedScores->first[2] = total1;

      alignedScores->second.resize(3);
      alignedScores->second[0] = matched2;
      alignedScores->second[1] = overlap2;
      alignedScores->second[2] = total2;
    }
    //cout << "shift = " << shift << ", shift_tol = " << shift_tol
    //    << ", score1 = " << score1 << ", score2 = " << score2 << "\n";
    //cout << "matched1/overlap1 = " <<  matched1/overlap1 << ", matched2/overlap2 = " <<  matched2/overlap2 << endl;

    //pair<float, float> scores(score1, score2);
    return pair<int, pair<float, float> > (numMP, pair<float, float> (score1,
                                                                      score2));
  }

  /**
   * Computes all possible shifts of spectrum2 in relation to spectrum1 and
   *   stores them in shift_mp and bestShifts
   * @param minMatchedPeaks only consider shifts with minimum number of
   *   matching peaks. Also, put any shifts with minMatchedPeaks in
   *   bestShifts
   * @return
   */
  void PRMAlignment::computeShifts(int minMatchedPeaks)
  {
    cleanup();

    list<pair<int, int> > mps;
    list<MZRange> shifts;
    pair<int, int> mp;
    pair<list<MZRange> , list<pair<int, int> > > shiftMP;

    pair<list<MZRange> , list<pair<int, int> > >* shiftMP_ptr;

    int maxJ = szSpec2 - minMatchedPeaks;
    int minJ = minMatchedPeaks - szSpec1 - 2;
    int jBoundMin, jBoundMax;
    int numMP;

    float tol1, tol2;

    MZRange new_shift;
    MZRange range1;
    MZRange range2;

    for (int i = 0; i < szSpec1; i++)
    {

      tol1 = spectrum1->getTolerance(i);

      range1.set((*spectrum1)[i][0], (*spectrum1)[i][1], tol1);

      maxJ++;
      minJ++;
      jBoundMin = max(0, minJ);
      jBoundMax = min(szSpec2, maxJ);

      for (int j = jBoundMin; j < jBoundMax; j++)
      {

        tol2 = spectrum2->getTolerance(j);

        range2.set((*spectrum2)[j][0], (*spectrum2)[j][1], tol2);

        new_shift = range1 - range2;

        mp.first = i;
        mp.second = j;

        if (shift_mp->count(new_shift) == 0)
        {
          mps.clear();
          mps.push_back(mp);
          shifts.clear();
          shifts.push_back(new_shift);
          shiftMP.first = shifts;
          shiftMP.second = mps;
          (*shift_mp)[new_shift] = shiftMP;
          numMP = 1;
        }
        else
        {
          shiftMP_ptr = &(*shift_mp)[new_shift];
          shiftMP_ptr->first.push_back(new_shift);
          shiftMP_ptr->second.push_back(mp);
          numMP = shiftMP_ptr->first.size();
        }

        if (numMP >= minMatchedPeaks && bestShifts->count(new_shift) == 0)
        {
          bestShifts->insert(new_shift);
        }
      }
    }
  }

  /**
   * Called internally to clear all shifts
   */
  void PRMAlignment::cleanup()
  {
    shift_mp->clear();
    bestShifts->clear();
  }

}

