// Header Includes
#include "Logger.h"
#include "SpectralPairs.h"

// External Module Includes
#include "alignment_scoring.h"
#include "label.h"
#include "msn.h"
#include "aminoacid.h"

// System Includes

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  //
  //  SplitSpectra - duplicates every i-th spectrum in specSet into
  ///    itself+reversed (positions 2*i/2*i+1 in specSetSplit)
  //
  // -------------------------------------------------------------------------
  void SplitSpectra(SpecSet &specSet, SpecSet &specSetSplit)
  {
    specSetSplit.resize(2 * specSet.size());
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specSetSplit[2 * i] = specSet[i];
      specSet[i].reverse(0, &specSetSplit[2 * i + 1]);
    }
  }

  // -------------------------------------------------------------------------
  void SplitPairs(SpecSet &specSet,
                  SpectrumPairSet &aligns,
                  float peakTol,
                  float pmTol,
                  int maxAAjump,
                  float maxModMass,
                  float penalty_sameVert,
                  float penalty_ptm,
                  vector<vector<TwoValues<int> > > &matches,
                  vector<bool> &specFlipped,
                  vector<float> &modPos,
                  unsigned int minMatchedPeaks,
                  unsigned int minEdgesToComponent,
                  bool forceSymmetry,
                  bool ignoreReversals,
                  vector<vector<float> > *alignStats,
                  vector<SpectrumPeakLabels> *labelsP)
  {

    vector<list<int> > alignsEntries(specSet.size()); // List of pairs in aligns that incide on the corresponding vertex
    vector<float> flipScores(specSet.size()); // Accumulated don't_flip/flip score balance per vertex
    SpecSet specSetRev;
    specSetRev.resize(specSet.size()); // Holds reversed versions of the spectra in the component
    vector<int> numPairsPerSpectrum(specSet.size()); // Number of pairs to from spectrum i to spectra included in the component
    list<int> toProcess; // Indices of spectra left to process
    TwoValues<vector<vector<TwoValues<int> > > > tmpMatches; // Temporary matches when a neighbor spectrum S1 in the i-th pair is matched directly/reversed (its pair S2 is assumed to already be in the component)
    vector<TwoValues<float> > tmpModPos(aligns.size()); // Temporary storage for the position of the modification
    vector<bool> matchComputed(aligns.size()); // Indicates whether the i-th entry in tmpAligns was already computed
    SpectrumPair curPair; // Used to access aligns for template class T = SpectrumPair or SpectrumPair
    bool isASP;
    unsigned int specsInComponent = 0; // Number of spectra already included in the component

    // Initializations
    matches.resize(aligns.size());
    modPos.resize(aligns.size());
    tmpMatches[0].resize(aligns.size());
    tmpMatches[1].resize(aligns.size());
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      if (specSet[aligns[i].spec1].size() == 0
          || specSet[aligns[i].spec1].size()
              != specSetRev[aligns[i].spec1].size())
      {
        specSet[aligns[i].spec1].reverse(0, &specSetRev[aligns[i].spec1]);
        toProcess.push_back(aligns[i].spec1);
      }
      if (specSet[aligns[i].spec2].size() == 0
          || specSet[aligns[i].spec2].size()
              != specSetRev[aligns[i].spec2].size())
      {
        specSet[aligns[i].spec2].reverse(0, &specSetRev[aligns[i].spec2]);
        toProcess.push_back(aligns[i].spec2);
      }
      matchComputed[i] = false;
      alignsEntries[aligns[i].spec1].push_back(i);
      alignsEntries[aligns[i].spec2].push_back(i);
      matches[i].resize(0);
    }

    list<int>::iterator pIter;
    for (pIter = toProcess.begin(); pIter != toProcess.end(); pIter++)
    {
      flipScores[*pIter] = 0;
    }

    // Use edge scores as percentages of spectrum score to select the edge/vertex to initialize the direction selection
    vector<float> specScores(specSet.size());
    float bestScore = 0, curScore;
    int bestScoreIdx = -1, specToAdd = -1; // specToAdd is index of the new spectrum being added to the component
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specScores[i] = 0;
      numPairsPerSpectrum[i] = 0;
      for (unsigned int j = 0; j < specSet[i].size(); j++)
      {
        specScores[i] += specSet[i][j][1];
      }
    }
    // Select starting spectrum by participation in highest-scoring alignment
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      curScore = min(aligns[i].score1 / specScores[aligns[i].spec1],
                     aligns[i].score2 / specScores[aligns[i].spec2]);
      numPairsPerSpectrum[aligns[i].spec1]++;
      numPairsPerSpectrum[aligns[i].spec2]++;
      if (bestScoreIdx < 0 || bestScore < curScore)
      {
        bestScoreIdx = i;
        bestScore = curScore;
        specToAdd = aligns[bestScoreIdx].spec1;
      }
    }
    // Select starting spectrum by highest number of pairs
    pIter = toProcess.begin();
    specToAdd = *pIter;
    bestScore = numPairsPerSpectrum[specToAdd];
    pIter++;
    for (; pIter != toProcess.end(); pIter++)
    {
      if (numPairsPerSpectrum[*pIter] > bestScore)
      {
        bestScore = numPairsPerSpectrum[*pIter];
        specToAdd = *pIter;
      }
      numPairsPerSpectrum[*pIter] = 0;
    }
    //	cout<<"  - (SplitPairs) Component initialized with spectrum "<<specToAdd<<" with "<<bestScore<<" neighbors\n";

    flipScores[specToAdd] = -1; // Arbitrarily decide that this spectrum is not flipped
    pIter = toProcess.begin();
    while (*pIter != specToAdd)
    {
      pIter++;
    }
    toProcess.erase(pIter);
    while (specToAdd >= 0)
    { // Process the inclusion of specToAdd to the oriented component
      // Flip/don't-flip specToAdd
      int flipDir = 0; // Spectrum not flipped (default)
      if (flipScores[specToAdd] > 0)
      {
        specSet[specToAdd] = specSetRev[specToAdd];
        flipDir = 1;
        specFlipped[specToAdd] = true;
        if (labelsP and labelsP->size()>specToAdd)
          (*labelsP)[specToAdd].reverse();
      }
      specsInComponent++;

      // Process matches between newly included specToAdd and its alignments to spectra not yet included in the oriented component
      for (pIter = alignsEntries[specToAdd].begin(); pIter
          != alignsEntries[specToAdd].end(); pIter++)
      {
        curPair = aligns[*pIter];
//        isASP = (abs(curPair.shift1) <= pmTol || abs(curPair.shift2) <= pmTol);
  			isASP = ((fabs(curPair.shift1)<=pmTol and fabs(curPair.shift2)<=maxModMass) or (fabs(curPair.shift1)<=maxModMass and fabs(curPair.shift2)<=pmTol));

//if(curPair.spec1==6876 or curPair.spec2==6876)
//	cerr<<" -> Pair ("<<curPair.spec1<<","<<curPair.spec2<<") with shifts "<<curPair.shift1<<"/"<<curPair.shift2<<", isASP = "<<isASP<<endl;

  			//DEBUG_VAR(isASP);
        if (!matchComputed[*pIter])
        {
          // otherSpec is not in the component yet. Compute how the addition of specToAdd changes otherSpec's orientation selection.
          int otherSpec, otherSpecPos; // otherSpec=index of current paired spectrum, otherSpecPos=position of otherSpec index in aligns
          if (curPair.spec1 == specToAdd)
          {
            otherSpec = curPair.spec2;
            otherSpecPos = 1;
          }
          else
          {
            otherSpec = curPair.spec1;
            otherSpecPos = 0;
          }

          // SANITY CHECK !
          if (specSet[specToAdd].size() == 0 || specSet[otherSpec].size() == 0) {
            continue;
          }

          TwoValues<float> matchScore1(0, 0), matchScore2(0, 0);
          vector<int> indices;
          vector<Spectrum> results(4);

          TwoValues<float> score1, score2;
          vector<int> matchA1, matchA2, matchB1, matchB2;
          bool gluesAccepted = false; // A pair's glues are accepted if the pair matches the minimum number of peaks (depending on orientation of the paired spectra)
          int extraMatchPeaks = 0; // ASP pairs always match either/both PRMs at 0/18 or/and PM-19/PM-1 so these have higher numbers of peaks that need to match

          //cerr<<"Spectrum "<<specToAdd<<":\n"; specSet[specToAdd].output(cerr);
          //cerr<<"Spectrum "<<otherSpec<<":\n"; specSet[otherSpec].output(cerr);

          if (isASP)
          {
            if (abs(curPair.shift1) <= peakTol)
            {
              extraMatchPeaks++; // matching 0/18 does not count towards achieving minMatchedPeaks
            }
            if (abs(curPair.shift2) <= peakTol)
            {
              extraMatchPeaks++; // matching PM-19/PM-1 does not count towards achieving minMatchedPeaks
            }

            tmpModPos[*pIter][0] = SpectrumAlignment(&specSet[specToAdd],
                                               &specSet[otherSpec],
                                               peakTol,
                                               &results[0],
                                               &results[1],
                                               maxAAjump,
                                               penalty_sameVert,
                                               penalty_ptm,
                                               forceSymmetry,
                                               false);
            /*if(specToAdd==9642) {
             cerr<<"  ==> ASP ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<results[0].size()<<" / "<<minMatchedPeaks+extraMatchPeaks<<"\n";
             if(otherSpec==9630)
             for(unsigned int i=0; i<results[0].size(); i++)
             cerr<<" ------ ("<<results[0][i][0]<<","<<results[0][i][1]<<")\n";
             }*/
            if (results[0].size() < minMatchedPeaks + extraMatchPeaks)
            {
              tmpMatches[0][*pIter].resize(0);
            }
            else
            {
              for (unsigned int i = 0; i < results[0].size(); i++)
              {
                matchScore1[0] += results[0][i][1];
                matchScore2[0] += results[1][i][1];
              }
              specSet[specToAdd].massesToIndices(results[0],
                                                 indices,
                                                 peakTol);
              tmpMatches[0][*pIter].resize(indices.size());
              for (unsigned int i = 0; i < indices.size(); i++)
              {
                tmpMatches[0][*pIter][i][1 - otherSpecPos] = indices[i];
              }
              specSet[otherSpec].massesToIndices(results[1],
                                                 indices,
                                                 peakTol);
              for (unsigned int i = 0; i < indices.size(); i++)
              {
                tmpMatches[0][*pIter][i][otherSpecPos] = indices[i];
                if (tmpMatches[0][*pIter][i][otherSpecPos] == -1)
                {
                  DEBUG_MSG("--------------------");
                  for (unsigned int d = 0; d < indices.size(); d++)
                  {
                    DEBUG_MSG(d << " " << results[1][d][0]);
                  }
                  DEBUG_MSG("--------------------");
                  DEBUG_VAR(curPair.shift1);
                  DEBUG_VAR(curPair.shift2);
                  DEBUG_MSG("ABORTING: tmpMatches[0][*pIter][i][otherSpecPos] = -1");
                  abort();
                }
              }
              gluesAccepted = true;
            }

            //if(specToAdd==111 and otherSpec==246)
            //	{ cerr<<"Direct match:\n";
            //	for(unsigned int i=0;i<results[0].size();i++) cerr<<"  ("<<results[0][i][0]<<","<<results[0][i][1]<<")/("<<results[1][i][0]<<","<<results[1][i][1]<<")\n"; }

            tmpModPos[*pIter][1] = SpectrumAlignment(&specSet[specToAdd],
                                               &specSetRev[otherSpec],
                                               peakTol,
                                               &results[2],
                                               &results[3],
                                               maxAAjump,
                                               penalty_sameVert,
                                               penalty_ptm,
                                               forceSymmetry,
                                               false);
            /*if(specToAdd==9642){
             cerr<<"  ==> ASP rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<results[2].size()<<" / "<<minMatchedPeaks+extraMatchPeaks<<"\n";
             if(otherSpec==9630)
             for(unsigned int i=0; i<results[0].size(); i++)
             cerr<<" ------ ("<<results[0][i][0]<<","<<results[0][i][1]<<")\n";
             }*/
            if (results[2].size() < minMatchedPeaks + extraMatchPeaks)
            {
              tmpMatches[1][*pIter].resize(0);
            }
            else
            {
              for (unsigned int i = 0; i < results[2].size(); i++)
              {
                matchScore1[1] += results[2][i][1];
                matchScore2[1] += results[3][i][1];
              }
              specSet[specToAdd].massesToIndices(results[2],
                                                 indices,
                                                 peakTol);
              tmpMatches[1][*pIter].resize(indices.size());
              for (unsigned int i = 0; i < indices.size(); i++)
              {
                tmpMatches[1][*pIter][i][1 - otherSpecPos] = indices[i];
              }
              specSetRev[otherSpec].massesToIndices(results[3],
                                                    indices,
                                                    peakTol);
              for (unsigned int i = 0; i < indices.size(); i++)
              {
                tmpMatches[1][*pIter][i][otherSpecPos] = indices[i];
              }
              gluesAccepted = true;
            }
            //if(specToAdd==111 and otherSpec==246)
            //	{ cerr<<"Reverse match:\n";
            //	for(unsigned int i=0;i<results[2].size();i++) cerr<<"  ("<<results[2][i][0]<<","<<results[2][i][1]<<")/("<<results[3][i][0]<<","<<results[3][i][1]<<")\n"; }
          }
          else
          {
            // Select shift for direct orientation
            float shift = 0;
            if (specToAdd == curPair.spec1)
            {
              shift = curPair.shift1;
            }
            else
            {
              shift = -curPair.shift1;
            }
            ScoreOverlap6(specSet[specToAdd],
                          specSet[otherSpec],
                          shift,
                          peakTol,
                          matchA1,
                          matchA2);
            //if(specToAdd==9642)
            //	cerr<<"  ==> PA A1 ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchA1.size()<<" / "<<minMatchedPeaks<<"\n";
            if (matchA1.size() < minMatchedPeaks)
            {
              matchA1.resize(0);
              matchA2.resize(0);
            }
            else
            {
              gluesAccepted = true;
            }
            score1.set(0, 0);
            for (unsigned int i = 0; i < matchA1.size(); i++)
            {
              score1[0] += specSet[specToAdd][matchA1[i]][1];
              score1[1] += specSet[otherSpec][matchA2[i]][1];
            }
            if (specToAdd == curPair.spec1)
            {
              shift = curPair.shift2;
            }
            else
            {
              shift = -curPair.shift2;
            }
            ScoreOverlap6(specSet[specToAdd],
                          specSet[otherSpec],
                          shift,
                          peakTol,
                          matchB1,
                          matchB2);
            //if(specToAdd==9642)
            //	cerr<<"  ==> PA B1 ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchB1.size()<<" / "<<minMatchedPeaks<<"\n";
            if (matchB1.size() < minMatchedPeaks)
            {
              matchB1.resize(0);
              matchB2.resize(0);
            }
            else
            {
              gluesAccepted = true;
            }
            score2.set(0, 0);
            for (unsigned int i = 0; i < matchB1.size(); i++)
            {
              score2[0] += specSet[specToAdd][matchB1[i]][1];
              score2[1] += specSet[otherSpec][matchB2[i]][1];
            }

            if ((score1[0] / specScores[specToAdd] + score1[1]
                / specScores[otherSpec]) > (score2[0] / specScores[specToAdd]
                + score2[1] / specScores[otherSpec]))
            {
              matchScore1[0] = score1[0];
              matchScore2[0] = score1[1];
              tmpMatches[0][*pIter].resize(matchA1.size());
              for (unsigned int i = 0; i < matchA1.size(); i++)
              {
                tmpMatches[0][*pIter][i][1 - otherSpecPos] = matchA1[i];
                tmpMatches[0][*pIter][i][otherSpecPos] = matchA2[i];
              }
            }
            else
            {
              matchScore1[0] = score2[0];
              matchScore2[0] = score2[1];
              tmpMatches[0][*pIter].resize(matchB1.size());
              for (unsigned int i = 0; i < matchB1.size(); i++)
              {
                tmpMatches[0][*pIter][i][1 - otherSpecPos] = matchB1[i];
                tmpMatches[0][*pIter][i][otherSpecPos] = matchB2[i];
              }
            }

            // Select shift for reversed orientation
            if (specToAdd == curPair.spec1)
            {
              shift = curPair.shift1;
            }
            else
            {
              shift = -curPair.shift1;
            }
            ScoreOverlap6(specSet[specToAdd],
                          specSetRev[otherSpec],
                          shift,
                          peakTol,
                          matchA1,
                          matchA2);

            if (ignoreReversals) {
              matchA1.resize(0);
            }
            //if(specToAdd==9642)
            //	cerr<<"  ==> PA A1 rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchA1.size()<<" / "<<minMatchedPeaks<<"\n";
            if (matchA1.size() < minMatchedPeaks)
            {
              matchA1.resize(0);
              matchA2.resize(0);
            }
            else
            {
              gluesAccepted = true;
            }
            score1.set(0, 0);
            for (unsigned int i = 0; i < matchA1.size(); i++)
            {
              score1[0] += specSet[specToAdd][matchA1[i]][1];
              score1[1] += specSetRev[otherSpec][matchA2[i]][1];
            }
            if (specToAdd == curPair.spec1)
            {
              shift = curPair.shift2;
            }
            else
            {
              shift = -curPair.shift2;
            }
            ScoreOverlap6(specSet[specToAdd],
                          specSetRev[otherSpec],
                          shift,
                          peakTol,
                          matchB1,
                          matchB2);

            if (ignoreReversals) {
              matchB1.resize(0);
            }
            //if(specToAdd==9642)
            //	cerr<<"  ==> PA B1 rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchB1.size()<<" / "<<minMatchedPeaks<<"\n";
            if (matchB1.size() < minMatchedPeaks)
            {
              matchB1.resize(0);
              matchB2.resize(0);
            }
            else
            {
              gluesAccepted = true;
            }
            score2.set(0, 0);
            for (unsigned int i = 0; i < matchB1.size(); i++)
            {
              score2[0] += specSet[specToAdd][matchB1[i]][1];
              score2[1] += specSetRev[otherSpec][matchB2[i]][1];
            }

            if ((score1[0] / specScores[specToAdd] + score1[1]
                / specScores[otherSpec]) > (score2[0] / specScores[specToAdd]
                + score2[1] / specScores[otherSpec]))
            {
              matchScore1[1] = score1[0];
              matchScore2[1] = score1[1];
              tmpMatches[1][*pIter].resize(matchA1.size());
              for (unsigned int i = 0; i < matchA1.size(); i++)
              {
                tmpMatches[1][*pIter][i][1 - otherSpecPos] = matchA1[i];
                tmpMatches[1][*pIter][i][otherSpecPos] = matchA2[i];
              }
            }
            else
            {
              matchScore1[1] = score2[0];
              matchScore2[1] = score2[1];
              tmpMatches[1][*pIter].resize(matchB1.size());
              for (unsigned int i = 0; i < matchB1.size(); i++)
              {
                tmpMatches[1][*pIter][i][1 - otherSpecPos] = matchB1[i];
                tmpMatches[1][*pIter][i][otherSpecPos] = matchB2[i];
              }
            }
          }

          //				flipScores[otherSpec] += (matchScore2[1]-matchScore2[0])/specScores[otherSpec];
          if (gluesAccepted)
          {
            flipScores[otherSpec] += (matchScore2[1] - matchScore2[0])
                / specScores[otherSpec] + (matchScore1[1] - matchScore1[0])
                / specScores[specToAdd];
            numPairsPerSpectrum[otherSpec]++;
          }
          matchComputed[*pIter] = true;
        }
        else
        {
          // Copy adequate entry (according to flipDir) from tmpMatch to matches
          //cerr<<" * copying match "<<curPair.spec1<<" vs "<<curPair.spec2<<endl;
          matches[*pIter].resize(tmpMatches[flipDir][*pIter].size());
          for (unsigned int i = 0; i < matches[*pIter].size(); i++)
          {
            matches[*pIter][i] = tmpMatches[flipDir][*pIter][i];
          }
          //				if(alignRatios and flipDir==1) { float tmp=(*alignRatios)[*pIter][0]; (*alignRatios)[*pIter][0]=(*alignRatios)[*pIter][1]; (*alignRatios)[*pIter][1]=tmp; }
          if (isASP)
          {
            modPos[*pIter] = tmpModPos[*pIter][flipDir];
          }
        }
      }

      if (toProcess.size() == 0)
      {
        specToAdd = -1;
      }
      else
      {
        //			list<int>::iterator bestNewSpec=toProcess.begin();  pIter=bestNewSpec;  pIter++;
        //			while(pIter!=toProcess.end()) { if(fabs(flipScores[*pIter])>fabs(flipScores[*bestNewSpec])) bestNewSpec=pIter; pIter++; }
        //			specToAdd = *bestNewSpec;   toProcess.erase(bestNewSpec);
        list<int>::iterator bestNewSpec = toProcess.begin();
        while (bestNewSpec != toProcess.end())
        {
          if (numPairsPerSpectrum[*bestNewSpec] >= min(minEdgesToComponent,
                                                       specsInComponent))
          {
            break;
          }
          bestNewSpec++;
        }
        if (bestNewSpec != toProcess.end())
        {
          pIter = bestNewSpec;
          pIter++;
          while (pIter != toProcess.end())
          {
            if (numPairsPerSpectrum[*pIter] >= min(minEdgesToComponent,
                                                   specsInComponent)
                && fabs(flipScores[*pIter]) > fabs(flipScores[*bestNewSpec]))
            {
              bestNewSpec = pIter;
            }
            pIter++;
          }
          specToAdd = *bestNewSpec;
          toProcess.erase(bestNewSpec);
        }
        else
        {
          specToAdd = -1;
        }
      }
    }
    //	cout<<"  - (SplitPairs) Component includes "<<specsInComponent<<" spectra.\n";

    if (alignStats)
    { // Cannot be computed correctly until all flipping decisions are made
      alignStats->resize(aligns.size());
      for (unsigned int pairIdx = 0; pairIdx < aligns.size(); pairIdx++)
      {
        float score1 = 0, score2 = 0;
        int s1 = aligns[pairIdx].spec1, s2 = aligns[pairIdx].spec2;
        (*alignStats)[pairIdx].resize(3);
        for (unsigned int i = 0; i < matches[pairIdx].size(); i++)
        {
          score1 += specSet[s1][matches[pairIdx][i][0]][1];
          score2 += specSet[s2][matches[pairIdx][i][1]][1];
        }
        (*alignStats)[pairIdx][0] = score1 / specScores[s1];
        (*alignStats)[pairIdx][1] = score2 / specScores[s2];
        (*alignStats)[pairIdx][2] = matches[pairIdx].size();
      }
    }
  }

  // -------------------------------------------------------------------------
  //
  //  SplitPairs - splits every (i,j) pair in aligns into either (i,j)/(rev(i),rev(j)) or (i,rev(j))/(rev(i),j)
  //
  // -------------------------------------------------------------------------
  void SplitPairs(SpecSet &specSet,
                  SpectrumPairSet &aligns,
                  float peakTol,
                  int maxAAjump,
                  float penalty_sameVert,
                  float penalty_ptm,
                  bool forceSymmetry,
                  SpecSet &specSetSplit,
                  SpectrumPairSet &alignsNew,
                  vector<vector<TwoValues<int> > > &matches,
                  vector<bool> &pairFlipped,
                  vector<vector<float> > *dbg_matchScores,
                  ofstream *debug)
  {
    SplitSpectra(specSet, specSetSplit);
    SplitAligns(specSetSplit,
                aligns,
                peakTol,
                maxAAjump,
                penalty_sameVert,
                penalty_ptm,
                forceSymmetry,
                alignsNew,
                matches,
                pairFlipped,
                dbg_matchScores,
                debug);
  }

  // -------------------------------------------------------------------------
  void SplitAligns(SpecSet &specSetSplit,
                   SpectrumPairSet &aligns,
                   float peakTol,
                   int maxAAjump,
                   float penalty_sameVert,
                   float penalty_ptm,
                   bool forceSymmetry,
                   SpectrumPairSet &alignsNew,
                   vector<vector<TwoValues<int> > > &matches,
                   vector<bool> &pairFlipped,
                   vector<vector<float> > *dbg_matchScores,
                   ofstream *debug)
  {
    TwoValues<float> matchScore(0, 0), matchScore1(0, 0), matchScore2(0, 0);
#ifdef DEBUG
    vector<vector<float> > modPos(aligns.size());
#endif

    /*	specSetSplit.resize(2*specSet.size());
     for(unsigned int i=0; i<specSet.size(); i++) {
     specSetSplit[2*i] = specSet[i];
     specSet[i].reverse(0, &specSetSplit[2*i+1]);
     }*/
    if (dbg_matchScores != 0)
    {
      (*dbg_matchScores).resize(aligns.size());
      for (unsigned int i = 0; i < aligns.size(); i++)
      {
        (*dbg_matchScores)[i].resize(4);
      }
    }

    int spec1, spec1rev, spec2, spec2rev, newIdx;
    vector<Spectrum> results(4);
    alignsNew.resize(2 * aligns.size());
    matches.resize(2 * aligns.size());
    pairFlipped.resize(aligns.size());
    for (unsigned int alignIdx = 0; alignIdx < aligns.size(); alignIdx++)
    {
      spec1 = 2 * aligns[alignIdx].spec1;
      spec1rev = 2 * aligns[alignIdx].spec1 + 1;
      spec2 = 2 * aligns[alignIdx].spec2;
      spec2rev = 2 * aligns[alignIdx].spec2 + 1;
      newIdx = 2 * alignIdx;
      matchScore.set(0, 0);
      matchScore1.set(0, 0);
      matchScore2.set(0, 0);

#ifdef DEBUG
      (*debug)<<"Entry pair: "<<aligns[alignIdx].spec1+1<<", "<<aligns[alignIdx].spec2+1<<endl;
      (*debug)<<"Direct matching results:\n";
#endif
#ifdef DEBUG
      modPos[alignIdx].resize(2);
      modPos[alignIdx][0]=SpectrumAlignment(&specSetSplit[spec1],&specSetSplit[spec2],peakTol,&results[0],&results[1],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,*debug);
#else
      SpectrumAlignment(&specSetSplit[spec1],
                  &specSetSplit[spec2],
                  peakTol,
                  &results[0],
                  &results[1],
                  maxAAjump,
                  penalty_sameVert,
                  penalty_ptm,
                  forceSymmetry);
#endif
      for (unsigned int j = 0; j < results[0].size(); j++)
      {
        //cerr<<"  ("<<results[0][j][0]<<","<<results[0][j][1]<<") / ("<<results[1][j][0]<<","<<results[1][j][1]<<")\n";
        matchScore[0] += results[0][j][1] + results[1][j][1];
        matchScore1[0] += results[0][j][1];
        matchScore2[0] += results[1][j][1];
      }

#ifdef DEBUG
      (*debug)<<"Reversed matching results:\n";
      modPos[alignIdx][1]=SpectrumAlignment(&specSetSplit[spec1rev],&specSetSplit[spec2],peakTol,&results[2],&results[3],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,*debug);
#else
      SpectrumAlignment(&specSetSplit[spec1rev],
                  &specSetSplit[spec2],
                  peakTol,
                  &results[2],
                  &results[3],
                  maxAAjump,
                  penalty_sameVert,
                  penalty_ptm,
                  forceSymmetry);
#endif
      for (unsigned int j = 0; j < results[2].size(); j++)
      {
        //cerr<<"  ("<<results[2][j][0]<<","<<results[2][j][1]<<") / ("<<results[3][j][0]<<","<<results[3][j][1]<<")\n";
        matchScore[1] += results[2][j][1] + results[3][j][1];
        matchScore1[1] += results[2][j][1];
        matchScore2[1] += results[3][j][1];
      }

      if (dbg_matchScores != 0)
      {
        (*dbg_matchScores)[alignIdx][0] = matchScore1[0];
        (*dbg_matchScores)[alignIdx][1] = matchScore1[1];
        (*dbg_matchScores)[alignIdx][2] = matchScore2[0];
        (*dbg_matchScores)[alignIdx][3] = matchScore2[1];
      }
#ifdef DEBUG
      (*debug)<<"matchScore[0] = ["<<matchScore[0]<<","<<matchScore[1]<<","<<max(matchScore[0],matchScore[1])/min(matchScore[0],matchScore[1])<<"], matchScore1 = ["<<matchScore1[0]<<","<<matchScore1[1]<<","<<max(matchScore1[0],matchScore1[1])/min(matchScore1[0],matchScore1[1])<<"], matchScore2 = ["<<matchScore2[0]<<","<<matchScore2[1]<<","<<max(matchScore2[0],matchScore2[1])/min(matchScore2[0],matchScore2[1])<<"]\n";
#endif

      if (matchScore[1] > matchScore[0])
      {
#ifdef DEBUG
        (*debug)<<"Additional matching of the first spectrum to the reversed second spectrum\n";
        SpectrumAlignment(&specSetSplit[spec1],&specSetSplit[spec2rev],peakTol,&results[0],&results[1],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,*debug);
#else
        SpectrumAlignment(&specSetSplit[spec1],
                    &specSetSplit[spec2rev],
                    peakTol,
                    &results[0],
                    &results[1],
                    maxAAjump,
                    penalty_sameVert,
                    penalty_ptm,
                    forceSymmetry);
#endif
        alignsNew[newIdx].spec1 = spec1;
        alignsNew[newIdx].spec2 = spec2rev;
        alignsNew[newIdx + 1].spec1 = spec1rev;
        alignsNew[newIdx + 1].spec2 = spec2;
        pairFlipped[alignIdx] = true;
      }
      else
      {
#ifdef DEBUG
        (*debug)<<"Additional matching of the reversed first spectrum to the reversed second spectrum\n";
        SpectrumAlignment(&specSetSplit[spec1rev],&specSetSplit[spec2rev],peakTol,&results[2],&results[3],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,*debug);
#else
        SpectrumAlignment(&specSetSplit[spec1rev],
                    &specSetSplit[spec2rev],
                    peakTol,
                    &results[2],
                    &results[3],
                    maxAAjump,
                    penalty_sameVert,
                    penalty_ptm,
                    forceSymmetry);
#endif
        alignsNew[newIdx].spec1 = spec1;
        alignsNew[newIdx].spec2 = spec2;
        alignsNew[newIdx + 1].spec1 = spec1rev;
        alignsNew[newIdx + 1].spec2 = spec2rev;
        pairFlipped[alignIdx] = false;
      }

      vector<int> indices;
      for (unsigned int j = 0; j < 2; j++)
      {
        specSetSplit[alignsNew[newIdx + j].spec1].massesToIndices(results[2 * j],
                                                                  indices,
                                                                  peakTol);
        matches[newIdx + j].resize(indices.size());
        for (unsigned int i = 0; i < indices.size(); i++)
        {
          matches[newIdx + j][i][0] = indices[i];
        }
        specSetSplit[alignsNew[newIdx + j].spec2].massesToIndices(results[2 * j
            + 1], indices, peakTol);
        for (unsigned int i = 0; i < indices.size(); i++)
        {
          matches[newIdx + j][i][1] = indices[i];
        }
      }
    }
#ifdef DEBUG
    Save_binArray("debug_modPos.bna", modPos);
#endif
  }

  // -------------------------------------------------------------------------
  void SplitLabels(vector<SpectrumPeakLabels> &labels, vector<
      SpectrumPeakLabels> &newLabels)
  {
    newLabels.resize(2 * labels.size());
    for (unsigned int i = 0; i < labels.size(); i++)
    {
      newLabels[2 * i] = labels[i];
      newLabels[2 * i + 1] = labels[i];
      newLabels[2 * i + 1].reverse();
    }
  }

  // -------------------------------------------------------------------------
  //
  //  SplitPairs2 - Decides on a consensus orientation for  a set of edges.
  //    Order of instantiation of orientations is: initialize orientations with the
  //    highest scoring unused edge; add vertex with highest scoring connection to
  //    component until no vertex is connected to the component.
  //
  //    NOTE: Flips spectra in specSet acoording to assigned orientations
  //
  // -------------------------------------------------------------------------
  void SplitPairs2(SpecSet &specSet,
                   SpectrumPairSet &aligns,
                   float peakTol,
                   int maxAAjump,
                   float penalty_sameVert,
                   float penalty_ptm,
                   vector<vector<TwoValues<int> > > &matches,
                   vector<bool> &specFlipped,
                   vector<float> &modPos,
                   bool forceSymmetry,
                   vector<SpectrumPeakLabels> *labelsP)
  {
    vector<list<int> > alignsEntries(specSet.size()); // List of pairs in aligns that incide on the corresponding vertex
    vector<float> flipScores(specSet.size()); // Accumulated don't_flip/flip score balance per vertex
    SpecSet specSetRev;
    specSetRev.resize(specSet.size()); // Holds reversed versions of the spectra in the component
    list<int> toProcess; // Indices of spectra left to process
    TwoValues<vector<vector<TwoValues<int> > > > tmpMatches; // Temporary matches when a neighbor spectrum S1 in the i-th pair is matched directly/reversed (its pair S2 is assumed to already be in the component)
    vector<TwoValues<float> > tmpModPos(aligns.size()); // Temporary storage for the position of the modification
    vector<bool> matchComputed(aligns.size()); // Indicates whether the i-th entry in tmpAligns was already computed
    matches.resize(aligns.size());
    modPos.resize(aligns.size());
    if (aligns.size() == 0)
      return;

    //cerr<<" * input aligns: ";
    //for(unsigned int i=0; i<aligns.size(); i++) cerr<<"["<<aligns[i].spec1<<","<<aligns[i].spec2<<"]"; cerr<<endl;

    // Initializations
    tmpMatches[0].resize(aligns.size());
    tmpMatches[1].resize(aligns.size());
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      if (specSet[aligns[i].spec1].size() != specSetRev[aligns[i].spec1].size())
      {
        specSet[aligns[i].spec1].reverse(0, &specSetRev[aligns[i].spec1]);
        toProcess.push_back(aligns[i].spec1);
      }
      if (specSet[aligns[i].spec2].size() != specSetRev[aligns[i].spec2].size())
      {
        specSet[aligns[i].spec2].reverse(0, &specSetRev[aligns[i].spec2]);
        toProcess.push_back(aligns[i].spec2);
      }
      matchComputed[i] = false;
      alignsEntries[aligns[i].spec1].push_back(i);
      alignsEntries[aligns[i].spec2].push_back(i);
    }
    list<int>::iterator pIter;
    for (pIter = toProcess.begin(); pIter != toProcess.end(); pIter++)
    {
      flipScores[*pIter] = 0;
    }

    // Use edge scores as percentages of spectrum score to select the edge/vertex to initialize the direction selection
    vector<float> specScores(specSet.size());
    float bestScore, curScore;
    int bestScoreIdx = -1;
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specScores[i] = 0;
      for (unsigned int j = 0; j < specSet[i].size(); j++)
      {
        specScores[i] += specSet[i][j][1];
      }
    }
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      curScore = min(aligns[i].score1 / specScores[aligns[i].spec1],
                     aligns[i].score2 / specScores[aligns[i].spec2]);
      if (bestScoreIdx < 0 || bestScore < curScore)
      {
        bestScoreIdx = i;
        bestScore = curScore;
      }
    }

    int specToAdd = aligns[bestScoreIdx].spec1; // Index of the new spectrum being added to the component
    //cerr<<" * initialized orientations component with vertex "<<specToAdd<<endl;
    flipScores[specToAdd] = -1; // Arbitrarily decide that this spectrum is not flipped
    pIter = toProcess.begin();
    while (*pIter != specToAdd)
    {
      pIter++;
    }
    toProcess.erase(pIter);
    while (specToAdd >= 0)
    {
      int flipDir = 0; // Spectrum not flipped (default)
      if (flipScores[specToAdd] > 0)
      {
        specSet[specToAdd] = specSetRev[specToAdd];
        flipDir = 1;
        specFlipped[specToAdd] = true;
        if (labelsP and labelsP->size()>specToAdd)
        {
          (*labelsP)[specToAdd].reverse();
        }
      }
      //cerr<<" * spectrum "<<specToAdd; if(flipDir) cerr<<" reversed"; else cerr<<" not reversed"; cerr<<", score = "<<flipScores[specToAdd]<<endl;

      for (pIter = alignsEntries[specToAdd].begin(); pIter
          != alignsEntries[specToAdd].end(); pIter++)
      {
        if (!matchComputed[*pIter])
        {
          // otherSpec is not in the component yet. Compute how the addition of specToAdd changes otherSpec's orientation selection.
          int otherSpec, otherSpecPos; // otherSpec=index of current paired spectrum, otherSpecPos=position of otherSpec index in aligns
          if (aligns[*pIter].spec1 == specToAdd)
          {
            otherSpec = aligns[*pIter].spec2;
            otherSpecPos = 1;
          }
          else
          {
            otherSpec = aligns[*pIter].spec1;
            otherSpecPos = 0;
          }

          TwoValues<float> matchScore1(0, 0), matchScore2(0, 0);
          vector<int> indices;
          vector<Spectrum> results(4);

          tmpModPos[*pIter][0] = SpectrumAlignment(&specSet[specToAdd],
                                             &specSet[otherSpec],
                                             peakTol,
                                             &results[0],
                                             &results[1],
                                             maxAAjump,
                                             penalty_sameVert,
                                             penalty_ptm,
                                             forceSymmetry,
                                             true);
          for (unsigned int i = 0; i < results[0].size(); i++)
          {
            matchScore1[0] += results[0][i][1];
            matchScore2[0] += results[1][i][1];
          }
          specSet[specToAdd].massesToIndices(results[0],
                                             indices,
                                             peakTol);
          tmpMatches[0][*pIter].resize(indices.size());
          for (unsigned int i = 0; i < indices.size(); i++)
          {
            if (indices[i] < 0)
            {
              cerr << "ERROR: Could not find mass "
                  << results[0][i][0] << " in spectrum 1 ("
                  << specToAdd << ", " << specSet[specToAdd].size() << ")!\n";
              specSet[specToAdd].output(cerr);
              exit(-1);
            }
            tmpMatches[0][*pIter][i][1 - otherSpecPos] = indices[i];
          }
          specSet[otherSpec].massesToIndices(results[1],
                                             indices,
                                             peakTol);
          for (unsigned int i = 0; i < indices.size(); i++)
          {
            if (indices[i] < 0)
            {
              cerr << "ERROR: Could not find mass "
                  << results[1][i][0] << " in spectrum 2 ("
                  << otherSpec << ", " << specSet[otherSpec].size() << ")!\n";
              specSet[otherSpec].output(cerr);
              exit(-1);
            }
            tmpMatches[0][*pIter][i][otherSpecPos] = indices[i];
          }

          //if(specToAdd==111 and otherSpec==246)
          //	{ cerr<<"Direct match:\n";
          //	for(unsigned int i=0;i<results[0].size();i++) cerr<<"  ("<<results[0][i][0]<<","<<results[0][i][1]<<")/("<<results[1][i][0]<<","<<results[1][i][1]<<")\n"; }

          tmpModPos[*pIter][1] = SpectrumAlignment(&specSet[specToAdd],
                                             &specSetRev[otherSpec],
                                             peakTol,
                                             &results[2],
                                             &results[3],
                                             maxAAjump,
                                             penalty_sameVert,
                                             penalty_ptm,
                                             forceSymmetry,
                                             true);
          for (unsigned int i = 0; i < results[2].size(); i++)
          {
            matchScore1[1] += results[2][i][1];
            matchScore2[1] += results[3][i][1];
          }
          specSet[specToAdd].massesToIndices(results[2],
                                             indices,
                                             peakTol);
          tmpMatches[1][*pIter].resize(indices.size());
          for (unsigned int i = 0; i < indices.size(); i++)
          {
            if (indices[i] < 0)
            {
              cerr << "ERROR: Could not find mass "
                  << results[2][i][0] << " in spectrum 1 ("
                  << specToAdd << ", " << specSet[specToAdd].size() << ")!\n";
              specSet[specToAdd].output(cerr);
              exit(-1);
            }
            tmpMatches[1][*pIter][i][1 - otherSpecPos] = indices[i];
          }
          specSetRev[otherSpec].massesToIndices(results[3],
                                                indices,
                                                peakTol);
          for (unsigned int i = 0; i < indices.size(); i++)
          {
            if (indices[i] < 0)
            {
              cerr << "ERROR: Could not find mass "
                  << results[3][i][0] << " in spectrum 2r ("
                  << otherSpec << ", " << specSetRev[otherSpec].size()
                  << ")!\n";
              specSetRev[otherSpec].output(cerr);
              specSet[otherSpec].output(cerr);
              exit(-1);
            }
            tmpMatches[1][*pIter][i][otherSpecPos] = indices[i];
          }

          //if(specToAdd==111 and otherSpec==246)
          //	{ cerr<<"Reverse match:\n";
          //	for(unsigned int i=0;i<results[2].size();i++) cerr<<"  ("<<results[2][i][0]<<","<<results[2][i][1]<<")/("<<results[3][i][0]<<","<<results[3][i][1]<<")\n"; }

          //				flipScores[otherSpec] += (matchScore2[1]-matchScore2[0])/specScores[otherSpec];
          flipScores[otherSpec] += (matchScore2[1] - matchScore2[0])
              / specScores[otherSpec] + (matchScore1[1] - matchScore1[0])
              / specScores[specToAdd];
          matchComputed[*pIter] = true;
          //cerr<<" * match "<<aligns[*pIter].spec1<<" vs "<<aligns[*pIter].spec2<<", matchScore1=("<<matchScore1[0]<<","<<matchScore1[1]<<")/"<<specScores[specToAdd]<<", matchScore2=("<<matchScore2[0]<<","<<matchScore2[1]<<")/"<<specScores[otherSpec]<<", flipScores["<<otherSpec<<"]="<<flipScores[otherSpec]<<endl;
        }
        else
        {
          // Simply copy adequate entry (acording to flipDir) from tmpMatch to matches
          //cerr<<" * copying match "<<aligns[*pIter].spec1<<" vs "<<aligns[*pIter].spec2<<endl;
          matches[*pIter].resize(tmpMatches[flipDir][*pIter].size());
          for (unsigned int i = 0; i < matches[*pIter].size(); i++)
            matches[*pIter][i] = tmpMatches[flipDir][*pIter][i];
          modPos[*pIter] = tmpModPos[*pIter][flipDir];
        }
      }

      if (toProcess.size() == 0)
        specToAdd = -1;
      else
      {
        list<int>::iterator bestNewSpec = toProcess.begin();
        pIter = bestNewSpec;
        pIter++;
        while (pIter != toProcess.end())
        {
          if (fabs(flipScores[*pIter]) > fabs(flipScores[*bestNewSpec]))
            bestNewSpec = pIter;
          pIter++;
        }
        specToAdd = *bestNewSpec;
        toProcess.erase(bestNewSpec);
      }
    }
  }

  // -------------------------------------------------------------------------
  //
  //  SplitPairs3 - Like SplitPairs2 but also processes Partial Overlap alignments (alignsPA)
  //    Decides on a consensus orientation for  a set of edges.
  //    Order of instantiation of orientations is: initialize orientations with the
  //    highest scoring unused edge; add vertex with highest scoring connection to
  //    component until no vertex is connected to the com ponent.
  //
  //    NOTE: Flips spectra in specSet acoording to assigned orientations
  //
  // -------------------------------------------------------------------------
  void SplitPairs3(SpecSet &specSet,
                   SpectrumPairSet &aligns,
                   SpectrumPairSet &alignsPA,
                   float peakTol,
                   int maxAAjump,
                   float penalty_sameVert,
                   float penalty_ptm,
                   vector<vector<TwoValues<int> > > &matches,
                   vector<vector<TwoValues<int> > > &matchesPA,
                   vector<bool> &specFlipped,
                   vector<float> &modPos,
                   bool forceSymmetry,
                   vector<SpectrumPeakLabels> *labelsP,
                   vector<TwoValues<float> > *alignRatios,
                   vector<TwoValues<float> > *alignRatiosPA)
  {
    //	vector<list<int> > alignsEntries(specSet.size());       // List of pairs in aligns that incide on the corresponding vertex
    vector<TwoValues<list<int> > > alignsEntries(specSet.size()); // List of pairs in aligns/alignsPA that incide on the corresponding vertex
    vector<float> flipScores(specSet.size()); // Accumulated don't_flip/flip score balance per vertex
    SpecSet specSetRev;
    specSetRev.resize(specSet.size()); // Holds reversed versions of the spectra in the component
    list<int> toProcess; // Indices of spectra left to process
    TwoValues<vector<vector<TwoValues<int> > > > tmpMatches; // Temporary matches when a neighbor spectrum S1 in the i-th pair is matched directly/reversed (its pair S2 is assumed to already be in the component)
    TwoValues<vector<vector<TwoValues<int> > > > tmpMatchesPA; // Temporary matches when a neighbor spectrum S1 in the i-th PA pair is matched directly/reversed (its pair S2 is assumed to already be in the component)
    vector<TwoValues<float> > tmpModPos(aligns.size()); // Temporary storage for the position of the modification
    vector<bool> matchComputed(aligns.size()); // Indicates whether the i-th entry in tmpAligns was already computed
    vector<bool> matchPAComputed(alignsPA.size()); // Indicates whether the i-th entry in tmpAlignsPA was already computed
    matches.resize(aligns.size());
    modPos.resize(aligns.size());
    matchesPA.resize(alignsPA.size());
    //	vector<int> idx1, idx2; idx1.reserve(500); idx2.reserve(500);

    //cerr<<" * input aligns: ";
    //for(unsigned int i=0; i<aligns.size(); i++) cerr<<"["<<aligns[i].spec1<<","<<aligns[i].spec2<<"]"; cerr<<endl;

    // Initializations
    tmpMatches[0].resize(aligns.size());
    tmpMatches[1].resize(aligns.size());
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      if (specSet[aligns[i].spec1].size() != specSetRev[aligns[i].spec1].size())
      {
        specSet[aligns[i].spec1].reverse(0, &specSetRev[aligns[i].spec1]);
        toProcess.push_back(aligns[i].spec1);
      }
      if (specSet[aligns[i].spec2].size() != specSetRev[aligns[i].spec2].size())
      {
        specSet[aligns[i].spec2].reverse(0, &specSetRev[aligns[i].spec2]);
        toProcess.push_back(aligns[i].spec2);
      }
      matchComputed[i] = false;
      alignsEntries[aligns[i].spec1][0].push_back(i);
      alignsEntries[aligns[i].spec2][0].push_back(i);
    }

    tmpMatchesPA[0].resize(alignsPA.size());
    tmpMatchesPA[1].resize(alignsPA.size());
    for (unsigned int i = 0; i < alignsPA.size(); i++)
    {
      if (specSet[alignsPA[i].spec1].size()
          != specSetRev[alignsPA[i].spec1].size())
      {
        specSet[alignsPA[i].spec1].reverse(0, &specSetRev[alignsPA[i].spec1]);
        toProcess.push_back(alignsPA[i].spec1);
      }
      if (specSet[alignsPA[i].spec2].size()
          != specSetRev[alignsPA[i].spec2].size())
      {
        specSet[alignsPA[i].spec2].reverse(0, &specSetRev[alignsPA[i].spec2]);
        toProcess.push_back(alignsPA[i].spec2);
      }
      matchPAComputed[i] = false;
      alignsEntries[alignsPA[i].spec1][1].push_back(i);
      alignsEntries[alignsPA[i].spec2][1].push_back(i);
    }
    list<int>::iterator pIter;
    for (pIter = toProcess.begin(); pIter != toProcess.end(); pIter++)
      flipScores[*pIter] = 0;

    // Use edge scores as percentages of spectrum score to select the edge/vertex to initialize the direction selection
    vector<float> specScores(specSet.size());
    float bestScore = 0, curScore;
    int bestScoreIdx = -1, specToAdd = -1; // specToAdd is index of the new spectrum being added to the component
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specScores[i] = 0;
      for (unsigned int j = 0; j < specSet[i].size(); j++)
        specScores[i] += specSet[i][j][1];
    }
    for (unsigned int i = 0; i < aligns.size(); i++)
    {
      curScore = min(aligns[i].score1 / specScores[aligns[i].spec1],
                     aligns[i].score2 / specScores[aligns[i].spec2]);
      if (bestScoreIdx < 0 or bestScore < curScore)
      {
        bestScoreIdx = i;
        bestScore = curScore;
        specToAdd = aligns[bestScoreIdx].spec1;
      }
    }
    for (unsigned int i = 0; i < alignsPA.size(); i++)
    {
      curScore = min(alignsPA[i].score1 / specScores[alignsPA[i].spec1],
                     alignsPA[i].score2 / specScores[alignsPA[i].spec2]);
      if (bestScoreIdx < 0 or bestScore < curScore)
      {
        bestScoreIdx = i;
        bestScore = curScore;
        specToAdd = alignsPA[bestScoreIdx].spec1;
      }
    }

    //cerr<<" * initialized orientations component with vertex "<<specToAdd<<endl;
    flipScores[specToAdd] = -1; // Arbitrarily decide that this spectrum is not flipped
    pIter = toProcess.begin();
    while (*pIter != specToAdd)
      pIter++;
    toProcess.erase(pIter);
    while (specToAdd >= 0)
    {
      int flipDir = 0; // Spectrum not flipped (default)
      if (flipScores[specToAdd] > 0)
      {
        specSet[specToAdd] = specSetRev[specToAdd];
        flipDir = 1;
        specFlipped[specToAdd] = true;
        if (labelsP and labelsP->size()>specToAdd)
          (*labelsP)[specToAdd].reverse();
      }
      //cerr<<" * spectrum "<<specToAdd; if(flipDir) cerr<<" reversed"; else cerr<<" not reversed"; cerr<<", score = "<<flipScores[specToAdd]<<endl;

      for (pIter = alignsEntries[specToAdd][0].begin(); pIter
          != alignsEntries[specToAdd][0].end(); pIter++)
      {
        if (!matchComputed[*pIter])
        {
          // otherSpec is not in the component yet. Compute how the addition of specToAdd changes otherSpec's orientation selection.
          int otherSpec, otherSpecPos; // otherSpec=index of current paired spectrum, otherSpecPos=position of otherSpec index in aligns
          if (aligns[*pIter].spec1 == specToAdd)
          {
            otherSpec = aligns[*pIter].spec2;
            otherSpecPos = 1;
          }
          else
          {
            otherSpec = aligns[*pIter].spec1;
            otherSpecPos = 0;
          }

          TwoValues<float> matchScore1(0, 0), matchScore2(0, 0);
          vector<int> indices;
          vector<Spectrum> results(4);

          tmpModPos[*pIter][0] = SpectrumAlignment(&specSet[specToAdd],
                                             &specSet[otherSpec],
                                             peakTol,
                                             &results[0],
                                             &results[1],
                                             maxAAjump,
                                             penalty_sameVert,
                                             penalty_ptm,
                                             forceSymmetry,
                                             true);
          for (unsigned int i = 0; i < results[0].size(); i++)
          {
            matchScore1[0] += results[0][i][1];
            matchScore2[0] += results[1][i][1];
          }
          specSet[specToAdd].massesToIndices(results[0],
                                             indices,
                                             peakTol);
          tmpMatches[0][*pIter].resize(indices.size());
          for (unsigned int i = 0; i < indices.size(); i++)
            tmpMatches[0][*pIter][i][1 - otherSpecPos] = indices[i];
          specSet[otherSpec].massesToIndices(results[1],
                                             indices,
                                             peakTol);
          for (unsigned int i = 0; i < indices.size(); i++)
            tmpMatches[0][*pIter][i][otherSpecPos] = indices[i];

          //if(specToAdd==111 and otherSpec==246)
          //	{ cerr<<"Direct match:\n";
          //	for(unsigned int i=0;i<results[0].size();i++) cerr<<"  ("<<results[0][i][0]<<","<<results[0][i][1]<<")/("<<results[1][i][0]<<","<<results[1][i][1]<<")\n"; }

          tmpModPos[*pIter][1] = SpectrumAlignment(&specSet[specToAdd],
                                             &specSetRev[otherSpec],
                                             peakTol,
                                             &results[2],
                                             &results[3],
                                             maxAAjump,
                                             penalty_sameVert,
                                             penalty_ptm,
                                             forceSymmetry,
                                             true);
          for (unsigned int i = 0; i < results[2].size(); i++)
          {
            matchScore1[1] += results[2][i][1];
            matchScore2[1] += results[3][i][1];
          }
          specSet[specToAdd].massesToIndices(results[2],
                                             indices,
                                             peakTol);
          tmpMatches[1][*pIter].resize(indices.size());
          for (unsigned int i = 0; i < indices.size(); i++)
            tmpMatches[1][*pIter][i][1 - otherSpecPos] = indices[i];
          specSetRev[otherSpec].massesToIndices(results[3],
                                                indices,
                                                peakTol);
          for (unsigned int i = 0; i < indices.size(); i++)
            tmpMatches[1][*pIter][i][otherSpecPos] = indices[i];

          //if(specToAdd==111 and otherSpec==246)
          //	{ cerr<<"Reverse match:\n";
          //	for(unsigned int i=0;i<results[2].size();i++) cerr<<"  ("<<results[2][i][0]<<","<<results[2][i][1]<<")/("<<results[3][i][0]<<","<<results[3][i][1]<<")\n"; }

          matchScore1[0] /= specScores[specToAdd];
          matchScore1[1] /= specScores[specToAdd];
          matchScore2[0] /= specScores[otherSpec];
          matchScore2[1] /= specScores[otherSpec];
          //				flipScores[otherSpec] += matchScore2[1]-matchScore2[0];
          flipScores[otherSpec] += matchScore2[1] - matchScore2[0]
              + matchScore1[1] - matchScore1[0];
          //				if(alignRatios) (*alignRatios)[*pIter][0]=min(matchScore1[0], matchScore2[0]);
          //				if(alignRatios) (*alignRatios)[*pIter][1]=min(matchScore1[1], matchScore2[1]);
          if (alignRatios)
          {
            if (min(matchScore1[0], matchScore2[0]) >= min(matchScore1[1],
                                                           matchScore2[1]))
            {
              (*alignRatios)[*pIter][1 - otherSpecPos] = matchScore1[0];
              (*alignRatios)[*pIter][otherSpecPos] = matchScore2[0];
            }
            else
            {
              (*alignRatios)[*pIter][1 - otherSpecPos] = matchScore1[1];
              (*alignRatios)[*pIter][otherSpecPos] = matchScore2[1];
            }
          }
          matchComputed[*pIter] = true;
          //cerr<<" * match "<<aligns[*pIter].spec1<<" vs "<<aligns[*pIter].spec2<<", matchScore1=("<<matchScore1[0]<<","<<matchScore1[1]<<")/"<<specScores[specToAdd]<<", matchScore2=("<<matchScore2[0]<<","<<matchScore2[1]<<")/"<<specScores[otherSpec]<<", flipScores["<<otherSpec<<"]="<<flipScores[otherSpec]<<endl;
        }
        else
        {
          // Simply copy adequate entry (acording to flipDir) from tmpMatch to matches
          //cerr<<" * copying ASP match "<<aligns[*pIter].spec1<<" vs "<<aligns[*pIter].spec2<<endl;
          matches[*pIter].resize(tmpMatches[flipDir][*pIter].size());
          for (unsigned int i = 0; i < matches[*pIter].size(); i++)
            matches[*pIter][i] = tmpMatches[flipDir][*pIter][i];
          modPos[*pIter] = tmpModPos[*pIter][flipDir];
          if (alignRatios && flipDir == 1)
          {
            float tmp = (*alignRatios)[*pIter][0];
            (*alignRatios)[*pIter][0] = (*alignRatios)[*pIter][1];
            (*alignRatios)[*pIter][1] = tmp;
          }
        }
      }

      // Repeat cycle for PA overlaps
      for (pIter = alignsEntries[specToAdd][1].begin(); pIter
          != alignsEntries[specToAdd][1].end(); pIter++)
      {
        if (!matchPAComputed[*pIter])
        {
          // otherSpec is not in the component yet. Compute how the addition of specToAdd changes otherSpec's orientation selection.
          int otherSpec, otherSpecPos; // otherSpec=index of current paired spectrum, otherSpecPos=position of otherSpec index in aligns
          if (alignsPA[*pIter].spec1 == specToAdd)
          {
            otherSpec = alignsPA[*pIter].spec2;
            otherSpecPos = 1;
          }
          else
          {
            otherSpec = alignsPA[*pIter].spec1;
            otherSpecPos = 0;
          }

          TwoValues<float> score1, score2, matchScore1, matchScore2;
          vector<int> matchA1, matchA2, matchB1, matchB2;

          // Select shift for direct orientation
          //	            FindMatchPeaksAll(specSet[specToAdd], specSet[otherSpec], alignsPA[*pIter].shift1, peakTol, idx1, idx2);
          ScoreOverlap6(specSet[specToAdd],
                        specSet[otherSpec],
                        alignsPA[*pIter].shift1,
                        peakTol,
                        matchA1,
                        matchA2);
          score1.set(0, 0);
          for (unsigned int i = 0; i < matchA1.size(); i++)
          {
            score1[0] += specSet[specToAdd][matchA1[i]][1];
            score1[1] += specSet[otherSpec][matchA2[i]][1];
          }
          //	            FindMatchPeaksAll(specSet[specToAdd], specSet[otherSpec], alignsPA[*pIter].shift2, peakTol, idx1, idx2);
          ScoreOverlap6(specSet[specToAdd],
                        specSet[otherSpec],
                        alignsPA[*pIter].shift2,
                        peakTol,
                        matchB1,
                        matchB2);
          score2.set(0, 0);
          for (unsigned int i = 0; i < matchB1.size(); i++)
          {
            score2[0] += specSet[specToAdd][matchB1[i]][1];
            score2[1] += specSet[otherSpec][matchB2[i]][1];
          }
          if ((score1[0] / specScores[specToAdd] + score1[1]
              / specScores[otherSpec]) > (score2[0] / specScores[specToAdd]
              + score2[1] / specScores[otherSpec]))
          {
            matchScore1[0] = score1[0];
            matchScore2[0] = score1[1];
            tmpMatchesPA[0][*pIter].resize(matchA1.size());
            for (unsigned int i = 0; i < matchA1.size(); i++)
            {
              tmpMatchesPA[0][*pIter][i][1 - otherSpecPos] = matchA1[i];
              tmpMatchesPA[0][*pIter][i][otherSpecPos] = matchA2[i];
            }
          }
          else
          {
            matchScore1[0] = score2[0];
            matchScore2[0] = score2[1];
            tmpMatchesPA[0][*pIter].resize(matchB1.size());
            for (unsigned int i = 0; i < matchB1.size(); i++)
            {
              tmpMatchesPA[0][*pIter][i][1 - otherSpecPos] = matchB1[i];
              tmpMatchesPA[0][*pIter][i][otherSpecPos] = matchB2[i];
            }
          }

          // Select shift for reversed orientation
          //	            FindMatchPeaksAll(specSet[specToAdd], specSetRev[otherSpec], alignsPA[*pIter].shift1, peakTol, idx1, idx2);
          ScoreOverlap6(specSet[specToAdd],
                        specSetRev[otherSpec],
                        alignsPA[*pIter].shift1,
                        peakTol,
                        matchA1,
                        matchA2);
          score1.set(0, 0);
          for (unsigned int i = 0; i < matchA1.size(); i++)
          {
            score1[0] += specSet[specToAdd][matchA1[i]][1];
            score1[1] += specSetRev[otherSpec][matchA2[i]][1];
          }
          //	            FindMatchPeaksAll(specSet[specToAdd], specSetRev[otherSpec], alignsPA[*pIter].shift2, peakTol, idx1, idx2);
          ScoreOverlap6(specSet[specToAdd],
                        specSetRev[otherSpec],
                        alignsPA[*pIter].shift2,
                        peakTol,
                        matchB1,
                        matchB2);
          score2.set(0, 0);
          for (unsigned int i = 0; i < matchB1.size(); i++)
          {
            score2[0] += specSet[specToAdd][matchB1[i]][1];
            score2[1] += specSetRev[otherSpec][matchB2[i]][1];
          }
          if ((score1[0] / specScores[specToAdd] + score1[1]
              / specScores[otherSpec]) > (score2[0] / specScores[specToAdd]
              + score2[1] / specScores[otherSpec]))
          {
            matchScore1[1] = score1[0];
            matchScore2[1] = score1[1];
            tmpMatchesPA[1][*pIter].resize(matchA1.size());
            for (unsigned int i = 0; i < matchA1.size(); i++)
            {
              tmpMatchesPA[1][*pIter][i][1 - otherSpecPos] = matchA1[i];
              tmpMatchesPA[1][*pIter][i][otherSpecPos] = matchA2[i];
            }
          }
          else
          {
            matchScore1[1] = score2[0];
            matchScore2[1] = score2[1];
            tmpMatchesPA[1][*pIter].resize(matchB1.size());
            for (unsigned int i = 0; i < matchB1.size(); i++)
            {
              tmpMatchesPA[1][*pIter][i][1 - otherSpecPos] = matchB1[i];
              tmpMatchesPA[1][*pIter][i][otherSpecPos] = matchB2[i];
            }
          }

          matchScore1[0] /= specScores[specToAdd];
          matchScore1[1] /= specScores[specToAdd];
          matchScore2[0] /= specScores[otherSpec];
          matchScore2[1] /= specScores[otherSpec];
          flipScores[otherSpec] += matchScore2[1] - matchScore2[0]
              + matchScore1[1] - matchScore1[0];
          //				flipScores[otherSpec] += (matchScore2[1]-matchScore2[0])/specScores[otherSpec] + (matchScore1[1]-matchScore1[0])/specScores[specToAdd];
          //              if(alignRatiosPA) (*alignRatiosPA)[*pIter][0]=min(matchScore1[0]/specScores[specToAdd], matchScore2[0]/specScores[otherSpec]);
          //              if(alignRatiosPA) (*alignRatiosPA)[*pIter][1]=min(matchScore1[1]/specScores[specToAdd], matchScore2[1]/specScores[otherSpec]);
          if (alignRatiosPA)
          {
            if (min(matchScore1[0], matchScore2[0]) >= min(matchScore1[1],
                                                           matchScore2[1]))
              matchScore1[1] = matchScore2[0];
            else
            {
              matchScore1[0] = matchScore1[1];
              matchScore1[1] = matchScore2[1];
            }
            //					if(otherSpecPos==1) (*alignRatiosPA)[*pIter]=matchScore1; else (*alignRatiosPA)[*pIter].set(matchScore1[1],matchScore1[0]);
            (*alignRatiosPA)[*pIter].set(matchScore1[1 - otherSpecPos],
                                         matchScore1[otherSpecPos]);
          }
          matchPAComputed[*pIter] = true;
        }
        else
        {
          // Simply copy adequate entry (acording to flipDir) from tmpMatch to matches
          //cerr<<" * copying PA match "<<aligns[*pIter].spec1<<" vs "<<aligns[*pIter].spec2<<endl;
          matchesPA[*pIter].resize(tmpMatchesPA[flipDir][*pIter].size());
          for (unsigned int i = 0; i < matchesPA[*pIter].size(); i++)
            matchesPA[*pIter][i] = tmpMatchesPA[flipDir][*pIter][i];
          //				if(alignRatiosPA and flipDir==1) { float tmp=(*alignRatiosPA)[*pIter][0]; (*alignRatiosPA)[*pIter][0]=(*alignRatiosPA)[*pIter][1]; (*alignRatiosPA)[*pIter][1]=tmp; }
        }
      }

      if (toProcess.size() == 0)
        specToAdd = -1;
      else
      {
        list<int>::iterator bestNewSpec = toProcess.begin();
        pIter = bestNewSpec;
        pIter++;
        while (pIter != toProcess.end())
        {
          if (fabs(flipScores[*pIter]) > fabs(flipScores[*bestNewSpec]))
            bestNewSpec = pIter;
          pIter++;
        }
        specToAdd = *bestNewSpec;
        toProcess.erase(bestNewSpec);
      }
    }

    if (alignRatios)
    { // Cannot be computed correctly until all flipping decisions are made
      alignRatios->resize(aligns.size());
      for (unsigned int pairIdx = 0; pairIdx < aligns.size(); pairIdx++)
      {
        float score1 = 0, score2 = 0;
        int s1 = aligns[pairIdx].spec1, s2 = aligns[pairIdx].spec2;
        for (unsigned int i = 0; i < matches[pairIdx].size(); i++)
        {
          score1 += specSet[s1][matches[pairIdx][i][0]][1];
          score2 += specSet[s2][matches[pairIdx][i][1]][1];
        }
        (*alignRatios)[pairIdx].set(score1 / specScores[s1], score2
            / specScores[s2]);
      }
    }

    if (alignRatiosPA)
    { // Cannot be computed correctly until all flipping decisions are made
      alignRatiosPA->resize(alignsPA.size());
      for (unsigned int pairIdx = 0; pairIdx < alignsPA.size(); pairIdx++)
      {
        float score1 = 0, score2 = 0;
        int s1 = alignsPA[pairIdx].spec1, s2 = alignsPA[pairIdx].spec2;
        for (unsigned int i = 0; i < matchesPA[pairIdx].size(); i++)
        {
          score1 += specSet[s1][matchesPA[pairIdx][i][0]][1];
          score2 += specSet[s2][matchesPA[pairIdx][i][1]][1];
        }
        (*alignRatiosPA)[pairIdx].set(score1 / specScores[s1], score2
            / specScores[s2]);
      }
    }
  }

  // -------------------------------------------------------------------------
  //
  //  ComputeSpectralStars - decides on consensus orientations for a set of
  //        spectral-pair spectra and uses it to compute the corresponding spectral star.
  //    Order of instantiation of orientations is:
  //      1) initialize orientations with the highest scoring spectrum (summed scores)
  //      2) add next spectrum with highest scoring match to consensus
  //      3) update consensus wit h additional oriented spectrum
  //      4) iterate until no more spectra left
  //
  // -------------------------------------------------------------------------
  void ComputeSpectralStars(SpecSet &specSet,
                            Spectrum &consensus,
                            float peakTol,
                            float resolution)
  {
    unsigned int specIdx, peakIdx, maxScoreIdx = 0, maxScoreDir, numSpecs =
        specSet.size();
    float curScore, curScoreRev, maxScore = -1000000.0;
    vector<char> done(specSet.size());
    SpecSet specSetRev(specSet.size());
    Spectrum tmpConsensus;

    // Set initial consensus to highest scoring spectrum
    for (specIdx = 0; specIdx < numSpecs; specIdx++) {
      done[specIdx] = 0;
      curScore = 0;
      for (peakIdx = 0; peakIdx < specSet[specIdx].size(); peakIdx++) {
        curScore += specSet[specIdx][peakIdx][1];
      }
      if (curScore > maxScore) {
        maxScore = curScore;
        maxScoreIdx = specIdx;
      }
      specSet[specIdx].setPeakTolerance(peakTol); // Set the per-peak tolerances
      specSet[specIdx].reverse(0, &specSetRev[specIdx]);
    }
    consensus = specSet[maxScoreIdx];
    done[maxScoreIdx] = 1;

    // Select spectrum with highest-scoring match to consensus and use it to update consensus
    vector<int> idx1, idx2; // Indices of matching peaks between spectra and consensus spectrum
    do
    {
      maxScoreIdx = numSpecs;
      for (specIdx = 0; specIdx < numSpecs; specIdx++) {
        if (not done[specIdx])
        {
          // Compute score of direct match
          FindMatchPeaks(consensus, specSet[specIdx], 0, peakTol, idx1, idx2);
          curScore = 0;
          for (peakIdx = 0; peakIdx < idx1.size(); peakIdx++)
            curScore += consensus[idx1[peakIdx]][1];
          for (peakIdx = 0; peakIdx < idx2.size(); peakIdx++)
            curScore += specSet[specIdx][idx2[peakIdx]][1];

          // Compute score of reversed match
          FindMatchPeaks(consensus, specSetRev[specIdx], 0, peakTol, idx1, idx2);
          curScoreRev = 0;
          for (peakIdx = 0; peakIdx < idx1.size(); peakIdx++)
            curScoreRev += consensus[idx1[peakIdx]][1];
          for (peakIdx = 0; peakIdx < idx2.size(); peakIdx++)
            curScoreRev += specSet[specIdx][idx2[peakIdx]][1];

          if (max(curScoreRev, curScore) > maxScore)
          {
            if (curScoreRev > curScore)
            {
              maxScore = curScoreRev;
              maxScoreDir = 1;
            }
            else
            {
              maxScore = curScoreRev;
              maxScoreDir = 0;
            }
            maxScoreIdx = specIdx;
          }
        }
      }

      // Update consensus
      if (maxScoreIdx < numSpecs)
      {
        if (maxScoreDir == 0) {
          consensus.mergeClosestPeaks(specSet[maxScoreIdx], 1);
        } else {
    	   consensus.mergeClosestPeaks(specSetRev[maxScoreIdx], 1);
        }
        done[maxScoreIdx] = 1;
      } // if (maxScoreIdx < numSpecs)

    } while (maxScoreIdx < numSpecs);
    consensus.copyNP(specSet[0]);
  }

  // -------------------------------------------------------------------------
  /*
   * projectSpectrum - specFrom onto specTo:
   *   projection of specFrom onto specTo: Let M=mass(specTo)-mass(specFrom).
   *   The projection of specFrom onto specTo increases the scores of
   *   peaks in specTo by the scores of peaks in specFrom whose masses
   *   match a modified version of specFrom (with mod mass M).
   *
   *@param specFrom   Spectrum to project from
   *@param specTo     Spectrum to project onto
   *@param peakTol    Tolerance for mass errors (in Daltons)
   *@param bestScore  Best projection score in specTo over all possible deltas
   *@param bestDeltas The locations of the delta resulting in the highest score (bestScore)
   *@param finalProj  Projected version of specTo after projecting from specFrom with bestDelta
   *@param idxMatched Indices of matched peaks in specFrom/specTo
   *@param minNumMatchedPeaks Minimum number of matched peaks to accept a projection
   */
  // -------------------------------------------------------------------------
  void ProjectSpectrum(Spectrum &specFrom,
                       Spectrum &specTo, // was specBase,
                       float peakTol,
                       float &bestScore,
                       float &bestDelta,
                       Spectrum *finalProj,
                       vector<TwoValues<int> > *idxMatched,
                       unsigned int minNumMatchedPeaks)
  {
    if (specFrom.size()==0 or specTo.size()==0)
      return; // Invalid function call - no projections with empty spectra
    Spectrum specLocal;
    specLocal = specFrom; // version of specFrom with increased scores for matched peaks (no delta)
    Spectrum specLocalDelta;
    specLocalDelta = specFrom; // version of specFrom with increased scores for matched peaks (with delta)
    unsigned int numPeaks = specLocal.size();
    float curDelta = specTo.parentMass - specFrom.parentMass; // Mass difference between From and To spectra

    vector<int> matchesLocal, matchesTo; // Sets of matched indices returned by ScoreOverlap6
    vector<int> matchedIdx(numPeaks), // Each specFrom entry contains the index of the corresponding matched peak in specTo (-1 if not matched)
        matchedIdxDelta(numPeaks); // Each specFrom entry contains the index of the corresponding matched peak in specTo (specTo.size() if not matched)
    vector<int> prevMatchedIdx(numPeaks), nextMatchedIdxDelta(numPeaks); // Each entry contains the index of the prev/next matched peak in specTo (-1 if not matched)
    for (unsigned int p = 0; p < numPeaks; p++)
    {
      matchedIdx[p] = -1;
      matchedIdxDelta[p] = specTo.size();
    }

    // Find specTo matched peaks with no delta (tracked in specLocal)
    ScoreOverlap6(specLocal,
    		      specTo,
                  0.0,
                  peakTol,
                  matchesLocal,
                  matchesTo);
//cerr<<" --> got "<<matchesLocal.size()<<" matches for shift 0: ";
    for (unsigned int p = 0; p < matchesLocal.size(); p++)
    {
      specLocal[matchesLocal[p]][1] += specTo[matchesTo[p]][1]; // Add the scores of the matched peaks
      matchedIdx[matchesLocal[p]] = matchesTo[p]; // Keep track of which peaks were matched in specTo/specFrom
//cerr<<"("<<matchesLocal[p]<<","<<matchesTo[p]<<") ";
    }
//cerr<<endl;

    // Find specTo matched peaks with delta (tracked in specLocalDelta)
    ScoreOverlap6(specLocalDelta,
    		      specTo,
                  -curDelta,
                  peakTol,
                  matchesLocal,
                  matchesTo);
//cerr<<" --> got "<<matchesLocal.size()<<" matches for shift "<<-curDelta<<": ";
    for (unsigned int p = 0; p < matchesLocal.size(); p++)
    {
      specLocalDelta[matchesLocal[p]][1] += specTo[matchesTo[p]][1]; // Add the scores of the matched peaks
      matchedIdxDelta[matchesLocal[p]] = matchesTo[p]; // Keep track of which peaks were matched in specTo/specFrom
//cerr<<"("<<matchesLocal[p]<<","<<matchesTo[p]<<") ";
    }
//cerr<<endl;

    // Cumulative registers of the previous/next matched peaks. Used to avoid matching the same peak twice when curDelta is negative (see nextPeak code below)
    prevMatchedIdx[0] = matchedIdx[0];
    nextMatchedIdxDelta[numPeaks - 1] = matchedIdxDelta[numPeaks - 1];
//cerr<<0<<"\t"<<specLocal[0][0]<<"\t"<<specLocal[0][1]<<"\t"<<matchedIdx[0]<<"\t"<<prevMatchedIdx[0]<<" -- \t -- "
//		<<specLocalDelta[numPeaks - 1][0]<<"\t"<<specLocalDelta[numPeaks - 1][1]<<"\t"<<matchedIdxDelta[numPeaks - 1]<<"\t"<<nextMatchedIdxDelta[numPeaks - 1]<<endl;
    for (unsigned int p = 1; p < numPeaks; p++)
    {
      prevMatchedIdx[p] = max(matchedIdx[p], prevMatchedIdx[p - 1]);
      nextMatchedIdxDelta[numPeaks - p - 1] = min(nextMatchedIdxDelta[numPeaks - p], matchedIdxDelta[numPeaks - p - 1]);

//cerr<<p<<"\t"<<specLocal[p][0]<<"\t"<<specLocal[p][1]<<"\t"<<matchedIdx[p]<<"\t"<<prevMatchedIdx[p]<<" -- \t -- "
//		<<specLocalDelta[numPeaks - p - 1][0]<<"\t"<<specLocalDelta[numPeaks - p - 1][1]<<"\t"<<matchedIdxDelta[numPeaks - p - 1]<<"\t"<<nextMatchedIdxDelta[numPeaks - p - 1]<<endl;
    }

    // compute match scores
      vector<float> cumFromLeft(numPeaks), cumFromRight(numPeaks); // Cumulative sums
      cumFromLeft[0] = specLocal[0][1];
      cumFromRight[numPeaks - 1] = specLocalDelta[numPeaks - 1][1];
      for (unsigned int p = 1; p < numPeaks; p++)
      {
    	  if(specLocal[p][0] > specTo.parentMass-AAJumps::massMH-peakTol)  // Peak is not overlapped with this alignment
    		  cumFromLeft[p] = -1000000;
    	  else
    	  {
    		  cumFromLeft[p] = specLocal[p][1];
    		  if(cumFromLeft[p - 1]>0) cumFromLeft[p] += cumFromLeft[p - 1];
    	  }

    	  if(specLocalDelta[numPeaks - p - 1][0] + curDelta < -peakTol)  // Peak is not overlapped with this alignment
    		  cumFromRight[numPeaks - p - 1] = -1000000;
    	  else
    	  {
    		  cumFromRight[numPeaks - p - 1] = specLocalDelta[numPeaks - p - 1][1];
    		  if(cumFromRight[numPeaks - p]>0) cumFromRight[numPeaks - p - 1] += cumFromRight[numPeaks - p];
    	  }

//cerr<<"cumFromLeft["<<p<<"] = "<<cumFromLeft[p]<<"\t cumFromRight["<<numPeaks-p-1<<"] = "<<cumFromRight[numPeaks-p-1]<<endl;
      }

      // Find local best placement of the modification
      unsigned int bestLocalDelta = numPeaks + 1,
    		       nextPeak,
    		       bestNextPeak; // Index of the first peak after the modification at the best location
      float bestLocalScore = -1, curScore;
      for (int deltaIdx = (int)numPeaks; deltaIdx >= 0; deltaIdx--)  // Iterator of the first peak after the modification
      {
        if (deltaIdx == (int)numPeaks) // placed on the C-term
        {
          curScore = cumFromLeft[numPeaks - 1];
          nextPeak = numPeaks;
        }
        else
        {
            nextPeak = deltaIdx;
            if (deltaIdx > 0) {
            	curScore = cumFromLeft[deltaIdx - 1];
            	while (nextPeak < numPeaks and
            			(prevMatchedIdx[deltaIdx - 1] >= nextMatchedIdxDelta[nextPeak]
            			 or specLocal[deltaIdx - 1][0] + AAJumps::minAAmass - 2*peakTol >= specLocalDelta[nextPeak][0] + curDelta))
            		nextPeak++; // Avoid matching the same peak twice on negative modifications
            } else {
                curScore = 0;
                while (nextPeak < numPeaks and
                		AAJumps::minAAmass - 2*peakTol >= specLocalDelta[nextPeak][0] + curDelta)
                  nextPeak++; // Avoid projecting to masses smaller than the minimum amino acid mass
            }

            if (nextPeak < numPeaks)
            	curScore += cumFromRight[nextPeak];
        }
        if (curScore > bestLocalScore)
        {
          bestLocalScore = curScore;
          bestLocalDelta = deltaIdx;
          bestNextPeak = nextPeak;
        }
      }

//cerr<<"bestLocalScore = "<<bestLocalScore<<", bestLocalDelta = "<<bestLocalDelta<<", bestNextPeak = "<<bestNextPeak<<endl;

      // Check if current solution outscores global best and replace the latter if it does.
      unsigned int matchedCount = 0;
      for (unsigned int i = 0; i < bestLocalDelta; i++)
        if (matchedIdx[i] != -1)
          matchedCount++;
      for (unsigned int i = bestNextPeak; i < numPeaks; i++)
        if (matchedIdxDelta[i] < (int)specFrom.size())
          matchedCount++;

      if (bestLocalScore > bestScore && matchedCount >= minNumMatchedPeaks)
      {
        bestScore = bestLocalScore;
        if (bestNextPeak < numPeaks)
          bestDelta = bestNextPeak;
        else
          bestDelta = -1; // -1 marks modification at the end

        if (finalProj)
        {
          (*finalProj) = specTo;
          finalProj->parentMass = specTo.parentMass;
          finalProj->parentCharge = specTo.parentCharge;
          finalProj->scan = specTo.scan;
          finalProj->resize(numPeaks);  // max size
          unsigned int curPeak = 0;
          for (curPeak = 0; curPeak < bestLocalDelta; curPeak++)
          {
            (*finalProj)[curPeak] = specLocal[curPeak];
          }
          for (unsigned int p = (unsigned int)bestNextPeak; p < numPeaks; p++)
          {
            (*finalProj)[curPeak] = specLocalDelta[p];
            (*finalProj)[curPeak++][0] += curDelta;
          }
          finalProj->resize(curPeak);  // used size
        }

        if (idxMatched)
        {
          idxMatched->resize(matchedCount);
          matchedCount = 0;
          for (unsigned int i = 0; i < bestLocalDelta; i++)
            if (matchedIdx[i] != -1)
              (*idxMatched)[matchedCount++].set(i, matchedIdx[i]);

          for (unsigned int i = bestNextPeak; i < numPeaks; i++)
            if (matchedIdxDelta[i] < (int)specFrom.size())
              (*idxMatched)[matchedCount++].set(i, matchedIdxDelta[i]);
        }
      }

  }

  // -------------------------------------------------------------------------
  /*
   * ProjectSpectrumOld - "projects" all specSet spectra with indices in specsToProcess onto specBase:
   *   projection of A onto B: Let M=mass(B)-mass(A). The projection of A onto B increases the scores of
   *   peaks in B by the scores of peaks in A whose masses match a modified version of A (with mod mass M).
   *
   *@param specSet    Set of all spectra
   *@param specBase   Spectrum to project onto
   *@param specsToProcess Indices of spectra (in specSet) to project from
   *@param curDeltas  Locations of the deltas (modifications) for the current subset of processed spectra (from specsToProcess, used in recursion)
   *@param bestScore  Best projection score in specBase over all possible deltas on all neighbors in specsToProcess
   *@param bestDeltas The locations of the deltas resulting in the highest score (bestScore)
   *@param peakTol    Tolerance for mass errors (in Daltons)
   *@param finalProj  Projected version of specBase after projecting all specsToProcess with bestDeltas
   *@param idxMatched Indices of matched p eaks in specBase/specsToProcess.front() (only when specsToProcess.size()==1)
   */
  // -------------------------------------------------------------------------
  void ProjectSpectrumOld(SpecSet &specSet,
                       const Spectrum &specBase,
                       list<int> &specsToProcess,
                       list<int> &curDeltas,
                       float &bestScore,
                       list<int> &bestDeltas,
                       float peakTol,
                       Spectrum *finalProj,
                       vector<TwoValues<int> > *idxMatched,
                       unsigned int minNumMatchedPeaks)
  {
    if (specsToProcess.empty() or specBase.size() == 0)
      return; // Invalid function call - specsToProcess should always have at least one element
    Spectrum specLocal;
    specLocal = specBase; // version of specBase with increased scores for matched peaks (no delta)
    Spectrum specLocalDelta;
    specLocalDelta = specBase; // version of specBase with increased scores for matched peaks (with delta)
    unsigned int projSpecIdx = specsToProcess.front();
    specsToProcess.pop_front(); // index of the spectrum to project from
    unsigned int numPeaks = specLocal.size();
    float curDelta = specSet[projSpecIdx].parentMass - specBase.parentMass; // Mass difference between From and To spectra

    vector<int> matchesLocal, matchesNew; // Sets of matched indices returned by FindMatchPeaksAll
    vector<int> matchedIdx(numPeaks), // Each specBase entry contains the index of the corresponding matched peak in projSpecIdx (-1 if not matched)
        matchedIdxDelta(numPeaks); // Each specBase entry contains the index of the corresponding matched peak in projSpecIdx (specSet[projSpecIdx].size() if not matched)
    vector<int> prevMatchedIdx(numPeaks), nextMatchedIdxDelta(numPeaks); // Each entry contains the index of the prev/next matched peak in projSpecIdx (-1 if not matched)
    for (unsigned int p = 0; p < numPeaks; p++)
    {
      matchedIdx[p] = -1;
      matchedIdxDelta[p] = specSet[projSpecIdx].size();
    }

    // Find specBase matched peaks with no delta (tracked in specLocal)
    //	FindMatchPeaksAll(specLocal, specSet[projSpecIdx], 0, peakTol, matchesLocal, matchesNew);
    ScoreOverlap6(specLocal,
                  specSet[projSpecIdx],
                  0.0,
                  peakTol,
                  matchesLocal,
                  matchesNew);
cerr<<" --> got "<<matchesLocal.size()<<" matches for shift 0: ";
    for (unsigned int p = 0; p < matchesLocal.size(); p++)
    {
      specLocal[matchesLocal[p]][1] += specSet[projSpecIdx][matchesNew[p]][1]; // Add the scores of the matched peaks
      matchedIdx[matchesLocal[p]] = matchesNew[p]; // Keep track of which peaks were matched in specBase/projSpecIdx
cerr<<"("<<matchesLocal[p]<<","<<matchesNew[p]<<") ";
    }
cerr<<endl;

    // Find specBase matched peaks with delta (tracked in specLocalDelta)
    //	FindMatchPeaksAll(specLocalDelta, specSet[projSpecIdx], -curDelta, peakTol, matchesLocal, matchesNew);
    ScoreOverlap6(specLocalDelta,
                  specSet[projSpecIdx],
                  -curDelta,
                  peakTol,
                  matchesLocal,
                  matchesNew);
cerr<<" --> got "<<matchesLocal.size()<<" matches for shift "<<-curDelta<<": ";
    for (unsigned int p = 0; p < matchesLocal.size(); p++)
    {
      specLocalDelta[matchesLocal[p]][1]
          += specSet[projSpecIdx][matchesNew[p]][1]; // Add the scores of the matched peaks
      matchedIdxDelta[matchesLocal[p]] = matchesNew[p]; // Keep track of which peaks were matched in specBase/projSpecIdx
cerr<<"("<<matchesLocal[p]<<","<<matchesNew[p]<<") ";
    }
cerr<<endl;

    // Cumulative registers of the previous/next matched peaks. Used to avoid matching the same peak twice when the mod is negative (see nextPeak code below)
    prevMatchedIdx[0] = matchedIdx[0];
    nextMatchedIdxDelta[numPeaks - 1] = matchedIdxDelta[numPeaks - 1];
    for (unsigned int p = 1; p < numPeaks; p++)
    {
      prevMatchedIdx[p] = max(matchedIdx[p], prevMatchedIdx[p - 1]);
      nextMatchedIdxDelta[numPeaks - p - 1] = min(nextMatchedIdxDelta[numPeaks - p], matchedIdxDelta[numPeaks - p - 1]);

cerr<<specLocal[p][0]<<"\t"<<specLocal[p][1]<<"\t"<<matchedIdx[p]<<"\t"<<prevMatchedIdx[p]<<" -- \t -- "
		<<specLocal[numPeaks - p - 1][0]<<"\t"<<specLocal[numPeaks - p - 1][1]<<"\t"<<matchedIdxDelta[numPeaks - p - 1]<<"\t"<<nextMatchedIdxDelta[numPeaks - p - 1]<<endl;
    }

    if (specsToProcess.size() == 0)
    { // This is the last spectrum being projected into - compute match scores
      vector<float> cumFromLeft(numPeaks), cumFromRight(numPeaks); // Cumulative sums
      cumFromLeft[0] = specLocal[0][1];
      cumFromRight[numPeaks - 1] = specLocalDelta[numPeaks - 1][1];
      for (unsigned int p = 1; p < numPeaks; p++)
      {
        cumFromLeft[p] = cumFromLeft[p - 1];
        if(specLocal[p][0] + curDelta) cumFromLeft[p] += specLocal[p][1];
        cumFromRight[numPeaks - p - 1] = cumFromRight[numPeaks - p]
            + specLocalDelta[numPeaks - p - 1][1];
cerr<<"cumFromLeft["<<p<<"] = "<<cumFromLeft[p]<<"\t cumFromRight["<<numPeaks-p-1<<"] = "<<cumFromRight[numPeaks-p-1]<<endl;
      }

      // Find local best placement of the modification
      unsigned int bestLocalDelta = numPeaks + 1,
    		       nextPeak,
    		       bestNextPeak; // Index of the first peak after the modification
      float bestLocalScore = -1, curScore;
      for (int deltaIdx = (int)numPeaks; deltaIdx >= 0; deltaIdx--)
      {
        if (deltaIdx == (int)numPeaks) // placed on the C-term
        {
          curScore = cumFromLeft[numPeaks - 1];
          nextPeak = numPeaks;
        }
        else
        {
          if (deltaIdx == 0) {
            curScore = cumFromRight[0];
            nextPeak = 0;
          }
          else
          {
            nextPeak = deltaIdx;
cerr<<"deltaIdx = "<<deltaIdx<<endl;
            while (nextPeak < numPeaks and
            		(prevMatchedIdx[deltaIdx - 1] >= nextMatchedIdxDelta[nextPeak]
                     or specLocal[deltaIdx - 1][0] + 56 >= specLocalDelta[nextPeak][0] + curDelta))
              nextPeak++; // Avoid matching the same peak twice on negative modifications
cerr<<"nextPeak = "<<nextPeak<<endl;
            if (nextPeak < numPeaks)
              curScore = cumFromLeft[deltaIdx - 1] + cumFromRight[nextPeak];
            else
              curScore = cumFromLeft[deltaIdx - 1];
          }
        }
        if (curScore > bestLocalScore)
        {
          bestLocalScore = curScore;
          bestLocalDelta = deltaIdx;
          bestNextPeak = nextPeak;
        }
      }

cerr<<"bestLocalScore = "<<bestLocalScore<<", bestLocalDelta = "<<bestLocalDelta<<", bestNextPeak = "<<bestNextPeak<<endl;

      // Check if current solution outscores global best and replace the latter if it does.
      unsigned int matchedCount = 0;
      for (unsigned int i = 0; i < bestLocalDelta; i++)
        if (matchedIdx[i] != -1)
          matchedCount++;
      for (unsigned int i = bestNextPeak; i < numPeaks; i++)
        if (matchedIdxDelta[i] < (int)specSet[projSpecIdx].size())
          matchedCount++;

      if (bestLocalScore > bestScore && matchedCount >= minNumMatchedPeaks)
      {
        bestScore = bestLocalScore;
        bestDeltas.clear();
        bestDeltas.assign(curDeltas.begin(), curDeltas.end());
        if (bestNextPeak < numPeaks)
          bestDeltas.push_back(bestNextPeak);
        else
          bestDeltas.push_back(-1); // -1 marks modification at the end

        if (finalProj)
        {
          (*finalProj) = specLocal;
          finalProj->parentMass = specSet[projSpecIdx].parentMass;
          finalProj->parentCharge = specSet[projSpecIdx].parentCharge;
          finalProj->scan = specSet[projSpecIdx].scan;
          for (unsigned int p = (unsigned int)bestNextPeak; p < numPeaks; p++)
          {
            (*finalProj)[p][1] = specLocalDelta[p][1];
            (*finalProj)[p][0] += curDelta;
          }
        }

        if (idxMatched)
        {
          idxMatched->resize(matchedCount);
          matchedCount = 0;
          for (unsigned int i = 0; i < bestLocalDelta; i++)
            if (matchedIdx[i] != -1)
              (*idxMatched)[matchedCount++].set(i, matchedIdx[i]);

          for (unsigned int i = bestNextPeak; i < numPeaks; i++)
            if (matchedIdxDelta[i] < (int)specSet[projSpecIdx].size())
              (*idxMatched)[matchedCount++].set(i, matchedIdxDelta[i]);
        }
      }

    }
    else
    { // Recurse to project onto remaining spectra
      specLocal.parentMass = specSet[projSpecIdx].parentMass; // Because specBase is being projected to specSet[projSpecIdx]
      specLocal.parentCharge = specSet[projSpecIdx].parentCharge;
      specLocal.scan = specSet[projSpecIdx].scan;
      for (int deltaIdx = (int)numPeaks; deltaIdx >= 0; deltaIdx--)
      { // Modification is placed right before the deltaIdx-th peak
        if (deltaIdx >= 0 && deltaIdx < (int)numPeaks)
        {
          specLocal[deltaIdx][1] = specLocalDelta[deltaIdx][1];
          specLocal[deltaIdx][0] += curDelta;
        }
        cerr
            << "(ERROR) Unfixed bug: negative mods of amino acid masses may cause the same peak to be matched twice in the spectrum with smaller parent mass!\n";
        curDeltas.push_back(deltaIdx);
        ProjectSpectrumOld(specSet,
                        specLocal,
                        specsToProcess,
                        curDeltas,
                        bestScore,
                        bestDeltas,
                        peakTol,
                        finalProj);
        curDeltas.pop_back();
      }
    }

    specsToProcess.push_front(projSpecIdx);
  }

} // namespace specnets

