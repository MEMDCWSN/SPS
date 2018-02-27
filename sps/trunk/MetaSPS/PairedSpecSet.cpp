/*
 * PairedSpecSet.cpp
 *
 *  Created on: Nov 11, 2011
 *      Author: aguthals
 */

#include "PairedSpecSet.h"
#include "math.h"
#include "Logger.h"

using namespace std;

namespace specnets
{

  PairedSpecSet::PairedSpecSet() :
    mergedSpectrum(0x0), mergedLabels(0x0), inputSpectra(0x0),
        reversedPeaks(0x0), usedPeaks(0x0), m_srmOffset(0),
        m_boostedSpectra(0x0)
  {

  }

  PairedSpecSet::PairedSpecSet(SpecSet* allSpectra, bool enforceDaTol, vector<
      vector<bool> >* reversedPeaks, float peakTol, float srmOffset) :
    mergedSpectrum(0x0), mergedLabels(0x0), inputSpectra(0x0),
        reversedPeaks(0x0), usedPeaks(0x0), enforceDaTolerance(enforceDaTol),
        m_peakTol(peakTol), m_srmOffset(srmOffset), m_boostedSpectra(0x0)
  {

    initialize(allSpectra, enforceDaTol, reversedPeaks, peakTol, srmOffset);
  }

  PairedSpecSet::~PairedSpecSet()
  {
    if (mergedSpectrum != 0)
    {
      delete mergedSpectrum;
    }

    if (mergedLabels != 0)
    {
      delete mergedLabels;
    }

    if (inputSpectra != 0)
    {
      delete inputSpectra;
    }

    if (reversedPeaks != 0)
    {
      delete reversedPeaks;
    }

    if (usedPeaks != 0)
    {
      delete usedPeaks;
    }

    if (m_boostedSpectra != 0)
    {
      delete m_boostedSpectra;
    }
  }

  void PairedSpecSet::initialize(SpecSet* allSpectra,
                                 bool enforceDaTol,
                                 vector<vector<bool> >* revPeaks,
                                 float peakTol,
                                 float srmOffset)
  {

    if (inputSpectra == 0)
    {
      inputSpectra = new SpecSet(allSpectra->size());
    }
    (*inputSpectra) = *allSpectra;

    if (revPeaks)
    {
      if (reversedPeaks == 0)
      {
        reversedPeaks = new vector<vector<bool> > (revPeaks->size());
      }
      (*reversedPeaks) = *revPeaks;
    }

    if (!mergedSpectrum)
    {
      mergedSpectrum = new Spectrum;
    }
    else
    {
      mergedSpectrum->resize(0);
    }

    if (!m_boostedSpectra)
    {
      m_boostedSpectra = new SpecSet(inputSpectra->size());
    }
    else
    {
      m_boostedSpectra->resize(inputSpectra->size());
    }
    (*m_boostedSpectra) = *inputSpectra;

    if (!mergedLabels)
    {
      mergedLabels = new vector<short> (0);
    }
    else
    {
      mergedLabels->resize(0);
    }

    if (!usedPeaks)
    {
      usedPeaks = new vector<vector<short> > (inputSpectra->size());
    }
    initializeUsedPeaks(usedPeaks, inputSpectra);

    enforceDaTolerance = enforceDaTol;
    m_peakTol = peakTol;
    m_srmOffset = srmOffset;

    if (enforceDaTolerance)
    {
      inputSpectra->setPeakTolerance(m_peakTol, false);
    }
    /*
     double totalMass = 0;
     double totalMassTol = 0;
     float totalCharges = 0;
     double numMasses = 0;
     */

    for (unsigned int i = 0; i < inputSpectra->size(); i++)
    {
      if ((*inputSpectra)[i].size() == 0)
      {
        continue;
      }
      float pmTolUse = (enforceDaTol) ? 0 : (*inputSpectra)[i].parentMassTol;

      MZRange massB0(0, 0, 0.0001);
      MZRange massBk((*inputSpectra)[i].parentMass - AAJumps::massMH,
                     0,
                     pmTolUse);
      MZRange massY0(AAJumps::massH2O, 0, 0.0001);
      MZRange massYk((*inputSpectra)[i].parentMass - AAJumps::massHion,
                     0,
                     pmTolUse);
      MZRange massZk((*inputSpectra)[i].parentMass - AAJumps::massMH
          - AAJumps::massNH, 0, pmTolUse);
      int idxFound;

      Spectrum* spec = &(*inputSpectra)[i];

      idxFound = spec->findPeaks(massB0);
      if (idxFound >= 0)
      {
        (*usedPeaks)[i][idxFound] = 3;
      }
      idxFound = spec->findPeaks(massBk);
      if (idxFound >= 0)
      {
        (*usedPeaks)[i][idxFound] = 3;
      }
      idxFound = spec->findPeaks(massY0);
      if (idxFound >= 0)
      {
        (*usedPeaks)[i][idxFound] = 3;
      }
      idxFound = spec->findPeaks(massYk);
      if (idxFound >= 0)
      {
        (*usedPeaks)[i][idxFound] = 3;
      }
      if (spec->msFragType == Spectrum::FragType_ETD)
      {
        idxFound = spec->findPeaks(massZk);
        if (idxFound >= 0)
        {
          (*usedPeaks)[i][idxFound] = 3;
        }
      }
      /*
       totalMass += spec->parentMass;
       totalMassTol += spec->parentMassTol;
       totalCharges += (float) spec->parentCharge;
       numMasses += 1.0;
       */
    }

    float endPtScore = 30.0;

    int insertIdx;
    mergedSpectrum->copyNP((*inputSpectra)[0]);
    //mergedSpectrum->parentMass = totalMass / numMasses;
    //mergedSpectrum->parentMassTol = totalMassTol / numMasses;
    //mergedSpectrum->setCharge((short) round(totalCharges / numMasses));
    float pmTolUse = (enforceDaTol) ? m_peakTol : mergedSpectrum->parentMassTol;
    mergedSpectrum->resize(0);
    float massB0 = 0;
    float massBk = mergedSpectrum->parentMass - AAJumps::massMH;
    float massY0 = AAJumps::massH2O;
    float massYk = mergedSpectrum->parentMass - AAJumps::massHion;
    insertIdx = mergedSpectrum->insertPeak(massB0, endPtScore, 0);
    insertIdx = mergedSpectrum->insertPeak(massBk, endPtScore, pmTolUse);
    insertIdx = mergedSpectrum->insertPeak(massY0, endPtScore, 0);
    insertIdx = mergedSpectrum->insertPeak(massYk, endPtScore, pmTolUse);
    (*mergedLabels).resize(0);
    (*mergedLabels).resize(4, (short)3);
  }

  void PairedSpecSet::mergePRMs(unsigned short finalStage)
  {
    if (finalStage >= 1)
    {
      mergePRMsStage1();
    }
    if (finalStage >= 2)
    {
      mergePRMsStage2();
    }
    if (finalStage >= 3)
    {
      mergePRMsStage3();
    }
    if (finalStage >= 4)
    {
      mergePRMsStage4();
    }
    if (finalStage >= 5)
    {
      mergePRMsStage5();
    }
  }

  void PairedSpecSet::boostPRMs(SpecSet& boostedSpecs)
  {
    list<int> ETDIdxs;
    list<int> CIDIdxs;

    for (int j = 0; j < inputSpectra->size(); j++)
    {
      if ((*inputSpectra)[j].size() == 0)
      {
        continue;
      }
      if ((*inputSpectra)[j].msFragType == Spectrum::FragType_ETD)
      {
        ETDIdxs.push_back(j);
      }
      else
      {
        CIDIdxs.push_back(j);
      }
    }

    for (int j = 0; j < inputSpectra->size(); j++)
    {
      if ((*inputSpectra)[j].size() == 0)
      {
        continue;
      }
      if ((*inputSpectra)[j].msFragType == Spectrum::FragType_ETD)
      {
        boostPRMsETD(j, CIDIdxs);
      }
      else
      {
        boostPRMsCID(j, ETDIdxs);
      }
    }

    boostedSpecs = *m_boostedSpectra;
  }

  void PairedSpecSet::mergePRMsStage1()
  {

    list<pair<int, int> > cidEtdPairs;
    getUniqueCIDETDPairs(cidEtdPairs);

    if (cidEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no CID/ETD pairs, skipping Stage 1 CID/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = cidEtdPairs.begin(); pairIt
          != cidEtdPairs.end(); pairIt++)
      {
        mergePRMsPair(pairIt->first, pairIt->second);
      }
    }

    list<pair<int, int> > hcdEtdPairs;
    getUniqueHCDETDPairs(hcdEtdPairs);

    if (hcdEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no HCD/ETD pairs, skipping Stage 1 HCD/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = hcdEtdPairs.begin(); pairIt
          != hcdEtdPairs.end(); pairIt++)
      {
        mergePRMsPair(pairIt->first, pairIt->second);
      }
    }
  }

  void PairedSpecSet::mergePRMsStage2()
  {
    list<pair<int, int> > cidEtdPairs;
    getUniqueCIDETDPairs(cidEtdPairs);

    if (cidEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no CID/ETD pairs, skipping Stage 2 CID/ETD merge");
    }
    else
    {
      list<int> disjointIdxs;
      for (list<pair<int, int> >::iterator pairIt = cidEtdPairs.begin(); pairIt
          != cidEtdPairs.end(); pairIt++)
      {
        getComplementIdxs(pairIt->first, pairIt->second, disjointIdxs);
        mergeSRMsPair(pairIt->first, pairIt->second, &disjointIdxs);
      }
    }

    list<pair<int, int> > hcdEtdPairs;
    getUniqueHCDETDPairs(hcdEtdPairs);

    if (hcdEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no HCD/ETD pairs, skipping Stage 2 HCD/ETD merge");
    }
    else
    {
      list<int> disjointIdxs;
      for (list<pair<int, int> >::iterator pairIt = hcdEtdPairs.begin(); pairIt
          != hcdEtdPairs.end(); pairIt++)
      {
        getComplementIdxs(pairIt->first, pairIt->second, disjointIdxs);
        mergeSRMsPair(pairIt->first, pairIt->second, &disjointIdxs);
      }
    }
  }

  void PairedSpecSet::mergePRMsStage3()
  {
    list<pair<int, int> > cidEtdPairs;
    getUniqueCIDETDPairs(cidEtdPairs);

    if (cidEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no CID/ETD pairs, skipping Stage 3 CID/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = cidEtdPairs.begin(); pairIt
          != cidEtdPairs.end(); pairIt++)
      {
        mergePRMsSRMsPair(pairIt->first, pairIt->second);
      }
    }

    list<pair<int, int> > hcdEtdPairs;
    getUniqueHCDETDPairs(hcdEtdPairs);

    if (hcdEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no HCD/ETD pairs, skipping Stage 3 HCD/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = hcdEtdPairs.begin(); pairIt
          != hcdEtdPairs.end(); pairIt++)
      {
        mergePRMsSRMsPair(pairIt->first, pairIt->second);
      }
    }
  }

  void PairedSpecSet::mergePRMsStage4()
  {

    list<pair<int, int> > cidEtdPairs;
    getUniqueCIDETDPairs(cidEtdPairs);

    if (cidEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no CID/ETD pairs, skipping Stage 4 CID/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = cidEtdPairs.begin(); pairIt
          != cidEtdPairs.end(); pairIt++)
      {
        mergeSRMsPair(pairIt->first, pairIt->second);
      }
    }

    list<pair<int, int> > hcdEtdPairs;
    getUniqueHCDETDPairs(hcdEtdPairs);

    if (hcdEtdPairs.size() == 0)
    {
      //DEBUG_MSG("Found no HCD/ETD pairs, skipping Stage 4 HCD/ETD merge");
    }
    else
    {
      for (list<pair<int, int> >::iterator pairIt = hcdEtdPairs.begin(); pairIt
          != hcdEtdPairs.end(); pairIt++)
      {
        mergeSRMsPair(pairIt->first, pairIt->second);
      }
    }
  }

  void PairedSpecSet::mergePRMsStage5()
  {
    for (int i = 0; i < inputSpectra->size(); i++)
    {
      mergeLeftovers(i);
    }
  }

  void PairedSpecSet::getMergedSpectrum(Spectrum* toSpec)
  {
    (*toSpec) = *mergedSpectrum;
  }

  void PairedSpecSet::initializeUsedPeaks(vector<vector<short> >* usedPeaks,
                                          SpecSet* spectra)
  {
    usedPeaks->resize(spectra->size());
    for (int i = 0; i < spectra->size(); i++)
    {
      (*usedPeaks)[i].resize((*spectra)[i].size());
      for (int j = 0; j < (*usedPeaks)[i].size(); j++)
      {
        (*usedPeaks)[i][j] = 0;
      }
    }
  }

  void PairedSpecSet::insertMergedLabel(vector<short>* peakLabels,
                                        int index,
                                        short label)
  {
    vector<short>::iterator labelIt = peakLabels->begin();
    labelIt += index;
    peakLabels->insert(labelIt, label);
  }

  void PairedSpecSet::boostPRMsCID(int CIDIdx, list<int>& ETDIdxs)
  {
    //DEBUG_VAR(CIDIdx);
    Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
    Spectrum* boostedSpec = &(*m_boostedSpectra)[CIDIdx];
    boostedSpec->operator =(*CIDspec);

    if (CIDspec->size() == 0 || boostedSpec->size() == 0)
    {
      return;
    }

    //DEBUG_VAR(CIDspec->size());

    vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];

    vector<bool> addedPeaksCID(CIDspec->size(), false);
    vector<set<int> > addedPeaksETD(inputSpectra->size());

    for (unsigned int i = 0; i < CIDlabels->size(); i++)
    {
      addedPeaksCID[i] = ((*CIDlabels)[i] == 3);
    }

    if (!enforceDaTolerance)
    {
      ERROR_MSG("Not yet supported!!");
      abort();
    }

    set<int> peaksToRemove;

    int numPeaks = CIDspec->size();

    //float pmTol = (enforceDaTolerance) ? 0 : CIDspec->parentMassTol;

    float massBk = CIDspec->parentMass - AAJumps::massMH;

    // Boost PRM/PRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3 || addedPeaksCID[j])
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = (*CIDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massCID) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      for (list<int>::iterator idxIt = ETDIdxs.begin(); idxIt != ETDIdxs.end(); idxIt++)
      {
        int ETDIdx = *idxIt;
        Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
        vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
        int closestETD = ETDspec->findClosest(massCID);
        float massETD = (*ETDspec)[closestETD][0];

        float intenETD = (*ETDspec)[closestETD][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massCID, massETD, m_peakTol)
            || (*ETDlabels)[closestETD] == 2 || (*ETDlabels)[closestETD] == 3
            || addedPeaksETD[ETDIdx].count(closestETD) > 0)
        {
          continue;
        }

        (*CIDlabels)[j] = 1;
        (*ETDlabels)[closestETD] = 1;
        (*boostedSpec)[j][1] += intenETD;
        addedPeaksETD[ETDIdx].insert(closestETD);
        addedPeaksCID[j] = true;

        float cPeakCID = massBk - massCID + AAJumps::massH2O;
        float cPeakETD = cPeakCID + m_srmOffset;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && !addedPeaksCID[cPeakIdxCID]
            && (*CIDlabels)[cPeakIdxCID] != 1)
        {
          (*CIDlabels)[cPeakIdxCID] = 2;
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
          peaksToRemove.insert(cPeakIdxCID);
          addedPeaksCID[cPeakIdxCID] = true;
        }

        if (cPeakIdxETD >= 0 && addedPeaksETD[ETDIdx].count(cPeakIdxETD) == 0
            && (*ETDlabels)[cPeakIdxETD] != 1)
        {
          (*ETDlabels)[cPeakIdxETD] = 2;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          addedPeaksETD[ETDIdx].insert(cPeakIdxETD);
        }
      }
    }

    // Boost SRM/SRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*CIDlabels)[j] == 1 || (*CIDlabels)[j] == 3 || addedPeaksCID[j])
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = (*CIDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massCID) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massETDSRM = massCID + m_srmOffset;

      for (list<int>::iterator idxIt = ETDIdxs.begin(); idxIt != ETDIdxs.end(); idxIt++)
      {
        int ETDIdx = *idxIt;
        Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
        vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
        int closestETD = ETDspec->findClosest(massETDSRM);
        float massETD = (*ETDspec)[closestETD][0];

        float intenETD = (*ETDspec)[closestETD][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massETDSRM, massETD, m_peakTol)
            || (*ETDlabels)[closestETD] == 1 || (*ETDlabels)[closestETD] == 3
            || addedPeaksETD[ETDIdx].count(closestETD) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massCID);
         DEBUG_VAR(massETDSRM);
         DEBUG_VAR(massETD);
         */

        (*CIDlabels)[j] = 2;
        (*ETDlabels)[closestETD] = 2;
        (*boostedSpec)[j][1] += intenETD;
        addedPeaksETD[ETDIdx].insert(closestETD);
        addedPeaksCID[j] = true;

        float cPeakCID = massBk - massCID + AAJumps::massH2O;
        float cPeakETD = cPeakCID;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] != 2
            && !addedPeaksCID[cPeakIdxCID])
        {
          (*CIDlabels)[cPeakIdxCID] = 1;
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
          addedPeaksCID[cPeakIdxCID] = true;
          peaksToRemove.insert(cPeakIdxCID);
        }
        else if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] != 2
            && addedPeaksETD[ETDIdx].count(cPeakIdxETD) == 0)
        {
          (*ETDlabels)[cPeakIdxETD] = 1;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          addedPeaksETD[ETDIdx].insert(cPeakIdxETD);
        }
      }
    }

    // Boost PRM/SRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3 || addedPeaksCID[j])
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = (*CIDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massCID) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massETDSRM = (massBk - massCID + AAJumps::massH2O) + m_srmOffset;

      for (list<int>::iterator idxIt = ETDIdxs.begin(); idxIt != ETDIdxs.end(); idxIt++)
      {
        int ETDIdx = *idxIt;
        Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
        vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
        int closestETD = ETDspec->findClosest(massETDSRM);
        float massETD = (*ETDspec)[closestETD][0];

        float intenETD = (*ETDspec)[closestETD][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massETDSRM, massETD, m_peakTol)
            || (*ETDlabels)[closestETD] == 1 || (*ETDlabels)[closestETD] == 3
            || addedPeaksETD[ETDIdx].count(closestETD) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massCID);
         DEBUG_VAR(massETDSRM);
         DEBUG_VAR(massETD);
         */

        (*CIDlabels)[j] = 1;
        (*ETDlabels)[closestETD] = 2;
        (*boostedSpec)[j][1] += intenETD;
        addedPeaksETD[ETDIdx].insert(closestETD);
        addedPeaksCID[j] = true;

        float cPeakCID = massBk - massCID + AAJumps::massH2O;
        float cPeakETD = massCID;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] != 1
            && !addedPeaksCID[cPeakIdxCID])
        {
          (*CIDlabels)[cPeakIdxCID] = 2;
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
          addedPeaksCID[cPeakIdxCID] = true;
          peaksToRemove.insert(cPeakIdxCID);
        }
        else if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] != 2
            && addedPeaksETD[ETDIdx].count(cPeakIdxETD) == 0)
        {
          (*ETDlabels)[cPeakIdxETD] = 1;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          addedPeaksETD[ETDIdx].insert(cPeakIdxETD);
        }
      }
    }

    // Boost SRM/PRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*CIDlabels)[j] == 1 || (*CIDlabels)[j] == 3 || addedPeaksCID[j])
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = (*CIDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massCID) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massETDPRM = massBk - massCID + AAJumps::massH2O;

      for (list<int>::iterator idxIt = ETDIdxs.begin(); idxIt != ETDIdxs.end(); idxIt++)
      {
        int ETDIdx = *idxIt;
        Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
        vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
        int closestETD = ETDspec->findClosest(massETDPRM);
        float massETD = (*ETDspec)[closestETD][0];

        float intenETD = (*ETDspec)[closestETD][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massETDPRM, massETD, m_peakTol)
            || (*ETDlabels)[closestETD] == 2 || (*ETDlabels)[closestETD] == 3
            || addedPeaksETD[ETDIdx].count(closestETD) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massCID);
         DEBUG_VAR(massETDPRM);
         DEBUG_VAR(massETD);
         */

        (*CIDlabels)[j] = 2;
        (*ETDlabels)[closestETD] = 1;
        (*boostedSpec)[j][1] += intenETD;
        addedPeaksETD[ETDIdx].insert(closestETD);
        addedPeaksCID[j] = true;

        float cPeakCID = massETDPRM;
        float cPeakETD = massCID + m_srmOffset;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && !addedPeaksCID[cPeakIdxCID]
            && (*CIDlabels)[cPeakIdxCID] != 2)
        {
          (*CIDlabels)[cPeakIdxCID] = 1;
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
          peaksToRemove.insert(cPeakIdxCID);
          addedPeaksCID[cPeakIdxCID] = true;
        }

        if (cPeakIdxETD >= 0 && addedPeaksETD[ETDIdx].count(cPeakIdxETD) == 0
            && (*ETDlabels)[cPeakIdxETD] != 2)
        {
          (*ETDlabels)[cPeakIdxETD] = 1;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          addedPeaksETD[ETDIdx].insert(cPeakIdxETD);
        }
      }
    }

    int idxShift = 0;
    for (int j = 0; j < numPeaks; j++)
    {
      int peakIdx = j - idxShift;
      if (peaksToRemove.count(j) > 0)
      {
        boostedSpec->removePeak(peakIdx);
        idxShift++;
        continue;
      }
      if ((*CIDlabels)[j] == 2)
      {
        (*boostedSpec)[peakIdx][0] = massBk - (*boostedSpec)[peakIdx][0]
            + AAJumps::massH2O;
      }
      if ((*CIDlabels)[j] != 0 && (*CIDlabels)[j] != 3)
      {
        (*boostedSpec)[peakIdx][1] *= 3;
      }
    }
    boostedSpec->sortPeaks();
  }

  void PairedSpecSet::boostPRMsETD(int ETDIdx, list<int>& CIDIdxs)
  {
    //DEBUG_VAR(ETDIdx);
    Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
    Spectrum* boostedSpec = &(*m_boostedSpectra)[ETDIdx];
    boostedSpec->operator =(*ETDspec);

    if (ETDspec->size() == 0 || boostedSpec->size() == 0)
    {
      return;
    }

    //DEBUG_VAR(ETDspec->size());

    //ETDspec->output(cerr);

    vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];

    vector<bool> addedPeaksETD(ETDspec->size(), false);
    vector<set<int> > addedPeaksCID(inputSpectra->size());

    for (unsigned int i = 0; i < ETDlabels->size(); i++)
    {
      addedPeaksETD[i] = ((*ETDlabels)[i] == 3);
    }

    if (!enforceDaTolerance)
    {
      ERROR_MSG("Not yet supported!!");
      abort();
    }

    set<int> peaksToRemove;

    int numPeaks = ETDspec->size();

    //float pmTol = (enforceDaTolerance) ? 0 : CIDspec->parentMassTol;

    float massBk = ETDspec->parentMass - AAJumps::massMH;

    // Boost PRM/PRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*ETDlabels)[j] == 2 || (*ETDlabels)[j] == 3 || addedPeaksETD[j])
      {
        continue;
      }

      float massETD = (*ETDspec)[j][0];
      float intenETD = (*ETDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massETD) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      for (list<int>::iterator idxIt = CIDIdxs.begin(); idxIt != CIDIdxs.end(); idxIt++)
      {
        int CIDIdx = *idxIt;
        Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
        vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
        int closestCID = CIDspec->findClosest(massETD);
        float massCID = (*CIDspec)[closestCID][0];

        float intenCID = (*CIDspec)[closestCID][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massCID, massETD, m_peakTol)
            || (*CIDlabels)[closestCID] == 2 || (*CIDlabels)[closestCID] == 3
            || addedPeaksCID[CIDIdx].count(closestCID) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massETD);
         DEBUG_VAR(massCID);
         */

        (*ETDlabels)[j] = 1;
        (*CIDlabels)[closestCID] = 1;
        (*boostedSpec)[j][1] += intenCID;
        addedPeaksCID[CIDIdx].insert(closestCID);
        addedPeaksETD[j] = true;

        float cPeakCID = massBk - massETD + AAJumps::massH2O;
        float cPeakETD = cPeakCID + m_srmOffset;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && addedPeaksCID[CIDIdx].count(cPeakIdxCID) == 0
            && (*CIDlabels)[cPeakIdxCID] != 1)
        {
          (*CIDlabels)[cPeakIdxCID] = 2;
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
          addedPeaksCID[CIDIdx].insert(cPeakIdxCID);
        }

        if (cPeakIdxETD >= 0 && !addedPeaksETD[cPeakIdxETD] == 0
            && (*ETDlabels)[cPeakIdxETD] != 1)
        {
          (*ETDlabels)[cPeakIdxETD] = 2;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          peaksToRemove.insert(cPeakIdxETD);
          addedPeaksETD[cPeakIdxETD] = true;
        }
      }
    }

    // Boost SRM/SRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*ETDlabels)[j] == 1 || (*ETDlabels)[j] == 3 || addedPeaksETD[j])
      {
        continue;
      }

      float massETD = (*ETDspec)[j][0];
      float intenETD = (*ETDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massETD) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massCIDSRM = massETD - m_srmOffset;

      for (list<int>::iterator idxIt = CIDIdxs.begin(); idxIt != CIDIdxs.end(); idxIt++)
      {
        int CIDIdx = *idxIt;
        Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
        vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
        int closestCID = CIDspec->findClosest(massCIDSRM);
        float massCID = (*CIDspec)[closestCID][0];

        float intenCID = (*CIDspec)[closestCID][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massCIDSRM, massCID, m_peakTol)
            || (*CIDlabels)[closestCID] == 1 || (*CIDlabels)[closestCID] == 3
            || addedPeaksCID[CIDIdx].count(closestCID) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massETD);
         DEBUG_VAR(massCIDSRM);
         DEBUG_VAR(massCID);
         */

        (*ETDlabels)[j] = 2;
        (*CIDlabels)[closestCID] = 2;
        (*boostedSpec)[j][1] += intenCID;
        addedPeaksCID[CIDIdx].insert(closestCID);
        addedPeaksETD[j] = true;

        float cPeakCID = massBk - massCID + AAJumps::massH2O;
        float cPeakETD = cPeakCID;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] != 2
            && addedPeaksCID[CIDIdx].count(cPeakIdxCID) == 0)
        {
          (*CIDlabels)[cPeakIdxCID] = 1;
          addedPeaksCID[CIDIdx].insert(cPeakIdxCID);
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
        }
        else if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] != 2
            && !addedPeaksETD[cPeakIdxETD])
        {
          (*ETDlabels)[cPeakIdxETD] = 1;
          addedPeaksETD[cPeakIdxETD] = true;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          peaksToRemove.insert(cPeakIdxETD);
        }
      }
    }

    // Boost PRM/SRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*ETDlabels)[j] == 2 || (*ETDlabels)[j] == 3 || addedPeaksETD[j])
      {
        continue;
      }

      float massETD = (*ETDspec)[j][0];
      float intenETD = (*ETDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massETD) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massCIDSRM = (massBk - massETD) + AAJumps::massH2O;

      for (list<int>::iterator idxIt = CIDIdxs.begin(); idxIt != CIDIdxs.end(); idxIt++)
      {
        int CIDIdx = *idxIt;
        Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
        vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
        int closestCID = CIDspec->findClosest(massCIDSRM);
        float massCID = (*CIDspec)[closestCID][0];

        float intenCID = (*CIDspec)[closestCID][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massCIDSRM, massCID, m_peakTol)
            || (*CIDlabels)[closestCID] == 1 || (*CIDlabels)[closestCID] == 3
            || addedPeaksCID[CIDIdx].count(closestCID) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massETD);
         DEBUG_VAR(massCIDSRM);
         DEBUG_VAR(massCID);
         */

        (*ETDlabels)[j] = 1;
        (*CIDlabels)[closestCID] = 2;
        (*boostedSpec)[j][1] += intenCID;
        addedPeaksCID[CIDIdx].insert(closestCID);
        addedPeaksETD[j] = true;

        float cPeakETD = massCIDSRM + m_srmOffset;
        float cPeakCID = massETD;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] != 2
            && addedPeaksCID[CIDIdx].count(cPeakIdxCID) == 0)
        {
          (*CIDlabels)[cPeakIdxCID] = 1;
          addedPeaksCID[CIDIdx].insert(cPeakIdxCID);
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
        }
        else if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] != 1
            && !addedPeaksETD[cPeakIdxETD])
        {
          (*ETDlabels)[cPeakIdxETD] = 2;
          addedPeaksETD[cPeakIdxETD] = true;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          peaksToRemove.insert(cPeakIdxETD);
        }
      }
    }

    // Boost SRM/PRM pairs
    for (int j = 0; j < numPeaks; j++)
    {
      if ((*ETDlabels)[j] == 1 || (*ETDlabels)[j] == 3 || addedPeaksETD[j])
      {
        continue;
      }

      float massETD = (*ETDspec)[j][0];
      float intenETD = (*ETDspec)[j][1];
      //float massTolCID = CIDspec->getTolerance(j);

      if (boostedSpec->findPeaks(massETD) < 0)
      {
        ERROR_MSG("Cannot find peak!");
        abort();
      }

      float massCIDPRM = massBk - (massETD - m_srmOffset) + AAJumps::massH2O;

      for (list<int>::iterator idxIt = CIDIdxs.begin(); idxIt != CIDIdxs.end(); idxIt++)
      {
        int CIDIdx = *idxIt;
        Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
        vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
        int closestCID = CIDspec->findClosest(massCIDPRM);
        float massCID = (*CIDspec)[closestCID][0];

        float intenCID = (*CIDspec)[closestCID][1];
        //float massTolETD = ETDspec->getTolerance(closestETD);

        if (!MZRange::EqualWithinRange(massCIDPRM, massCID, m_peakTol)
            || (*CIDlabels)[closestCID] == 2 || (*CIDlabels)[closestCID] == 3
            || addedPeaksCID[CIDIdx].count(closestCID) > 0)
        {
          continue;
        }
        /*
         DEBUG_VAR(massETD);
         DEBUG_VAR(massCIDPRM);
         DEBUG_VAR(massCID);
         */

        (*ETDlabels)[j] = 2;
        (*CIDlabels)[closestCID] = 1;
        (*boostedSpec)[j][1] += intenCID;
        addedPeaksCID[CIDIdx].insert(closestCID);
        addedPeaksETD[j] = true;

        float cPeakETD = massCIDPRM;
        float cPeakCID = massETD - m_srmOffset;

        int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
        int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

        if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] != 1
            && addedPeaksCID[CIDIdx].count(cPeakIdxCID) == 0)
        {
          (*CIDlabels)[cPeakIdxCID] = 2;
          addedPeaksCID[CIDIdx].insert(cPeakIdxCID);
          (*boostedSpec)[j][1] += (*CIDspec)[cPeakIdxCID][1];
        }
        else if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] != 2
            && !addedPeaksETD[cPeakIdxETD])
        {
          (*ETDlabels)[cPeakIdxETD] = 1;
          addedPeaksETD[cPeakIdxETD] = true;
          (*boostedSpec)[j][1] += (*ETDspec)[cPeakIdxETD][1];
          peaksToRemove.insert(cPeakIdxETD);
        }
      }
    }

    int idxShift = 0;
    for (int j = 0; j < numPeaks; j++)
    {
      int peakIdx = j - idxShift;
      if (peaksToRemove.count(j) > 0)
      {
        boostedSpec->removePeak(peakIdx);
        idxShift++;
        continue;
      }
      if ((*ETDlabels)[j] == 2)
      {
        (*boostedSpec)[peakIdx][0] = massBk - ((*boostedSpec)[peakIdx][0]
            - m_srmOffset - AAJumps::massH2O);
      }

      if ((*ETDlabels)[j] != 0 && (*ETDlabels)[j] != 3)
      {
        (*boostedSpec)[peakIdx][1] *= 3;
      }
    }

    boostedSpec->sortPeaks();
    /*
     boostedSpec->output(cerr);
     abort();
     */
  }

  void PairedSpecSet::mergePRMsPair(int CIDIdx, int ETDIdx)
  {

    Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
    Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
    Spectrum* mergedSpec = mergedSpectrum;

    if (CIDspec->size() == 0 || ETDspec->size() == 0)
    {
      return;
    }
    /*
     bool debug = false;
     if (mergedIdx == 2207) {
     debug = true;
     }
     */
    vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
    vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
    vector<short>* mergedLocLabels = mergedLabels;

    MZRange rangeCID;
    MZRange rangeETD;
    MZRange cPeakCID;
    MZRange cPeakETD;
    MZRange mergeRange;
    list<MZRange> peaksToAdd;

    int numPeaks = CIDspec->size();

    float pmTol = (enforceDaTolerance) ? 0 : mergedSpec->parentMassTol;

    float massBk = mergedSpec->parentMass - AAJumps::massMH;

    for (int j = 0; j < numPeaks; j++)
    {

      if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3)
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
      float massTolCID = CIDspec->getTolerance(j);
      rangeCID.set(massCID, intenCID, massTolCID);

      int closestETD = ETDspec->findClosest(massCID);
      float massETD = (*ETDspec)[closestETD][0];

      float intenETD = ((*ETDlabels)[closestETD] == 0)
          ? (*ETDspec)[closestETD][1] : 0;
      float massTolETD = ETDspec->getTolerance(closestETD);
      rangeETD.set(massETD, intenETD, massTolETD);

      if (rangeCID != rangeETD || (*ETDlabels)[closestETD] == 2
          || (*ETDlabels)[closestETD] == 3)
      {
        continue;
      }

      if ((*CIDlabels)[j] == 1 && (*ETDlabels)[closestETD] == 1)
      {
        continue;
      }
      else
      {
        (*CIDlabels)[j] = 1;
        (*ETDlabels)[closestETD] = 1;
      }

      mergeRange = (massTolCID < massTolETD) ? rangeCID : rangeETD;
      mergeRange.setIntensity(intenCID + intenETD);

      cPeakCID.set(massBk - mergeRange.getMass() + AAJumps::massH2O, 0, pmTol
          + mergeRange.getTolerance());
      cPeakETD.set(cPeakCID.getMass() + m_srmOffset, 0, pmTol
          + mergeRange.getTolerance());

      if (enforceDaTolerance)
      {
        cPeakCID.setTolerance(0);
        cPeakETD.setTolerance(0);
      }
      int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
      int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

      if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] == 0)
      {
        mergeRange.setIntensity(mergeRange.getIntensity()
            + (*CIDspec)[cPeakIdxCID][1]);
        (*CIDlabels)[cPeakIdxCID] = 2;
      }

      if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] == 0)
      {
        mergeRange.setIntensity(mergeRange.getIntensity()
            + (*ETDspec)[cPeakIdxETD][1]);
        (*ETDlabels)[cPeakIdxETD] = 2;
      }

      mergeRange.setIntensity(mergeRange.getIntensity() * 2);
      peaksToAdd.push_back(mergeRange);
    }
    addNewPeaks(&peaksToAdd, 1);
  }

  void PairedSpecSet::mergeSRMsPair(int CIDIdx, int ETDIdx, list<int>* checkPRM)
  {

    Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
    Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
    Spectrum* mergedSpec = mergedSpectrum;

    if (CIDspec->size() == 0 || ETDspec->size() == 0)
    {
      return;
    }

    vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
    vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
    vector<short>* mergedLocLabels = mergedLabels;

    int numPeaks = CIDspec->size();

    float pmTol = (enforceDaTolerance) ? 0 : mergedSpec->parentMassTol;

    float massBk = mergedSpec->parentMass - AAJumps::massMH;

    MZRange rangeETD;
    MZRange rangeCID;
    MZRange mergeRange;
    list<MZRange> peaksToAdd;

    for (int j = 0; j < numPeaks; j++)
    {

      if ((*CIDlabels)[j] == 1 || (*CIDlabels)[j] == 3)
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
      float massTolCID = CIDspec->getTolerance(j);
      rangeCID.set(massCID + m_srmOffset, intenCID, massTolCID);

      int closestETD = ETDspec->findClosest(rangeCID.getMass());

      float massETD = (*ETDspec)[closestETD][0];
      float intenETD = ((*ETDlabels)[closestETD] == 0)
          ? (*ETDspec)[closestETD][1] : 0;
      float massTolETD = ETDspec->getTolerance(closestETD);
      rangeETD.set(massETD, intenETD, massTolETD);

      if (rangeCID != rangeETD || (*ETDlabels)[closestETD] == 1
          || (*ETDlabels)[closestETD] == 3)
      {
        continue;
      }

      if ((*CIDlabels)[j] == 2 && (*ETDlabels)[closestETD] == 2)
      {
        continue;
      }

      mergeRange = (massTolCID < massTolETD) ? rangeCID : rangeETD;
      mergeRange.setMass(massBk - (mergeRange.getMass() - m_srmOffset
          - AAJumps::massH2O));
      mergeRange.setTolerance(pmTol + mergeRange.getTolerance());
      mergeRange.setIntensity(intenCID + intenETD);

      int prmIdx = -1;
      if (checkPRM)
      {
        for (list<int>::iterator idxIt = checkPRM->begin(); idxIt
            != checkPRM->end(); idxIt++)
        {
          int i = *idxIt;
          Spectrum* checkSpec = &(*inputSpectra)[i];
          vector<short>* checkLabels = &(*usedPeaks)[i];
          int prmIdxTemp = -1;
          if (enforceDaTolerance)
          {
            prmIdxTemp = checkSpec->findPeaks(mergeRange.getMass());
          }
          else
          {
            prmIdxTemp = checkSpec->findPeaks(mergeRange);
          }
          if (prmIdxTemp >= 0 && (*checkLabels)[prmIdxTemp] != 2
              && (*checkLabels)[prmIdxTemp] != 3)
          {
            prmIdx = prmIdxTemp;
            if ((*checkLabels)[prmIdx] == 0)
            {
              mergeRange.setIntensity(mergeRange.getIntensity()
                  + (*checkSpec)[prmIdx][1]);
            }
            (*checkLabels)[prmIdx] = 1;
            if (checkSpec->getTolerance(prmIdx) < mergeRange.getTolerance())
            {
              mergeRange.setMass((*checkSpec)[prmIdx][0]);
              mergeRange.setTolerance(checkSpec->getTolerance(prmIdx));
            }
          }
        }

        if (prmIdx < 0)
        {
          continue;
        }
      }

      (*CIDlabels)[j] = 2;
      (*ETDlabels)[closestETD] = 2;

      mergeRange.setIntensity(mergeRange.getIntensity() * 2);
      peaksToAdd.push_back(mergeRange);
    }
    addNewPeaks(&peaksToAdd, 1);
  }

  void PairedSpecSet::mergePRMsSRMsPair(int CIDIdx, int ETDIdx)
  {

    Spectrum* CIDspec = &(*inputSpectra)[CIDIdx];
    Spectrum* ETDspec = &(*inputSpectra)[ETDIdx];
    Spectrum* mergedSpec = mergedSpectrum;

    if (CIDspec->size() == 0 || ETDspec->size() == 0)
    {
      return;
    }

    vector<short>* CIDlabels = &(*usedPeaks)[CIDIdx];
    vector<short>* ETDlabels = &(*usedPeaks)[ETDIdx];
    vector<short>* mergedLocLabels = mergedLabels;

    MZRange rangeCID;
    MZRange rangeETD;
    MZRange cPeakETD;
    MZRange mergeRange;
    list<MZRange> peaksToAdd;

    int numPeaks = CIDspec->size();

    float pmTol = (enforceDaTolerance) ? 0 : mergedSpec->parentMassTol;

    float massBk = mergedSpec->parentMass - AAJumps::massMH;

    for (int j = 0; j < numPeaks; j++)
    {

      if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3)
      {
        continue;
      }

      float massCID = (*CIDspec)[j][0];
      float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
      float massTolCID = CIDspec->getTolerance(j);
      rangeCID.set(massCID, intenCID, massTolCID);
      cPeakETD.set(massBk - massCID + AAJumps::massH2O + m_srmOffset,
                   intenCID,
                   massTolCID + pmTol);

      int closestETD = ETDspec->findClosest(cPeakETD.getMass());

      float massETD = (*ETDspec)[closestETD][0];
      float intenETD = ((*ETDlabels)[closestETD] == 0)
          ? (*ETDspec)[closestETD][1] : 0;
      float massTolETD = ETDspec->getTolerance(closestETD);
      rangeETD.set(massETD, intenETD, massTolETD);

      if (rangeETD != cPeakETD || (*ETDlabels)[closestETD] == 1
          || (*ETDlabels)[closestETD] == 3)
      {
        continue;
      }

      if ((*CIDlabels)[j] == 1 && (*ETDlabels)[closestETD] == 2)
      {
        continue;
      }
      else
      {
        (*CIDlabels)[j] = 1;
        (*ETDlabels)[closestETD] = 2;
      }

      mergeRange = rangeCID;
      mergeRange.setIntensity(intenCID + intenETD);

      mergeRange.setIntensity(mergeRange.getIntensity() * 2);
      peaksToAdd.push_back(mergeRange);
    }

    MZRange cPeakCID;
    MZRange closestETD;
    numPeaks = ETDspec->size();

    for (int j = 0; j < numPeaks; j++)
    {

      if ((*ETDlabels)[j] == 2 || (*ETDlabels)[j] == 3)
      {
        continue;
      }

      float massETD = (*ETDspec)[j][0];
      float intenETD = ((*ETDlabels)[j] == 0) ? (*ETDspec)[j][1] : 0;
      float massTolETD = ETDspec->getTolerance(j);
      rangeETD.set(massETD, intenETD, massTolETD);
      cPeakCID.set(massBk - massETD + AAJumps::massH2O, intenETD, massTolETD
          + pmTol);

      int closestCID = CIDspec->findClosest(cPeakCID.getMass());

      float massCID = (*CIDspec)[closestCID][0];
      float intenCID = ((*CIDlabels)[closestCID] == 0)
          ? (*CIDspec)[closestCID][1] : 0;
      float massTolCID = CIDspec->getTolerance(closestCID);
      rangeCID.set(massCID, intenCID, massTolCID);

      if (rangeCID != cPeakCID || (*CIDlabels)[closestCID] == 1
          || (*ETDlabels)[j] == 3)
      {
        continue;
      }

      if ((*ETDlabels)[j] == 1 && (*CIDlabels)[closestCID] == 2)
      {
        continue;
      }
      else
      {
        (*CIDlabels)[closestCID] = 2;
        (*ETDlabels)[j] = 1;
      }

      mergeRange = rangeETD;
      mergeRange.setIntensity(intenCID + intenETD);

      mergeRange.setIntensity(mergeRange.getIntensity() * 2);
      peaksToAdd.push_back(mergeRange);
    }
    addNewPeaks(&peaksToAdd, 1);
  }

  void PairedSpecSet::mergeLeftovers(int childIdx)
  {

    Spectrum* childSpec = &(*inputSpectra)[childIdx];
    Spectrum* mergedSpec = mergedSpectrum;

    vector<short>* childLabels = &(*usedPeaks)[childIdx];

    vector<bool>* reversedLabels = (reversedPeaks)
        ? &(*reversedPeaks)[childIdx] : 0;
    vector<short>* mergedLocLabels = mergedLabels;

    MZRange mergeRange;
    list<MZRange> peaksToAdd;

    float pmTol = (enforceDaTolerance) ? 0 : mergedSpec->parentMassTol;

    float massBk = mergedSpec->parentMass - AAJumps::massMH;

    int numPeaks = childSpec->size();

    for (int j = 0; j < numPeaks; j++)
    {

      float mass = (*childSpec)[j][0];
      float intensity = (*childSpec)[j][1];
      float tolerance = childSpec->getTolerance(j);

      // These peaks were called PRMs by PepNovo but we reversed them to reduce mass error
      if (reversedLabels && (*reversedLabels)[j])
      {
        float srmOffset = (childSpec->msFragType == Spectrum::FragType_ETD)
            ? 0.0 - m_srmOffset - AAJumps::massH2O : 0.0 - AAJumps::massH2O;

        mass = massBk - (mass + srmOffset);
        tolerance += pmTol;
      }

      if ((*childLabels)[j] != 0)
      {
        continue;
      }

      mergeRange.set(mass, intensity, tolerance);

      peaksToAdd.push_back(mergeRange);
    }
    addNewPeaks(&peaksToAdd, 0);
  }

  void PairedSpecSet::addNewPeaks(list<MZRange>* newPeaks, short label)
  {
    Spectrum* mergedSpec = mergedSpectrum;
    vector<short>* mergedLocLabels = mergedLabels;

    for (list<MZRange>::iterator peakIt = newPeaks->begin(); peakIt
        != newPeaks->end(); peakIt++)
    {
      int mergeIdx = mergedSpec->findPeaks(peakIt->getMass());
      if (mergeIdx >= 0)
      {
        (*mergedSpec)[mergeIdx][1] += peakIt->getIntensity();
      }
      else
      {
        int insertIdx = mergedSpec->insertPeak(&(*peakIt));
        if (enforceDaTolerance)
        {
          mergedSpec->setTolerance(insertIdx, m_peakTol);
        }

        insertMergedLabel(mergedLocLabels, insertIdx, label);
      }
    }
  }

  void PairedSpecSet::getUniqueCIDETDPairs(list<pair<int, int> >& outputPairs)
  {
    outputPairs.clear();

    for (int i = 0; i < inputSpectra->size(); i++)
    {
      Spectrum* spec1 = &(*inputSpectra)[i];
      for (int j = i + 1; j < inputSpectra->size(); j++)
      {

        Spectrum* spec2 = &(*inputSpectra)[j];
        if (spec1->msFragType == Spectrum::FragType_CID && spec2->msFragType
            == Spectrum::FragType_ETD)
        {
          outputPairs.push_back(pair<int, int> (i, j));
        }
        else if (spec2->msFragType == Spectrum::FragType_CID
            && spec1->msFragType == Spectrum::FragType_ETD)
        {
          outputPairs.push_back(pair<int, int> (j, i));
        }

      }
    }
  }

  void PairedSpecSet::getUniqueHCDETDPairs(list<pair<int, int> >& outputPairs)
  {
    outputPairs.clear();

    for (int i = 0; i < inputSpectra->size(); i++)
    {
      Spectrum* spec1 = &(*inputSpectra)[i];
      for (int j = i + 1; j < inputSpectra->size(); j++)
      {
        Spectrum* spec2 = &(*inputSpectra)[j];
        if (spec1->msFragType == Spectrum::FragType_HCD && spec2->msFragType
            == Spectrum::FragType_ETD)
        {
          outputPairs.push_back(pair<int, int> (i, j));
        }
        else if (spec2->msFragType == Spectrum::FragType_HCD
            && spec1->msFragType == Spectrum::FragType_ETD)
        {
          outputPairs.push_back(pair<int, int> (j, i));
        }

      }
    }
  }

  void PairedSpecSet::getComplementIdxs(int idx1,
                                        int idx2,
                                        list<int>& outputIdxs)
  {
    outputIdxs.clear();

    for (int i = 0; i < inputSpectra->size(); i++)
    {
      if (i != idx1 && i != idx2)
      {
        outputIdxs.push_back(i);
      }
    }
  }
}

