/*
 * MappedContig.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#include "MappedContig.h"

namespace specnets
{
  const char* MappedVertex::Labels[] = { "un-annotated", "incorrect",
                                         "chimeric", "correct" };

  const char* MappedGap::Labels[] = { "un-annotated", "incorrect", "correct" };

  /**
   * Annotates this mapped contig by setting all class variables
   * @param parentContigs SpecSet containing contig
   * @param idx index of contig
   * @param parentComponents ab_info of contigs
   * @param overlaps for each mapped contig (spectrum idx), this matches
   *   contig peak indicies (mass) to protein residue indicies (intensity).
   * @param protIdx for each contig, this details the mapped protein (0) and
   *   whether the contig is reversed (2)
   * @param mappedSpecs mapped star spectra the contigs assemble
   * @param protein_spectra all target proteins as cummulative masses
   * @param peakTol Da peak tolerance
   * @return
   */
  void MappedContig::mapProt(SpecSet* parentContigs,
                             int idx,
                             abinfo_t* parentComponents,
                             SpecSet* overlaps,
                             vector<vector<int> >* protIdx,
                             vector<MappedSpectrum>* mappedSpecs,
                             SpecSet* protein_spectra,
                             float peakTol)
  {
    //DEBUG_MSG("received contig " << idx);

    Spectrum::operator =((*parentContigs)[idx]);
    this->setPeakTolerance(peakTol);
    /*
     DEBUG_MSG("copied contig " << idx);

     DEBUG_VAR(peakList.size());

     output(cerr);
     */
    reset();
    index = idx;

    if (size() == 0)
    {
      return;
    }

    int prot_idx = (*protIdx)[index][0];
    reversed = (*protIdx)[index][2] == 1;
    vertMapped = prot_idx >= 0 && (*overlaps)[index].size() > 0;
    vertProtIdx = prot_idx;
    length = size();
    chimeric = false;

    vector<int>* spectrum_idxs = &(*parentComponents)[index].first.first;
    numSpecs = spectrum_idxs->size();
    mappedSpectra.resize(numSpecs);
    MZRange parentMassRange;
    set<MZRange> massRanges;

    map<int, int> prot_count;

    int maxProt;
    int protCount = -1;

    for (int i = 0; i < numSpecs; i++)
    {
      int specIdx = (*spectrum_idxs)[i];
      mappedSpectra[i] = &(*mappedSpecs)[specIdx];
      parentMassRange.set(mappedSpectra[i]->parentMass, 1.0, 2.5);
      massRanges.insert(parentMassRange);

      if (mappedSpectra[i]->mapped)
      {
        for (map<int, list<int> >::iterator protIt =
            mappedSpectra[i]->residueIdxs.begin(); protIt
            != mappedSpectra[i]->residueIdxs.end(); protIt++)
        {
          if (prot_count.count(protIt->first) == 0)
          {
            prot_count[protIt->first] = 1;
          }
          else
          {
            prot_count[protIt->first]++;
          }

          if (prot_count[protIt->first] > protCount)
          {
            maxProt = protIt->first;
            protCount = prot_count[protIt->first];
          }
        }
      }
    }

    starMapped = protCount >= 0;
    starProtIdx = (starMapped) ? maxProt : -1;
    numPeptides = massRanges.size();

    if (starMapped && prot_count.count(vertProtIdx) > 0)
    {
      starProtIdx = vertProtIdx;
    }

    if (index == 15)
    {
      DEBUG_VAR(vertProtIdx);
      DEBUG_VAR(maxProt);
    }

    abruijnVerts.resize(length);
    abruijnGaps.resize(length - 1);

    vector<pair<vector<int> , vector<double> > > *abruijn_verts =
        &(*parentComponents)[index].second;
    vector<int>* spectrumIdxs;
    vector<double>* peakMasses;
    set<int> mappedProteinIndicies;
    set<int> mappedResidueIndicies;
    int spectrumIdx, peakIdx, proteinIdx, residueIdx;
    float peakMass;

    int overlap_idx = -1;
    int mappedPeak = -2, mappedResidue = -2;

    map<int, int> specToProt;
    for (int i = 0; i < (*overlaps)[index].size(); i++)
    {
      mappedPeak = floatToInt((*overlaps)[index][i][0]);
      mappedPeak = (reversed) ? length - 1 - mappedPeak : mappedPeak;
      mappedResidue = floatToInt((*overlaps)[index][i][1]);

      specToProt[mappedPeak] = mappedResidue;
    }

    /*
     DEBUG_VAR(vertMapped);
     DEBUG_VAR(overlaps->size());
     DEBUG_VAR(vertProtIdx);
     DEBUG_VAR((*overlaps)[index].size());
     */

    if (vertMapped)
    {
      int firstMapped = floatToInt((*overlaps)[index][0][1]);
      int lastMapped = floatToInt((*overlaps)[index][(*overlaps)[index].size()
          - 1][1]);
      firstResidue = min(firstMapped, lastMapped);
      lastResidue = max(firstMapped, lastMapped) - 1;
    }
    else
    {
      firstResidue = -1;
      lastResidue = -1;
    }

    map<int, pair<bool, bool> > prevAnnotPeaks;
    map<int, pair<bool, bool> > nextAnnotPeaks;
    map<int, pair<bool, bool> > pathAnnotPeaks;
    pair<bool, bool> pairAnnot;
    //check every contig peak


    for (int i = 0; i < abruijn_verts->size(); i++)
    {
      spectrumIdxs = &(*abruijn_verts)[i].first;
      peakMasses = &(*abruijn_verts)[i].second;

      mappedProteinIndicies.clear();
      mappedResidueIndicies.clear();

      abruijnVerts[i].mapped = specToProt.count(i) > 0;
      if (abruijnVerts[i].mapped)
      {
        abruijnVerts[i].residueIdx = specToProt[i];
      }
      else
      {
        abruijnVerts[i].residueIdx = -1;
      }

      abruijnVerts[i].set(peakList[i][0], peakList[i][1], peakTol);
      MZRange range = (MZRange)abruijnVerts[i];

      abruijnVerts[i].endpt = (range == parentMass - AAJumps::massHion)
          || (range == 0) || (range == AAJumps::massH2O) || (range
          == parentMass - AAJumps::massMH);

      abruijnVerts[i].numPeaks = spectrumIdxs->size();
      abruijnVerts[i].numBPeaks = 0;
      abruijnVerts[i].numYPeaks = 0;
      abruijnVerts[i].numIncorPeaks = 0;
      abruijnVerts[i].starPeaks.resize(abruijnVerts[i].numPeaks);

      map<int, int> resIdxs;
      set<int> checkedRes;
      int maxCountRes = 0;
      int numMapped = 0;
      nextAnnotPeaks.clear();

      if (i > 0)
      {
        abruijnGaps[i - 1].sameSpec = false;
        abruijnGaps[i - 1].label = 0;
      }

      //check all derived spectrum peaks
      for (int j = 0; j < spectrumIdxs->size(); j++)
      {
        spectrumIdx = (*spectrumIdxs)[j];
        MappedSpectrum* spec = &(*mappedSpecs)[spectrumIdx];
        peakMass = (*peakMasses)[j];
        MZRange peakMassR(peakMass, 1.0, 0.01);
        if (spec->massToPeak.count(peakMassR) == 0)
        {
          ERROR_MSG("Inconsistent peaks mass " << peakMass << " between abinfo and MappedSpectrum " << spectrumIdx << " (" << spec->reversed << ")");
          //spec->output(cout);
          continue;
        }

        peakIdx = spec->massToPeak[peakMassR];

        MappedPeak* pk = &spec->peaks[peakIdx];
        abruijnVerts[i].starPeaks[j] = pk;

        if (spec->identified)
        {
          abruijnVerts[i].annotated = true;
        }

        pairAnnot.first = pk->annotation.find("b") != string::npos;
        pairAnnot.second = pk->annotation.find("y") != string::npos;
        nextAnnotPeaks[spectrumIdx] = pairAnnot;

        if (i > 0 && prevAnnotPeaks.count(spectrumIdx) > 0)
        {
          abruijnGaps[i - 1].sameSpec = true;
          if (pk->mapped && ((pairAnnot.first
              && prevAnnotPeaks[spectrumIdx].first) || (pairAnnot.second
              && prevAnnotPeaks[spectrumIdx].second))
              && spec->residueIdxs.count(starProtIdx) > 0)
          {
            abruijnGaps[i - 1].label = 2;
          }
          else if (abruijnGaps[i - 1].label != 2 && spec->identified)
          {
            abruijnGaps[i - 1].label = 1;
          }
        }

        if (pk->mapped)
        {
          if (spec->residueIdxs.size() > 0)
          {
            numMapped++;
          }
          abruijnVerts[i].numBPeaks += (pairAnnot.first) ? 1 : 0;
          abruijnVerts[i].numYPeaks += (pairAnnot.second) ? 1 : 0;
          checkedRes.clear();
          if (spec->residueIdxs.count(starProtIdx) > 0)
          {
            for (list<int>::iterator resIt =
                spec->residueIdxs[starProtIdx].begin(); resIt
                != spec->residueIdxs[starProtIdx].end(); resIt++)
            {

              int idxCheck[] = { pk->BpeptideIdx, pk->YpeptideIdx };
              for (int z = 0; z < 2; z++)
              {
                if (idxCheck[z] < 0)
                {
                  continue;
                }
                int globalResIdx = *resIt + idxCheck[z];

                if (checkedRes.count(globalResIdx) > 0)
                {
                  continue;
                }
                checkedRes.insert(globalResIdx);

                if (resIdxs.count(globalResIdx) == 0)
                {
                  resIdxs[globalResIdx] = 1;
                  maxCountRes = max(maxCountRes, 1);
                }
                else if (resIdxs.count(globalResIdx) > 0)
                {
                  resIdxs[globalResIdx]++;
                  maxCountRes = max(maxCountRes, resIdxs[globalResIdx]);
                }
              }
            }
          }
        }
        else if (spec->identified)
        {
          abruijnVerts[i].numIncorPeaks++;
        }

      }
      //check all derived spectrum peaks
      for (int j = 0; j < spectrumIdxs->size(); j++)
      {
        MappedPeak* pk = abruijnVerts[i].starPeaks[j];
        spectrumIdx = pk->specIdx;
        MappedSpectrum* spec = &(*mappedSpecs)[pk->specIdx];
        pairAnnot.first = pk->annotation.find("b") != string::npos;
        pairAnnot.second = pk->annotation.find("y") != string::npos;

        if (abruijnGaps[i - 1].label != 2 && i > 0 && pk->mapped
            && pathAnnotPeaks.count(spectrumIdx) > 0
            && spec->residueIdxs.count(starProtIdx) > 0)
        {
          if ((pairAnnot.first && pathAnnotPeaks[spectrumIdx].first)
              || (pairAnnot.second && pathAnnotPeaks[spectrumIdx].second))
          {
            abruijnGaps[i - 1].label = 2;
            //DEBUG_MSG("Caught correct gap in contig " << idx << " between vertices " << i -1 << " and " << i);
          }
        }
      }

      abruijnVerts[i].starMapped = numMapped > 0;
      if (abruijnVerts[i].numBPeaks > 0)
      {
        abruijnVerts[i].annotation = "b";
      }
      else if (abruijnVerts[i].numYPeaks > 0)
      {
        abruijnVerts[i].annotation = "y";
      }
      else
      {
        abruijnVerts[i].annotation = "";
      }

      /*if (idx == 18 && i == 4)
       {
       DEBUG_MSG("found " << maxCountRes << " maxCountRes and " << numMapped << " numMapped");
       } */

      if (!abruijnVerts[i].annotated)
      {
        abruijnVerts[i].label = 0;
      }
      else if (maxCountRes > 0 && maxCountRes == numMapped)
      {
        abruijnVerts[i].label = 3;
      }
      else if (maxCountRes == 0)
      {
        abruijnVerts[i].label = 1;
      }
      else
      {
        abruijnVerts[i].label = 2;
        chimeric = true;
      }

      if (i > 0)
      {
        MZRange gapM = ((MZRange)abruijnVerts[i]) - ((MZRange)abruijnVerts[i
            - 1]);
        abruijnGaps[i - 1].set(gapM);

        abruijnGaps[i - 1].vertMapped = abruijnVerts[i].mapped
            && abruijnVerts[i - 1].mapped;

        abruijnGaps[i - 1].starMapped = abruijnVerts[i].starMapped
            && abruijnVerts[i - 1].starMapped;

        abruijnGaps[i - 1].annotated = abruijnGaps[i - 1].label > 0;

        if (abruijnGaps[i - 1].vertMapped)
        {
          abruijnGaps[i - 1].mappedMass
              = abs((*protein_spectra)[vertProtIdx][specToProt[i]][0]
                  - (*protein_spectra)[vertProtIdx][specToProt[i - 1]][0]);
        }
        else
        {
          abruijnGaps[i - 1].mappedMass = -1.0;
        }

        if (abruijnGaps[i - 1].label == 2)
        {
          for (map<int, pair<bool, bool> >::iterator annotIt =
              prevAnnotPeaks.begin(); annotIt != prevAnnotPeaks.end(); annotIt++)
          {
            pathAnnotPeaks[annotIt->first] = annotIt->second;
          }
        }
        else
        {
          pathAnnotPeaks.clear();
        }
      }
      prevAnnotPeaks = nextAnnotPeaks;
    }

    pathAnnotPeaks.clear();
    prevAnnotPeaks.clear();
    nextAnnotPeaks.clear();
    for (int i = abruijn_verts->size() - 1; i >= 0; i--)
    {

      //check all derived spectrum peaks
      for (int j = 0; j < abruijnVerts[i].numPeaks; j++)
      {
        MappedPeak* pk = abruijnVerts[i].starPeaks[j];
        spectrumIdx = pk->specIdx;
        MappedSpectrum* spec = &(*mappedSpecs)[pk->specIdx];
        pairAnnot.first = pk->annotation.find("b") != string::npos;
        pairAnnot.second = pk->annotation.find("y") != string::npos;
        nextAnnotPeaks[spectrumIdx] = pairAnnot;

        if (i < abruijn_verts->size() - 1 && abruijnGaps[i].label != 2
            && pk->mapped && pathAnnotPeaks.count(spectrumIdx) > 0
            && spec->residueIdxs.count(starProtIdx) > 0)
        {
          if ((pairAnnot.first && pathAnnotPeaks[spectrumIdx].first)
              || (pairAnnot.second && pathAnnotPeaks[spectrumIdx].second))
          {
            abruijnGaps[i].label = 2;
            //DEBUG_MSG("Caught correct gap in contig " << idx << " between vertices " << i << " and " << i + 1);
          }
        }
      }

      if (i < abruijn_verts->size() - 1)
      {
        abruijnGaps[i].annotated = abruijnGaps[i].label > 0;
        if (abruijnGaps[i].label == 2)
        {
          for (map<int, pair<bool, bool> >::iterator annotIt =
              prevAnnotPeaks.begin(); annotIt != prevAnnotPeaks.end(); annotIt++)
          {
            pathAnnotPeaks[annotIt->first] = annotIt->second;
          }
        }
        else
        {
          pathAnnotPeaks.clear();
        }
      }
      prevAnnotPeaks = nextAnnotPeaks;
    }
  }

  /**
   * @return number of assembled annotated spectra
   */
  int MappedContig::getNumAnnotSpec(void)
  {
    int count = 0;
    for (int i = 0; i < mappedSpectra.size(); i++)
    {
      count += (mappedSpectra[i]->identified) ? 1 : 0;
    }
    return count;
  }

  /**
   * @return number of assembled mapped spectra
   */
  int MappedContig::getNumMappedSpec(void)
  {
    int count = 0;
    for (int i = 0; i < mappedSpectra.size(); i++)
    {
      count += (mappedSpectra[i]->mapped) ? 1 : 0;
    }
    return count;
  }

  /**
   * @return % of annotated vertices from b ions
   */
  float MappedContig::getPercBVerts(void)
  {
    float numerator = 0;
    float denominator = 0;

    for (int i = 0; i < length; i++)
    {
      if (abruijnVerts[i].annotated)
      {
        denominator += 1.0;
        numerator += (abruijnVerts[i].annotation.find("b") != string::npos)
            ? 1.0 : 0;
      }
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @return % of annotated vertices from y ions
   */
  float MappedContig::getPercYVerts(void)
  {
    float numerator = 0;
    float denominator = 0;

    for (int i = 0; i < length; i++)
    {
      if (abruijnVerts[i].annotated)
      {
        denominator += 1.0;
        numerator += (abruijnVerts[i].annotation.find("y") != string::npos)
            ? 1.0 : 0;
      }
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @return % of annotated vertices from b or y ions
   */
  float MappedContig::getPercBYVerts(void)
  {
    return getPercYVerts() + getPercBVerts();
  }

  /**
   * @return % of vertices that match label. if _label == 0, denominator is
   *   # of all vertices. if _label > 0, denominator is # of annotated vertices
   */
  float MappedContig::getPercVerts(short _label)
  {
    float numerator = 0;
    float denominator = 0;

    for (int i = 0; i < length; i++)
    {
      if (_label > 0)
      {
        if (abruijnVerts[i].annotated)
        {
          denominator += 1.0;
          numerator += (abruijnVerts[i].label == _label) ? 1.0 : 0;
        }
      }
      else
      {
        denominator += 1.0;
        numerator += (abruijnVerts[i].label == _label) ? 1.0 : 0;
      }
    }
    /*
     if (_label == 2 && (numerator / denominator) > 0.3) {
     DEBUG_MSG("Check contig " << index << " with " << (numerator / denominator) * 100.0 << "% chimeric vertices");
     }
     if (_label == 1 && (numerator / denominator) > 0.3) {
     DEBUG_MSG("Check contig " << index << " with " << (numerator / denominator) * 100.0 << "% incorrect vertices");
     }
     */
    return (numerator / denominator) * 100.0;
  }

  /**
   * @return % of gaps that match accuracy label. if _label == 0, denominator is
   *   # of all gaps. if _label > 0, denominator is # of annotated gaps
   */
  pair<float, float> MappedContig::getPercCallsAcc(short _label)
  {
    float numerator = 0;
    float denominator = 0;

    for (int i = 0; i < length - 1; i++)
    {
      if (_label > 0)
      {
        if (abruijnGaps[i].annotated)
        {
          denominator += 1.0;
          numerator += (abruijnGaps[i].label == _label) ? 1.0 : 0;
        }
      }
      else
      {
        denominator += 1.0;
        numerator += (!abruijnGaps[i].annotated) ? 1.0 : 0;
      }
    }
    return pair<float, float> (numerator, denominator);
  }

  /**
   * @return % of regions between mapped vertices that match length label:
   *   0 : region spans one AA in the database with no no un-mapped vertices in between
   *   1 : region spans multiple AA in the database with no un-mapped vertices in between
   *   2 : region spans one or more AA in the database with un-mapped vertices in between
   */
  pair<float, float> MappedContig::getPercCallsLen(short _label)
  {
    float numerator = 0;
    float denominator = 0;

    if (abruijnVerts.size() <= 1)
    {
      return pair<float, float> (numerator, denominator);
    }

    for (int i = 0; i < abruijnVerts.size() - 1; i++)
    {
      if (!abruijnVerts[i].mapped)
      {
        continue;
      }
      bool foundBetween = false;

      for (int j = i + 1; j < abruijnVerts.size(); j++)
      {
        if (!abruijnVerts[j].mapped)
        {
          foundBetween = true;
          continue;
        }

        int gapSpan = abs(abruijnVerts[j].residueIdx - abruijnVerts[i].residueIdx);

        denominator += 1.0;

        switch (_label)
        {
        case 0:
          numerator += (gapSpan == 1 && !foundBetween) ? 1.0 : 0.0;
          break;
        case 1:
          numerator += (gapSpan > 1 && !foundBetween) ? 1.0 : 0.0;
          break;
        case 2:
          numerator += (foundBetween) ? 1.0 : 0.0;
          break;
        default:
          ERROR_MSG("Un-specified label " << _label << "!!!")
          ;
          abort();
          break;
        }

        break;
      }
    }

    return pair<float, float> (numerator, denominator);
  }
}
