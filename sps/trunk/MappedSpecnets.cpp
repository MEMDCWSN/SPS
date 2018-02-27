/*
 * MappedSpecnets.cpp
 *
 *  Created on: Mar 2, 2011
 *      Author: aguthals
 */

#include "MappedSpecnets.h"

namespace specnets
{

  MappedSpecnets::MappedSpecnets(void)
  {
    contigs = new vector<MappedContig>;
    spectra = new vector<MappedSpectrum>;
  }

  MappedSpecnets::~MappedSpecnets(void)
  {
    delete contigs;
    delete spectra;
  }

  /**
   * Maps this SPS project to proteins
   * @param sps_contigs specnets contigs
   * @param sps_components contig component information detailing which
   *   spectra and which peaks were included in each contig
   * @param star_spectra PRM spectra assembled into contigs
   * @param matchma_overlaps matchma output detailing where mapped
   *   contigs overlap on their protein
   * @param matchma_prot_idx matchma output detailing which contigs are
   *   mapped to which proteins
   * @param _proteins sequence of target proteins in same order as matchma
   *   saw them
   * @param protein_spectra target proteins as PRM spectra
   * @param spectrum_ids sequence of PRM spectra annotations in same order
   *   as in_star_spectra. Must be in specnets format
   * @param input_ion_types loaded MS2ScoringModel with b and y ion types
   *   specified
   * @param in_peak_tol peak tolerance in Da
   * @return
   */
  void MappedSpecnets::mapProt(SpecSet* sps_contigs,
                               abinfo_t* sps_components,
                               SpecSet* star_spectra,
                               SpecSet* matchma_overlaps,
                               vector<vector<int> >* matchma_prot_idx,
                               vector<string>* _proteins,
                               SpecSet* protein_spectra,
                               vector<string>* spectrum_ids,
                               MS2ScoringModel& input_ion_types,
                               float peak_tol)
  {

    contigs->resize(sps_contigs->size());
    spectra->resize(star_spectra->size());
    proteins = _proteins;
    peptides = spectrum_ids;

    star_spectra->addZPMpeaks(peak_tol, 0, true);

    DEBUG_TRACE;

    vector<bool> stars_reversed(star_spectra->size(), false);
    for (abinfo_t::iterator seqIt = sps_components->begin();
        seqIt != sps_components->end(); seqIt++)
    {
      if (seqIt->second.second.size() == 0)
        continue;

      for (int j = 0; j < seqIt->second.second.size(); j++)
      {
        pair<vector<int>, vector<double> >* vert = &seqIt->second.second[j];
        for (int k = 0; k < vert->first.size(); k++)
        {
          int specIdx = (*vert).first[k];
          float mass = (*vert).second[k];

          /*if (specIdx == 16299) {
           DEBUG_MSG("Contig " << seqIt->first << " peak=" << j << " spectrum " << specIdx << " w/ mass " << mass);
           }*/
          int peakIdx = (*star_spectra)[specIdx].findClosest(mass);

          if (!MZRange::EqualWithinRange(mass,
                                         (*star_spectra)[specIdx][peakIdx][0],
                                         0.01))
          {
            /*if (specIdx == 16299) {
             DEBUG_MSG("Peak mass " << mass << " causing reversal of spectrum " << 11366 << " in contig " << seqIt->first);
             }*/
            stars_reversed[specIdx] = true;
          }
        }
      }
    }

    DEBUG_TRACE;

    for (int i = 0; i < spectra->size(); i++)
    {

      (*spectra)[i].mapProt(star_spectra,
                            i,
                            (*spectrum_ids)[i],
                            proteins,
                            stars_reversed[i],
                            input_ion_types,
                            peak_tol);

    }

    DEBUG_TRACE;

    for (int i = 0; i < contigs->size(); i++)
    {
      //DEBUG_VAR(i);
      (*contigs)[i].mapProt(sps_contigs,
                            i,
                            sps_components,
                            matchma_overlaps,
                            matchma_prot_idx,
                            spectra,
                            protein_spectra,
                            peak_tol);
    }

    DEBUG_TRACE;

    for (int p = 0; p < proteins->size(); p++)
    {
      if ((*proteins)[p].length() > 0)
      {
        target_prots.insert(p);
      }
    }

    DEBUG_TRACE;

  }

  /**
   * @return number of contigs w/ > 0 peaks
   */
  int MappedSpecnets::getNumContigs(void)
  {
    int count = 0;
    for (int i = 0; i < contigs->size(); i++)
    {
      count += ((*contigs)[i].length > 0) ? 1 : 0;
    }
    return count;
  }

  /**
   * @return number of spectra w/ > 0 peaks
   */
  int MappedSpecnets::getNumSpectra(void)
  {
    int count = 0;
    for (int i = 0; i < spectra->size(); i++)
    {
      count += ((*spectra)[i].size() > 0) ? 1 : 0;
    }
    return count;
  }

  /**
   * @return number of identified spectra
   */
  int MappedSpecnets::getNumSpecIdent(void)
  {
    int count = 0;

    for (int i = 0; i < spectra->size(); i++)
    {
      count += ((*spectra)[i].identified) ? 1 : 0;
    }
    return count;
  }

  /**
   * @return number of spectra mapped to target proteins
   */
  int MappedSpecnets::getNumSpecMapped(int prot_idx)
  {
    int count = 0;

    for (int i = 0; i < spectra->size(); i++)
    {
      if (prot_idx < 0 && (*spectra)[i].mapped)
      {
        count++;
      }
      else if (prot_idx >= 0 && (*spectra)[i].residueIdxs.count(prot_idx) > 0)
      {
        count++;
      }
    }
    return count;
  }

  /**
   * @return number of modified spectra
   */
  int MappedSpecnets::getNumSpecModified(void)
  {
    int count = 0;

    for (int i = 0; i < spectra->size(); i++)
    {
      count += ((*spectra)[i].modified) ? 1 : 0;
    }
    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return number of contigs mapped to prot_idx by matchma or spectrum IDs
   */
  int MappedSpecnets::getNumContigsMapped(int prot_idx)
  {

    int count = 0;

    int num6mix = 0;
    int numBoth = 0;
    int numVenom = 0;

    for (int i = 0; i < contigs->size(); i++)
    {

      if ((*contigs)[i].length == 0)
      {
        continue;
      }

      int found6Mix = 0;
      int foundVenom = 0;

      for (int j = 0; j < (*contigs)[i].mappedSpectra.size(); j++)
      {
        int specIdx = (*contigs)[i].mappedSpectra[j]->specIdx;

        if (specIdx < 11010 || (specIdx >= 44879 && specIdx < 58919))
        {
          found6Mix++;
        }
        else
        {
          foundVenom++;
        }
      }

      if (found6Mix > 3 && foundVenom > 3)
      {
        numBoth++;
      }
      else if (found6Mix > 3)
      {
        num6mix++;
      }
      else if (foundVenom > 3)
      {
        numVenom++;
      }

      if (prot_idx < 0
          && (target_prots.count((*contigs)[i].vertProtIdx) > 0
              || target_prots.count((*contigs)[i].starProtIdx) > 0))
      {
        count++;
      }
      else if (prot_idx >= 0
          && ((*contigs)[i].vertProtIdx == prot_idx
              || (*contigs)[i].starProtIdx == prot_idx))
      {
        count++;
      }
    }

    //DEBUG_VAR(numBoth);
    //DEBUG_VAR(num6mix);
    //DEBUG_VAR(numVenom);
    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return number of contigs mapped to prot_idx by matchma
   */
  int MappedSpecnets::getNumContigsVertMapped(int prot_idx)
  {

    int count = 0;

    for (int i = 0; i < contigs->size(); i++)
    {

      if ((*contigs)[i].length == 0)
      {
        continue;
      }

      if (prot_idx < 0 && target_prots.count((*contigs)[i].vertProtIdx) > 0)
      {
        count++;
      }
      else if (prot_idx >= 0 && (*contigs)[i].vertProtIdx == prot_idx)
      {
        count++;
      }
    }
    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return number of contigs mapped to prot_idx by spectrum IDs
   */
  int MappedSpecnets::getNumContigsSpecMapped(int prot_idx)
  {

    int count = 0;

    for (int i = 0; i < contigs->size(); i++)
    {

      if ((*contigs)[i].length == 0)
      {
        continue;
      }

      if (prot_idx < 0 && target_prots.count((*contigs)[i].starProtIdx) > 0)
      {
        count++;
      }
      else if (prot_idx >= 0 && (*contigs)[i].starProtIdx == prot_idx)
      {
        count++;
      }
    }
    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return number of contigs mapped to prot_idx by matchma and spectrum IDs
   */
  int MappedSpecnets::getNumContigsVertSpecMapped(int prot_idx)
  {

    int count = 0;

    for (int i = 0; i < contigs->size(); i++)
    {

      if ((*contigs)[i].length == 0)
      {
        continue;
      }

      if (prot_idx < 0
          && (target_prots.count((*contigs)[i].vertProtIdx) > 0
              && target_prots.count((*contigs)[i].starProtIdx) > 0))
      {
        count++;
      }
      else if (prot_idx >= 0
          && ((*contigs)[i].vertProtIdx == prot_idx
              && (*contigs)[i].starProtIdx == prot_idx))
      {
        count++;
      }
    }
    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return percent of protein covered by identified spectra
   */
  float MappedSpecnets::getPercSpecCov(int prot_idx)
  {
    vector<bool> deflt_vec(0);
    vector<vector<bool> > covered_bins(proteins->size(), deflt_vec);

    float numerator = 0;
    float denominator = 0;

    MappedSpectrum* spec;

    for (int p = 0; p < proteins->size(); p++)
    {
      if ((*proteins)[p].length() == 0 || (prot_idx >= 0 && p != prot_idx))
      {
        continue;
      }
      covered_bins[p].resize((*proteins)[p].length(), false);
      denominator += (float)covered_bins[p].size();

      for (int i = 0; i < spectra->size(); i++)
      {

        spec = &(*spectra)[i];

        if ((!spec->mapped) || (spec->residueIdxs.count(p) == 0))
        {
          continue;
        }

        for (list<int>::iterator startIt = spec->residueIdxs[p].begin();
            startIt != spec->residueIdxs[p].end(); startIt++)
        {
          int start_idx = *startIt;
          int end_idx = start_idx + spec->AAPeptideSeq.length();
          for (int r = start_idx; r < end_idx; r++)
          {
            if (!covered_bins[p][r])
            {
              covered_bins[p][r] = true;
              numerator += 1.0;
            }
          }
        }
      }
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return percent of protein covered by matchma mapped contigs
   */
  float MappedSpecnets::getPercSeqCov(int prot_idx)
  {
    vector<bool> deflt_vec(0);
    vector<vector<bool> > covered_bins(proteins->size(), deflt_vec);

    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;
    set<int> prots_count;
    for (int p = 0; p < proteins->size(); p++)
    {
      if ((*proteins)[p].length() == 0 || (prot_idx >= 0 && p != prot_idx))
      {
        continue;
      }
      covered_bins[p].resize((*proteins)[p].length(), false);
      prots_count.insert(p);
      //DEBUG_MSG("prot " << contig->vertProtIdx << " has size " << covered_bins[p].size());
      denominator += (float)covered_bins[p].size();
    }

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (contig->size() == 0 || (!contig->vertMapped)
          || (prots_count.count(contig->vertProtIdx) == 0))
      {
        continue;
      }

      for (int r = contig->firstResidue; r <= contig->lastResidue; r++)
      {
        //if (r < 0 || r > covered_bins[contig->vertProtIdx].size()) {
        //  DEBUG_MSG("invalid index " << r << " for prot " << contig->vertProtIdx);
        //}
        if (!covered_bins[contig->vertProtIdx][r])
        {
          covered_bins[contig->vertProtIdx][r] = true;
          numerator += 1.0;
        }
      }

    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return average number of contigs covering each covered residue
   */
  float MappedSpecnets::getCovRedundancy(int prot_idx)
  {
    vector<bool> deflt_vec(0);
    vector<vector<bool> > covered_bins(proteins->size(), deflt_vec);

    float numerator = 0;
    float denominator = 0;
    set<int> prots_count;
    MappedContig* contig;
    for (int p = 0; p < proteins->size(); p++)
    {
      if ((*proteins)[p].length() == 0 || (prot_idx >= 0 && p != prot_idx))
      {
        continue;
      }
      covered_bins[p].resize((*proteins)[p].length(), false);
      prots_count.insert(p);
    }

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if ((!contig->vertMapped)
          || (prots_count.count(contig->vertProtIdx) == 0))
      {
        continue;
      }

      for (int r = contig->firstResidue; r <= contig->lastResidue; r++)
      {

        if (!covered_bins[contig->vertProtIdx][r])
        {
          covered_bins[contig->vertProtIdx][r] = true;
          denominator += 1.0;
          numerator += 1.0;
        }
        else
        {
          numerator += 1.0;
        }
      }
    }
    return numerator / denominator;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return number of spectra assembled into contigs mapping to protein by matchma or spectrum IDs
   */
  int MappedSpecnets::getNumAssembledSpec(int prot_idx)
  {
    int count = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];

      if (target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      count += contig->numSpecs;
    }

    return count;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return average number of spectra per contig mapping to protein by matchma or spectrum IDs
   */
  float MappedSpecnets::getSpecPerContig(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      numerator += (float)contig->numSpecs;
      denominator += 1.0;
    }

    return numerator / denominator;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return average number of peptides per contig mapping to protein by matchma or spectrum IDs
   */
  float MappedSpecnets::getPepPerContig(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      numerator += (float)contig->numPeptides;
      denominator += 1.0;
    }

    return numerator / denominator;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return average Da of each contig mapping to protein by matchma or spectrum IDs
   */
  float MappedSpecnets::getDaLengthPerContig(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (prot_idx < 0 && target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      numerator += contig->abruijnVerts[contig->length - 1].getMass()
          - contig->abruijnVerts[0].getMass();
      denominator += 1.0;
    }

    return numerator / denominator;

  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return average AA of each contig mapping to protein by matchma
   */
  float MappedSpecnets::getAALengthPerContig(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->vertProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx))
      {
        continue;
      }

      numerator += (float)(contig->lastResidue - contig->firstResidue + 1);
      denominator += 1.0;
    }

    return numerator / denominator;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return pair of longest contig (second) mapping to protein by matchma and its length (first)
   */
  pair<int, int> MappedSpecnets::getLongestAAContig(int prot_idx)
  {
    int maxLength = -1;
    int maxIdx = -1;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->vertProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx))
      {
        continue;
      }

      int length = contig->lastResidue - contig->firstResidue + 1;

      if (length > maxLength)
      {
        maxLength = length;
        maxIdx = i;
      }
    }

    return pair<int, int>(maxLength, maxIdx);
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @param label 0, 1, 2, or 3
   * @return if label != 0, percent of abruijn vertices that have label out of
   *  those that are annotated. if label == 0, percent of abruijn vertices
   *  that have label out of those that are not annotated
   */
  float MappedSpecnets::getPercVerts(int prot_idx, short label)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      int checkC = 0;

      for (int j = 0; j < contig->length; j++)
      {
        if (label > 0 && contig->abruijnVerts[j].annotated)
        {
          denominator += 1.0;
          if (contig->abruijnVerts[j].label == label)
          {
            numerator += 1.0;
            checkC++;
          }
        }
        else if (label == 0)
        {
          denominator += 1.0;
          numerator += (contig->abruijnVerts[j].annotated) ? 0.0 : 1.0;
        }
      }
      /*
       if (label == 1 && checkC > 5) {
       DEBUG_MSG("Check contig " << i << " with " << checkC << " incorrect verts");
       }
       */
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return pair of percent b ions (first) and percent y ions (second) out
   *   of all annotated vertices
   */
  pair<float, float> MappedSpecnets::getPercBYVerts(int prot_idx)
  {
    float numeratorB = 0;
    float numeratorY = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      for (int j = 0; j < contig->length; j++)
      {
        if (contig->abruijnVerts[j].annotated)
        {
          denominator += 1.0;
          if (contig->abruijnVerts[j].annotation.find("b") != string::npos)
          {
            numeratorB += 1.0;
          }
          else if (contig->abruijnVerts[j].annotation.find("y") != string::npos)
          {
            numeratorY += 1.0;
          }
        }
      }

    }
    return pair<float, float>((numeratorB / denominator) * 100.0,
                              (numeratorY / denominator) * 100.0);
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @param label 0, 1, or 2
   * @return if label != 0, percent of abruijn gaps that have accuracy label out of
   *  those that are annotated. if label == 0, percent of abruijn gaps
   *  that are annotated out of all
   */
  float MappedSpecnets::getPercCallsAcc(int prot_idx, short label)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      pair<float, float> res = contig->getPercCallsAcc(label);
      denominator += res.second;
      numerator += res.first;
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @param label 0, 1, or 2
   * @return % of regions between mapped vertices that match length label:
   *   0 : region spans one AA in the database
   *   1 : region spans multiple AA in the database with no un-mapped vertices in between
   *   2 : region spans multiple AA in the database with un-mapped vertices in between
   */
  float MappedSpecnets::getPercCallsLen(int prot_idx, short label)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      pair<float, float> res = contig->getPercCallsLen(label);
      denominator += res.second;
      numerator += res.first;
    }
    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return if prot_idx < 0, percent of all assembled spectra that are
   *   identified. Otherwise, the percent of all assembled spectra in a contig
   *   mapped to protein prot_idx that are identified.
   */
  float MappedSpecnets::getPercAssemSpecIdent(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];

      if (target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      denominator += contig->numSpecs;

      for (int s = 0; s < contig->mappedSpectra.size(); s++)
      {
        numerator += (contig->mappedSpectra[s]->identified) ? 1.0 : 0;
      }
    }

    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return if prot_idx < 0, percent of all assembled identified spectra that are
   *   identified for contigs mapping to prot_idx (-1 for all proteins)
   */
  float MappedSpecnets::getPercAssemIdentSpecMod(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];

      if (target_prots.count(contig->vertProtIdx) == 0
          && target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx)
          && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      for (int s = 0; s < contig->mappedSpectra.size(); s++)
      {
        denominator += (contig->mappedSpectra[s]->identified) ? 1.0 : 0;
        numerator += (contig->mappedSpectra[s]->modified) ? 1.0 : 0;
      }
    }

    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @return of all contigs mapped to prot_idx, the percent of all matchma
   *   and inspect mapped abruijn vertices that assemble at least one
   *   identified MS/MS b/y peak mapping to the same matchma residue index.
   */
  float MappedSpecnets::getMatchmaAccuracy(int prot_idx)
  {
    float numerator = 0;
    float denominator = 0;

    MappedContig* contig;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];

      if (!contig->vertMapped)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->vertProtIdx != prot_idx))
      {
        continue;
      }
      int prot_idx_use = (prot_idx < 0) ? contig->vertProtIdx : prot_idx;

      if (target_prots.count(prot_idx_use) == 0)
      {
        continue;
      }

      for (int s = 0; s < contig->size(); s++)
      {
        if ((!contig->abruijnVerts[s].mapped)
            || (!contig->abruijnVerts[s].annotated))
        {
          continue;
        }
        denominator += 1.0;
        int resIdx = contig->abruijnVerts[s].residueIdx;
        for (int v = 0; v < contig->abruijnVerts[s].starPeaks.size(); v++)
        {
          if (!contig->abruijnVerts[s].starPeaks[v]->mapped)
          {
            continue;
          }

          int specIdx = contig->abruijnVerts[s].starPeaks[v]->specIdx;
          int BpepIdx = contig->abruijnVerts[s].starPeaks[v]->BpeptideIdx;
          int YpepIdx = contig->abruijnVerts[s].starPeaks[v]->YpeptideIdx;

          if ((*spectra)[specIdx].residueIdxs.count(prot_idx_use) == 0)
          {
            continue;
          }
          bool broke = false;
          for (list<int>::iterator resIt =
              (*spectra)[specIdx].residueIdxs[prot_idx_use].begin();
              resIt != (*spectra)[specIdx].residueIdxs[prot_idx_use].end();
              resIt++)
          {

            if (BpepIdx >= 0 && (*resIt) + BpepIdx == resIdx)
            {
              numerator += 1.0;
              broke = true;
              break;
            }
            if (YpepIdx >= 0 && (*resIt) + YpepIdx == resIdx)
            {
              numerator += 1.0;
              broke = true;
              break;
            }
          }
          if (broke)
          {
            break;
          }
        }
      }
    }

    return (numerator / denominator) * 100.0;
  }

  /**
   * @param prot_idx index of target protein or -1 to specify all proteins
   * @param outValues output vector that, over all annotated contig gaps,
   *   will contain a pair of integers for each index i where the first is
   *   the number of annotated gaps at position i and the second is the
   *   number of correct annotated gaps at position i. "position" is the
   *   number of gaps between a gap and the closest end of the contig
   *   sequence.
   * @param outputTable optional output data structure that, if specified,
   *   will mirror 'outValues' but in a table format suitable for OutputTable.cpp
   * @return
   */
  void MappedSpecnets::getGapAccVsPos(int prot_idx,
                                      vector<pair<int, int> >& outValues,
                                      vector<vector<pair<string, bool> > >* outputTable)
  {
    MappedContig* contig;
    outValues.resize(100, pair<int, int>(0, 0));
    int largestPos = 0;

    for (int i = 0; i < contigs->size(); i++)
    {
      contig = &(*contigs)[i];
      if (target_prots.count(contig->starProtIdx) == 0)
      {
        continue;
      }
      if (prot_idx >= 0 && (contig->starProtIdx != prot_idx))
      {
        continue;
      }

      for (int j = 0; j < contig->length - 1; j++)
      {
        if (contig->abruijnGaps[j].annotated)
        {
          int pos = min(j, contig->length - 2 - j);
          if (pos >= outValues.size())
          {
            outValues.resize(pos + 1, pair<int, int>(0, 0));
          }
          largestPos = max(largestPos, pos);

          outValues[pos].first += 1;
          if (contig->abruijnGaps[j].label == 2)
          {
            outValues[pos].second += 1;
          }
        }
      }
    }
    outValues.resize(largestPos + 1);

    if (outputTable != 0)
    {
      outputTable->resize(outValues.size() + 1);
      pair<string, bool> cell;

      cell.first = "Position";
      cell.second = true;
      (*outputTable)[0].resize(3);
      (*outputTable)[0][0] = cell;
      cell.first = "Annotated";
      (*outputTable)[0][1] = cell;
      cell.first = "Accurate";
      (*outputTable)[0][2] = cell;

      cell.second = false;

      for (int i = 0; i < outValues.size(); i++)
      {
        (*outputTable)[i + 1].resize(3);
        cell.first = parseInt(i + 1);
        (*outputTable)[i + 1][0] = cell;
        cell.first = parseInt(outValues[i].first);
        (*outputTable)[i + 1][1] = cell;
        cell.first = parseInt(outValues[i].second);
        (*outputTable)[i + 1][2] = cell;
      }
    }
  }
}
