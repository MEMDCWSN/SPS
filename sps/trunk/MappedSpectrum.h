/*
 * MappedSpectrum.h
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#ifndef MAPPEDSPECTRUM_H_
#define MAPPEDSPECTRUM_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "aminoacid.h"
#include "mzrange.h"
#include "Logger.h"
#include "spectrum.h"
#include "spectrum_scoring.h"
#include "SpecSet.h"

using namespace std;

namespace specnets
{
  /*
   * Spectrum peak mapped to protein by database search
   */
  class MappedPeak: public MZRange
  {
  public:

    // whether this peak is mapped to a peptide residue
    bool mapped;

    // b/y annotation
    string annotation;

    // peptide residue idx mapped to if b ion
    int BpeptideIdx;

    // peptide residue idx mapped to if y ion
    int YpeptideIdx;

    // whether the peak is an end point
    bool endpt;

    // index of peak in parent spectrum
    int peakIdx;

    // index of parent spectrum
    int specIdx;

    MappedPeak() :
      MZRange()
    {
      reset();

    }

    MappedPeak(MZRange& other) :
      MZRange(other)
    {
      reset();
    }

    void reset(void)
    {
      mapped = false;
      annotation = "";
      BpeptideIdx = -1;
      YpeptideIdx = -1;
      endpt = false;
      peakIdx = -1;
      specIdx = -1;
      set(0, 0, 0);
    }

  };

  /*
   * Spectrum mapped to a protein by database search
   */
  class MappedSpectrum: public Spectrum
  {
  public:
    // whether this spectrum is mapped to a protein or peptide
    bool mapped;

    // whether this spectrum is identified to a modified peptide sequence
    bool modified;

    // whether this spectrum is identified
    bool identified;

    // whether this star spectrum needs to be reversed to get assembled masses
    bool reversed;

    // peptide annotation
    string peptide;

    // peptide annotation without mods or special characters
    string AAPeptideSeq;

    // protein idxs -> residue idx of first mapped b0/y0 peak
    map<int, list<int> > residueIdxs;

    // index of spectrum
    int specIdx;

    // number of annotated PRM masses
    int numPRMs;

    // number of annotated SRM masses
    int numSRMs;

    // all mapped/un-mapped peaks
    vector<MappedPeak> peaks;

    // maps peak mass to peak index (abinfo stores masses, not peak indices)
    map<MZRange, int> massToPeak;

    MappedSpectrum(void) :
      Spectrum()
    {
      reset();
    }

    void reset(void)
    {
      mapped = false;
      modified = false;
      identified = false;
      reversed = false;
      peptide = "";
      AAPeptideSeq = "";
      residueIdxs.clear();
      specIdx = -1;
      numPRMs = 0;
      numSRMs = 0;
      peaks.resize(0);
      massToPeak.clear();
    }

    /**
     * Maps spectrum to peptide and/or proteins by filling in all MappedSpectrum and
     *   MappedPeak class variables
     * @param parentSpecs SpecSet containing spectrum
     * @param idx index of spectrum
     * @param _peptide specnets-style peptide annotation ("" if no annotation)
     * @param proteins sequence of target proteins
     * @param _reverse whether this spectrum should be reversed before annotating
     * @param inputIonTypes needed for Spectrum.annotate
     * @param peakTol Da peak tolerance
     * @return
     */
    void mapProt(SpecSet* parentSpecs,
                 int idx,
                 string _peptide,
                 vector<string>* proteins,
                 bool _reverse,
                 MS2ScoringModel& inputIonTypes,
                 float peakTol);
  };
}

#endif /* MAPPEDSPECTRUM_H_ */
