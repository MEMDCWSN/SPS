/*
 * prm_alignment.h
 *
 *  Created on: Feb 17, 2011
 *      Author: aguthals
 */

#ifndef PRM_ALIGNMENT_H_
#define PRM_ALIGNMENT_H_

#include "aminoacid.h"
#include "mzrange.h"
#include "spectrum.h"
#include "twovalues.h"
#include "SpectrumPair.h"
#include "utils.h"

namespace specnets
{

  class PRMAlignment: public SpectrumPair
  {
  public:

    /**
     * Computes the alignment score of 2 PRM spectra
     * @param matchedPeaks number of peaks matching in each spectrum
     * @param matchedIntensity summed intensity of all matching peaks in spectrum
     * @param overlapIntensity summed intensity all peaks in spectrum that overlap
     *   with second spectrum in pair
     * @return alignment score of spectrum
     */
    static float overlapScore(int matchedPeaks,
                              float matchedIntensity,
                              float overlapIntensity)
    {
      return ((float) matchedPeaks) * (matchedIntensity / overlapIntensity);
    }

    /**
     * Computes the alignment score of 2 PRM spectra
     * @param matchedIntensity summed intensity of all matching peaks in spectrum
     * @param overlapIntensity summed intensity all peaks in spectrum
     * @return alignment score of spectrum
     */
    static float totalScore(float matchedIntensity, float totalIntensity)
    {
      return matchedIntensity / totalIntensity;
    }

    // maps shifts to list of overlapping shifts and matched peaks
    map<MZRange, pair<list<MZRange>, list<pair<int, int> > > >* shift_mp;

    // all shifts with minimum number of matching peaks
    set<MZRange>* bestShifts;

    /**
     * Default constructor
     */
    PRMAlignment(void) :
      SpectrumPair()
    {
      shift_mp = new map<MZRange, pair<list<MZRange>, list<pair<int, int> > > > ;
      bestShifts = new set<MZRange> ;
      spec2R = new Spectrum;
    }

    /**
     * Constructor
     * @param _spec1
     * @param _spec2
     * @param _peak_tol peak tolerance
     * @param _tol_type how to compute mass tolerance
     *   0: _peak_tol taken as Da tolerance
     *   1: _peak_tol taken as PPM tolerance
     */
    PRMAlignment(Spectrum* _spec1,
                 Spectrum* _spec2) :
      SpectrumPair()
    {
      shift_mp = new map<MZRange, pair<list<MZRange>, list<pair<int, int> > > > ;
      bestShifts = new set<MZRange> ;
      spec2R = new Spectrum;
      setSpec1(_spec1);
      setSpec2(_spec2);
    }

    /**
     * De-constructor
     */
    ~PRMAlignment(void)
    {
      delete shift_mp;
      delete bestShifts;
      delete spec2R;
    }

    /**
     * Set first spectrum in alignment
     * @param setSpec1
     */
    void setSpec1(Spectrum* _spec1)
    {
      spectrum1 = _spec1;
      szSpec1 = spectrum1->size();
      cleanup();
    }

    /**
     * Set second spectrum in alignment
     * @param setSpec2
     */
    void setSpec2(Spectrum* _spec2)
    {
      spec2F = _spec2;
      spectrum2 = spec2F;

      szSpec2 = spectrum2->size();

      spec2F->reverse(0.0 - AAJumps::massH2O, spec2R);
      cleanup();
    }

    void reverseSpec2(bool rev)
    {
      if (rev) {
        spectrum2 = spec2R;
      }
      else {
        spectrum2 = spec2F;
      }
    }

    /**
     * Computes the highest-scoring shift of spectrum2 in relation to spectrum1
     *   and reversed(spectrum2) in relation to spectrum1. This alignment is
     *   stored in the SpectrumPair class variables
     * @param minMatchedPeaks minimum number of matching peaks allowed
     * @param minScore minimum alignment score allowed
     * @param score_type see method computShiftScore
     * @return true if alignment was found, false if none met criteria
     */
    bool computeShiftNoModFR(int minMatchedPeaks,
                             float minScore,
                             short score_type);

    /**
     * Computes the highest-scoring shift of spectrum2 in relation to spectrum1
     * @param minMatchedPeaks minimum number of matching peaks allowed
     * @param minScore minimum alignment score allowed
     * @param score_type see method computShiftScore
     * @param alignment output alignment
     * @return true if alignment was found, false if none met criteria
     */
    bool getShiftNoMod(int minMatchedPeaks,
                       float minScore,
                       short score_type,
                       SpectrumPair& alignment);

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
    pair<int, pair<float, float> >
    getShiftScore(float shift, float shift_tol, short score_type, pair<vector<
        float> , vector<float> >* alignedScores = 0) const;

    /**
     * Computes all possible shifts of spectrum2 in relation to spectrum1 and
     *   stores them in shift_mp and bestShifts
     * @param minMatchedPeaks only consider shifts with minimum number of
     *   matching peaks. Also, put any shifts with minMatchedPeaks in
     *   bestShifts
     * @return
     */
    void computeShifts(int minMatchedPeaks);

  private:

    /**
     * Called internally to clear all shifts
     */
    void cleanup(void);

    Spectrum* spectrum1;
    Spectrum* spectrum2;

    Spectrum* spec2F;
    Spectrum* spec2R;

    int szSpec1;
    int szSpec2;
  };

}

#endif /* PRM_ALIGNMENT_H_ */
