#ifndef SPECTRAL_PAIRS_H
#define SPECTRAL_PAIRS_H
// Module Includes
#include "SpectrumPairSet.h"

// External Module Includes
#include "spectrum.h"
#include "SpecSet.h"

// System Includes
#include <fstream>

namespace specnets
{
  /**
   *@deprecated use SplitSpectra/SplitPairs instead
   */
  void SplitPairs(SpecSet &specSet,
                  SpectrumPairSet &aligns,
                  float peakTol,
                  int maxAAjump,
                  float penalty_sameVert,
                  float penalty_ptm,
                  bool forceSymmetry,
                  SpecSet &specSetNew,
                  SpectrumPairSet &alignsNew,
                  vector<vector<TwoValues<int> > > &matches,
                  vector<bool> &pairFlipped,
                  vector<vector<float> > *dbg_matchScores = 0,
                  ofstream *debug = 0);

  /**
   * Decides on a consensus orientation for  a set of edges. Flips spectra
   * according to assigned orientation.
   *
   *@param specSet - star spectra
   *@param aligns  - pairwise alignments
   *@param peakTol - peak mass tolerance
   *@param pmTol   - parent mass tolerance
   *@param maxAAjump - maximum amino-acid mass jump. Consecutive matched peaks with mass
   *                     difference in this range must match an amino acid mass.
   *@param penalty_sameVert
   *@param penalty_ptm
   *@param matches     - position [i] contains set of peak indices for matched peaks in aligns[i].
   *@param specFlipped - indicates whether spectra were reversed when added to the ABruijn graph
   *@param modPos      - location of the mass difference when the spectrum was added to the ABruijn graph
   *@param minMatchedPeaks - Minimum number of matched peaks to keep matches (ABruijn glues) between two spectra
   *@param minEdgesToComponent - Minimum number of edges to spectra in the component for additional spectra to be added
   *@param forceSymmetry
   *@param alignStats - position [i] %matched score in spec1/spec2 (pos 0/1) and
   #matched peaks (pos 2) for aligns[i]
   *@param labelsP
   *@return
   */
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
                  unsigned int minMatchedPeaks = 0,
                  unsigned int minEdgesToComponent = 0,
                  bool forceSymmetry = false,
                  bool ignoreReversals = false,
                  vector<vector<float> > *alignStats = NULL,
                  vector<SpectrumPeakLabels> *labelsP = NULL);

  /**
   * Decides on a consensus orientation for  a set of edges. Flips spectra
   * according to assigned orientation.
   *
   *@param specSet
   *@param aligns
   *@param peakTol
   *@param maxAAjump
   *@param penalty_sameVert
   *@param penalty_ptm
   *@param matches
   *@param specFlipped
   *@param modPos
   *@param forceSymmetry
   *@param labelsP
   */
  void SplitPairs2(SpecSet &specSet,
                   SpectrumPairSet &aligns,
                   float peakTol,
                   int maxAAjump,
                   float penalty_sameVert,
                   float penalty_ptm,
                   vector<vector<TwoValues<int> > > &matches,
                   vector<bool> &specFlipped,
                   vector<float> &modPos,
                   bool forceSymmetry = false,
                   vector<SpectrumPeakLabels> *labelsP = NULL);

  /**
   * Like SplitPairs2 but also processes Partial Overlap alignments (alignsPA).
   *
   *@param specSet
   *@param aligns
   *@pa ram alignsPA*@param peakTol
   *@param maxAAjump
   *@param penalty_sameVert
   *@param penalty_ptm
   *@param matches
   *@param matchesPA
   *@param specFlipped
   *@param modPos
   *@param forceSymmetry
   *@param labelsP
   *@param alignRatios
   *@param alignRatiosPA
   */
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
                   bool forceSymmetry = false,
                   vector<SpectrumPeakLabels> *labelsP = NULL,
                   vector<TwoValues<float> > *alignRatios = NULL,
                   vector<TwoValues<float> > *alignRatiosPA = NULL);
  /**
   * TODO: add description
   *
   *@param specSet
   *@param specSetSplit
   */
  void SplitSpectra(SpecSet &specSet, SpecSet &specSetSplit);

  /**
   * TODO: add description
   *
   *@param specSetSplit
   *@param aligns
   *@param peakTol
   *@param maxAAjump
   *@param penalty_sameVert
   *@param penalty_ptm
   *@param forceSymmetry
   *@param alignsNew
   *@param matches
   *@param pairFlipped
   *@param dbg_matchScores
   *@param debug
   */
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
                   vector<vector<float> > *dbg_matchScores = 0,
                   ofstream *debug = 0);

  /**
   * TODO: add description
   *
   *@param labels
   *@param newLabels
   */
  void SplitLabels(vector<SpectrumPeakLabels> &labels, vector<
      SpectrumPeakLabels> &newLabels);

  /**
   * TODO: add description
   *
   *@param specSet
   *@param consen sus
   *@param peakTol
   *@param resolution
   */
  void ComputeSpectralStars(SpecSet &specSet,
                            Spectrum &consensus,
                            float peakTol,
                            float resolution);

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
                       unsigned int minNumMatchedPeaks);

  /*
   * projectSpectrumOld - "projects" all specSet spectra with indices in specsToProcess onto specBase:
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
   *@param idxMatched Indices of matched peaks in specBase/specsToProcess.front() (only when specsToProcess.size()==1)
   */
  void ProjectSpectrumOld(SpecSet &specSet,
                       const Spectrum &specBase,
                       list<int> &specsToProcess,
                       list<int> &curDeltas,
                       float &bestScore,
                       list<int> &bestDeltas,
                       float peakTol,
                       Spectrum *finalProj = NULL,
                       vector<TwoValues<int> > *idxMatched = NULL,
                       unsigned int minNumMatchedPeaks = 0);

} // namespace specnets

#endif
