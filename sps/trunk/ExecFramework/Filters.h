#ifndef _Filters_H_
#define _Filters_H_

// Module Includes
#include "SpectrumPairSet.h"

// External Includes
#include "spectrum.h"
#include "SpecSet.h"
//#include "batch.h"

// System Includes
#include <cmath>
#include <iostream>
#include <vector>

namespace specnets
{
  /**
   * TODO: add description
   *
   *@param edge1
   *@param edge2
   *@param edge3
   *@param pmTol
   *@return
   */
  bool validateTriangle(SpectrumPair &edge1,
                        SpectrumPair &edge2,
                        SpectrumPair &edge3,
                        float pmTol);

  /**
   * TODO: add description
   *
   *@param aligns
   *@param idxStart
   *@param idxEnd
   *@param pmTol
   *@param selectedIdx
   *@return
   */
  /*
   * FilterTriangles - Retains only edges forming a triangle (see ValidateTriangle
   *   for an illustration of a triangle).
   *
   * Note: Assumes that:
   *         - aligns is sorted by (.spec1,.spec2)
   *         - contains no repeated edges
   *         - for every aligns[i], aligns[i].spec1 < aligns[i].spec2
   */
  unsigned int filterTriangles(SpectrumPairSet &aligns,
                               unsigned int idxStart,
                               unsigned int idxEnd,
                               float pmTol,
                               std::vector<unsigned int> &selectedIdx);

  /**
   * TODO: add description
   *
   *@param aligns
   *@param idxKept
   *@param pvalues
   *@param gaussianParams
   *@param ratios
   *@param minPValue
   *@param minRatio
   *@param pmTol
   *@param filterTrigs
   *@return
   */
  unsigned int filterAligns(SpectrumPairSet &aligns,
                            std::vector<unsigned int> &idxKept,
                            std::vector<TwoValues<float> > &pvalues,
                            std::vector<TwoValues<float> > &means,
                            std::vector<float> &variances,
                            std::vector<TwoValues<float> > &ratios,
                            float minPValue,
                            float minRatio,
                            float pmTol,
                            bool filterTrigs);
  /**
   * TODO: add description
   *
   *@param specSet
   *@param results
   *@param ratioType
   *@param ratios
   */
  void filterMatchRatioASP(SpecSet &specSet,
                           SpectrumPairSet &results,
                           short ratioType,
                           std::vector<TwoValues<float> > &ratios);

  /**
   * TODO: add description
   *
   *@param specSet
   *@param results
   *@param scoreMeans is mean (col.0) and sample size (col.1)
   */
  void computeScoreMeans(SpecSet &specSet,
                         SpectrumPairSet &results,
                         std::vector<TwoValues<float> > &scoreMeans);

} //namespace specnets

#endif
