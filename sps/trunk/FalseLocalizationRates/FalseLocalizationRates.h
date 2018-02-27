/*
 * FalseLocalizationRates.h
 *
 *  Created on: May 14, 2011
 *      Author: jsnedecor
 */

#ifndef FALSELOCALIZATIONRATES_H_
#define FALSELOCALIZATIONRATES_H_

//system includes
#include <pthread.h>
#include <limits.h>

// Module Includes

// External Includes
#include "spectrum.h"
#include "PeptideSpectrumMatch.h"
#include "SpectrumAnnotStatistics.h"
#include "LinearEquation.h"
#include "aminoacid.h"
#include "Logger.h"
#include "MathUtils.h"
#include "Timer.h"

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_set>
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_set>
#  include <unordered_map>
#endif

//System includes
namespace specnets
{
  class FalseLocalizationRates
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief constructor for PSM class
     *
     */
    FalseLocalizationRates();
    //@}

    FalseLocalizationRates(AAJumps jumps, MS2ScoringModel model, float peakTol);

    //! \name DESTRUCTOR
    //@{
    virtual ~FalseLocalizationRates();
    //@}

    static bool groupHasValidVariant(vector<string> &validSites,
                                     const vector<float> &modifications,
                                     vector<string> &variants,
                                     vector<string> &validVariants);

    bool generateTheoreticalSpectrum(string &unmodifiedPeptide,
                                     Spectrum &unmodifiedSpectrum,
                                     vector<string> &variants,
                                     vector<float> &variantQuantities,
                                     Spectrum &allPossiblePeaks,
                                     Spectrum &outputSpectrum);

    bool generateTheoreticalSpectrum(string &unmodifiedPeptide,
                                     Spectrum &unmodifiedSpectrum,
                                     vector<vector<string> > &variants,
                                     vector<float> &variantQuantities,
                                     Spectrum &allPossiblePeaks,
                                     Spectrum &outputSpectrum);

    static void generateVariantSequences(string &unmodifiedPeptide,
                                         const vector<float> &modifications,
                                         vector<string> &outputVariants);

    static void
        generateVariantSequences(string &unmodifiedPeptide,
                                 const vector<vector<float> > &modifications,
                                 vector<vector<string> > &outputVariants);

    void generateAllIonMasses(const string &peptide,
                              const vector<float> &modifications,
                              vector<float> &outputMasses);

    void generateAllIonMasses(const string &peptide,
                              const vector<vector<float> > &massShifts,
                              vector<float> &outputMasses);

    /* Use AVERAGE distance between groups instead of max distance.*/
    void groupVariants4(Spectrum &modifiedSpectrum,
                        Spectrum &allMassesSpectrum,
                        vector<vector<string> > &variantSeq,
                        const vector<float> &quantities,
                        const vector<float> &massShifts,
                        vector<vector<string> > &variantGroups,
                        vector<float> &variantQuantities,
                        vector<float> &minDistinguishingIntensity,
                        float minimumGroupingPercentIntensity);

    /* Use AVERAGE distance between groups instead of max distance.*/
    void groupVariants4(Spectrum &modifiedSpectrum,
                        Spectrum &allMassesSpectrum,
                        vector<string> &variantSeq,
                        const vector<float> &quantities,
                        const vector<float> &massShifts,
                        vector<vector<string> > &variantGroups,
                        vector<float> &variantQuantities,
                        vector<float> &minDistinguishingIntensity,
                        float minimumGroupingPercentIntensity);

    /* Use maximum distance between groups, hierarchical clustering */
    void groupVariants3(Spectrum &modifiedSpectrum,
                        Spectrum &allMassesSpectrum,
                        vector<string> &variantSeq,
                        const vector<float> &quantities,
                        const vector<float> &massShifts,
                        vector<vector<string> > &variantGroups,
                        vector<float> &variantQuantities,
                        vector<float> &minDistinguishingIntensity,
                        float minimumGroupingPercentIntensity);

    void groupVariants2(Spectrum &modifiedSpectrum,
                        vector<string> &variantSeq,
                        const vector<float> &quanitites,
                        const vector<float> &massShifts,
                        vector<vector<string> > &variantGroups,
                        vector<float> &variantQuantities,
                        vector<float> &minDistinguishingIntensity,
                        float minimumGroupingPercentIntensity);

    void groupVariants(Spectrum &modifiedSpectrum,
                       Spectrum &allMassesSpectrum,
                       vector<string> &variantSeq,
                       const vector<float> &quantities,
                       const vector<float> &massShifts,
                       vector<vector<string> > &variantGroups,
                       vector<float> &variantQuantities,
                       vector<float> &minDistinguishingIntensity,
                       float minimumGroupingPercentIntensity);

    void generateFlrLP(Spectrum &unmodifiedSpectrum,
                       Spectrum &modifiedSpectrum,
                       const string &unmodifiedPeptide,
                       const vector<float> &massShifts,
                       const char * outputCplexFile,
                       const char * outputAllPeaksSpectrum,
                       const char * outputVariantMaps);

    void generateFlrLP(Spectrum &unmodifiedSpectrum,
                       Spectrum &modifiedSpectrum,
                       const string &unmodifiedPeptide,
                       const vector<float> &modifications,
                       const vector<vector<float> > &massShifts,
                       const char * outputCplexFile,
                       const char * outputAllPeaksSpectrum,
                       const char * outputVariantMaps);

    void generateFlrLPPhospho(Spectrum &unmodifiedSpectrum,
                              Spectrum &modifiedSpectrum,
                              const string &unmodifiedPeptide,
                              const vector<float> &modifications,
                              const vector<vector<float> > &massShifts,
                              const char * outputCplexFile,
                              const char * outputAllPeaksSpectrum,
                              const char * outputVariantMaps);

    AAJumps m_jumps; //! Amino acid masses
    MS2ScoringModel m_model; //! ions to include
    float m_peakTol; //! peak tolerance in Da
  protected:
    float getTotalIonCurrent(Spectrum &modifiedSpectrum,
                             Spectrum &allMassesSpectrum);

    float getTotalIonCurrent(PeptideSpectrumMatch &psm1,
                             PeptideSpectrumMatch &psm2);
  };
}
#endif /* FALSELOCALIZATIONRATES_H_ */

