/*
 * PeptideSpectrumMatchNetwork.h
 *
 *  Created on: Sep 2011
 *      Author: cboucher
 */

#ifndef PEPTIDESPECTRUMMATCHNETWORK_H_
#define PEPTIDESPECTRUMMATCHNETWORK_H_

//Module includes
#include "utils.h"
#include "spectrum.h"
#include "spectrum_scoring.h"
#include "aminoacid.h"
#include "PeptideSpectrumMatch.h"

//System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <tr1/unordered_map>

namespace specnets
{
  /**
   * @see spectrum.h
   */
  class Spectrum;

  /*! \brief Peptide spectrum match class. A single spectrum is matched to a single annotation
   *
   * It is possible to have more than one annotation per spectrum, however, they will show up
   * as different PSMs.
   */
  class PeptideSpectrumMatchNetwork : public PeptideSpectrumMatch
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief constructor for PSM class
     *
     */
    PeptideSpectrumMatchNetwork();
    //@}

    //! \name DESTRUCTOR
    //@{
    ~PeptideSpectrumMatchNetwork();
    //@}

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    PeptideSpectrumMatchNetwork(const PeptideSpectrumMatchNetwork &other);

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    PeptideSpectrumMatchNetwork &operator=(const PeptideSpectrumMatchNetwork &other);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    void
        setAnnotationToMatches(vector<int> &matches,
                               vector<pair<const ftIonFragment*, short> > &annotation,
                               int ionIdx,
                               const ftIonFragment* currIonFrag);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector, taking only the higher intensity match in case of duplicates
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    void
        setAnnotationToMatchesNoDups(vector<int> &matches,
                                     vector<pair<const ftIonFragment*, short> > &annotation,
                                     int ionIdx,
                                     const ftIonFragment* currIonFrag);

    /**
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     * @param jumps - AAJumps indicating amino acid masses
     */
    bool annotate(string &peptide,
                  string &ionNamesInclude,
                  MS2ScoringModel &inputIonTypes,
                  float prmOffset,
                  float srmOffset,
                  float peakTol,
                  AAJumps &jumps,
                  bool removeDuplicates = true,
                  bool ignoreParentCharge = false);

    /**
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     */
    bool annotate(string &peptide,
                  string &ionNamesInclude,
                  MS2ScoringModel &inputIonTypes,
                  float prmOffset,
                  float srmOffset,
                  float peakTol,
                  bool removeDuplicates = true,
                  bool ignoreParentCharge = false);

    /** Maps peakList indices to ion names.
     *
     *
     * @param outputMatches output for matched indices
     * @param ionMap mapping between ion name + ion index keyed to
     * index of outputMatches
     */
    void mapIons(vector<vector<int> > &outputMatches, std::tr1::unordered_map<
        string, int> &ionMap);

    vector<int> m_parentMasses;
    vector<int> m_spectrumIndices;
    vector<string> m_partialAnnotations;
    vector<int> m_partialPeptideScores;
    vector<double> m_partialSpecProbabilities;
    int m_individualScore;
    bool m_isTuple;
    vector< Spectrum *> m_spectra; //! associated spectrum match
  };
}
#endif /* PEPTIDESPECTRUMMATCHNETWORK_H_ */
