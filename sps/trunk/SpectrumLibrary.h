/*
 * SpectrumLibrary.h
 *
 *  Created on: Apr 21, 2011
 *      Author: jsnedecor
 */

#ifndef SPECTRUMLIBRARY_H_
#define SPECTRUMLIBRARY_H_

// Module Includes
#include "PeptideSpectrumMatch.h"
#include "PeptideSpectrumMatchSet.h"
#include "spectrum.h"
#include "Logger.h"
#include "cosine_similarity.h"

// System Includes
#include <vector>

using namespace std;

namespace specnets
{
  /*! \brief Class to build spectrum library

   This class takes a single peptide spectrum match set and
   spectrum set.
   */
  class SpectrumLibrary
  {
  public:
    /*! \brief Return associated  SpectrumLibraryMatch for
     that vector position.
     */
    psmPtr & operator[](unsigned int i);

    /*! \brief Return associated  SpectrumLibraryMatch
     for that vector position.
     */
    const psmPtr & operator[](unsigned int i) const;

    /*! \brief Set value of one SpectrumLibraryMatchSet to another

     */
    SpectrumLibrary & operator=(SpectrumLibrary &other);

    /*! \brief Returns size of m_library vector

     */
    unsigned int size();

    /*! \brief Resizes m_library vector

     @param newSize the new size of the parameter vector
     */
    unsigned int resize(unsigned int newSize);

    /*! \brief adds new values to peptide spectrum library.

     @param spectrum candidate spectrum for current psm
     @param psm candidate spectrum psm
     @param score used to compare to existing spectrum library.
     @param scoreAscending indicates whether we want a higher or lower score.
     True = p-value (higher values worse), false = MQScore (higher values better)
     @return true if this spectrum has been added to the library, false if not.
     */
    bool addToLibrary(Spectrum &spectrum,
                      PeptideSpectrumMatch &psm,
                      float score,
                      bool scoreAscending = 0);

    /*! \brief adds new values to peptide spectrum library with key as parameter.

     @param spectrum candidate spectrum for current psm
     @param psm candidate spectrum psm
     @param score used to compare to existing spectrum library.
     @param key indicates the key we wish to store this peptide under
     @param scoreAscending indicates whether we want a higher or lower score.
     True = p-value (higher values worse), false = MQScore (higher values better)
     @return true if this spectrum has been added to the library, false if not.
     */
    bool addToLibrary(Spectrum &spectrum,
                      PeptideSpectrumMatch &psm,
                      float score,
                      string &key,
                      bool scoreAscending = 0);

    /*! \brief adds new values to peptide spectrum library regardless of score

     @param spectrum candidate spectrum for current psm
     @param psm candidate spectrum psm
     */
    bool addToLibrary(Spectrum &spectrum, PeptideSpectrumMatch &psm);

    /*! \brief adds new values to peptide spectrum library regardless of score
     * with key as parameter.

     @param spectrum candidate spectrum for current psm
     @param psm candidate spectrum psm
     @param key indicates the key we wish to store this peptide under
     */
    bool addToLibrary(Spectrum &spectrum,
                      PeptideSpectrumMatch &psm,
                      string &key);

    /*! Takes best match of inputPsm to psmsToMatch and returns it into outputPsm
         * If there is only psm, returns that psm.
         *
         * @param inputPsm - Input psm to compare to
         * @param psmsToMatch - Input psms to compare to inputPsm
         * @param iontToExtract - list of ions to use in cosine
         * @param model - AA model to use for cosine
         * @param outputPsm - output psm for best match to input psm.
         */
    static bool pickBestHit(psmPtr &inputPsm,
                            vector<psmPtr> &psmsToMatch,
                            vector<string> &ionsToExtract,
                            MS2ScoringModel &model,
                            psmPtr &outputPsm);


    /*! \brief returns spectrum library match for annotation and charge.
     * If there is more than one psm associated with key, returns first psm.
     *
     * Returns matching spectrum if current annotation and charge is defined.
     * Returns false if annotation and charge do not have candidate spectrum.
     * @param annotation - Input peptide annotation
     * @param charge - Input charge specification.
     */
    bool getSpectrumLibraryMatch(string &annotation,
                                 int charge,
                                 psmPtr &outputPsm);

    /*! brief returns all spectrum library matches for annotation and charge.
     *
     * Returns matching spectrum if current annotation and charge is defined.
     * Returns false if annotation and charge do not have candidate spectrum.
     * @param annotation - Input peptide annotation
     * @param charge - Input charge specification.
     */
    bool getSpectrumLibraryMatches(string &annotation, int charge, vector<
        psmPtr> &outputPsms);

    PeptideSpectrumMatchSet m_library; //! vector to contain top match spectra for library.
  protected:
    /*! \brief build m_map from m_library
     *
     */
    void buildKeyMap(void);

    /*! \brief helper function for addToLibrary
     *
     */
    void addSpectrumMatch(Spectrum &spectrum,
                          PeptideSpectrumMatch &psm,
                          string &key,
                          float score);

    vector<float> m_libraryScores; //! vector of score associated with m_library
    std::tr1::unordered_multimap<string, unsigned int> m_map; //! associates annotation_charge with the index in m_library.
  };
}

#endif /* SPECTRUMLIBRARY_H_ */
