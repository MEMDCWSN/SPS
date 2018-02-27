/*
 * MappedSpecnets.h
 *
 *  Created on: Mar 2, 2011
 *      Author: aguthals
 */

#ifndef MAPPEDSPECNETS_H_
#define MAPPEDSPECNETS_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "MappedContig.h"
#include "Logger.h"
#include "mzrange.h"
#include "utils.h"

namespace specnets
{
  class MappedSpecnets
  {
  public:

    // mapped contigs
    vector<MappedContig>* contigs;

    // mapped star spectra
    vector<MappedSpectrum>* spectra;

    // target proteins
    vector<string>* proteins;

    // identified peptides
    vector<string>* peptides;

    set<int> target_prots;

    MappedSpecnets(void);

    ~MappedSpecnets(void);

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
    void mapProt(SpecSet* sps_contigs,
                 abinfo_t* sps_components,
                 SpecSet* star_spectra,
                 SpecSet* matchma_overlaps,
                 vector<vector<int> >* matchma_prot_idx,
                 vector<string>* _proteins,
                 SpecSet* protein_spectra,
                 vector<string>* spectrum_ids,
                 MS2ScoringModel& input_ion_types,
                 float peak_tol);

    /**
     * @return number of contigs w/ > 0 peaks
     */
    int getNumContigs(void);

    /**
     * @return number of spectra w/ > 0 peaks
     */
    int getNumSpectra(void);

    /**
     * @return number of identified spectra
     */
    int getNumSpecIdent(void);

    /**
     * @return number of spectra mapped to target proteins
     */
    int getNumSpecMapped(int prot_idx);

    /**
     * @return number of modified spectra
     */
    int getNumSpecModified(void);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return number of contigs mapped to prot_idx by matchma or spectrum IDs
     */
    int getNumContigsMapped(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return number of contigs mapped to prot_idx by matchma
     */
    int getNumContigsVertMapped(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return number of contigs mapped to prot_idx by spectrum IDs
     */
    int getNumContigsSpecMapped(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return number of contigs mapped to prot_idx by matchma and spectrum IDs
     */
    int getNumContigsVertSpecMapped(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return percent of protein covered by identified spectra
     */
    float getPercSpecCov(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return percent of protein covered by matchma mapped contigs
     */
    float getPercSeqCov(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return average number of contigs covering each covered residue
     */
    float getCovRedundancy(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return number of spectra assembled into contigs mapping to protein by matchma or spectrum IDs
     */
    int getNumAssembledSpec(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return average number of spectra per contig mapping to protein by matchma or spectrum IDs
     */
    float getSpecPerContig(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return average number of peptides per contig mapping to protein by matchma or spectrum IDs
     */
    float getPepPerContig(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return average Da of each contig mapping to protein by matchma or spectrum IDs
     */
    float getDaLengthPerContig(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return average AA of each contig mapping to protein by matchma
     */
    float getAALengthPerContig(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return pair of longest contig (second) mapping to protein by matchma and its length (first)
     */
    pair<int, int> getLongestAAContig(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @param label 0, 1, 2, or 3
     * @return if label != 0, percent of abruijn vertices that have label out of
     *  those that are annotated. if label == 0, percent of abruijn vertices
     *  that have label out of those that are not annotated
     */
    float getPercVerts(int prot_idx, short label);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return pair of percent b ions (first) and percent y ions (second) out
     *   of all annotated vertices
     */
    pair<float, float> getPercBYVerts(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @param label 0, 1, or 2
     * @return if label != 0, percent of abruijn gaps that have accuracy label out of
     *  those that are annotated. if label == 0, percent of abruijn gaps
     *  that are annotated out of all
     */
    float getPercCallsAcc(int prot_idx, short label);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @param label 0, 1, or 2
     * @return % of regions between mapped vertices that match length label:
     *   0 : region spans one AA in the database
     *   1 : region spans multiple AA in the database with no un-mapped vertices in between
     *   2 : region spans multiple AA in the database with un-mapped vertices in between
     */
    float getPercCallsLen(int prot_idx, short label);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return if prot_idx < 0, percent of all assembled spectra that are
     *   identified. Otherwise, the percent of all assembled spectra in a contig
     *   mapped to protein prot_idx that are identified.
     */
    float getPercAssemSpecIdent(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return if prot_idx < 0, percent of all assembled identified spectra that are
     *   identified for contigs mapping to prot_idx (-1 for all proteins)
     */
    float getPercAssemIdentSpecMod(int prot_idx);

    /**
     * @param prot_idx index of target protein or -1 to specify all proteins
     * @return of all contigs mapped to prot_idx, the percent of all matchma
     *   and inspect mapped abruijn vertices that assemble at least one
     *   identified MS/MS b/y peak mapping to the same matchma residue index.
     */
    float getMatchmaAccuracy(int prot_idx);

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
    void getGapAccVsPos(int prot_idx,
                        vector<pair<int, int> >& outValues,
                        vector<vector<pair<string, bool> > >* outputTable);
  };
}

#endif /* MAPPEDSPECNETS_H_ */
