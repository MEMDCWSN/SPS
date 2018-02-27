/*
 * Contig.h
 *
 *  Created on: Jun 23, 2011
 *      Author: aguthals
 */

#ifndef CONTIG_H_
#define CONTIG_H_

#include <map>

#include "abruijn.h"
#include "aminoacid.h"
#include "spectrum.h"
#include "clusters.h"

#include "PRMSpecEdge.h"

namespace specnets
{
  class Contig : public Spectrum
  {
  public:

    /**
     * Unique index. By default, this is the index of the root assembled spectrum
     */
    int index;

    /**
     * true = this contig is reversed from its native orientation
     */
    bool reversed;

    /**
     * Da distance between left-most assembled spectrum and beginning of contig sequence.
     *   endGaps.first = distance when this contig is in its forward orientation
     *   endGaps.second = distance when this contig is in its reversed orientation
     */
    pair<float, float> endGaps;

    /*
     * Abinfo detailing which peaks were assembled into which abruijn vertices (at key 0)
     */
    abinfo_t* assembledStars;

    /**
     * Alignments used by SPS to sequence consensus
     */
    list<PRMSpecEdge>* innerEdges;

    /**
     * Da distance from root assembled spectrum to every other assembled spectrum
     *   key = index of assembled spectrum
     *   value.first = distance when this contig is in its forward orientation
     *   value.second = distance when this contig is in its reversed orientation
     */
    map<int, pair<float, float> >* rootRef;

    /**
     * All child spectra and their orientations
     *   key = index of assembled spectrum
     *   value.first = child spectrum index in childSpecSet
     *   value.second = orientation of spectrum (true = reversed, false = native)
     */
    map<int, pair<int, bool> >* childSpectra;

    SpecSet* childSpecSet;

    Contig();

    Contig(int idx, SpecSet* contigs, abinfo_t* _assembledStars = 0);

    Contig(const Contig& other);

    ~Contig();

    Contig& operator=(const Contig &other);

    virtual void initialize(int idx,
                            SpecSet* contigs,
                            abinfo_t* _assembledStars = 0);

    virtual void copy(const Contig& other);

    virtual void reverse(void);

    virtual unsigned int assembleConsensus(int minNumMatchedPeaks,
                                           float pkTol,
                                           float pmTol);

    virtual bool merge(PRMSpecEdge* edge,
                       Contig* other,
                       int minNumMatchedPeaks,
                       float pktol,
                       float pmTol);

  protected:
    void create(void);

  };
}

#endif /* CONTIG_H_ */
