/*
 * MappedContig.h
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#ifndef MAPPEDCONTIG_H_
#define MAPPEDCONTIG_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "abruijn.h"
#include "mzrange.h"
#include "Logger.h"
#include "MappedSpectrum.h"
#include "spectrum.h"
#include "utils.h"

using namespace std;

namespace specnets
{
  /*
   * Abruijn vertex mapped to a protein residue by matchma and/or database search
   */
  class MappedVertex : public MZRange
  {
  public:
    static const char* Labels[];

    // whether this vertex is mapped to a protein residue
    bool mapped;

    // whether this vertex assembles at least one peak mapped to a protein residue
    bool starMapped;

    // whether this vertex assembles at least one peak from an identified spectrum
    bool annotated;

    // b/y annotation
    string annotation;

    // best matching protein residue index
    int residueIdx;

    // whether the vertex is an end point
    bool endpt;

    /**
     * un-annotated (0), incorrect (1), chimeric (2), or correct (3)
     * un-annotated: All derived peaks come from un-identified spectra
     * incorrect: All derived peaks from identified spectra were not b or y ions
     * chimeric: All identified b/y ions were mapped to different
     *   proteins or different residues
     * correct: All identified b/y ions were mapped to the same protein
     *   and the same residue
     */
    short label;

    // number of assembled peaks
    int numPeaks;

    // number of assembled bions
    int numBPeaks;

    // number of assembled yions
    int numYPeaks;

    // number of assembled peaks that came from identified spectra but were not b or y ions
    int numIncorPeaks;

    // assembled mapped/un-mapped peaks
    vector<MappedPeak*> starPeaks;

    MappedVertex() :
      MZRange()
    {
      reset();
    }

    MappedVertex(MZRange& other) :
      MZRange(other)
    {
      reset();
    }

    void reset(void)
    {
      mapped = false;
      starMapped = false;
      annotated = false;
      annotation = "";
      residueIdx = -1;
      endpt = false;
      label = 0;
      numPeaks = 0;
      numBPeaks = 0;
      numYPeaks = 0;
      numIncorPeaks = 0;
      starPeaks.resize(0);
      set(0, 0, 0);
    }
  };

  /*
   * Abruijn vertex difference mapped to a protein residue by matchma and/or database search
   */
  class MappedGap : public MZRange
  {
  public:
    static const char* Labels[];

    // whether both flanking vertices are mapped to a protein residue
    bool vertMapped;

    /** A gap is annotated iff (1) or (2):
     *   (1) Both flanking vertices assemble peaks from the same annotated spectrum
     *   (2) There exists a path of correct gaps following the same ion type such that
     *         the path starts and ends at annotated ions in the same spectrum and this
     *         gap is flanked by one of the starting or ending points.
     */
    bool annotated;

    // whether both flanking vertices each assemble at least one peak mapped to a protein residue
    bool starMapped;

    // if vertMapped, this is the mass of the gap in the protein
    float mappedMass;

    // whether both flanking vertices assemble peaks from the same spectrum
    float sameSpec;

    /** (0) un-annotated, (1) incorrect, or (2) correct
     * un-annotated - not annotated
     * incorrect - annotated but all masses in flanking vertices are in different ion series or
     *   they do not map to contig->starProtIdx
     * correct - annotated and at least one pair of masses in flanking vertices are in the same
     *   ion series and map to contig->starProtIdx. or (2) from above
     */
    short label;

    MappedGap() :
      MZRange()
    {
      reset();
    }

    MappedGap(MZRange& other) :
      MZRange(other)
    {
      reset();
    }

    void reset(void)
    {
      vertMapped = false;
      annotated = false;
      starMapped = false;
      mappedMass = -1.0;
      sameSpec = false;
      label = 0;
      set(0, 0, 0);
    }
  };

  /*
   * Contig mapped to a protein residue by matchma and/or database search
   */
  class MappedContig : public Spectrum
  {
  public:

    // index of contig
    int index;

    // whether abjuijn vertices are mapped to a protein residue
    bool vertMapped;

    // whether at least one assembled star spectrum is mapped to a protein
    bool starMapped;

    // best protein matched to by abjuijn vertices
    int vertProtIdx;

    // protein matched by most mapped assembled spectra
    int starProtIdx;

    // whether this contig has at least one chimeric vertex
    bool chimeric;

    // whether should be reversed when mapped to reference protein
    bool reversed;

    // length of contig
    int length;

    // number of assembled peptides
    int numPeptides;

    // number of assembled spectra
    int numSpecs;

    // index of first mapped AA residue if vertMapped
    int firstResidue;

    // index of last mapped AA residue if vertMapped
    int lastResidue;

    // abruijn vertices
    vector<MappedVertex> abruijnVerts;

    // abruijn gaps
    vector<MappedGap> abruijnGaps;

    // assembled star spectra
    vector<MappedSpectrum*> mappedSpectra;

    MappedContig()
    {
      reset();
    }

    void reset(void)
    {
      index = -1;
      vertMapped = false;
      starMapped = false;
      vertProtIdx = -1;
      starProtIdx = -1;
      chimeric = false;
      reversed = false;
      length = 0;
      numPeptides = 0;
      numSpecs = 0;
      firstResidue = -1;
      lastResidue = -1;
      abruijnVerts.resize(0);
      abruijnGaps.resize(0);
      mappedSpectra.resize(0);
    }

    /**
     * Annotates this mapped contig by setting all class variables
     * @param parentContigs SpecSet containing contig
     * @param idx index of contig
     * @param parentComponents ab_info of contigs
     * @param overlaps for each mapped contig (spectrum idx), this matches
     *   contig peak indicies (mass) to protein residue indicies (intensity).
     * @param protIdx for each contig, this details the mapped protein (0) and
     *   whether the contig is reversed (2)
     * @param mappedSpecs mapped star spectra the contigs assemble
     * @param protein_spectra all target proteins as cummulative masses
     * @param peakTol Da peak tolerance
     * @return
     */
    void mapProt(SpecSet* parentContigs,
                 int idx,
                 abinfo_t* parentComponents,
                 SpecSet* overlaps,
                 vector<vector<int> >* protIdx,
                 vector<MappedSpectrum>* mappedSpecs,
                 SpecSet* protein_spectra,
                 float peakTol);

    /**
     * @return number of assembled annotated spectra
     */
    int getNumAnnotSpec(void);

    /**
     * @return number of assembled mapped spectra
     */
    int getNumMappedSpec(void);

    /**
     * @return % of annotated vertices from b ions
     */
    float getPercBVerts(void);

    /**
     * @return % of annotated vertices from y ions
     */
    float getPercYVerts(void);

    /**
     * @return % of annotated vertices from b or y ions
     */
    float getPercBYVerts(void);

    /**
     * @return % of vertices that match label. if _label == 0, denominator is
     *   # of all vertices. if _label > 0, denominator is # of annotated vertices
     */
    float getPercVerts(short _label);

    /**
     * @return % of calls that match accuracy label. if _label == 0, denominator is
     *   # of all gaps. if _label > 0, denominator is # of annotated gaps
     */
    pair<float, float> getPercCallsAcc(short _label);

    /**
     * @return % of regions between mapped vertices that match length label:
     *   0 : region spans one AA in the database
     *   1 : region spans multiple AA in the database with no un-mapped vertices in between
     *   2 : region spans multiple AA in the database with un-mapped vertices in between
     */
    pair<float, float> getPercCallsLen(short _label);

  };
}

#endif /* MAPPEDCONTIG_H_ */
