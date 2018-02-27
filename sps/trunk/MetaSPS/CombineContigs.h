/*
 * CombineContigs.h
 *
 *  Created on: Jun 11, 2012
 *      Author: aguthals
 */

#ifndef COMBINECONTIGS_H_
#define COMBINECONTIGS_H_

#include <ctype.h>
#include <cstdio>
#include <cmath>
#include <list>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <set>
#include <map>
#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>

#include "abruijn.h"
#include "ion.h"
#include "label.h"
#include "range.h"
#include "spectrum.h"
#include "aminoacid.h"
#include "twovalues.h"
#include "db_fasta.h"
#include "utils.h"
#include "Logger.h"
#include "SpectrumPairSet.h"
#include "prm_alignment.h"

using namespace std;

namespace specnets
{

  struct CombineContigsParams
  {
    //required contig-contig merging parameters
    SpecSet* contigs;
    SpecSet* star_spectra;
    SpecSet consensus;
    abinfo_t* in_abinfo;
    abinfo_t out_abinfo;
    vector<int> out_reversed;
    list<SpectrumPair>* contig_alignments;

    float peak_tol;
    float parent_mass_tol;
    float resolution;
    float contig_overlap_score;
    int min_matched_peaks_contig_align;

    //required contig-matepair merging parameters
    //SpecSet* matepairs;
    //list<SpectrumPair>* matepair_alignments;

    //short min_matepair_charge;
    //float min_ratio_matepair_align;
   // float min_ratio_matepair_contig_align;
    //int min_matched_peaks_matepair_align;
    //int min_num_matepairs;
    //int min_edges_to_component;
    int min_component_size;

    //required matchma parameters
    vector<vector<int> >* contig_prot_idx;
    SpecSet* contig_peak_idx;
    DB_fasta* proteins;
    set<int>* protein_idx_ignore;
  };

  void copy(CombineContigsParams& from, CombineContigsParams& to);

  class CombineContigs
  {
  public:

    CombineContigsParams* merging_params;

    CombineContigs();
    ~CombineContigs();

    CombineContigs(CombineContigsParams* params);

    void construct(CombineContigsParams* params);

    // mergeType = 0 (combine w/ contig shifts), 1 (combine w/ connector shifts)
    void combineEdges(short mergeType);

    bool saveComponents(const char* outcomponents);

    bool getContigShiftStats(const char* outfile,
                             int startMinNumMatchedPeaks = 6,
                             int endMinNumMatchedPeaks = 9,
                             int stepMinNumMatchedPeaks = 1,
                             float startMinRatio = 0.2,
                             float endMinRatio = 0.8,
                             float stepMinRatio = 0.01,
                             float startScore = 0,
                             float endScore = 6,
                             float stepScore = 0.1);
    /*
     void getConnectorShiftStats(const char* outfile,
     int precCharge,
     int startMinNumMatchedPeaks = 6,
     int endMinNumMatchedPeaks = 10,
     int stepMinNumMatchedPeaks = 1,
     float startMinRatio = 0.2,
     float endMinRatio = 0.6,
     float stepMinRatio = 0.1);
     */
  private:
    SpecSet connectors;
    SpecSet contigs;
    abinfo_t contig_abinfo;
    SpecSet oriented;
    SpecSet consensus;
    abinfo_t consensus_abinfo;
    abinfo_t parent_abinfo;

    map<int, map<int, int> > graph;
    vector<vector<float> > edges;
    vector<set<int> > components;
    vector<list<TwoValues<int> > > root_alignments;
    vector<bool> reversed;

    vector<TwoValues<float> > extraShifts;
    list<SpectrumPair> contig_alignments;
    list<SpectrumPair> connector_alignments;
    SpecSet overlaps;
    DB_fasta fasta;
    vector<vector<int> > prot_match;
    set<int> idxNotUsed;

    bool haveRes;
    bool haveMatePairs;
    float peakTol;
    float parentMassTol;
    int intPeakTol;
    int intParentMassTol;
    vector<int> edgeNum;
    int numCorrectPairs;
    map<int, set<int> > correctPairs;
    map<int, set<int> > seenBef;
    map<int, set<int> > possiblePairs;
    vector<float> edge;
    TwoValues<float> scoredEdge;
    map<int, int> edgeMap;

    void init();
    /*
     void
     getAlignmentGraph(map<unsigned int, map<unsigned int,
     list<unsigned int> > >* graphCont,
     map<unsigned int, map<unsigned int, SpectrumPair> >* connContRes,
     short precCharge,
     unsigned int minNumMatchedPeaks,
     float minRatio,
     float minContigOverlapIntensity);
     */

    void outputGraph(ostream& output,
                     map<int, map<int, int> >& _graph,
                     vector<vector<float> >& _edges,
                     set<int>* erasedEdges = 0,
                     vector<set<int> >* _components = 0);

    float getComponentOverlapArea(int idx1,
                                  int idx2,
                                  float shift_use,
                                  bool rev,
                                  map<int, map<int, int> >& _graph,
                                  vector<vector<float> >& _edges,
                                  vector<set<int> >& _components);

    void outputComponents(FILE* output,
                          map<int, map<int, int> >& _graph,
                          vector<vector<float> >& _edges,
                          vector<set<int> >& _components,
                          SpecSet& orientedContigs,
                          vector<bool>& _reversed,
                          bool outputStats = false,
                          bool printComps = true,
                          set<int>* countNodes = 0,
                          bool justContigs = false);

    void recomputeEdges(int idx,
                        map<int, map<int, int> >& _graph,
                        vector<vector<float> >& _edges,
                        vector<set<int> >& _components,
                        set<int>& erasedEdges,
                        vector<int>& edgeRef,
                        vector<list<TwoValues<int> > >& _root_alignments,
                        list<TwoValues<float> >& scoredEdges);

    TwoValues<float> getOverlapScore(int idx1,
                                     int idx2,
                                     float fShift,
                                     float rShift,
                                     map<int, map<int, int> >& _graph,
                                     vector<vector<float> >& _edges,
                                     vector<set<int> >& _components);

    TwoValues<float> getConsensusShifts(int idx1,
                                        int idx2,
                                        float fShift,
                                        float rShift,
                                        map<int, map<int, int> >& _graph,
                                        vector<vector<float> >& _edges,
                                        vector<set<int> >& _components,
                                        bool debug = false);

    TwoValues<float>
    getConsensusOverlap(int idx1,
                        int idx2,
                        float fShift,
                        float rShift,
                        map<int, map<int, int> >& _graph,
                        vector<vector<float> >& _edges,
                        vector<set<int> >& _components,
                        vector<list<TwoValues<int> > >& _root_alignments,
                        bool debug = false);

    bool tryMergeContigs(int idx1,
                         int idx2,
                         int edgeIdx,
                         map<int, map<int, int> >& _graph,
                         vector<vector<float> >& _edges,
                         vector<set<int> >& _components,
                         vector<list<TwoValues<int> > >& _root_alignments,
                         bool rev,
                         bool debug);

    TwoValues<float>
    mergeContigs(int idx, map<int, map<int, int> >& _graph, vector<
        vector<float> >& _edges, vector<set<int> >& _components, vector<list<
        TwoValues<int> > >& _root_alignments, Spectrum& putSpec, pair<pair<
        vector<int> , vector<int> > , vector<
        pair<vector<int> , vector<double> > > >& putAb, bool debug);

    void reverseNode(int idx,
                     map<int, map<int, int> >& _graph,
                     vector<vector<float> >& _edges,
                     vector<set<int> >& combinations,
                     vector<bool>& _reversed);

    void getCondensedContigAlignmentGraph(map<int, map<int, int> >& _graph,
                                          vector<vector<float> >& _edges,
                                          list<TwoValues<float> >& scoredEdges,
                                          float minCombScore,
                                          int minMatchedPeaks);
    /*
     void
     getCondensedAlignmentGraph(vector<set<int> >& _components,
     map<int, map<int, int> >& _graph,
     vector<vector<float> >& _edges,
     vector<bool>& _reversed,
     list<TwoValues<float> >& scoredEdges,
     vector<list<TwoValues<int> > >& _root_alignments,
     short precCharge,
     float minRatioConnContig,
     float minRatio,
     int minNumMatchedPeaks,
     int minNumConnectors,
     int minEdgesToComponent);
     */
    bool validShift(float FFShift, int idx1, int idx2, bool isReversed = false);

    bool
    validShiftMod(float FFShift, int idx1, int idx2, bool isReversed = false);

    float getContigShift(int index, bool reverse = false);

    void getContigDistances(map<int, map<int, float> >& contigContigShifts);
    /*
     void outputConnectorOverlap(int connIdx,
     int contig1,
     int contig2,
     SpectrumPair& res1,
     SpectrumPair& res2);
     */
  };

}

#endif /* COMBINECONTIGS_H_ */
