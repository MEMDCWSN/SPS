#ifndef __ABRUIJN__H__
#define __ABRUIJN__H__

#include <map>
#include <vector>
#include <cstdio>
#include <numeric>
#include <stdint.h>

#include "spectrum.h"
#include "graph.h"
//#include "batch.h"
#include "SpectrumPair.h"
#include "tuple.h"
#include "AbruijnGraph.h"

namespace specnets
{

  /**
   * Match and MatchEdge are used when constructing the A-Bruijn graph - the composition
   * of the final vertices is not known until the whole graph is constructed
   * because vertices get merged along the way. Keeping this additional Match
   * and MatchEdge information along the way avoids more complicated merging
   * of A-Bruijn vertices and facilitates splitting of composite vertices.
   *
   * Match - stores a pair of peaks matched between a pair of spectra (one peak from each spectrum)
   * MatchEdge - stores consecutive peak matches (represented by Match objects) between two
   *             aligned spectra
   *
   */
  struct Match
  {

    /**
     * Peak matched in the first spectrum: spectrum index in specPeak1[0]
     * and peak index in specPeak1[1]
     */
    TwoValues<int> specPeak1;

    /**
     * Peak matched in the second spectrum: spectrum index in specPeak2[0]
     * and peak index in specPeak2[1]
     */
    TwoValues<int> specPeak2;

    /**
     * Default constructor; sets specPeak1/specPeak2 to (-1,-1).
     */
    Match()
    {
      specPeak1.set(-1, -1);
      specPeak2.set(-1, -1);
    }

    /**
     * Sets match to p1/p2
     *
     *@param p1 Matched peak in the first spectrum (formatted as specPeak1)
     *@param p2 Matched peak in the second spectrum (formatted as specPeak2)
     */
    Match(TwoValues<int> &p1, TwoValues<int> &p2)
    {
      specPeak1 = p1;
      specPeak2 = p2;
    }

    /**
     * Sets match to p1/p2
     *
     *@param p1 Matched peak in the first spectrum (formatted as specPeak1)
     *@param p2 Matched peak in the second spectrum (formatted as specPeak2)
     *@return Reference to self
     */
    Match &set(TwoValues<int> &p1, TwoValues<int> &p2)
    {
      specPeak1 = p1;
      specPeak2 = p2;
      return *this;
    }

    /**
     * Sets specPeak1 to (s1,p1) and specPeak2 to (s2,p2)
     *
     *@param s1 Index of the first spectrum
     *@param p1 Index of the matched peak in the first spectrum
     *@param s2 Index of the second spectrum
     *@param p2 Index of the matched peak in the second spectrum
     *@return Reference to self
     */
    Match &set(int s1, int p1, int s2, int p2)
    {
      specPeak1[0] = s1;
      specPeak1[1] = p1;
      specPeak2[0] = s2;
      specPeak2[1] = p2;
      return *this;
    }

    /**
     * Sets specPeak1/specPeak2 as in other (parameter)
     *
     *@param other Source match
     *@return Reference to self
     */
    Match &operator=(Match other)
    {
      specPeak1 = other.specPeak1;
      specPeak2 = other.specPeak2;
      return *this;
    }
  };

  /**
   * Match and MatchEdge are used when constructing the A-Bruijn graph - the composition
   * of the final vertices is not known until the whole graph is constructed
   * because vertices get merged along the way. Keeping this additional Match
   * and MatchEdge information along the way avoids more complicated merging
   * of A-Bruijn vertices and facilitates splitting of composite vertices.
   *
   * Match - stores a pair of peaks matched between a pair of spectra (one peak from each spectrum)
   * MatchEdge - stores consecutive peak matches (represented by Match objects) between two
   *             aligned spectra
   *
   */
  struct MatchEdge
  {

    /**
     * Mass assigned to the edge
     *
     * TODO: in which spectrum??
     *
     */
    float edgeMass;

    /**
     * First pair of matched peaks
     */
    Match from;

    /**
     * Second pair of matched peaks
     */
    Match to;

    /**
     * Constructor without info on matched peaks
     *
     *@param mass Mass assigned to the edge
     */
    MatchEdge(float mass = 0)
    {
      edgeMass = mass;
    }

    /**
     * Constructor with info on matched peaks
     *
     *@param mass Mass assigned to the edge
     *@param m1 First pair of matched peaks
     *@param m2 Second pair of matched peaks
     */
    MatchEdge(float mass, Match m1, Match m2)
    {
      edgeMass = mass;
      from = m1;
      to = m2;
    }
  };

  /**
   * Edge between peaks in the same spectrum
   */
  struct SimpleEdge
  {

    /**
     * Mass assigned to the edge
     */
    float edgeMass;

    /**
     * Score assigned to the edge
     */
    float edgeScore;

    /**
     * Index of the spectrum where the edge occurs
     */
    int specIdx;

    /**
     * Index of destination vertex
     *
     * ONLY SET AFTER CALLING VertexSet::addParallelPaths
     */
    int vertexIdx;

    /**
     * Indices of the peaks connected by the edge  (in spectrum specIdx)
     */
    TwoValues<int> peaks;

    /**
     * Constructor with edge mass
     *
     *@param mass Mass assigned to the edge
     */
    SimpleEdge(float mass = 0)
    {
      edgeMass = mass;
      vertexIdx = -1;
    }

    /**
     * Constructor with all info
     *
     *@param mass Mass assigned to the edge
     *@param score Score assigned to the edge
     *@param in_specIdx Index of the spectrum where the edge occurs
     *@param peakIdx1 Index of the source peak connected by the edge (in spectrum specIdx)
     *@param peakIdx2 Index of the sink peak connected by the edge (in spectrum specIdx)
     */
    SimpleEdge(float mass,
               float score,
               int in_specIdx,
               int peakIdx1,
               int peakIdx2,
               int vertIdx = -1)
    {
      set(mass, score, in_specIdx, peakIdx1, peakIdx2, vertIdx);
    }

    /**
     * Sets values as specified
     *
     *@param mass Mass assigned to the edge
     *@param score Score assigned to the edge
     *@param in_specIdx Index of the spectrum where the edge occurs
     *@param peakIdx1 Index of the source peak connected by the edge
     *@param peakIdx2 Index of the sink peak connected by the edge
     *@return Reference to self
     */
    SimpleEdge &set(float mass,
                    float score,
                    int in_specIdx,
                    int peakIdx1,
                    int peakIdx2,
                    int vertIdx = -1)
    {
      edgeMass = mass;
      edgeScore = score;
      specIdx = in_specIdx;
      peaks.set(peakIdx1, peakIdx2);
      vertexIdx = vertIdx;
      return *this;
    }

    /**
     * Tests whether two SimpleEdge objects indicate the same connected peaks in
     * the same spectrum; edge score and mass do not affect the outcome.
     *
     *@param other SimpleEdge object to compare against
     *@return true if self equals other
     */
    bool operator==(SimpleEdge &other)
    {
      return specIdx == other.specIdx and peaks == other.peaks and vertexIdx == other.vertexIdx;
    }

    /**
     * Sets values as specified in other (parameter)
     *
     *@param other Source object
     *@return Reference to self
     */
    SimpleEdge &operator=(const SimpleEdge &other)
    {
      specIdx = other.specIdx;
      peaks = other.peaks;
      edgeMass = other.edgeMass;
      edgeScore = other.edgeScore;
      vertexIdx = other.vertexIdx;
      return *this;
    }

    /**
     * self precedes other if at least one of the following is true:
     *   1) specIdx < other.specIdx
     *   2) specIdx==other.specIdx and peaks[0]<other.peaks[0]
     *   3) specIdx==other.specIdx and peaks[0]==other.peaks.values[0] and peaks[1]<other.peaks[1]
     *
     *@param other SimpleEdge object to compare against
     *@return true if self precedes other
     */
    bool operator<(const SimpleEdge &other)
    {
      if (specIdx < other.specIdx or (specIdx == other.specIdx
          and peaks[0] < other.peaks[0]) or (specIdx
          == other.specIdx and peaks[0] == other.peaks[0]
          and peaks[1] < other.peaks[1]))
        return true;
      return false;
    }
  };

  /**
   * Vertex in the A-Bruijn graph, defined as a set of matched spectrum peaks.
   */
  class Vertex
  {
  public:
    /**
     * Sorted list of (specIdx,peakIdx) spectrum peak members.
     */
    list<TwoValues<int> > specPeaks;

    /**
     * List of outgoing edges; MatchEdge objects indicate consecutive matched peaks
     *  between two aligned spectra.
     */
    list<MatchEdge> outMatchEdges;

    /**
     * List of outgoing edges; SimpleEdge objects indicate pairs of connected peaks in
     *  the same spectrum.
     */
    list<SimpleEdge> outEdges;

    /**
     * True if Vertex contains two different peaks from the same spectrum (in specPeaks).
     */
    bool compositeVertex;

    /**
     * Default constructor (see reset()).
     */
    Vertex()
    {
      reset();
    }

    /**
     * Copy constructor.
     */
    Vertex(const Vertex &other);

    /**
     * Adds a new spectrum peak to Vertex.
     *
     *@param specPeak Pair of values: spectrum index, peak index
     */
    void addPeak(TwoValues<int> specPeak);

    /**
     * Adds a new outgoing edge from consecutive aligned peaks (currently NOT USED).
     *
     *@param edgeMass Mass assigned to the edge
     *@param from First pair of matched peaks
     *@param to Second pair of matched peaks
     */
    void addMatchEdge(float edgeMass, Match from, Match to);

    /**
     * Adds a new outgoing edge.
     *
     *@param edgeMass Mass assigned to the edge
     *@param edgeScore Score assigned to the edge
     *@param specIdx Index of the spectrum where the edge occurs
     *@param peakIdx1 Index of the source peak (in spectrum specIdx)
     *@param peakIdx2 Index of the sink peak (in spectrum specIdx)
     */
    void addEdge(float edgeMass,
                 float edgeScore,
                 int specIdx,
                 int peakIdx1,
                 int peakIdx2);

    /**
     * Merges two ABruijn vertices (spectrum/peak indices and outgoing edges).
     *
     *@param withVertex Vertex to merge with
     */
    void merge(Vertex &withVertex);

    /**
     * Resets the vertex back to empty - zero spectrum peaks and outgoing edges.
     */
    void reset();

    /**
     * Size as number of spectrum peaks grouped in the vertex.
     *
     *@return Number of spectrum peaks in the vertex
     */
    unsigned int size()
    {
      return specPeaks.size();
    }

    /**
     * Modifies an outgoing edge from a spectrum peak in the vertex. This function is
     * used by VertexSet::consolidatePaths().
     *
     *@param toSpecIdx Index of the spectrum containing the peak
     *@param toPeakIdx Index of the current edge sink peak
     *@param newPeakIdx Index of the new edge sink peak
     *@param addEdgeMass Mass difference in relation to current edge mass
     *@param addEdgeScore Score difference in relation to current score mass
     */
    bool replaceEdge(int toSpecIdx,
                     int toPeakIdx,
                     int newPeakIdx,
                     float addEdgeMass,
                     float addEdgeScore);
  };

  /**
   * Helper structure for the construction of an A-Bruijn graph. Keeps
   * track of existing vertices, their composition and takes care of
   * merging (when new spectrum peak matches result in merged ABruijn
   * vertices).
   */
  class VertexSet
  {
    /**
     * Number of new entries to allocate if all vertices are
     *   being used and more are needed (in vertices/freeVertices)
     */
    static const int SZEXPAND = 1024;

    /**
     * Simply-linked list of free vertex cells in vertices; avoids need for frequent
     *   allocations/deletions on Vertex objects.
     */
    vector<int> freeVertices;

    /**
     * Index of the first free vertex cell in vertices.
     */
    int firstFreeVertex;

    /**
     * Returns a Vertex in vertices to the pool of available Vertex objects
     *
     *@param idx Index of the vertex to release
     */
    void addFreeVertex(int idx)
    {
      freeVertices[idx] = firstFreeVertex;
      firstFreeVertex = idx;
    }

    /**
     * Gets an available Vertex object from the pool of free Vertex objects.
     *
     *@return Index of the Vertex object
     */
    int getFreeVertex()
    {
      if (firstFreeVertex == -1)
      { // Current structure is full, need to allocate more space
        unsigned int oldSize = vertices.size();
        vertices.resize(vertices.size() + SZEXPAND);
        freeVertices.resize(freeVertices.size() + SZEXPAND);
        for (unsigned int i = oldSize; i < freeVertices.size() - 1; i++)
          freeVertices[i] = i + 1;
        freeVertices[freeVertices.size() - 1] = -1;
        firstFreeVertex = oldSize;
      }
      int t = firstFreeVertex;
      firstFreeVertex = freeVertices[t];
      return t;
    }

    /**
     * Merges two Vertex objects using a pair of matched spectrum peaks
     *
     *@param m Pair of consecutive matched spectrum peaks
     *@return Index of the merged Vertex object (in vertices)
     */
    int mergeVertices(Match m);

    /**
     * Creates a new Vertex object with a single spectrum peak
     *
     *@param specIdx Index of the spectrum containing the peak
     *@param peakIdx Index of the peak in specIdx
     *@return Index of the Vertex object (in vertices)
     */
    int addVertex(int specIdx, int peakIdx);

    /**
     * Resets a Vertex in vertices and returns it to the pool of available Vertex objects
     *
     *@param vertex Index of the vertex to release
     *@return Updated number of used Vertex objects
     */
    int releaseVertex(int vertex)
    {
      vertices[vertex].reset();
      addFreeVertex(vertex);
      return --numUsedVertices;
    }

    /**
     * Indices of all composite vertices left to process
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<int> scv_compVerts;

    /**
     * List of predecessor edges (per vertex).
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<list<SimpleEdge> > scv_predEdges;

    /**
     * List of successor edges (per vertex).
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<list<SimpleEdge> > scv_succEdges;

    /**
     * List of internal edges (per vertex).
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<list<SimpleEdge> > scv_intnEdges;

    /**
     * List of composite/composite edges (per vertex).
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<list<SimpleEdge> > scv_compEdges;

    /**
     * Lists of in/out split edges per composite vertex.
     *   (used when splitting composite vertices with splitComposite)
     */
    vector<TwoValues<list<SimpleEdge> > > scv_splitEdges;

    /**
     * Adds edges to the pred/succs/intn/comp lists of adjacent edges
     *   (used when splitting composite vertices with splitComposite)
     *
     *@param edges
     *@param vertexAdjChanged
     */
    void scv_addEdges(list<SimpleEdge> &edges,
                      vector<TwoValues<bool> > &vertexAdjChanged);

    /**
     * Finds the best split edge for a given composite vertex cvIdx and a
     * pred/succ direction (updates scv_splitEdges)
     *   (used when splitting composite vertices with splitComposite)
     *
     *@param cvIdx Index of the composite vertex (in scv_compVerts)
     *@param direction Set to 0/1 to consider in/out edges, respectively
     *@param peakTol Allowed tolerance for mass errors (in Da)
     */
    void scv_findSplitEdge(int cvIdx, int direction, float peakTol);

    /**
     * Uses a split edge to break a composite vertex into two vertices: a new
     * vertex which is guaranteed not to be composite and a leftover vertex with
     * the remaining glued spectrum peaks that may still be composite.
     *   (used when splitting composite vertices with splitComposite)
     *
     *@param cvIdx Index of the composite vertex (in scv_compVerts)
     *@param direction Set to 0/1 to consider in/out edges, respectively
     *@param vertexAdjChanged
     */
    void scv_useSplitEdge(int cvIdx,
                          int direction,
                          vector<TwoValues<bool> > &vertexAdjChanged);

  public:

    /**
     * Edge Score Type EST_EDGE_MULT: Edge multiplicity.
     */
    static const short EST_EDGE_MULT = 0;

    /**
     * Edge Score Type EST_EDGE_SCORES: each edge gets its score from the
     *   destination peak in the original spectrum.
     */
    static const short EST_EDGE_SCORES = 1;

    /**
     * Edge Score Type EST_ABVERTEX_SCORES: each edge gets its score from the
     * destination A-Bruijn vertex. Edge multiplicity is added to vertex scores to
     * help distinguish between edges
     */
    static const short EST_ABVERTEX_SCORES = 2;

    /**
     * peakToVertex[i][j] - index of the Vertex object containing peak j from spectrum i
     */
    vector<vector<int> > peakToVertex;

    /**
     * Set of ABruijn vertices
     */
    vector<Vertex> vertices;

    /**
     * Number of used vertices
     */
    int numUsedVertices;

    /**
     * Pointer to the set containing the spectra used to construct the ABruijn graph
     */
    SpecSet *specSet;

    /**
     * Constructor from set of spectra and number of vertices
     *
     *@param specs Set containing the spectra used to construct the ABruijn graph
     *@param numVertices Number of Vertex objects to pre-allocate
     */
    VertexSet(SpecSet &specs, int numVertices)
    {
      peakToVertex.resize(specs.size());
      for (unsigned int i = 0; i < peakToVertex.size(); i++)
      {
        peakToVertex[i].resize(specs[i].size());
        for (unsigned int j = 0; j < peakToVertex[i].size(); j++)
          peakToVertex[i][j] = -1;
      }
      specSet = &specs;
      vertices.resize(numVertices);
      freeVertices.resize(numVertices);
      numUsedVertices = 0;
      for (unsigned int i = 0; i < freeVertices.size() - 1; i++)
        freeVertices[i] = i + 1;
      freeVertices[freeVertices.size() - 1] = -1;
      firstFreeVertex = 0;
    }

    /**
     * Copy constructor.
     */
    VertexSet(const VertexSet &other);

    /**
     * Glues/merges ABruijn vertices matched by pairwise alignment, as indicated in matches.
     *   If spectrumGraphs is NULL then also adds ABruijn edges for every pair of consecutive
     *   matched peaks (edge scores are set to zero).
     *
     *@param specIdx1 Index of the first aligned spectrum
     *@param specIdx2 Index of the second aligned spectrum
     *@param matches Peaks matched by spectral alignment
     *@param spectrumGraphs If not NULL then ABruijn edges are created for all consecutive matched peaks
     *@param addMissingEdges If true, spectrumGraphs are extended to include edges between consecutive matched peaks (not implemented)
     */
    void addGlues(int specIdx1,
                  int specIdx2,
                  vector<TwoValues<int> > &matches,
                  vector<MSGraph> *spectrumGraphs = 0,
                  bool addMissingEdges = false);

    /**
     * Creates edges between ABruijn vertices using edges from each spectrum's
     *   spectrum graph - two ABruijn vertices A, A' become connected if A contains (i,j)
     *   and A' contains (i,k) such that peaks j and k are connected by an edge in spectrum i.
     *
     *@param spectrumGraphs Set of spectrum graphs, one per spectrum
     *@param usedSpectra If specified, it's updated to indicate which spectra have at
     *                      least one peak in an ABruijn vertex.
     */
    void addEdges(vector<MSGraph> &spectrumGraphs, vector<bool> *usedSpectra =
        0);

    /**
     * Adds edges between endpoints (start/end) of paired spectra to account
     *   for the placement of the modification IFF the mass difference between aligned
     *   spectra is contained in jumps (i.e., masses larger than those in jumps do not
     *   generate edges).
     *
     * WARNING: Assumes that the first two peaks in each spectrum are b0/y0 and that
     *            the last two peaks are b(n-1)/y(n-1)
     *
     *@param aligns Set of pairwise alignments
     *@param matches Sets of matched peaks for each pairwise alignment
     *@param modPos Mass-location for the modification mass between each pair of aligned spectra
     *@param jumps Set of allowed modification masses
     *@param peakTol Allowed tolerance for mass errors (in Da)
     */
    void addEndpointEdges(SpectrumPairSet &aligns,
                          vector<vector<TwoValues<int> > > &matches,
                          vector<float> &modPos,
                          AAJumps &jumps,
                          float maxModMass,
                          float peakTol,
                          float pmTol);

    /**
     * Uses edges between composite and non-composite vertices to separate composite
     * vertices into two or more non-composite vertices (similar to RepeatGluer).
     * In brief, composite vertices are split by the multiplicity of incoming/outgoing
     * edges.
     *
     *@param spectrumGraphs Sets of edges used to split composite vertices
     *@param peakTol Allowed tolerance for mass errors (in Da)
     *@param usedSpectra If specified, it's updated to indicate which spectra have at
     *                      least one peak in an ABruijn vertex.
     */
    void splitComposite(vector<MSGraph> &spectrumGraphs, float peakTol, vector<
        bool> *usedSpectra = 0);

    /**
     * TODO: add description
     *
     *@param typeB
     *@param peakTol
     */
    void removeEndpoints(bool typeB, float peakTol);

    /**
     * Converts VertexSet into a graph structure by eliminating all references to
     *   spectrum peaks and merging edges by their amino acid masses as defined by
     *   jumps (edges are merged if they're closest to the same entry in jumps, within
     *   peakTol).
     *
     *@param g ABruijn graph output as an MSGraph
     *@param jumps Sets of amino acid masses
     *@param peakTol Allowed tolerance for mass errors (in Da)
     *@param vSet_index (by reference) vSet_index[i] is set to the index of the i-th vertex of g in this VertexSet.
     *@param edgeScoreType One of EST_EDGE_MULT, EST_EDGE_SCORES, EST_ABVERTEX_SCORES (as defined in this class)
     *@param labels
     */
    void buildGraph(MSGraph &g,
                    AAJumps &jumps,
                    float peakTol,
                    vector<int> &vSet_index,
                    vector<SpectrumPeakLabels> &labels,
                    short edgeScoreType = EST_EDGE_SCORES);

    /**
     * Converts VertexSet to an AbruijnGraph class structure for computing ParallelPaths and inserts expanded paths
     */
    //void addParallelPaths(float peakTol);

    int addAbruijnNode(abruijn::AbruijnNode* abNode, map<string, unsigned int>& specIdToIdx);

    /**
     * TODO: add description
     *
     */
    void consolidatePaths();

    /**
     * TODO: add description
     *
     *@param g
     *@param vSet_index
     *@param enumSpecs
     *@param enumSpecsIdx
     *@return
     */
    void enumerateHeaviestPath(MSGraph &g,
                               vector<int> &vSet_index,
                               SpecSet &enumSpecs,
                               vector<int> &enumSpecsIdx);

    /**
     * Selects a subset of vertices from the ABruijn graph.
     *
     *@param matchedVertices Set of vertices to select from the ABruijn graph; returns all vertices if matchedVertices.size()==0
     *@param putHere One vertex per vector position (putHere[i]), each with a list of (spectrum index, peak index) pairs
     */
    void getMatchedPeaks(vector<int> &matchedVertices, vector<list<TwoValues<
        int> > > &putHere);

    /**
     * TODO: add description
     *
     *@param labels
     *@return
     */
    void outputGraph(vector<SpectrumPeakLabels> &labels);

    /**
     * Outputs an A-Bruijn multiple alignment in graphviz format.
     *
     *@param filename
     *@param matchedVertices
     *@return
     */
    void output_graphviz_ma(char *filename, vector<int> &matchedVertices);
  };

  // File format for abinfo / component_info.bin (in masab.cpp) is:
  //     numUsedSpectra [int]
  //     numUsedSpectra-by-1 unsigned int array of all used spectrum indices
  //     numUsedSpectra-by-2 unsigned short array of [component index,specFlipped]
  //     numComponents [int]
  //     numComponents-by-1 unsigned short array of number of ABruijn vertices per component
  //     numABVertices-by-1 unsigned short array of number of spectrum peaks per ABruijn vertex
  //     totNumSpecPeaks-by-1 unsigned int array of spectrum index per peak
  // OLD    totNumSpecPeaks-by-1 unsigned short array of peak index per spectrum peak
  //     totNumSpecPeaks-by-1 float array of peak masses per spectrum peak

  /**
   * TODO: add description
   *
   *@param filename
   *@param specSet
   *@param cSpectra
   *@param specFlipped
   *@param abVertices
   */
  void Save_abinfo(const char *filename,
                   SpecSet &specSet,
                   vector<list<int> > &cSpectra,
                   vector<bool> &specFlipped,
                   vector<vector<list<TwoValues<int> > > > &abVertices);

  typedef std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > > abinfo_t;

  /**
   * Computes the shift of each assembled spectrum wrt mass 0 in its contig spectrum
   * @param inputContigSpectra set of contig spectra parallel to inputContigAbinfo
   * @param inputContigAbinfo set of abruijn components parallel to inputContigSpectra
   * @param outputAssembledShifts output data structure:
   *     1st index = 0-based index of contig spectrum (matches inputContigSpectra and inputContigAbinfo)
   *          2nd index = local index to lookup shift for star spectrum
   *              tuple.m0 = global 0-based index of star spectrum (taken from inputContigAbinfo)
   *              tuple.m1 = Da shift of mass 0 in star spectrum from mass 0 in contig spectrum
   *              tuple.m2 = true if star spectrum is reversed in reference to consensus, false if not
   * @return
   */
  void
      getAssembledShifts(SpecSet& inputContigSpectra,
                         abinfo_t& inputContigAbinfo,
                         vector<vector<sps::tuple<unsigned int, float, bool> > >& outputAssembledShifts);

  void Merge_abinfo_v1_0(SpecSet& specSet,
                         vector<list<int> >& cSpectra,
                         vector<bool>& specFlipped,
                         vector<vector<list<TwoValues<int> > > >& abVertices,
                         abinfo_t& merged_abinfo);

  /**
   * Appends multiple abinfo_t from parallel SPS projects
   *@inabinfo order of abinfo_t to append
   *@specCount parallel order of spectral count for SPS projects
   *@outabinfo output of appended abinfo_t
   *@return
   */
  void Append_abinfo_v1_0(list<abinfo_t>& inabinfo,
                          list<int>& specCount,
                          abinfo_t& outabinfo);

  /**
   * Combines abinfo when contigs have been merged to meta-contigs
   * @param childAbinfo abinfo of original contigs
   * @param contig_rev whether original contigs have been reversed
   * @param star_spectra spectra assembled into original contigs
   * @param contigsMerged whether each original contig was merged
   * @param contigs original contigs
   * @param parentAbinfo abinfo of parent contigs (no original contigs
   *   should be marked as reversed here)
   * @param outabinfo abinfo output mapping meta-contigs and their abruijn
   *   vertices to star spectrum indices
   * @return
   */
  void Combine_abinfo_v1_0(abinfo_t& childAbinfo,
                           vector<bool>& contig_rev,
                           SpecSet& star_spectra,
                           vector<bool>& contigsMerged,
                           SpecSet& contigs,
                           abinfo_t& parentAbinfo,
                           abinfo_t& outabinfo);

  /**
   * Physically merges 2 abinfo's into a new one.
   * ASSUMPTIONS:
   * 1.  The two abinfos were two interpretations of the same set of spectra, (1-N), but now we want one of the abinfos to be on spectra 1-N
   * and the other abinfo to be on spectra N+1 - 2N (spectrum 1 = spectrum N+1)
   */
  bool merge_abinfo(abinfo_t & info1,
                    abinfo_t & info2,
                    const char * filename,
                    unsigned int totalSpectra);

  void Copy_abinfo(abinfo_t& from, abinfo_t& to);

  bool Save_abinfo_v1_0(const char *filename, abinfo_t & abinfo);

  int Load_abinfo(const char *filename, abinfo_t & abinfo);



  int dumpAbInfo(const char * filename, abinfo_t & abinfo);


  //////////////////////////////////////////////////////////////////////////////
  // Data structures needed to return the abruijn comparison results
  //

  // container for a single node entry
  typedef struct {
    int     starId1, starId2;   // star ids
    double  val1, val2;         // peak values
  } abDiffPeak_t;

  // contains por the peaks in a node
  typedef vector<abDiffPeak_t>
  abDiffPeaks_t;

  // Container for the node info
  typedef struct {
    int size1i, size1d, size2i, size2d;
    int           idx;            // node index
    abDiffPeaks_t peaks;          // node data
  } abDiffNode_t;

  // abruijn nodes for the star
  typedef struct {
    int size1, size2;
    vector<abDiffNode_t> nodeData;
  } abDiffNodes_t;

  // IDs of start with different reversal flag
  typedef vector<int>
  abDiffRever_t;

  // contig container
  typedef struct {
    int             id1, id2;             // ids of the contigs at this position
    int             size1, size2;         // number of stars
    abDiffRever_t   reversedList;         // the diff for reversed list
    abDiffNodes_t   nodesList;            // the diff for the node list
  } abDiffContig_t;

  // all contigs. unsigned is contig index
  typedef vector<abDiffContig_t>
  abDiffContigs_t;

  // abruijn data container
  typedef struct {
    bool            different;      // are abruijns different?
    int             size1, size2;   // what are the abruijns sizes?
    abDiffContigs_t contigs;        // the abruijns diff data
  } abDiff_t;

  // Compares 2 abruijn graphs and returns difference data
  int compareAbInfo(abinfo_t & abinfo1, abinfo_t & abinfo2, abDiff_t &diffData);

}
#endif
