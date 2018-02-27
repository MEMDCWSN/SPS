#ifndef ALIGNMENT_MOTMUT

/** */
#define ALIGNMENT_MOTMUT

#include "spectrum.h"
#include "db_fasta.h"
#include "clusters.h"
#include "SpectrumPairSet.h"
#include <vector>
#include <deque>
#include <stack>
namespace specnets
{
  using namespace std;

  class AMM_match; // Defined below but used in AMM_peak_match

  /**
   *  Used inside scoreOverlapAMM/scoreOverlapAMMme to keep track of peak matches
   *  (current score + predecessor) while running the DP spectral alignment
   *  algorithm.
   */
  class AMM_peak_match
  {
  public:

    /**
     * Score of the best match up to the current cell.
     */
    float score;

    /**
     * Previous index for number of mods/muts.
     */
    int predMods;

    /**
     * Previous minimum number of matched peaks.
     */
    int predMPeaks;

    /**
     * Index of the previous PRM in the spectrum.
     */
    int predSpecIdx;

    /**
     * Index of the previous peak in the db sequence.
     */
    int predDbSpecIdx;

    /**
     *  Best match to a previous exon (if an immediate predecessor of the current match).
     */
    AMM_match *predExonMatch;

    //    AMM_peak_match() { score=-1; predMods=-1; predSpecIdx=-1; predDbSpecIdx=-1; }
    //    AMM_peak_match() { score=0; predMods=0; predSpecIdx=0; predDbSpecIdx=0; }

    /**
     * Default constructor - all values set to zero except score (-infinity) and predMPeaks (-1)
     *
     */
    AMM_peak_match()
    {
      score = -3.4e37;
      predMods = -1;
      predMPeaks = -1;
      predSpecIdx = -1;
      predDbSpecIdx = -1;
      predExonMatch = (AMM_match *) 0;
    }

    /**
     * Replaces field values with parameter values
     *
     *@param newScore
     *@param newPredMods
     *@param newPredSpecIdx
     *@param newPredDbSpecIdx
     *@param newPredMPeaks
     *@param newPredExonMatch
     */
    void set(float newScore,
             int newPredMods,
             int newPredSpecIdx,
             int newPredDbSpecIdx,
             int newPredMPeaks = -1,
             AMM_match *newPredExonMatch = 0)
    {
      score = newScore;
      predMods = newPredMods;
      predMPeaks = newPredMPeaks;
      predSpecIdx = newPredSpecIdx;
      predDbSpecIdx = newPredDbSpecIdx;
      predExonMatch = newPredExonMatch;
    }

    /**
     * Replaces field values with parameter values if newScore is equal to or higher than current score
     *
     *@param newScore
     *@param newPredMods
     *@param newPredSpecIdx
     *@param newPredDbSpecIdx
     *@param newPredExonMatch
     */
    void replaceWithMax(float newScore,
                        int newPredMods,
                        int newPredSpecIdx,
                        int newPredDbSpecIdx,
                        int newPredMPeaks = -1,
                        AMM_match *newPredExonMatch = 0)
    {
      if (newScore >= score)
        set(newScore,
            newPredMods,
            newPredSpecIdx,
            newPredDbSpecIdx,
            newPredMPeaks,
            newPredExonMatch);
    }

    /**
     * Replaces field values with parameter values if newScore is equal to or higher than current score
     *
     *@param other
     */
    void replaceWithMax(AMM_peak_match &other)
    {
      if (other.score >= score)
        set(other.score,
            other.predMods,
            other.predSpecIdx,
            other.predDbSpecIdx,
            other.predMPeaks,
            other.predExonMatch);
    }
  };

  /**
   * Used to describe a match between one spectrum and a specific protein sequence.
   */
  class AMM_match
  {
  public:

    /**
     * Index of the protein that the spectrum was matched to (in a companion DB_fasta object)
     *   (if matched to a single protein in non-multi-exon mode)
     */
    int proteinIdx;

    /**
     * Index of the (exon type, exon allele) that the spectrum was matched to (in companion DB_fasta objects)
     *   (if matched in multi-exon mode)
     */
    TwoValues<int> exonAllele;

    /**
     * Pointer to previous matches to a different (exon type, exon allele)
     *   (if matched in multi-exon mode)
     */
    AMM_match *predExon;

    /**
     * Index of the first matched amino acid (in the protein sequence in a companion DB_fasta object)
     *   (if matched to a single protein in non-multi-exon mode)
     */
    int aaStart;

    /**
     * Index of the last matched amino acid (in the protein sequence in a companion DB_fasta object)
     *   (if matched to a single protein in non-multi-exon mode)
     */
    int aaEnd;

    /**
     * Score of the spectrum/protein match
     */
    float matchScore;

    /**
     * Mass of the modification (smaller mod sizes are considered better).
     */
    float modSize;

    /**
     * Parent mass of the matched experimental spectrum.
     * Used to check whether the match completely explains the spectrum.
     */
    float specParentMass;

    /**
     * 0 if PRMs are interpreted as b-ions, 1 if PRMs are interpreted as y ions
     */
    short orientationPRMs;

    /**
     * List of matched masses between two spectra
     *
     * col 0 for experimental spec, col 1 for dbSpec
     */
    deque<TwoValues<float> > matchedMasses;

    /**
     * List of matched peak indices between two spectra
     *
     */
    deque<TwoValues<int> > matchedIndices;

    /**
     * Default constructor
     */
    AMM_match()
    {
      reset();
    }

    /**
     *
     * Copy constructor.
     * **Note: does NOT copy matchedMasses.
     *
     *@param other
     */
    AMM_match(const AMM_match &other);

    /**
     * reset - returns the object to default values
     */
    void reset()
    {
      matchScore = -3.4e37;
      proteinIdx = -1;
      aaStart = -1;
      aaEnd = -1;
      orientationPRMs = 0;
      exonAllele.set(-1, -1);

      specParentMass = -1;
      matchedMasses.clear();
      matchedIndices.clear();
      predExon = (AMM_match *) 0;
    }

    //    ~AMM_match() { for(int i=0; i<matchedMasses.size(); i++) matchedMasses.pop_back(); }

    //    void setSequence(char *newSeq) { if(sequence!=(char *)0) delete[] sequence;  sequence = new char[strlen(newSeq)+1];   strcpy(sequence,newSeq); }
    //    void setProteinID(char *newID) { if(proteinId!=(char *)0) delete[] proteinId;  proteinId = new char[strlen(newID)+1];   strcpy(proteinId,newID); }


    /**
     * TODO: add description
     *
     *@param spec
     *@param dbSpec
     *@param start
     *@param matchMatrix
     *@param curModSize
     */
    void setMatch(Spectrum spec, Spectrum dbSpec, AMM_peak_match start, vector<
        vector<vector<AMM_peak_match> > > &matchMatrix, float curModSize);

    /**
     * TODO: add description
     *
     *@param spec
     *@param dbSpec
     *@param start
     *@param matchMatrix
     *@param curModSize
     */
    void
        setMatch(Spectrum spec,
                 Spectrum dbSpec,
                 AMM_peak_match start,
                 vector<vector<vector<vector<AMM_peak_match> > > > &matchMatrix,
                 float curModSize);

    /**
     * TODO: add description
     *
     */
    void toString();

    /**
     * TODO: add description
     *
     *@param output
     *@param db
     *@param tolerance
     *@param separator
     *@param outputPRMs
     *@param specIdx
     *@param numMods
     */
    AMM_match *output_csv(ostream &output,
                          DB_fasta &db,
                          float tolerance,
                          char separator,
                          short outputPRMs = 0,
                          int specIdx = -1,
                          int numMods = -1);
  };

  /**
   * Used to keep track of set of best matches of one spectrum on all protein
   * sequences (per number of mod/muts).
   */
  class AMM_match_spec
  {

  public:

    /**
     * Number of matches to keep per number of mod/muts.
     */
    short numMatchesToKeep;

    /**
     * Lowest match score per number of mod/muts.
     */
    vector<float> lowestScore;

    /**
     * Top k matches per number of mod/muts.
     */
    vector<deque<AMM_match> > matches;

    /**
     * Initialize data structures to keep track of the best K matches
     *   with up to maxNumMods modifications/mutations per match
     */
    void init(short k, short maxNumMods)
    {
      lowestScore.resize(maxNumMods + 1);
      for (int i = 0; i < maxNumMods + 1; i++)
        lowestScore[i] = 0;
      matches.resize(maxNumMods + 1);
      numMatchesToKeep = k;
    }

    /**
     * Insert a new match with numMods mods/muts
     */
    void addMatch(short numMods, AMM_match match);

    /**
     * Output all matches to output
     */
    void output_csv(ostream &output,
                    DB_fasta &db,
                    float tolerance,
                    char separator = ';',
                    short outputPRMs = 0,
                    int specIdx = -1);
  };

  class AMM_match_simple
  { // Simplified version of AMM_match above
  public:

    /**
     * Index of the spectrum matched to the protein.
     */
    unsigned int index;

    /**
     * Maximum index matched on the protein.
     */
    unsigned int endAA;

    /**
     * Set of matches (spectrum index, protein mass index).
     */
    Spectrum matchedPeaksIdx;
  };

  /**
   * Keeps track of the results for a whole experiment - all the matches of all spectra
   * onto all protein sequences. Used to merge multiple matches into a single long
   * protein match.
   */
  class AMM_match_set
  {
  public:

    /**
     * Database where the matches point to.
     */
    DB_fasta *db;

    /**
     * Sets of matched peaks between spectra and proteins.
     */
    SpecSet matchedPeaks;

    /**
     * matches[i] - set of spectra matched to the i-th protein. Each pair contains
     *              (index of first matched protein mass, spectrum idx), list is sorted
     *              by increasing matched protein mass (position 0).
     */
    vector<list<TwoValues<int> > > matches;

    /**
     * Matched protein idx, #mods, b/y-match (cols 0-2, per matched spectrum).
     */
    vector<vector<int> > specProtMatches;

    /**
     * Constructor
     *
     *@param db
     */
    AMM_match_set(DB_fasta &db, unsigned int szSpecSet);

    /**
     * Resizes for different numbers of matched spectra (matchedPeaks and specProtMatches)
     *
     *@param szSpecSet  Number of matched spectra
     */
    unsigned int resize(unsigned int szSpecSet);

    /**
     * Set entries for the spectrum with index specIdx to the given spectrum/protein match.
     *
     *@param specIdx Index of the matched spectrum
     *@param match Details of the spectrum/protein match
     *@return true if the assignment was successful, failure if any indices were out of bounds (spectrum becomes unmatched)
     */
    bool set(unsigned int specIdx, AMM_match &match);

    /**
     * TODO: add description
     *
     */
    bool SetMatches(vector<vector<int> > &in_matchedProts,
                    SpecSet &in_matchedPeaks);
                    
    bool SetMatches(SpecSet * spectra);

    /**
     * TODO: add description
     *
     */
    bool LoadMatches(const char *matchedProts, const char *matchedPeaksIdx);

    /**
     * Defines spectrum overlaps based on spectrum-protein matches:
     *   two spectra are defined to overlap if both match to overlapping
     *   regions on the same protein and have at least one common matched
     *   mass on the matched protein sequence.
     */
    void GenerateGlues(SpectrumPairSet &pairsPA,
                       vector<vector<TwoValues<int> > > &vPairMatchedPeaks);

    /**
     * TODO: add description
     *
     */
    void MergeIntoReference(int refIdx,
                            int otherIdx,
                            vector<TwoValues<int> > &matchedIndices);

    /**
     * Outputs a string representation of a spectrum/protein match to an output stream.
     *
     *@param specIdx
     *@param output
     *@param db
     *@param tolerance
     *@param separator
     *@param outputPRMs
     */
    void output_csv(unsigned int specIdx,
                    ostream &output,
                    float tolerance,
                    char separator,
                    short outputPRMs = 0)
    {
    }

    /**
     * Outputs the set of matched spectrum/protein masses (per spectrum)
     *
     *@param filename Name of the output file
     *@return True if the file saved correctly, false otherwise
     */
    bool output_midx(const char *filename)
    {
      return false;
    }

    /**
     * Outputs matched-protein-index and #mods (best match per spectrum)
     *
     *@param filename Name of the output file
     *@return True if the file saved correctly, false otherwise
     */
    bool output_mp(const char *filename)
    {
      return false;
    }
  };

  //vector<vector<AMM_match> > scoreOverlapAMM(SpecSet spectra, DB_fasta db, Clusters clst, short maxNumMods, short keepTopK, float pmTolerance, float tolerance);
  //vector<AMM_match> scoreOverlapAMM(Spectrum spec, Spectrum dbSpec, short maxNumMods, short keepTopK, float pmTolerance, float tolerance);

  /**
   * scoreOverlapAMM - scores overlaps between two prefix-only PRM spectra allowing
   * up to maxNumMods Modifications/Mutations (mass offsets) between the two spectra.
   *
   *@param spec - first spectrum, usually from a contig sequence
   *@param dbSpec - second spectrum, usually from a protein sequence
   *@param maxNumMods - maximum number of allowed mass offsets between the two spectra
   *@param topKmatches - number of top-scoring solutions to be returned
   *@param curMatch - information on current match: orientation of PRMs/endpoints, proteinIdx
   *@param pmTolerance - parent mass tolerance
   *@param tolerance - matched-mass (e.g., amino acid) mass tolerance
   *@param minSpecDist - minimum distance between matched peaks in spec
   *@param maxDbSpecMod - maximum modification mass
   *@param minDbSpecMod - minimum modification mass
   *@param enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks (consuming available modifications as necessary)
   */
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       short maxNumMods,
                       AMM_match_spec &topKmatches,
                       AMM_match &curMatch,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks = false);

  /**
   * scoreOverlapAMM - scores overlaps between two prefix-only PRM spectra allowing
   * up to maxNumMods Modifications/Mutations (mass offsets) between the two spectra and
   * requiring a minimum number of matched peaks.
   *
   *@param spec - first spectrum, usually from a contig sequence
   *@param alleleSpecs - set of spectra, usually from a set of alleles for the same type of exons
   *@param maxNumMods - maximum number of allowed mass offsets
   *@param minMatchedPeaks - minimum acceptable number of matched peaks in spec
   *@param topKmatches - number of top-scoring solutions to be returned
   *@param curMatch - information on current match: orientation of PRMs/endpoints, proteinIdx
   *@param pmTolerance - parent mass tolerance
   *@param tolerance - matched-mass (e.g., amino acid) mass tolerance
   *@param minSpecDist - minimum distance between matched peaks in spec
   *@param maxDbSpecMod - maximum modification mass
   *@param minDbSpecMod - minimum modification mass
   *@param enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks (consuming available modifications as necessary)
   */
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       short maxNumMods,
                       int minMatchedPeaks,
                       AMM_match_spec &topKmatches,
                       AMM_match &curMatch,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks = false);

  /**
   * scoreOverlapAMM - scores overlaps between two prefix-only PRM spectra allowing
   * up to maxNumMods Modifications/Mutations (mass offsets) between the two spectra and
   * requiring a minimum number of matched peaks.
   *
   *@param spec - first spectrum, usually from a contig sequence
   *@param dbSpec - database spectra
   *@param dbIndex - index of spectrum in database
   *@param matchOrientation - orientation of spec (0 for forward, 1 for reverse)
   *@param startRange - set of possible starting masses
   *@param maxNumMods - maximum number of allowed mass offsets
   *@param maxNumMods - maximum number of allowed mass offsets
   *@param minMatchedPeaks - minimum acceptable number of matched peaks in spec
   *@param pmTolerance - parent mass tolerance
   *@param tolerance - matched-mass (e.g., amino acid) mass tolerance
   *@param minSpecDist - minimum distance between matched peaks in spec
   *@param maxDbSpecMod - maximum modification mass
   *@param minDbSpecMod - minimum modification mass
   *@param enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks (consuming available modifications as necessary)
   */
  void scoreOverlapAMM(Spectrum &spec,
                       Spectrum &dbSpec,
                       int dbIndex,
                       int matchOrientation,
                       set<float> & startRange,
                       short maxNumMods,
                       int minMatchedPeaks,
                       float pmTolerance,
                       float tolerance,
                       float minSpecDist,
                       float maxDbSpecMod,
                       float minDbSpecMod,
                       bool enforceEndpeaks = false);

  /**
   * scoreOverlapAMMme - scores overlaps between a prefix-only PRM spectrum and a set of
   * exon types, each with one or more possible alleles.
   *
   *@param spec - first spectrum, usually from a contig sequence
   *@param alleleSpecs - set of spectra, usually from a set of alleles for the same type of exons
   *@param maxNumMods - maximum number of allowed mass offsets
   *@param minMatchedPeaks - minimum acceptable number of matched peaks in spec
   *@param topKmatches - number of top-scoring solutions to be returned
   *@param curMatch - information on current match: orientation of PRMs/endpoints, proteinIdx
   *@param pmTolerance - parent mass tolerance
   *@param tolerance - matched-mass (e.g., amino acid) mass tolerance
   *@param minSpecDist - minimum distance between matched peaks in spec
   *@param maxDbSpecMod - maximum modification mass
   *@param minDbSpecMod - minimum modification mass
   *@param enforceEndpeaks - if true then peaks 0 and spec.size()-1 must match database peaks (consuming available modifications as necessary)
   */
  void scoreOverlapAMMme(Spectrum &spec,
                         Spectrum &dbSpec,
                         short maxNumMods,
                         int minMatchedPeaks,
                         AMM_match_spec &topKmatches,
                         AMM_match &curMatch,
                         vector<vector<AMM_match> > *prevMaxMatch,
                         vector<vector<AMM_match> > *curMaxMatch,
                         float pmTolerance,
                         float tolerance,
                         float minSpecDist,
                         float maxDbSpecMod,
                         float minDbSpecMod,
                         bool enforceEndpeaks = false);
}

#endif

