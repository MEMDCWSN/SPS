#ifndef ALIGNMENT_PENALTY_BASED2
#define ALIGNMENT_PENALTY_BASED2

#include "spectrum.h"
#include "PenaltyMatrix.h"
#include <set>

namespace specnets
{
  using namespace std;

  class AlignmentPenaltyBased
  {
  public:

    AlignmentPenaltyBased(int             maxSpecGap,
                          PenaltyMatrix * modPenaltyMatrix,
                          PenaltyMatrix * blossumPenaltyMatrix,
                          float           penaltyAlpha,
                          float           penaltyBeta,
                          float           maxMod,
                          float           minMod);
    ~AlignmentPenaltyBased();

    void clearCache(void);
    bool isCached(string & aaString);
    void cacheProteinStrings(string & proteinString, int maxLength);
    void sortStringLetters(string & toSort);

    void computeAlignment(Spectrum &	spec,
                          Spectrum &	dbSpec,
                          char * 		dbSeq,
                          int 		dbIndex,
                          int 		matchOrientation,
                          set<float> & 	startRange,
                          map<int, float> & scanSpecificPenalties,
                          int 		minMatchedPeaks,
                          int 		maxGapSize,
                          float 		pmTolerance,
                          float 		tolerance,
                          float           avgPeakIntensity,
                          bool 		enforceEndpeaks = false);

    void computeAllGapAnnotations(string & stringAnnotationIn,
                                  string & stringAnnotationOut);

    float computeCleavagePenalty(string dbSeq,
                                 int firstDbIndex,
                                 int lastDbIndex,
                                 bool isPepStart,
                                 bool isPepEnd,
                                 bool reverse,
                                 bool debugFlag);

    void getGapAnnotation(float    specGapLengthFloat,
                          string & aaString,
                          string & stringAnnotation,
                          int      ntermOrCterm,
                          bool     knownModsOnly);

    float getGapPenalty(int      specGapLength,
                        string   aaString,
                        float    avgPeakIntensity,
                        int      ntermOrCterm,
                        bool     knownModsOnly);

    void setDebugFlags(bool debugRange,
                       bool debugRange2,
                       bool debugAas,
                       bool debugSpecs,
                       bool debugAlign,
                       bool debugAlign1,
                       bool debugAlign2,
                       bool debugAlign3,
                       bool debugAlign4,
                       bool debugAlignExact,
                       bool debugAlignNterm,
                       bool debugAlignCterm,
                       bool debugGap,
                       bool debugGapAnno,
                       bool debugCache,
                       bool debugLarge);

    int 			      m_cacheHits;
    int 			      m_cacheMisses;
    int 			      m_cacheHitsTotal;
    int 			      m_cacheMissesTotal;
  private:
    int 			      m_kMer;
    int 			      m_maxSpecGap;
    PenaltyMatrix * m_modPenaltyMatrix;
    PenaltyMatrix * m_blossumPenaltyMatrix;
    float           m_penaltyAlpha;
    float           m_penaltyBeta;
    float           m_maxMod;
    float           m_minMod;
    map<string, vector<float> > m_mapGapVectors;
    map<string, vector<float> > m_mapCachedGapVectors;

    void createGapVectors();

    void combineVectors(vector<float> & scores1,
                        vector<float> & scores2,
                        vector<float> & outputScores);
			
    float combineVectorsSingleValue(vector<float> & scores1,
                                    vector<float> & scores2,
                                    int             specGapLength);
                        

    void combineVectorsWithAnnotation(vector<float>  & scores1,
                                      vector<float>  & scores2,
                                      vector<float>  & outputScores,
                                      vector<string> & annotationsOld,
                                      vector<string> & annotationsNew,
                                      string         & newAA,
                                      int            & newAAMass);

    void fillVectorSingle(vector<float> & scores,
                          string          strAA,
                          int             ntermOrCterm,
                          bool            knownModsOnly);

    float getPeptideMass(string & stringPeptide);

    void initializeGapVector(vector<float> & scores);

    string makeAnnotation(string & strAA, int mass);

    void makeDbString(char * dbGapStringOut, char * dbSeq, int start, int end);

    string massToString(float mass);

    void setStartingFlagArray(Spectrum &              spec,
                              Spectrum &              dbSpec,
                              set<float> &            startRange,
                              vector<vector<char> > & startFlags,
                              int                   & firstValidDbIndex,
                              float                   tolerance);
  };

}

#endif

