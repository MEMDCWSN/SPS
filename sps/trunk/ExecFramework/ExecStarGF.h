/*
 * ExecStarGF.h
 *
 *  Created on: Nov 15, 2013
 *      Author: aguthals
 */

#ifndef EXECSTARGF_H_
#define EXECSTARGF_H_

#include "ExecBase.h"
#include "SpecSet.h"
#include "PeptideSpectrumMatchSet.h"
#include "ClusterSet.h"
#include "SpectrumPairSet.h"
#include "db_fasta.h"

using namespace std;

namespace specnets
{
  class ExecStarGF : public ExecBase
  {
  public:

    static const int DEBUG_SCAN1;
    static const int DEBUG_SCAN2;

    struct PeptideStarMatch
    {
      int specIdx;
      int scan;
      int useYEndpts;
      int pepScore;
      int charge;
      // int dictScore;
      int bestIntersectPVIdx;
      // int bestIntersectDictIdx;
      double specProb;
      // double optimalSpecProb;
      double starProb;
      double fdr;
      double pepFdr;
      double pairSpecProb;
      double pairStarProb;
      double msgfSpecProb;
      float alignGFProb;
      int isDecoy;
      int pairIsDecoy;
      int rank;
      int pairRank;
      int pairMP;
      // double approxStarProb;
      // double approxStarProb2;
      // long long dictionarySz;
      // long long intersectDictSz;
      string filename;
      string peptide;
      string bestPairPeptide;
      string peptideOvlp;
      string protein;
      float precursor;
      int bestPairScan;
      string bestPairFile;

      PeptideStarMatch() :
          specIdx(-1), pepScore(-1), bestIntersectPVIdx(-1), peptide(""),
              protein(""), charge(0), peptideOvlp(""), useYEndpts(0),
              specProb(1), starProb(1), filename(""), scan(-1), precursor(-1),
              bestPairPeptide(""), bestPairScan(-1), bestPairFile(""),
              fdr(-1.0), pepFdr(-1.0), isDecoy(-1), pairIsDecoy(-1),
              pairSpecProb(-1), pairStarProb(-1), msgfSpecProb(1), rank(-1), pairRank(-1), pairMP(-1), alignGFProb(1)
      //dictScore(0), bestIntersectDictIdx(-1), optimalSpecProb(1), approxStarProb2(1), approxStarProb(1), dictionarySz(0), intersectDictSz(0)
      {
      }

      void setToPSM(PeptideSpectrumMatch &other)
      {
        other.m_annotation = peptide;
        other.m_scanNum = scan;
        other.m_protein = protein;
        other.m_fdr = fdr;
        other.m_pepFdr = pepFdr;
        other.m_charge = charge;
        other.m_pValue = starProb;
        other.m_score = pepScore;
        other.m_isDecoy = isDecoy == 1;
      }

      PeptideStarMatch &operator=(const PeptideStarMatch &other)
      {
        if (this == &other)
        {
          return *this;
        }
        specIdx = other.specIdx;
        peptide = other.peptide;
        protein = other.protein;
        filename = other.filename;
        scan = other.scan;
        precursor = other.precursor;
        charge = other.charge;
        isDecoy = other.isDecoy;
        pairIsDecoy = other.pairIsDecoy;
        //dictScore = other.dictScore;
        //approxStarProb2 = other.approxStarProb2;
        peptideOvlp = other.peptideOvlp;
        useYEndpts = other.useYEndpts;
        pepScore = other.pepScore;
        bestIntersectPVIdx = other.bestIntersectPVIdx;
        pairSpecProb = other.pairSpecProb;
        pairStarProb = other.pairStarProb;
        //bestIntersectDictIdx = other.bestIntersectDictIdx;
        specProb = other.specProb;
        starProb = other.starProb;
        //optimalSpecProb = other.optimalSpecProb;
        bestPairPeptide = other.bestPairPeptide;
        //bestPairIdx = other.bestPairIdx;
        bestPairScan = other.bestPairScan;
        bestPairFile = other.bestPairFile;
        fdr = other.fdr;
        pepFdr = other.pepFdr;
        msgfSpecProb = other.msgfSpecProb;
        rank = other.rank;
        pairRank = other.pairRank;
        pairMP = other.pairMP;
        alignGFProb = other.alignGFProb;
        //approxStarProb = other.approxStarProb;
        //dictionarySz = other.dictionarySz;
        //intersectDictSz = other.intersectDictSz;
        return *this;
      }
    };

    static bool SavePeptideStarMatches(const string &filename,
                                       const vector<PeptideStarMatch> &outputPSMs);

    ExecStarGF(void);

    ExecStarGF(const ParameterList & inputParams);

    ExecStarGF(const ParameterList & inputParams,
               SpecSet * inputSpectra,
               PeptideSpectrumMatchSet * inputPSMs,
               ClusterSet *inputClusters,
               SpectrumPairSet * inputPairs,
               SpectrumPairSet * inputAlignedPairs,
               vector<vector<double> > *inputProbs,
               DB_fasta *inputDB);

    ExecStarGF(const ParameterList & inputParams,
               SpecSet * inputSpectra,
               PeptideSpectrumMatchSet * inputPSMs,
               ClusterSet *inputClusters,
               SpectrumPairSet * inputPairs,
               SpectrumPairSet * inputAlignedPairs,
               vector<vector<double> > *inputProbs,
               DB_fasta *inputDB,
               SpecSet *outputSpectra,
               vector<PeptideStarMatch> *outputPSMs);

    virtual ~ExecStarGF(void);

    virtual ExecBase * clone(const ParameterList & inputParams) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

  protected:

    bool ownInput;
    bool ownOutput;

    SpecSet *m_inputSpectra;
    PeptideSpectrumMatchSet *m_inputPSMs;
    SpectrumPairSet *m_inputPairs;
    SpectrumPairSet *m_inputAlignedPairs;
    ClusterSet *m_inputClusters;
    vector<vector<double> > *m_inputProbs;
    DB_fasta *m_inputDB;

    SpecSet *m_outputSpectra;
    vector<PeptideStarMatch> *m_outputPSMs;

    bool saveUnidentifiedMS2(const std::set<int> & identified_psm_scans);
  };
}

#endif /* EXECSTARGF_H_ */
