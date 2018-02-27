#include "PeptideSpectrumMatchNetwork.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchNetwork::PeptideSpectrumMatchNetwork()
  {
    m_spectrumFile = "";
    m_scanNum = -1;
    m_annotation = "";
    m_protein = "";
    m_charge = -1;
    m_score = 0;
    m_pValue = 0;
    m_isDecoy = false;
    m_spectra.resize(0);
    m_partialAnnotations.resize(0);
    m_partialPeptideScores.resize(0);
    m_partialSpecProbabilities.resize(0);
    m_spectrumIndices.resize(0);
    m_individualScore = -9;
    m_isTuple = true;
    m_parentMasses.resize(0);

  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchNetwork::~PeptideSpectrumMatchNetwork()
  {
    //EMPTY
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchNetwork::PeptideSpectrumMatchNetwork(const PeptideSpectrumMatchNetwork &other)
  {
    m_spectrumFile = other.m_spectrumFile;
    m_annotation = other.m_annotation;
    m_protein = other.m_protein;
    m_charge = other.m_charge;
    m_score = other.m_score;
    m_pValue = other.m_pValue;
    m_isDecoy = other.m_isDecoy;
    m_individualScore = other.m_individualScore;
    m_isTuple = other.m_isTuple;


    // cboucher: changed to copy vectors
    m_spectra = other.m_spectra;
    m_partialAnnotations = other.m_partialAnnotations;
    m_partialPeptideScores = other.m_partialPeptideScores;
    m_partialSpecProbabilities = other.m_partialSpecProbabilities;
    m_spectrumIndices = other.m_spectrumIndices;
    m_parentMasses = other.m_parentMasses;

  }

  // -------------------------------------------------------------------------
  PeptideSpectrumMatchNetwork &PeptideSpectrumMatchNetwork::operator=(const PeptideSpectrumMatchNetwork &other)
  {
    m_spectrumFile = other.m_spectrumFile;
    m_annotation = other.m_annotation;
    m_protein = other.m_protein;
    m_charge = other.m_charge;
    m_score = other.m_score;
    m_pValue = other.m_pValue;
    m_isDecoy = other.m_isDecoy;
    m_individualScore = other.m_individualScore;
       m_isTuple = other.m_isTuple;

    
    // cboucher: changed to copy vectors
    m_spectra = other.m_spectra;
    m_partialAnnotations = other.m_partialAnnotations;
      m_partialPeptideScores = other.m_partialPeptideScores;
      m_partialSpecProbabilities = other.m_partialSpecProbabilities;
      m_spectrumIndices = other.m_spectrumIndices;
      m_parentMasses = other.m_parentMasses;


    return (*this);
  }


  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchNetwork::setAnnotationToMatches(vector<int> &matches,
                                                    vector<pair<
                                                        const ftIonFragment*,
                                                        short> > &annotation,
                                                    int ionIdx,
                                                    const ftIonFragment* currIonFrag)
  {
	  DEBUG_MSG("Not implemented for this class::should not use");
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchNetwork::setAnnotationToMatchesNoDups(vector<int> &matches,
                                                          vector<
                                                              pair<
                                                                  const ftIonFragment*,
                                                                  short> > &annotation,
                                                          int ionIdx,
                                                          const ftIonFragment* currIonFrag)
  {
	  DEBUG_MSG("Not implemented for this class::should not use");
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchNetwork::annotate(string &peptide,
                                      string &ionNamesInclude,
                                      MS2ScoringModel &inputIonTypes,
                                      float prmOffset,
                                      float srmOffset,
                                      float peakTol,
                                      bool removeDuplicates,
                                      bool ignoreParentCharge)
  {
	  DEBUG_MSG("Not implemented for this class::should not use");
	  return false;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchNetwork::annotate(string &peptide,
                                      string &ionNamesInclude,
                                      MS2ScoringModel &inputIonTypes,
                                      float prmOffset,
                                      float srmOffset,
                                      float peakTol,
                                      AAJumps &jumps,
                                      bool removeDuplicates,
                                      bool ignoreParentCharge)
  {
	DEBUG_MSG("Not implemented for this class::should not use");
    return false;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchNetwork::mapIons(vector<vector<int> > &outputMatches,
                                     std::tr1::unordered_map<string, int> &ionMap)
  {
	  DEBUG_MSG("Not implemented for this class::should not use");
  }

}
