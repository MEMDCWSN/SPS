// Header Include
#include "ExecFdrPeptide.h"

static bool DEBUG_FDR_PEPTIDE = true;

using namespace std;

namespace specnets
{

  // -------------------------------------------------------------------------
  bool comparePValue(psmPtr i, psmPtr j)
  {
    return i->m_pValue < j->m_pValue;
  }

  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(void) :
    m_inputPeptides(0x0), m_inputTargets(0x0), m_inputDecoys(0x0),
        m_fdrResults(0x0), ownInput(true), ownOutput(true)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
    m_scalingFactor = 1.0;
  }

  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputPeptides(0x0), m_inputTargets(0x0),
        m_inputDecoys(0x0), m_fdrResults(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
    m_scalingFactor = 1.0;
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams,
                                 PeptideSpectrumMatchSet * inputPeptides,
                                 PeptideSpectrumMatchSet * fdrResults) :
    ExecBase(inputParams), m_inputPeptides(inputPeptides),
        m_inputTargets(0x0), m_inputDecoys(0x0), m_fdrResults(fdrResults),
        ownInput(false), ownOutput(false)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
    m_scalingFactor = 1.0;
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams,
                                 PeptideSpectrumMatchSet * inputTargets,
                                 PeptideSpectrumMatchSet * inputDecoys,
                                 PeptideSpectrumMatchSet * fdrResults) :
    ExecBase(inputParams), m_inputTargets(inputTargets),
        m_inputPeptides(0x0), m_inputDecoys(inputDecoys),
        m_fdrResults(fdrResults), ownInput(false), ownOutput(false)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
    m_scalingFactor = 1.0;
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::~ExecFdrPeptide(void)
  {
    if (ownInput)
    {
      if (m_inputPeptides)
      {
        delete m_inputPeptides;
        m_inputPeptides = 0x0;
      }
      if (m_inputTargets)
      {
        delete m_inputTargets;
        m_inputTargets = 0x0;
      }
      if (m_inputDecoys)
      {
        delete m_inputDecoys;
        m_inputDecoys = 0x0;
      }
    }
    if (ownOutput)
    {
      if (m_fdrResults)
      {
        delete m_fdrResults;
        m_fdrResults = 0x0;
      }
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFdrPeptide::clone(const ParameterList & inputParams) const
  {
    return new ExecFdrPeptide(inputParams);
  }
  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::invoke(void)
  {
    DEBUG_TRACE;
    if (m_fdrResults == 0x0)
    {
      ownOutput = true;
      m_fdrResults = new PeptideSpectrumMatchSet;
    }

    DEBUG_FDR_PEPTIDE = m_params.getValueBool("DEBUG_FDR_PEPTIDE");
    DEBUG_VAR(DEBUG_FDR_PEPTIDE);
    
    m_tdaType = m_params.getValue("TDA_TYPE");
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_tdaType);
    if (m_tdaType != "concatenated" && 
        m_tdaType != "separate" && 
        m_tdaType != "split" && m_tdaType != "split_mods" && 
        m_tdaType != "split_charge" && m_tdaType != "split_tryptic" &&
        m_tdaType != "split_moda")
    {
      ERROR_MSG("Unknown TDA_TYPE! [" << m_tdaType << "] Valid options are: ");
      ERROR_MSG("    'concatenated', 'separate', 'split', 'split_mods', 'split_charge', 'split_tryptic', 'split_moda'.");
      return false;
    }  
    
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("TRUE_FALSE_DATABASE_RATIO"));
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("PEPTIDE_FDR_AUTOSCALE"));
    if (m_params.exists("TRUE_FALSE_DATABASE_RATIO")) {
      m_scalingFactor = m_params.getValueDouble("TRUE_FALSE_DATABASE_RATIO");
    } else if (m_params.exists("PEPTIDE_FDR_AUTOSCALE")) {
      if (m_inputDecoys->size() != 0) {
        m_scalingFactor = (float)m_inputTargets->size() / (float) m_inputDecoys->size();
      }
    } 
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_scalingFactor);

    /*
    if ((m_inputPeptides != 0x0 && m_inputPeptides->size() == 0) ||
        (m_inputTargets != 0x0 && m_inputTargets->size() == 0 &&
         m_inputDecoys != 0x0 && m_inputDecoys->size() == 0)) {
      DEBUG_MSG("No input peptides");


    } else
*/
      if (m_tdaType == "split" || m_tdaType == "split_mods" ||
               m_tdaType == "split_charge" || m_tdaType == "split_tryptic" ||
               m_tdaType == "split_moda") {
      if (DEBUG_FDR_PEPTIDE) DEBUG_MSG("Doing Split FDR");
      if (!performSplitFdr()) {
        return false;
      }
    } else if (m_inputPeptides != 0x0) {
      if (DEBUG_FDR_PEPTIDE) DEBUG_MSG("Doing Single List FDR");
      if (!performFdrSingleList(*m_inputPeptides,
                                *m_fdrResults,
                                m_tdaType,
                                m_scalingFactor)) {
          return false;
       }
    } else {
      if (DEBUG_FDR_PEPTIDE) DEBUG_MSG("Doing Two List FDR");
      if (!performFdrTwoLists(*m_inputTargets,
                              *m_inputDecoys,
                              *m_fdrResults,
                              m_tdaType,
                              m_scalingFactor)) {
          return false;
       }
    } 
    if (DEBUG_FDR_PEPTIDE) m_fdrResults->saveToFile("debug_fdr_results.txt", true, true);

    if (m_tdaType != "split" && m_tdaType != "split_mods" && 
        m_tdaType != "split_charge" && m_tdaType != "split_tryptic" &&
        m_tdaType != "split_moda") {
      if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("PEPTIDE_FDR_CUTOFF"));
      if (m_params.exists("PEPTIDE_FDR_CUTOFF")) {
        if (DEBUG_FDR_PEPTIDE) DEBUG_TRACE;
        double cutoff = m_params.getValueDouble("PEPTIDE_FDR_CUTOFF");
        if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(cutoff);
        if (!FdrPeptide::filterByPValue(*m_fdrResults, cutoff)) {
          return false;
        }
        if (DEBUG_FDR_PEPTIDE) m_fdrResults->saveToFile("debug_fdr_cutoff.txt", true, true);
      }
    }

    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("PEPTIDE_FDR_REMOVE_DECOYS"));
    if (m_params.exists("PEPTIDE_FDR_REMOVE_DECOYS"))
    {
      if (!FdrPeptide::removeDecoys(*m_fdrResults))
      {
        return false;
      }
      if (DEBUG_FDR_PEPTIDE) m_fdrResults->saveToFile("debug_fdr_nodecoy.txt", true, true);
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::performFdrSingleList(PeptideSpectrumMatchSet & inputPeptides,
                                            PeptideSpectrumMatchSet & outputPeptides,
                                            string                  & tdaType,
                                            double                    scalingFactor)
  {
    if (DEBUG_FDR_PEPTIDE) inputPeptides.saveToFile("debug_fdr_input.txt", true, true);

    // No input peptides is not strictly an error.
    if (inputPeptides.size() == 0) {
      DEBUG_MSG("No input peptides to performFdrSingleList");
      return true;
    }

    int fdrType = m_params.getValueInt("FDR_TYPE", FDR_TYPE_SPECTRUM);
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(fdrType);

    bool usePvalueSort = m_params.getValueBool("PEPTIDE_FDR_USE_PVALUE_SORT");
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(usePvalueSort);

    bool usePvalueReplace = m_params.getValueBool("PEPTIDE_FDR_USE_PVALUE_REPLACE");
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(usePvalueReplace);

    if (usePvalueReplace) {
      fdrType |= FDR_REPLACE_PVALUE;
    } else {
      fdrType |= FDR_REPLACE_SCORE;
    }

    if (usePvalueSort) {
      fdrType |= FDR_SORT_PVALUE;
    } else {
      fdrType |= FDR_SORT_SCORE;
    }
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(fdrType);

    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("TOLERANCE_PM"));
    float pmTol = m_params.getValueFloat("TOLERANCE_PM");
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(pmTol);

    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(tdaType);
    if (DEBUG_FDR_PEPTIDE) DEBUG_MSG("fdrType = " << std::hex << fdrType);

    bool proteinFilter = m_params.exists("FILTER_BY_PROTEIN_FDR");
    DEBUG_VAR(proteinFilter);

    if (proteinFilter) {
      filterByProteinFdr(inputPeptides,
                         tdaType,
                         fdrType,
                         scalingFactor);
    }

    if (tdaType == "concatenated") {
      if (!FdrPeptide::concatenatedTargetDecoy(inputPeptides,
                                               outputPeptides,
                                               scalingFactor,
                                               fdrType)) {
        ERROR_MSG("FdrPeptide::concatenatedTargetDecoy failed!");
        return false;
      }
    }
    else if (tdaType == "separate")
    {
      if (!FdrPeptide::separateTargetDecoy(inputPeptides,
                                           outputPeptides,
                                           scalingFactor,
                                           fdrType)) {
        ERROR_MSG("FdrPeptide::separateTargetDecoy failed!");
        return false;
      }
    }
    if (DEBUG_FDR_PEPTIDE) outputPeptides.saveToFile("debug_fdr_output.txt", true, true);
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::performFdrTwoLists(PeptideSpectrumMatchSet & targetPeptides,
                                          PeptideSpectrumMatchSet & decoyPeptides,
                                          PeptideSpectrumMatchSet & outputPeptides,
                                          string                  & tdaType,
                                          double                    scalingFactor)
  {
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(targetPeptides.size());
    if (DEBUG_FDR_PEPTIDE) targetPeptides.saveToFile("debug_fdr_targets.txt", true, true);
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(decoyPeptides.size());
    if (DEBUG_FDR_PEPTIDE) decoyPeptides.saveToFile("debug_fdr_decoys.txt", true, true);

    // Make sure the isDecoy fields are set on the decoy PSMs
    for (int iDecoy = 0; iDecoy < decoyPeptides.size(); iDecoy++) {
      decoyPeptides[iDecoy]->m_isDecoy = true;
    }      

    PeptideSpectrumMatchSet mergedPeptides;
    if (DEBUG_FDR_PEPTIDE) DEBUG_TRACE;
    for (int i = 0; i < targetPeptides.m_psmSet.size(); i++) {
      mergedPeptides.m_psmSet.push_back(targetPeptides.m_psmSet[i]);
    }
    for (int i = 0; i < decoyPeptides.m_psmSet.size(); i++) {
      mergedPeptides.m_psmSet.push_back(decoyPeptides.m_psmSet[i]);
    }
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(mergedPeptides.size());
    if (DEBUG_FDR_PEPTIDE) mergedPeptides.saveToFile("debug_fdr_merged.txt", true, true);

    if (!performFdrSingleList(mergedPeptides,
                              outputPeptides,
                              tdaType,
                              scalingFactor)) {
      return false;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::filterByProteinFdr(PeptideSpectrumMatchSet & inputPeptides,
                                          string                  & tdaType,
                                          int                       fdrType,
                                          double                    scalingFactor)
  {
    DEBUG_VAR(inputPeptides.size());

    int fdrProteinType = fdrType & ~FDR_TYPE_MASK;
    fdrProteinType |= FDR_TYPE_PROTEIN;
    fdrProteinType = fdrProteinType & ~FDR_REPLACE_MASK;
    fdrProteinType |= FDR_REPLACE_PVALUE;
    fdrProteinType = fdrProteinType & ~FDR_SORT_MASK;
    fdrProteinType |= FDR_SORT_PVALUE;
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(fdrType);

    PeptideSpectrumMatchSet unmodedPeptides;
    for (int i = 0; i < inputPeptides.m_psmSet.size(); i++) {
      vector<float> modifications;
      inputPeptides.m_psmSet[i]->getModifications(modifications);
      int numMods = modifications.size();
      if (numMods == 0) {
        unmodedPeptides.push_back(inputPeptides.m_psmSet[i]);
      }
    }
    if (DEBUG_FDR_PEPTIDE) unmodedPeptides.saveToFile("debug_unmoded.txt", true, true);
    DEBUG_VAR(unmodedPeptides.size());

    PeptideSpectrumMatchSet proteinPeptides;
    if (!FdrPeptide::concatenatedTargetDecoy(unmodedPeptides,
                                             proteinPeptides,
                                             scalingFactor,
                                             fdrProteinType)) {
      ERROR_MSG("FdrPeptide::concatenatedTargetDecoy for Proteins failed!");
      return false;
    }

    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("PEPTIDE_FDR_CUTOFF"));
    if (m_params.exists("PEPTIDE_FDR_CUTOFF")) {
      if (DEBUG_FDR_PEPTIDE) DEBUG_TRACE;
      double cutoff = m_params.getValueDouble("PEPTIDE_FDR_CUTOFF");
      if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(cutoff);
      if (!FdrPeptide::filterByPValue(proteinPeptides, cutoff)) {
        return false;
      }
      if (DEBUG_FDR_PEPTIDE) proteinPeptides.saveToFile("debug_protein_fdr.txt", true, true);
    }
    DEBUG_VAR(proteinPeptides.size());

    set<string> proteins;
    for (int i= 0; i < proteinPeptides.m_psmSet.size(); i++) {
      if (!proteinPeptides[i]->m_isDecoy) {
        proteins.insert(proteinPeptides[i]->m_protein);
        string decoyProtein = "XXX" + proteinPeptides[i]->m_protein;
        proteins.insert(decoyProtein);
      } else {
        proteins.insert(proteinPeptides[i]->m_protein);
      }
    }
/*
    DEBUG_VAR(proteins.size());
    for (int i= 0; i < unmodedPeptides.m_psmSet.size(); i++) {
      if (unmodedPeptides[i]->m_isDecoy) {
        proteins.insert(unmodedPeptides[i]->m_protein);
      }
    }
    DEBUG_VAR(proteins.size());
*/

    PeptideSpectrumMatchSet filteredSet;
    for (int i= 0; i < inputPeptides.m_psmSet.size(); i++) {
      if (proteins.find(inputPeptides[i]->m_protein) !=  proteins.end()) {
        filteredSet.push_back(inputPeptides[i]);
      }
    }
    if (DEBUG_FDR_PEPTIDE) filteredSet.saveToFile("debug_filter_by_protein.txt", true, true);

    inputPeptides.m_psmSet.clear();
    for (int i= 0; i < filteredSet.m_psmSet.size(); i++) {
      inputPeptides.push_back(filteredSet[i]);
    }
    if (DEBUG_FDR_PEPTIDE) filteredSet.saveToFile("debug_filter_by_protein_copy.txt", true, true);

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::performSplitFdr(void)
  {
    map<int, PeptideSpectrumMatchSet > mapPsmTargetSets;
    map<int, PeptideSpectrumMatchSet > mapPsmDecoySets;

    DEBUG_TRACE;

    if (m_inputPeptides != 0x0 && m_inputPeptides->size() != 0) {
      DEBUG_VAR(m_inputPeptides->m_psmSet.size());
      for (int i= 0; i < m_inputPeptides->m_psmSet.size(); i++) {

        int iSplit = 0;
        if (m_tdaType == "split" || m_tdaType == "split_mods") {
          vector<float> modifications;
          m_inputPeptides->m_psmSet[i]->getModifications(modifications);
          int numMods = modifications.size();
          if (numMods > 0) {
            iSplit = 1;
          }
        } else if (m_tdaType == "split_charge") {
          iSplit = m_inputPeptides->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
        } else if (m_tdaType == "split_tryptic") {

          string anno = m_inputPeptides->m_psmSet[i]->m_annotation;
          //DEBUG_VAR(anno);
          //DEBUG_VAR(m_inputPeptides->m_psmSet[i]->m_precedingAA);
          if (m_inputPeptides->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputPeptides->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputPeptides->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          if (m_inputPeptides->m_psmSet[i]->m_annotation.size() != 0) {
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(anno, annoClean);
            //DEBUG_VAR(annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }
          //DEBUG_VAR(iSplit);

        } else if (m_tdaType == "split_moda") {

          iSplit = m_inputPeptides->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
          iSplit *= 3;
          if (m_inputPeptides->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputPeptides->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputPeptides->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          int length = m_inputPeptides->m_psmSet[i]->m_annotation.size();
          if (length != 0) {
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(
                            m_inputPeptides->m_psmSet[i]->m_annotation, 
                            annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }

        }
        
        if (m_inputPeptides->m_psmSet[i]->m_isDecoy) {
          mapPsmDecoySets[iSplit].m_psmSet.push_back(m_inputPeptides->m_psmSet[i]);
        } else {
          mapPsmTargetSets[iSplit].m_psmSet.push_back(m_inputPeptides->m_psmSet[i]);
        }

      }
      
      DEBUG_VAR(mapPsmTargetSets.size());
      for (int i = 0; i < mapPsmTargetSets.size(); i++) {
        DEBUG_VAR(mapPsmTargetSets[i].size());
      }
      int totalDebug = 0;
      DEBUG_VAR(mapPsmDecoySets.size());
      for (int i = 0; i < mapPsmDecoySets.size(); i++) {
        DEBUG_VAR(mapPsmDecoySets[i].size());
        totalDebug += mapPsmDecoySets[i].size();
      }
      
      if (totalDebug == 0) {
        WARN_MSG("No decoys found");
      }
      
    } else {

      DEBUG_VAR(m_tdaType);
      DEBUG_VAR(m_inputTargets->m_psmSet.size());
      for (int i= 0; i < m_inputTargets->m_psmSet.size(); i++) {

        int iSplit = 0;
        if (m_tdaType == "split" || m_tdaType == "split_mods") {
          vector<float> modifications;
          m_inputTargets->m_psmSet[i]->getModifications(modifications);
          int numMods = modifications.size();
          if (numMods > 0) {
            iSplit = 1;
          }
        } else if (m_tdaType == "split_charge") {
          iSplit = m_inputTargets->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
        } else if (m_tdaType == "split_tryptic") {
          if (m_inputTargets->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputTargets->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputTargets->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          if (m_inputTargets->m_psmSet[i]->m_annotation.size() != 0) {
            string anno = m_inputTargets->m_psmSet[i]->m_annotation;
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(anno, annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }
        } else if (m_tdaType == "split_moda") {
          iSplit = m_inputTargets->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
          iSplit *= 3;
          if (m_inputTargets->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputTargets->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputTargets->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          if (m_inputTargets->m_psmSet[i]->m_annotation.size() != 0) {
            string anno = m_inputTargets->m_psmSet[i]->m_annotation;
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(anno, annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }
        }
        
        mapPsmTargetSets[iSplit].m_psmSet.push_back(m_inputTargets->m_psmSet[i]);
      }
      DEBUG_VAR(mapPsmTargetSets.size());

      for (int i= 0; i < m_inputDecoys->m_psmSet.size(); i++) {

        int iSplit = 0;
        if (m_tdaType == "split" || m_tdaType == "split_mods") {
          vector<float> modifications;
          m_inputDecoys->m_psmSet[i]->getModifications(modifications);
          int numMods = modifications.size();
          if (numMods > 0) {
            iSplit = 1;
          }
        } else if (m_tdaType == "split_charge") {
          iSplit = m_inputDecoys->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
        } else if (m_tdaType == "split_tryptic") {
          if (m_inputDecoys->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputDecoys->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputDecoys->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          int length = m_inputDecoys->m_psmSet[i]->m_annotation.size();
          if (length != 0) {
            string anno = m_inputDecoys->m_psmSet[i]->m_annotation;
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(anno, annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }
        } else if (m_tdaType == "split_moda") {
          iSplit = m_inputDecoys->m_psmSet[i]->m_charge;
          if (iSplit > 4) iSplit = 4; // Max bin at charge state 4
          iSplit *= 3;
          if (m_inputDecoys->m_psmSet[i]->m_precedingAA == '-' ||
              m_inputDecoys->m_psmSet[i]->m_precedingAA == 'K' ||
              m_inputDecoys->m_psmSet[i]->m_precedingAA == 'R') {
            iSplit++;
          }
          int length = m_inputDecoys->m_psmSet[i]->m_annotation.size();
          if (length != 0) {
            string anno = m_inputDecoys->m_psmSet[i]->m_annotation;
            string annoClean;
            PeptideSpectrumMatch::getUnmodifiedPeptide(anno, annoClean);
            if (annoClean[annoClean.size() - 1] == 'K' ||
                annoClean[annoClean.size() - 1] == 'R') {
              iSplit++;
            }
          }
        }

        mapPsmDecoySets[iSplit].m_psmSet.push_back(m_inputDecoys->m_psmSet[i]);
      }

      DEBUG_VAR(mapPsmDecoySets.size());
      for (int i = 0; i < mapPsmTargetSets.size(); i++) {
        DEBUG_VAR(mapPsmTargetSets[i].size());
      }
      for (int i = 0; i < mapPsmDecoySets.size(); i++) {
        DEBUG_VAR(mapPsmDecoySets[i].size());
      }
      
    }

    map<int, PeptideSpectrumMatchSet >::iterator itrMap = mapPsmTargetSets.begin();
    map<int, PeptideSpectrumMatchSet >::iterator itrMapEnd = mapPsmTargetSets.end();
    for ( ; itrMap != itrMapEnd; itrMap++) {
      PeptideSpectrumMatchSet & targets = itrMap->second;
      PeptideSpectrumMatchSet & decoys = mapPsmDecoySets[itrMap->first];

      stringstream sst;
      stringstream ssd;
      sst << "debug_fdr_target_" << itrMap->first << ".txt";
      ssd << "debug_fdr_decoy_" << itrMap->first << ".txt";
      if (DEBUG_FDR_PEPTIDE) targets.saveToFile(sst.str().c_str());
      if (DEBUG_FDR_PEPTIDE) decoys.saveToFile(ssd.str().c_str());

      DEBUG_MSG(itrMap->first << "  " << targets.size() <<  "  " << decoys.size());
      if (decoys.size() == 0) {
        // Do nothing because we can't be sure of targets if we have no decoys
      } else {

        // Make sure the isDecoy fields are set on the decoy PSMs
        for (int iDecoy = 0; iDecoy < decoys.size(); iDecoy++) {
          decoys[iDecoy]->m_isDecoy = true;
        }

        float scalingFactor = m_scalingFactor;
        if (m_params.exists("PEPTIDE_FDR_AUTOSCALE")) {
          scalingFactor = (float)targets.size() / (float)decoys.size();
        } 
        PeptideSpectrumMatchSet partialResult;
        string tdaType("concatenated");
        if (!performFdrTwoLists(targets,
                                decoys,
                                partialResult,
                                tdaType,
                                scalingFactor)) {
          ERROR_MSG("performFdrTwoLists failed for split [" << itrMap->first << "]");
        }

        if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(m_params.getValue("PEPTIDE_FDR_CUTOFF"));
        if (m_params.exists("PEPTIDE_FDR_CUTOFF")) {
          double cutoff = m_params.getValueDouble("PEPTIDE_FDR_CUTOFF");
          if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(cutoff);
          if (!FdrPeptide::filterByPValue(partialResult, cutoff)) {
            ERROR_MSG("filterByPValue failed for split [" << itrMap->first << "]");
            return false;
          }
          if (DEBUG_FDR_PEPTIDE) m_fdrResults->saveToFile("debug_fdr_cutoff.txt", true, true);
        }

        DEBUG_VAR(partialResult.m_psmSet.size());
        for (int i= 0; i < partialResult.m_psmSet.size(); i++) {
          m_fdrResults->m_psmSet.push_back(partialResult.m_psmSet[i]);
        }
      }
    }    
    DEBUG_VAR(m_fdrResults->m_psmSet.size());
    if (DEBUG_FDR_PEPTIDE) m_fdrResults->saveToFile("debug_fdr_allsplit.txt", true, true);

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::loadInputData(void)
  {
    // Load in statistics if they haven't been passed in by another function
    // LEB: I'm leaving this as it was because I don't know the effect of
    //      changing it on all code that exists.
    if (ownInput) {
      if (!m_inputTargets)
        m_inputTargets = new PeptideSpectrumMatchSet;
      if (!m_inputDecoys)
        m_inputDecoys = new PeptideSpectrumMatchSet;
    }

    if (ownOutput) {
      if (!m_fdrResults) {
        m_fdrResults = new PeptideSpectrumMatchSet;
      }
    }
    
    //Load in inspect / specnets results
    if (m_params.exists("INPUT_RESULTS"))
    {
      if (ownInput) {
        if (!m_inputPeptides)
          m_inputPeptides = new PeptideSpectrumMatchSet;
      }
      DEBUG_MSG("Input results: " << m_params.getValue("INPUT_RESULTS"));

      if (m_params.getValue("INPUT_RESULTS_TYPE").compare("inspect") == 0)
      {
        if (!m_inputPeptides->loadInspectResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                      m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
      else
      {
        if (!m_inputPeptides->loadSpecnetsResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                       m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
      DEBUG_VAR(m_inputPeptides->size());
    }

    if (m_params.exists("INPUT_TARGET")) {
      if (ownInput) {
        if (!m_inputTargets)
          m_inputTargets = new PeptideSpectrumMatchSet;
      }
      DEBUG_MSG("Input Target: " << m_params.getValue("INPUT_TARGET"));
      if (!m_inputTargets->loadFromFile(m_params.getValue("INPUT_TARGET").c_str())) {
        DEBUG_MSG("Could not load" << m_params.getValue("INPUT_TARGET"));
        return false;
      }
      DEBUG_VAR(m_inputTargets->size());
    }
    if (m_params.exists("INPUT_DECOY")) {
      if (ownInput) {
        if (!m_inputDecoys)
          m_inputDecoys = new PeptideSpectrumMatchSet;
      }
      DEBUG_MSG("Input Decoy: " << m_params.getValue("INPUT_DECOY"));
      if (!m_inputDecoys->loadFromFile(m_params.getValue("INPUT_DECOY").c_str())) {
        DEBUG_MSG("Could not load" << m_params.getValue("INPUT_DECOY"));
        return false;
      }
      DEBUG_VAR(m_inputDecoys->size());
    }

    return true;
    DEBUG_TRACE;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_FDR_RESULTS")) {
      string fileName = m_params.getValue("OUTPUT_FDR_RESULTS");
      if (m_params.exists("OUTPUT_FDR_RESULTS_DIR")) {
        string outDir = m_params.getValue("OUTPUT_FDR_RESULTS_DIR");
        fileName = outDir + "/" + fileName;
      }
      DEBUG_MSG("Saving FDR results to [" << fileName << "]");
      if (!m_fdrResults->saveToFile(fileName.c_str(), true)) {
        return false;
      }
    }
    
    if (m_params.exists("OUTPUT_PSM_MOD_FREQS")) {
      string fileName = m_params.getValue("OUTPUT_PSM_MOD_FREQS");
      if (m_params.exists("OUTPUT_PSM_MODS_DIR")) {
        string outDir = m_params.getValue("OUTPUT_PSM_MODS_DIR");
        fileName = outDir + "/" + fileName;
      }
      DEBUG_MSG("Saving PTM frequencies to [" << fileName << "]");
      if (!m_fdrResults->saveModMatrix(fileName.c_str(), false)) {
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::loadOutputData(void)
  {
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFdrPeptide::split(int numSplit)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::merge(void)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::validateParams(std::string & error)
  {
    m_isValid = false;
    VALIDATE_PARAM_EXIST("TDA_TYPE");
    string tdaType = m_params.getValue("TDA_TYPE");
    if (DEBUG_FDR_PEPTIDE) DEBUG_VAR(tdaType);
    if (tdaType != "concatenated" && 
        tdaType != "separate" && 
        tdaType != "split")
    {
      ERROR_MSG("Unknown TDA_TYPE! [" << m_tdaType << "] Valid options are 'concatenated', 'separate' and 'split'.");
      return false;
    }  
    m_isValid = true;
    return true;
  }

}
