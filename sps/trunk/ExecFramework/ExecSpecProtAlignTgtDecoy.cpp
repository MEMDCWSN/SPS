#include "ExecSpecProtAlignTgtDecoy.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "ExecFdrPeptide.h"
#include "FileUtils.h"
#include "Logger.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParallelThreadedExecution.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>
#include <time.h>


static bool DEBUG_SPECPROTALIGN = false;
static bool DEBUG_SPECPROTALIGN_ANNO = false;
static bool DEBUG_SPECPROTALIGN_SPECPROB = false;
static bool DEBUG_SPECPROTALIGN_TIME = false;
static bool DEBUG_SPECPROTALIGN_SPRINKLE = false;

using namespace std;
using namespace specnets;

namespace specnets
{
  void ExecSpecProtAlignTgtDecoy::loadScanSpecificPenalties(
                                  string & filename,
                                  map<int, map<int, float> > & penalties)
  {
    vector<string> header;
    vector<vector<string> > lines;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;
    DelimitedTextReader::loadDelimitedFile(filename.c_str(),
                              "\t",
                              "",
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
    for (int i = 0; i < lines.size(); i++) {
      //DEBUG_VAR(lines[i][0]);
      string scanString = lines[i][0];
      string massString = lines[i][1];
      string scoreString = lines[i][2];
      int scan = atoi(scanString.c_str());
      int mass = atoi(massString.c_str());
      float score = atof(scoreString.c_str());
      penalties[scan][mass] = score;
    }
    return;
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(void) :
        m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), m_dbDecoy(0x0), 
        m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), 
        m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0),
        m_scanSpecificPenalties(0x0),
        m_contigAbinfo(m_contigAbinfo), m_filterPsmSet(m_filterPsmSet),
        ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    abinfo_t * m_contigAbinfo;
    PeptideSpectrumMatchSet * m_filterPsmSet;
  
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_prmSpectra(0x0), 
        m_db(0x0), m_dbDecoy(0x0), m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), 
        m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0),
        m_scanSpecificPenalties(0x0),
        m_contigAbinfo(0x0), m_filterPsmSet(0x0),
        ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       DB_fasta * dbDecoy,
                                       PeptideSpectrumMatchSet * psmSetTag,
                                       PeptideSpectrumMatchSet * psmSetTagDecoy,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,
                                       map<int, map<int, float> > * scanSpecificPenalties,
                                       abinfo_t * contigAbinfo,
                                       PeptideSpectrumMatchSet * filterPsmSet,
                                       SpecSet * matchedSpectraAll,
                                       SpecSet * matchedSpectra,
                                       PeptideSpectrumMatchSet * psmSet,
                                       PeptideSpectrumMatchSet * psmSetDecoy,
                                       PeptideSpectrumMatchSet * psmSetFdr) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra),
        m_db(db), m_dbDecoy(dbDecoy), m_psmSetTag(psmSetTag), m_psmSetTagDecoy(psmSetTagDecoy),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum), 
        m_penaltyMatrixMods(penaltyMatrixMods),
        m_scanSpecificPenalties(scanSpecificPenalties), ownInput(false),
        m_contigAbinfo(contigAbinfo), m_filterPsmSet(filterPsmSet),
        m_matchedSpectraAll(matchedSpectraAll), m_matchedSpectra(matchedSpectra),
        m_psmSet(psmSet), m_psmSetDecoy(psmSetDecoy), m_psmSetFdr(psmSetFdr), ownOutput(false)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }
                			

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::~ExecSpecProtAlignTgtDecoy(void)
  {
    if (ownInput) {
      if (m_inputSpectra)
        delete m_inputSpectra;
      if (m_prmSpectra)
        delete m_prmSpectra;
      if (m_db)
        delete m_db;
      if (m_dbDecoy)
        delete m_dbDecoy;
      if (m_psmSetTag)
        delete m_psmSetTag;
      if (m_psmSetTagDecoy)
        delete m_psmSetTagDecoy;
      if (m_penaltyMatrixBlosum)
        delete m_penaltyMatrixBlosum;
      if (m_penaltyMatrixMods)
        delete m_penaltyMatrixMods;
      if (m_scanSpecificPenalties)
        delete m_scanSpecificPenalties;
    }
    if (ownOutput) {
      if (m_matchedSpectraAll)
        delete m_matchedSpectraAll;
      if (m_matchedSpectra)
        delete m_matchedSpectra;
      if (m_psmSet)
        delete m_psmSet;
      if (m_psmSetDecoy)
        delete m_psmSetDecoy;
      if (m_psmSetFdr)
        delete m_psmSetFdr;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecProtAlignTgtDecoy::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlignTgtDecoy(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::invoke(void)
  {
    DEBUG_SPECPROTALIGN = m_params.getValueBool("DEBUG_SPECPROTALIGN");
    DEBUG_VAR(DEBUG_SPECPROTALIGN);
    DEBUG_SPECPROTALIGN_SPECPROB = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPECPROB");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPECPROB);
    DEBUG_SPECPROTALIGN_TIME = m_params.getValueBool("DEBUG_SPECPROTALIGN_TIME");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_TIME);
    DEBUG_SPECPROTALIGN_SPRINKLE = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPRINKLE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPRINKLE);

    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll or !m_matchedSpectra) {
      ERROR_MSG("ERROR: empty returned data pointer");
      return false;
    }

    if (m_filterPsmSet) DEBUG_VAR(m_filterPsmSet->size());
    if (m_contigAbinfo) DEBUG_VAR(m_contigAbinfo->size());

    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(peakTol);
    bool gridExecutionFlag = m_params.getValueBool("GRID_EXECUTION_FLAG", false);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(gridExecutionFlag);
    bool resume = m_params.getValueBool("GRID_RESUME_FLAG", false);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(resume);

    // Set up all the grid processing directories (on Proteosafe these need to be separate)
    m_params.setValue("GRID_DATA_DIR_PARAMS", m_params.getValue("GRID_DATA_DIR_TARGET"));
    m_params.setValue("GRID_DATA_DIR_IN", m_params.getValue("GRID_DATA_DIR_TARGET"));
    m_params.setValue("GRID_DATA_DIR_INTERMEDIATE", m_params.getValue("GRID_DATA_DIR_TARGET"));
    m_params.setValue("GRID_DATA_DIR_OUT", m_params.getValue("GRID_DATA_DIR_TARGET"));

    m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_TGT"));
    m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_TGT"));
    m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_TGT"));
    m_params.setValue("INPUT_TAGS_FOR_SPLITTING", m_params.getValue("INPUT_TAGS_TARGET"));
    m_params.setValue("INPUT_PSM", m_params.getValue("INPUT_TAGS_TARGET"));

    bool returnStatus;
    // -- TARGET PROTEIN ALIGNMENT

    m_inputSpectra->clearPsms(); // Clear all the PSMs
    if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
      if (m_params.exists("INPUT_TAGS_TARGET")) {
        PeptideSpectrumMatchSet psmSetTag;
        psmSetTag.loadFromFile(m_params.getValue("INPUT_TAGS_TARGET").c_str());
        m_inputSpectra->clearPsms();
        psmSetTag.addSpectra(m_inputSpectra);
      }
    }

    ExecSpecProtAlign moduleSpecProtAlign(m_params,
                                          m_inputSpectra,
                                          m_prmSpectra,
                                          m_db,
                                          m_penaltyMatrixBlosum,
                                          m_penaltyMatrixMods,
                                          m_scanSpecificPenalties,
                                          m_contigAbinfo,
                                          m_filterPsmSet,
                                          m_matchedSpectraAll);

    if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
      returnStatus = moduleSpecProtAlign.invoke();
    } else {
      if (DEBUG_SPECPROTALIGN) DEBUG_TRACE;
      int numNodes = m_params.getValueInt("GRID_NUMNODES");

      string gridType = m_params.getValue("GRID_TYPE");
      if (gridType == "pbs") {
        ParallelPbsExecution exec(&moduleSpecProtAlign,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      } else if (gridType == "sge") {
        ParallelSgeExecution exec(&moduleSpecProtAlign,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      } else if (gridType == "debug") {
        ParallelSgeExecution exec(&moduleSpecProtAlign,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                  resume,
                                  true);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "threaded")
      {
	      ParallelThreadedExecution exec(&moduleSpecProtAlign);
	      returnStatus = exec.invoke(numNodes);
      }
    }

    returnStatus = moduleSpecProtAlign.saveOutputData();
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(returnStatus);

    float fdrThreshold = m_params.getValueFloat("FDR_THRESHOLD", -1.0);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(fdrThreshold);

    // -- DECOY PROTEIN ALIGNMENT
    if (fdrThreshold >= 0.0 || !m_params.getValueInt("ENFORCE_ENDPEAKS", 0)) {

      // Set up all the grid processing directories (on Proteosafe these need to be separate)
      m_params.setValue("GRID_DATA_DIR_PARAMS", m_params.getValue("GRID_DATA_DIR_DECOY"));
      m_params.setValue("GRID_DATA_DIR_IN", m_params.getValue("GRID_DATA_DIR_DECOY"));
      m_params.setValue("GRID_DATA_DIR_INTERMEDIATE", m_params.getValue("GRID_DATA_DIR_DECOY"));
      m_params.setValue("GRID_DATA_DIR_OUT", m_params.getValue("GRID_DATA_DIR_DECOY"));

      m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_DEC"));
      m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_DEC"));
      m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_DEC"));
      m_params.setValue("INPUT_TAGS_FOR_SPLITTING", m_params.getValue("INPUT_TAGS_DECOY"));
      m_params.setValue("INPUT_PSM", m_params.getValue("INPUT_TAGS_DECOY"));
      m_params.setValue("USE_DECOY_DATABASE", "1");

      m_inputSpectra->clearPsms(); // Clear all the PSMs
      if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
        if (m_params.exists("INPUT_TAGS_DECOY")) {
          PeptideSpectrumMatchSet psmSetTag;
          psmSetTag.loadFromFile(m_params.getValue("INPUT_TAGS_DECOY").c_str());
          psmSetTag.addSpectra(m_inputSpectra);
        }
      }
      
      ExecSpecProtAlign moduleSpecProtAlign2(m_params,
                                            m_inputSpectra,
                                            m_prmSpectra,
                                            m_dbDecoy,
                                            m_penaltyMatrixBlosum,
                                            m_penaltyMatrixMods,
                                            m_scanSpecificPenalties,
                                            m_contigAbinfo,
                                            m_filterPsmSet,
                                            &m_matchedSpectraAllDecoy);

      if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
        returnStatus = moduleSpecProtAlign2.invoke();
      } else {
        if (DEBUG_SPECPROTALIGN) DEBUG_TRACE;
        int numNodes = m_params.getValueInt("GRID_NUMNODES");

        string gridType = m_params.getValue("GRID_TYPE");
        if (gridType == "pbs") {
          ParallelPbsExecution exec(&moduleSpecProtAlign2,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                    resume);
          returnStatus = exec.invoke(numNodes);
        } else if (gridType == "sge") {
          ParallelSgeExecution exec(&moduleSpecProtAlign2,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                    resume);
          returnStatus = exec.invoke(numNodes);
        } else if (gridType == "debug") {
          ParallelSgeExecution exec(&moduleSpecProtAlign2,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                    resume,
                                    true);
          returnStatus = exec.invoke(numNodes);
        }
	      else if (gridType == "threaded")
	      {
	        ParallelThreadedExecution exec(&moduleSpecProtAlign2);
	        returnStatus = exec.invoke(numNodes);
	      }
      }

      if ((!m_params.exists("GRID_NUMNODES") || m_params.getValueInt("GRID_NUMNODES") > 0)
          && (!resume && !gridExecutionFlag)) {
        return true;
      }
      
      returnStatus = moduleSpecProtAlign2.saveOutputData();
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(returnStatus);

      m_psmSet->getPSMSet(m_matchedSpectraAll);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSet->size());
      if (DEBUG_SPECPROTALIGN) m_psmSet->saveToFile("debug_psm_target.txt");
      
      m_psmSetDecoy->getPSMSet(&m_matchedSpectraAllDecoy);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSetDecoy->size());

      // Make sure the isDecoy fields are set on the decoy PSMs
      for (int iDecoy = 0; iDecoy < m_psmSetDecoy->size(); iDecoy++) {
        (*m_psmSetDecoy)[iDecoy]->m_isDecoy = true;
      }
      if (DEBUG_SPECPROTALIGN) m_psmSetDecoy->saveToFile("debug_psm_decoy.txt");

      //m_params.setValue("PEPTIDE_FDR_AUTOSCALE", "1");
      m_params.setValue("TDA_TYPE", m_params.getValue("PEPTIDE_FDR_TDA_TYPE", "concatenated"));
      m_params.setValue("PEPTIDE_FDR_CUTOFF", m_params.getValue("FDR_THRESHOLD", "0.01"));
      m_params.setValue("PEPTIDE_FDR_USE_PVALUE_COMPARE", "1");
//      m_params.setValue("PEPTIDE_FDR_COMPUTE_VARIANT_GROUPS", "1");

      PeptideSpectrumMatchSet psmSetCombined;
      ExecFdrPeptide moduleFdrPeptide(m_params,
                                      m_psmSet,
                                      m_psmSetDecoy,
                                      m_psmSetFdr);
      returnStatus = moduleFdrPeptide.invoke();

      if (DEBUG_SPECPROTALIGN) m_psmSetFdr->saveToFile("debug_psm_cutoff.txt");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSetFdr->size());
    }

    if ((!m_params.exists("GRID_NUMNODES") || m_params.getValueInt("GRID_NUMNODES") > 0)
        && (!resume && !gridExecutionFlag)) {
      return true;
    }
      

    // Resize to make enough room
    m_matchedSpectraIndices.resize(m_matchedSpectraAll->size());
    m_matchedSpectra->resize(m_matchedSpectraAll->size());

    int keepIdx = 0;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(i);
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].size());
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
        m_matchedSpectraIndices[keepIdx] = i;
        (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
      }
    }

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(keepIdx);
    m_matchedSpectra->resize(keepIdx);
    m_matchedSpectraIndices.resize(keepIdx);

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectra->size());
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraIndices.size());

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_prmSpectra)
        m_prmSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
      if (!m_dbDecoy)
        m_dbDecoy = new DB_fasta;
      if (!m_psmSetTag)
        m_psmSetTag = new PeptideSpectrumMatchSet;
      if (!m_psmSetTagDecoy)
        m_psmSetTagDecoy = new PeptideSpectrumMatchSet;
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
      if (!m_psmSet)
        m_psmSet = new PeptideSpectrumMatchSet;
      if (!m_psmSetDecoy)
        m_psmSetDecoy = new PeptideSpectrumMatchSet;
      if (!m_psmSetFdr)
        m_psmSetFdr = new PeptideSpectrumMatchSet;
    }
    m_matchedSpectraAll->resize(0);

    if (m_params.exists("INPUT_SPECS_PKLBIN")) {
      DEBUG_MSG("Loading: " << m_params.getValue("INPUT_SPECS_PKLBIN"));
      if (m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (m_params.exists("INPUT_PRM_PKLBIN"))
    {
      DEBUG_MSG("Loading: " << m_params.getValue("INPUT_PRM_PKLBIN"));
      if (m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_PKLBIN").c_str())
          < 0)
      {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    DEBUG_MSG("Loading: " << m_params.getValue("FASTA_DATABASE"));
    if (!m_params.exists("FASTA_DATABASE")) {
      ERROR_MSG("Parameters are incomplete. FASTA_DATABASE is missing.");
      return false;
    }
    else if (m_db->Load(m_params.getValue("FASTA_DATABASE").c_str()) <= 0) {
      ERROR_MSG("Error reading database sequences from "
          << m_params.getValue("FASTA_DATABASE"));
      return false;
    }

    if (!m_params.exists("FASTA_DATABASE_DECOY")) {
      if (m_dbDecoy->Load(m_params.getValue("FASTA_DATABASE_DECOY").c_str()) <= 0) {
        ERROR_MSG("Error reading decoy database sequences from "
            << m_params.getValue("FASTA_DATABASE_DECOY"));
        return false;
      }
    }

    //---------------------------------------------------------------------------
    // Load penalty matrices for new alignment
    //---------------------------------------------------------------------------

    bool penaltyAlign = m_params.exists("PENALTY_ALIGNMENT") ?
    		m_params.getValueInt("PENALTY_ALIGNMENT", false) : false;

    if (penaltyAlign) {

      // Load amino acid masses
      AAJumps jumps(1);
      if (!m_params.exists("AMINO_ACID_MASSES")) {
        DEBUG_MSG("Loading: " << m_params.getValue("AMINO_ACID_MASSES"));
        jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
      }

      float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
      float knownModPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_KNOWN_PENALTY", 0.01);
      float unknownPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY", 1.0);
      float unknownMultiplier = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 1.5);
      float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100.0);
      float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", -100.0);

      if (!m_penaltyMatrixBlosum)
        m_penaltyMatrixBlosum = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);
      if (!m_penaltyMatrixMods)
        m_penaltyMatrixMods = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);

      if (m_params.exists("BLOSUM_PENALTY_FILE")) {
        string blosumFileName = m_params.getValue("BLOSUM_PENALTY_FILE"); 
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE"); 
        string cleavagePenaltiesFileName = m_params.getValue("CLEAVAGE_PENALTY_FILE");
        DEBUG_MSG("Loading: " << m_params.getValue("BLOSUM_PENALTY_FILE"));
        if (!m_penaltyMatrixBlosum->load(blosumFileName,
                                         knowmModsFileName,
                                         cleavagePenaltiesFileName)) {
          ERROR_MSG("Error loading blosum penalties from " << blosumFileName);
          return false;
        }
      }

      if (!m_params.exists("MODS_PENALTY_FILE")) {
        ERROR_MSG("Parameters are incomplete. MODS_PENALTY_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string modFileName = m_params.getValue("MODS_PENALTY_FILE"); 
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE"); 
        string cleavagePenaltiesFileName = m_params.getValue("CLEAVAGE_PENALTY_FILE");
        DEBUG_MSG("Loading: " << m_params.getValue("MODS_PENALTY_FILE"));
        if (!m_penaltyMatrixMods->load(modFileName,
                                       knowmModsFileName,
                                       cleavagePenaltiesFileName)) {
          ERROR_MSG("Error loading mod penalties from " << modFileName);
          return false;
        }
      }

      if (m_params.exists("SCAN_SPECIFIC_PENALTIES_FILE")) {
        string filename = m_params.getValue("SCAN_SPECIFIC_PENALTIES_FILE");
        DEBUG_MSG("Loading: " << filename);
        if (!m_scanSpecificPenalties)
          m_scanSpecificPenalties = new map<int, map<int, float> >();
        loadScanSpecificPenalties(filename,
                                  *m_scanSpecificPenalties);
      }

    } // if (penaltyAlign)

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveOutputData(void)
  {
    bool penaltyAlign = m_params.exists("PENALTY_ALIGNMENT") ?
  	  	m_params.getValueBool("PENALTY_ALIGNMENT", false) : false;
    DEBUG_VAR(penaltyAlign);

    string outDir = m_params.getValue("GRID_DATA_DIR");
    if (outDir.empty()) {
      outDir = ".";
    }
    outDir  = outDir + "/";
    DEBUG_VAR(outDir);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX")) {
      int keepIdx = 0;
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
        }
      }

      DEBUG_MSG("Saving matched peaks...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectra->size());
      for (int i = 0; i < m_matchedSpectra->size(); i++) {
        if ((*m_matchedSpectra)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr = (*m_matchedSpectra)[i].psmList.begin();
          int peakListSize = (*litr)->m_matchedPeaks.size();
          //DEBUG_VAR(peakListSize);
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            //DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      string fileName = outDir + m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX").c_str();
      tempMatchedPeaks.savePklBin(fileName.c_str());
    }

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectra->size());
    if (m_matchedSpectra and m_params.exists("OUTPUT_MATCHED_SPECS")) {

      string specsName = outDir + m_params.getValue("OUTPUT_MATCHED_SPECS").c_str();
      string psmName = outDir + m_params.getValue("OUTPUT_PSM").c_str();
      string peaksName = outDir + m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX").c_str();
      DEBUG_VAR(specsName);
      DEBUG_VAR(psmName);
      DEBUG_VAR(peaksName);

      if (m_params.exists("OUTPUT_PSM") && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX")) {
        DEBUG_MSG("Saving matched specs (3 files)...");
        m_matchedSpectra->savePklBin(specsName.c_str(),
                                     psmName.c_str(),
                                     peaksName.c_str());

      } else if (m_params.exists("OUTPUT_PSM")) {
        DEBUG_MSG("Saving matched specs (2 files)...");
        m_matchedSpectra->savePklBin(specsName.c_str(),
                                     psmName.c_str());

      } else {
        DEBUG_MSG("Saving matched specs (1 file)...");
        m_matchedSpectra->savePklBin(peaksName.c_str());
      }
    }

    //=================================================
    // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
    //=================================================
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraAll->size());
    if (m_params.exists("OUTPUT_MATCHED_PROTS_ALL")) {
      DEBUG_MSG("Saving matched prots all...");
      string fileName = outDir + m_params.getValue("OUTPUT_MATCHED_PROTS_ALL").c_str();
      DEBUG_MSG(fileName);
      m_matchedSpectraAll->saveMatchedProts(fileName.c_str());
    }
    if (m_params.exists("OUTPUT_MATCHED_PROTS")) {
      DEBUG_MSG("Saving matched prots...");
      string fileName = outDir + m_params.getValue("OUTPUT_MATCHED_PROTS").c_str();
      DEBUG_MSG(fileName);
      m_matchedSpectra->saveMatchedProts(fileName.c_str());
    }
    if (m_params.exists("OUTPUT_MATCHED_SPECS_IDX")) {
      DEBUG_MSG("Saving matched specs idx...");
      string fileName = outDir + m_params.getValue("OUTPUT_MATCHED_SPECS_IDX").c_str();
      DEBUG_MSG(fileName);
      Save_binArray(fileName.c_str(), m_matchedSpectraIndices);
    }
    //=================================================

    if (m_params.exists("OUTPUT_PSM_MOD_FREQS")) {
      string fileName = outDir + m_params.getValue("OUTPUT_PSM_MOD_FREQS").c_str();
      DEBUG_MSG(fileName);
      m_psmSetFdr->saveModMatrix(fileName.c_str(), false);
    }

    if (m_params.exists("OUTPUT_PSM_FDR")) {
      string fileName = outDir + m_params.getValue("OUTPUT_PSM_FDR").c_str();
      DEBUG_MSG(fileName);
      m_psmSetFdr->saveToFile(fileName.c_str());
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlignTgtDecoy::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::merge(void)
  {
    return false;
  }


  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

}

