//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_perform_specprotalign.h"

// Module Includes
#include "AlignmentPenaltyBased.h"
#include "Logger.h"
#include "ExecPenaltyGen.h"
#include "ExecSpecProtAlignTgtDecoy.h"
#include "ExecSpecTagGen.h"
#include "FileUtils.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
#include "utils.h"  // stringSplit
#include "clusters.h"
#include "ClusterData.h"
#include "abruijn.h"
#include "copyright.h"
#include "SpecSet.h"
#include "Specific.h"
#include "StatusFile.h"

// System Includes
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

using namespace specnets;
using namespace std;

#define DEBUG_SPEC_ALIGN 1

namespace specnets
{
  //---------------------------------------------------------------------------
  bool performPenaltyGen(ParameterList & ip,
                         SpectrumPairSet & filteredPairs,
                         SpecSet & prmSpectra,
                         PenaltyMatrix & penaltyMatrixMods,
                         map<int, map<int, float> > & scanSpecificPenalties)
  {

    ParameterList specPenaltyGenParams;
    
    specPenaltyGenParams.addIfExists(ip, "ALIGNMENT_RESOLUTION");
    specPenaltyGenParams.addIfDoesntExist("ALIGNMENT_RESOLUTION", "1.0");
    specPenaltyGenParams.addIfExists(ip, "MAX_PEAK_EQUIVALENTS");
    specPenaltyGenParams.addIfDoesntExist("MAX_PEAK_EQUIVALENTS", "1.5");
    specPenaltyGenParams.addIfExists(ip, "MIN_PEAK_EQUIVALENTS");
    specPenaltyGenParams.addIfDoesntExist("MIN_PEAK_EQUIVALENTS", "1.0");
    specPenaltyGenParams.addIfExists(ip, "MIN_PENALTY_FREQUENCY");
    specPenaltyGenParams.addIfDoesntExist("MIN_PENALTY_FREQUENCY", "0.005");
    specPenaltyGenParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_PENALTY");
    specPenaltyGenParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_PENALTY", "1.0");
    specPenaltyGenParams.addIfExists(ip, "PENALTY_ALIGNMENT_KNOWN_PENALTY");
    specPenaltyGenParams.addIfDoesntExist("PENALTY_ALIGNMENT_KNOWN_PENALTY", "0.01");
    specPenaltyGenParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER");
    specPenaltyGenParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", "1.5");

    specPenaltyGenParams.addIfExists(ip, "MIN_MOD_MASS");
    specPenaltyGenParams.addIfDoesntExist("MIN_MOD_MASS", "-100.0");
    specPenaltyGenParams.addIfExists(ip, "MAX_MOD_MASS");
    specPenaltyGenParams.addIfDoesntExist("MAX_MOD_MASS", "100.0");

    specPenaltyGenParams.addIfExists(ip, "CLEAVAGE_PENALTY_FILE");
    specPenaltyGenParams.addIfExists(ip, "KNOWN_MODS_FILE");

#if 0
    specTagGenParams.setValue("AMINO_ACID_MASSES",
                              getProjPath(ip,
                                          "homology/specprotalign_aa_masses.txt"));
    specTagGenParams.setValue("INPUT_BLOSUM_FILE", 
                              ip.getValue("BLOSUM_PENALTY_FILE"));
    specTagGenParams.setValue("INPUT_KNOWN_MODS", 
                              ip.getValue("KNOWN_MODS_FILE"));
    specTagGenParams.setValue("INPUT_CLEAVAGE_PENALTIES", 
                              ip.getValue("CLEAVAGE_PENALTY_FILE"));
#endif
                              
    specPenaltyGenParams.setValue("OUTPUT_MOD_PENALTIES",
                                  getProjPath(ip, "homology/specprotalign_mod_penal.txt"));
    specPenaltyGenParams.setValue("OUTPUT_KNOWN_MODS",
                                  getProjPath(ip, "homology/specprotalign_known_mods.txt"));
    specPenaltyGenParams.setValue("OUTPUT_CLEAVAGE_PEN",
                                  getProjPath(ip, "homology/specprotalign_cleave_pen.txt"));
    specPenaltyGenParams.setValue("OUTPUT_SCAN_SPECIFIC_PENALTIES",
                                  getProjPath(ip, "homology/specprotalign_scan_specific_pen.txt"));
                                  

    specPenaltyGenParams.addIfExists(ip, "DEBUG_PENALTY_GEN");
    specPenaltyGenParams.addIfExists(ip, "DEBUG_SCAN_SPECIFIC");

    DEBUG_TRACE;
    ExecPenaltyGen modulePenaltyGen(specPenaltyGenParams,
                                    &filteredPairs,
                                    &prmSpectra,
                                    &penaltyMatrixMods,
                                    &scanSpecificPenalties);
    DEBUG_TRACE;

    string errorString;
    bool isValid = modulePenaltyGen.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = modulePenaltyGen.invoke();

    // Test for return status
    TEST_RETURN_STATUS("modulePenaltyGen");

    // No need - nothing to do here (its all done in invoke)
    //returnStatus = modulePenaltyGen.saveOutputData();
    DEBUG_VAR(returnStatus);
    
    TEST_SAVE_OUPUT_DATA("modulePenaltyGen");

    return true;
  }


  //---------------------------------------------------------------------------
  bool performContigProtAlign(ParameterList & ip,
                              SpecSet & contigSpectra,
                              DB_fasta & db,
                              DB_fasta & dbDecoy,
                              PenaltyMatrix & penaltyMatrixBlosum,
                              PenaltyMatrix & penaltyMatrixMods,
                              map<int, map<int, float> > & scanSpecificPenalties,
                              abinfo_t & contigAbinfo,
                              PeptideSpectrumMatchSet & filterPsmSet,
                              SpecSet & matchedSpectraAll,
                              SpecSet & matchedSpectra,
                              PeptideSpectrumMatchSet & psmSet,
                              PeptideSpectrumMatchSet & psmSetDecoy,
                              PeptideSpectrumMatchSet & psmSetFdr,
                              bool gridExecutionFlag,
                              bool resume)
  {
    ParameterList specProtAlignParams;
    specProtAlignParams.addIfExists(ip, "ALIGNMENT_RESOLUTION");
    specProtAlignParams.addIfExists(ip, "TOLERANCE_PEAK");
    specProtAlignParams.addIfExists(ip, "TOLERANCE_PM");

    specProtAlignParams.setValue("GRID_EXECUTION_FLAG",
                                 gridExecutionFlag ? "1" : "0");
    specProtAlignParams.setValue("GRID_RESUME_FLAG", resume ? "1" : "0");
    specProtAlignParams.addIfExists(ip, "GRID_TYPE");

    //If we have an argument called "GRID_NUMNODES_SPECPROTALIGN", then replace GRID NUMNODES with that
    if (ip.exists("GRID_NUMNODES_SPECPROTALIGN"))
    {
      string val = ip.getValue("GRID_NUMNODES_SPECPROTALIGN");
      specProtAlignParams.setValue("GRID_NUMNODES", val);
    }
    else
      specProtAlignParams.addIfExists(ip, "GRID_NUMNODES");

    specProtAlignParams.addIfExists(ip, "GRID_NUMCPUS");
    specProtAlignParams.addIfExists(ip, "GRID_EXE_DIR");
    specProtAlignParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
    specProtAlignParams.addIfExists(ip, "GRID_PARAMS");
    specProtAlignParams.setValue("GRID_DATA_DIR", 
                                 getProjPath(ip, "", true));
    specProtAlignParams.setValue("GRID_DATA_DIR_TARGET",
                                 getProjPath(ip, "homology/grid_t", true));
    specProtAlignParams.setValue("GRID_DATA_DIR_DECOY",
                                 getProjPath(ip, "homology/grid_d", true));
    specProtAlignParams.addIfExists(ip, "RELATIVE_DIR");

    specProtAlignParams.addIfExists(ip, "PROJECT_DIR");
    specProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    if (ip.exists("MIN_MATCHED_PEAKS_DB"))
      specProtAlignParams.setValue("MIN_MATCHED_PEAKS_DB",
                                   ip.getValue("MIN_MATCHED_PEAKS_DB"));
    else
      specProtAlignParams.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7");
    specProtAlignParams.addIfExists(ip, "MAX_NUM_MODS");
    specProtAlignParams.addIfDoesntExist("MAX_NUM_MODS", "2");
    specProtAlignParams.addIfExists(ip, "MIN_MOD_MASS");
    specProtAlignParams.addIfDoesntExist("MIN_MOD_MASS", "-100");
    specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
    specProtAlignParams.addIfDoesntExist("MAX_MOD_MASS", "100");
    specProtAlignParams.addIfExists(ip, "MAX_PARSIMONY");
    specProtAlignParams.addIfDoesntExist("MAX_PARSIMONY", "1");
    specProtAlignParams.addIfExists(ip, "ALIGNMENT_SCORE_THRESHOLD");
    specProtAlignParams.addIfDoesntExist("ALIGNMENT_SCORE_THRESHOLD", "0.75");

    specProtAlignParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
    specProtAlignParams.addIfDoesntExist("MAX_ALIGN_DB_GAP_AAS", "8");
    specProtAlignParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");
    specProtAlignParams.addIfDoesntExist("MAX_ALIGN_SPECTRUM_GAP_DALTONS", "1500");

    if (ip.exists("FDR_CONTIG_THRESHOLD"))
    {
      specProtAlignParams.setValue("FDR_THRESHOLD",
                                   ip.getValue("FDR_CONTIG_THRESHOLD"));
    }

    if (ip.exists("CONTIGPROTALIGN_START_IDX"))
    {
      specProtAlignParams.setValue("IDX_START",
                                   ip.getValue("CONTIGPROTALIGN_START_IDX"));
    }
    if (ip.exists("CONTIGPROTALIGN_END_IDX"))
    {
      specProtAlignParams.setValue("IDX_END",
                                   ip.getValue("CONTIGPROTALIGN_END_IDX"));
    }

    specProtAlignParams.addIfExists(ip, "FASTA_DATABASE");
    specProtAlignParams.addIfExists(ip, "AMINO_ACID_MASSES");
    specProtAlignParams.addIfExists(ip, "BLOSUM_PENALTY_FILE");

    specProtAlignParams.addIfExists(ip, "INPUT_FILTER_PSMS");
    specProtAlignParams.addIfExists(ip, "INPUT_FILTER_MSGFDB_PSMS");

    specProtAlignParams.setValue("INPUT_ABINFO", 
        getProjPath(ip, "assembly/component_info.bin", true));
    specProtAlignParams.setValue("MODS_PENALTY_FILE", 
        getProjPath(ip, "homology/specprotalign_mod_penal.txt", true));
    specProtAlignParams.setValue("KNOWN_MODS_FILE", 
        getProjPath(ip, "homology/specprotalign_known_mods.txt", true));
    specProtAlignParams.setValue("CLEAVAGE_PENALTY_FILE", 
        getProjPath(ip, "homology/specprotalign_cleave_pen.txt", true));
    specProtAlignParams.setValue("SCAN_SPECIFIC_PENALTIES_FILE",
        getProjPath(ip, "homology/specprotalign_scan_specific_pen.txt", true));

    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_TGT", "homology/contigs_midx_tgt.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_TGT", "homology/contigs_tgt.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM_TGT", "homology/contigs_psm_tgt.txt");
    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_DEC", "homology/contigs_midx_dec.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_DEC", "homology/contigs_dec.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM_DEC", "homology/contigs_psm_dec.txt");

    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", "homology/contigs_midx_all.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX", "homology/contigs_midx.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS", "homology/contigs.pklbin");
    specProtAlignParams.setValue("OUTPUT_SPECTRA_SPRINKLED", "homology/contigs_sprinkled.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM", "homology/contigs_psm.txt");

    specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS", "homology/contigs_mp.bin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS_ALL", "homology/contigs_mp_all.bin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_IDX", "homology/contigs_indices.bin");

    specProtAlignParams.setValue("OUTPUT_PSM_MOD_FREQS", "homology/contigs_psm_mod_freqs.txt");
    specProtAlignParams.setValue("OUTPUT_PSM_FDR", "homology/contigs_psm_fdr.txt");

    specProtAlignParams.setValue("ENFORCE_ENDPEAKS", "0"); // Do not enforce endpeaks for contigs

    specProtAlignParams.addIfExists(ip, "ALIGNMENT_SCORE_THRESHOLD");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_PENALTY");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_KNOWN_PENALTY");
    specProtAlignParams.addIfExists(ip, "SPECPROTALIGN_AVERAGE_TYPE");
    specProtAlignParams.addIfExists(ip, "SPECPROTALIGN_USE_TAG_SEEDING");

    specProtAlignParams.addIfExists(ip, "PEPTIDE_FDR_TDA_TYPE");

    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECPROB");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECTRA");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TIME");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_ANNO");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPRINKLE");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SINGLESPECTRUM");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_MODSPECTRUM");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPLIT");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_RANGE");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TAG_SEEDING");
    
    specProtAlignParams.addIfExists(ip, "DEBUG_TAGSEARCH_SINGLE");

    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_RANGE");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_RANGE2");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_AAS");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_SPECS");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN1");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN2");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN3");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN4");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_EXACT");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_NTERM");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_CTERM");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP_ANNO");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP_CACHE");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_LARGE_GAP");

    specProtAlignParams.setValue("SPEC_TYPE_MSMS", "0");
    specProtAlignParams.setValue("MIN_RATIO", "0.4");

    // Parameters for Tag Matching
    specProtAlignParams.addIfExists(ip, "RESOLUTION");
    specProtAlignParams.addIfExists(ip, "TAG_LEN");
    specProtAlignParams.addIfDoesntExist("TAG_LEN", "6");
    specProtAlignParams.addIfExists(ip, "DOUBLE_AA_JUMPS");
    specProtAlignParams.addIfDoesntExist("DOUBLE_AA_JUMPS", "1");
    specProtAlignParams.addIfExists(ip, "MAX_NUM_TAGS");
    specProtAlignParams.addIfDoesntExist("MAX_NUM_TAGS", "0");
    specProtAlignParams.addIfExists(ip, "MATCH_TAG_FLANKING_MASSES");
    specProtAlignParams.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0");
    specProtAlignParams.addIfExists(ip, "TAG_MATCH_TOP_SCORING_ONLY");
    specProtAlignParams.addIfExists(ip, "TAG_PEAK_SKIP_PENALTY");
    specProtAlignParams.addIfExists(ip, "DEBUG_TAGSEARCH_CORRECT");
    specProtAlignParams.addIfExists(ip, "DEBUG_TAGSEARCH_SINGLE");

    DEBUG_VAR(contigSpectra.size());

    specProtAlignParams.writeToFile(getProjPath(ip,
                                                "debug_contigprotalign.params"));

    PeptideSpectrumMatchSet psmTag; // Empty psm sets (tags generated inside module now
    PeptideSpectrumMatchSet psmTagDecoy;
    ExecSpecProtAlignTgtDecoy moduleSpecProtAlignTgtDecoy(specProtAlignParams,
                                                          &contigSpectra,
                                                          0x0, // No PRM spectra to "sprinkle" into contig spectra
                                                          &db,
                                                          &dbDecoy,
                                                          &psmTag,
                                                          &psmTagDecoy,
                                                          &penaltyMatrixBlosum,
                                                          &penaltyMatrixMods,
                                                          &scanSpecificPenalties,
                                                          &contigAbinfo,
                                                          &filterPsmSet,
                                                          &matchedSpectraAll,
                                                          &matchedSpectra,
                                                          &psmSet,
                                                          &psmSetDecoy,
                                                          &psmSetFdr);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleSpecProtAlignTgtDecoy.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleSpecProtAlignTgtDecoy.invoke();

    // Test for return status
    TEST_RETURN_STATUS("ExecSpecProtAlignTgtDecoy");

    DEBUG_VAR(matchedSpectra.size());

    returnStatus = moduleSpecProtAlignTgtDecoy.saveOutputData();
    DEBUG_VAR(returnStatus);
    TEST_SAVE_OUPUT_DATA("ExecSpecProtAlignTgtDecoy");

    return true;
  }



  //---------------------------------------------------------------------------
  bool performSpecTagGen(ParameterList & ip,
                         abinfo_t & contigAbinfo,
                         SpecSet & contigSpectra,
                         SpecSet & starSpectra,
                         SpecSet & matchedContigs,
                         PeptideSpectrumMatchSet & psmContigTags,
                         PeptideSpectrumMatchSet & psmContigTagsDecoy)
  {

    ParameterList specTagGenParams;
    specTagGenParams.setValue("INPUT_ABINFO",
                              getProjPath(ip, "assembly/component_info.bin"));
                                             
    specTagGenParams.setValue("INPUT_CONTIGS",
                              getProjPath(ip, "homology/grid_t/contigs_tgt.pklbin"));

    specTagGenParams.setValue("INPUT_STARS_PKLBIN",
                              getProjPath(ip, "spectra/stars.pklbin"));
    
    specTagGenParams.setValue("INPUT_MATCHED_CONTIGS_PKLBIN",
                              getProjPath(ip, "homology/grid_t/contigs_tgt.pklbin").c_str());
    
    specTagGenParams.setValue("INPUT_CONTIG_TAGS_TGT",
                              getProjPath(ip, "homology/grid_t/contigs_psm_tgt.txt").c_str());
    
    specTagGenParams.setValue("INPUT_CONTIG_TAGS_DEC",
                              getProjPath(ip, "homology/grid_d/contigs_psm_dec.txt").c_str());
    
    specTagGenParams.setValue("OUTPUT_CONTIG_STAR_TAGS_TGT",
                              getProjPath(ip, "homology/contig_star_tags_tgt.txt").c_str());
    specTagGenParams.setValue("OUTPUT_CONTIG_STAR_TAGS_DEC",
                              getProjPath(ip, "homology/contig_star_tags_dec.txt").c_str());

                              
    specTagGenParams.writeToFile(getProjPath(ip,
                                             "debug_spectaggen.params"));

    specTagGenParams.addIfExists(ip, "DEBUG_SINGLECONTIG");

    DEBUG_TRACE;
    ExecSpecTagGen moduleSpecTagGen(specTagGenParams,
                                    &contigAbinfo,
                                    &contigSpectra,
                                    &starSpectra,
                                    &matchedContigs,
                                    &psmContigTags,
                                    &psmContigTagsDecoy);
                                  
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleSpecTagGen.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleSpecTagGen.invoke();

    // Test for return status
    TEST_RETURN_STATUS("moduleSpecTagGen");

    // No need - nothing to do here (its all done in invoke)
    //returnStatus = moduleSpecTagGen.saveOutputData();
    DEBUG_VAR(returnStatus);
    
    TEST_SAVE_OUPUT_DATA("moduleSpecTagGen");

    return true;
  }

  //---------------------------------------------------------------------------
  bool performSpecProtAlign(ParameterList & ip,
                            SpecSet & contigSpectra,
                            SpecSet & prmSpectra,
                            DB_fasta & db,
                            DB_fasta & dbDecoy,
                            PeptideSpectrumMatchSet & psmTag,
                            PeptideSpectrumMatchSet & psmTagDecoy,
                            PenaltyMatrix & penaltyMatrixBlosum,
                            PenaltyMatrix & penaltyMatrixMods,
                            map<int, map<int, float> > & scanSpecificPenalties,
                            abinfo_t & contigAbinfo,
                            PeptideSpectrumMatchSet & filterPsmSet,
                            SpecSet & matchedSpectraAll,
                            SpecSet & matchedSpectra,
                            PeptideSpectrumMatchSet & psmSet,
                            PeptideSpectrumMatchSet & psmSetDecoy,
                            PeptideSpectrumMatchSet & psmSetFdr,
                            bool gridExecutionFlag,
                            bool resume)
  {
    ParameterList specProtAlignParams;
    specProtAlignParams.addIfExists(ip, "ALIGNMENT_RESOLUTION");
    specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
    specProtAlignParams.addIfExists(ip, "TOLERANCE_PEAK");
    specProtAlignParams.addIfExists(ip, "TOLERANCE_PM");

    // We need to reset endpoints for spectrum alignment
    specProtAlignParams.setValue("RESET_ENDPOINTS", "1");
    specProtAlignParams.addIfExists(ip, "SPEC_TYPE_MSMS");

    specProtAlignParams.setValue("GRID_EXECUTION_FLAG",
                                 gridExecutionFlag ? "1" : "0");
    specProtAlignParams.setValue("GRID_RESUME_FLAG", resume ? "1" : "0");
    specProtAlignParams.addIfExists(ip, "GRID_TYPE");
    specProtAlignParams.addIfExists(ip, "GRID_NUMNODES");
    specProtAlignParams.addIfExists(ip, "GRID_NUMCPUS");
    specProtAlignParams.addIfExists(ip, "GRID_EXE_DIR");
    specProtAlignParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
    specProtAlignParams.addIfExists(ip, "GRID_PARAMS");
    specProtAlignParams.setValue("GRID_DATA_DIR", 
                                 getProjPath(ip, "", true));
    specProtAlignParams.setValue("GRID_DATA_DIR_TARGET",
                                 getProjPath(ip, "./homology/grid2_t", true));
    specProtAlignParams.setValue("GRID_DATA_DIR_DECOY",
                                 getProjPath(ip, "./homology/grid2_d", true));
    specProtAlignParams.addIfExists(ip, "RELATIVE_DIR");

    specProtAlignParams.addIfExists(ip, "PROJECT_DIR");
    specProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    if (ip.exists("MIN_MATCHED_PEAKS_DB"))
      specProtAlignParams.setValue("MIN_MATCHED_PEAKS_DB",
                                   ip.getValue("MIN_MATCHED_PEAKS_DB"));
    else
      specProtAlignParams.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7");
    specProtAlignParams.addIfExists(ip, "MAX_NUM_MODS");
    specProtAlignParams.addIfDoesntExist("MAX_NUM_MODS", "2");
    specProtAlignParams.addIfExists(ip, "MIN_MOD_MASS");
    specProtAlignParams.addIfDoesntExist("MIN_MOD_MASS", "-100");
    specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
    specProtAlignParams.addIfDoesntExist("MAX_MOD_MASS", "100");
    specProtAlignParams.addIfExists(ip, "MAX_PARSIMONY");
    specProtAlignParams.addIfDoesntExist("MAX_PARSIMONY", "1");

    specProtAlignParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
    specProtAlignParams.addIfDoesntExist("MAX_ALIGN_DB_GAP_AAS", "8");
    specProtAlignParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");
    specProtAlignParams.addIfDoesntExist("MAX_ALIGN_SPECTRUM_GAP_DALTONS",
                                         "1500");

    specProtAlignParams.setValue("ENFORCE_ENDPEAKS", "1"); // Always enforce end peaks for spectrum matches

    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_PENALTY");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_KNOWN_PENALTY");
    specProtAlignParams.addIfExists(ip, "USE_ALIGNMENT_SCORE_FOR_PROBABILITY");
    specProtAlignParams.addIfExists(ip,
                                    "SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER");
    specProtAlignParams.addIfExists(ip,
                                    "SPECPROTALIGN_PRM_CONTRIBUTION_RANK");
    specProtAlignParams.addIfExists(ip, "SPECPROTALIGN_ROUND_ANNOTATION_MAX");
    specProtAlignParams.addIfExists(ip, "SPECPROTALIGN_AVERAGE_TYPE");

    specProtAlignParams.addIfExists(ip, "PEPTIDE_FDR_TDA_TYPE");
    specProtAlignParams.addIfExists(ip, "FDR_USE_VARIANTS");

    // This is really for ProteoSAFe (but have it here so it can be tested)
    specProtAlignParams.addIfExists(ip, "CHANGE_GAP_ANNOTATIONS_TO_SINGLE");

    DEBUG_VAR(ip.getValue("FDR_SPECTRUM_THRESHOLD"));
    if (ip.exists("FDR_SPECTRUM_THRESHOLD"))
    {
      specProtAlignParams.setValue("FDR_THRESHOLD",
                                   ip.getValue("FDR_SPECTRUM_THRESHOLD"));
    }

    if (ip.exists("SPECPROTALIGN_START_IDX"))
    {
      specProtAlignParams.setValue("IDX_START",
                                   ip.getValue("SPECPROTALIGN_START_IDX"));
    }
    if (ip.exists("SPECPROTALIGN_END_IDX"))
    {
      specProtAlignParams.setValue("IDX_END",
                                   ip.getValue("SPECPROTALIGN_END_IDX"));
    }

    specProtAlignParams.setValue("CLUSTERED_SPECTRUM_FILENAME", SPECS_MS_MGF_FILE);

    AAJumps jumps(1); // Amino acid masses
    jumps.saveJumps(getProjPath(ip, "homology/specprotalign_aa_masses.txt").c_str());

    specProtAlignParams.addIfExists(ip, "FASTA_DATABASE");
    specProtAlignParams.addIfExists(ip, "AMINO_ACID_MASSES");
    specProtAlignParams.addIfExists(ip, "BLOSUM_PENALTY_FILE");

    specProtAlignParams.setValue("INPUT_ABINFO", 
        getProjPath(ip, "assembly/component_info.bin", true));
    specProtAlignParams.setValue("MODS_PENALTY_FILE", 
        getProjPath(ip, "homology/specprotalign_mod_penal.txt", true));
    specProtAlignParams.setValue("KNOWN_MODS_FILE", 
        getProjPath(ip, "homology/specprotalign_known_mods.txt", true));
    specProtAlignParams.setValue("CLEAVAGE_PENALTY_FILE", 
        getProjPath(ip, "homology/specprotalign_cleave_pen.txt", true));
    specProtAlignParams.setValue("SCAN_SPECIFIC_PENALTIES_FILE",
        getProjPath(ip, "homology/specprotalign_scan_specific_pen.txt", true));
    
    specProtAlignParams.setValue("INPUT_TAGS_TARGET", getProjPath(ip, "homology/contig_star_tags_tgt.txt", true));
    specProtAlignParams.setValue("INPUT_TAGS_DECOY", getProjPath(ip, "homology/contig_star_tags_dec.txt", true));

    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_TGT","homology/spectrum_midx_tgt.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_TGT","homology/spectrum_tgt.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM_TGT","homology/spectrum_psm_tgt.txt");
    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_DEC","homology/spectrum_midx_dec.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_DEC","homology/spectrum_dec.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM_DEC","homology/spectrum_psm_dec.txt");

    specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX","homology/spectrum_midx.pklbin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS","homology/spectrum.pklbin");
    specProtAlignParams.setValue("OUTPUT_SPECTRA_SPRINKLED","homology/spectrum_sprinkled.pklbin");
    specProtAlignParams.setValue("OUTPUT_PSM","homology/spectrum_psm.txt");

    specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS","homology/spectrum_mp.bin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS_ALL","homology/spectrum_mp_all.bin");
    specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_IDX","homology/spectrum_indices.bin");

    specProtAlignParams.setValue("OUTPUT_PSM_MOD_FREQS","homology/spectrum_psm_mod_freqs.txt");
    specProtAlignParams.setValue("OUTPUT_PSM_FDR","homology/spectrum_psm_fdr.txt");

    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECPROB");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECTRA");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TIME");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_ANNO");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPRINKLE");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SINGLESPECTRUM");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_MODSPECTRUM");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPLIT");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_MERGE");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_RANGE");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS");
    specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TAG_SEEDING");
    specProtAlignParams.addIfExists(ip, "DEBUG_FDR_PEPTIDE");

    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_RANGE");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_RANGE2");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_AAS");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_SPECS");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN1");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN2");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN3");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN4");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_EXACT");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_NTERM");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP_ANNO");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_GAP_CACHE");
    specProtAlignParams.addIfExists(ip, "PENALTY_ALIGN_DEBUG_ALIGN_LARGE_GAP");

    specProtAlignParams.writeToFile(getProjPath(ip,
                                                "debug_specprotalign.params"));

    DEBUG_VAR(contigSpectra.size());

    DEBUG_TRACE;
    ExecSpecProtAlignTgtDecoy moduleSpecProtAlignTgtDecoy(specProtAlignParams,
                                                          &contigSpectra,
                                                          &prmSpectra,
                                                          &db,
                                                          &dbDecoy,
                                                          &psmTag,
                                                          &psmTagDecoy,
                                                          &penaltyMatrixBlosum,
                                                          &penaltyMatrixMods,
                                                          &scanSpecificPenalties,
                                                          0x0,
                                                          0x0,
                                                          &matchedSpectraAll,
                                                          &matchedSpectra,
                                                          &psmSet,
                                                          &psmSetDecoy,
                                                          &psmSetFdr);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleSpecProtAlignTgtDecoy.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleSpecProtAlignTgtDecoy.invoke();

    // Test for return status
    TEST_RETURN_STATUS("ExecSpecProtAlignTgtDecoy");

    DEBUG_VAR(matchedSpectra.size());

    returnStatus = moduleSpecProtAlignTgtDecoy.saveOutputData();
    DEBUG_VAR(returnStatus);
    TEST_SAVE_OUPUT_DATA("ExecSpecProtAlignTgtDecoy");

    return true;
  }

} // end name specnets

