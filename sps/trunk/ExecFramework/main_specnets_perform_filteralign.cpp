//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_perform_filteralign.h"

// Module Includes
#include "AlignmentPenaltyBased.h"
#include "Logger.h"
#include "ExecAlignment.h"
#include "ExecFilterAligns.h"
#include "ExecFilterContigPairs.h"
#include "ExecFilterPairs.h"
#include "ExecFilterStarPairs.h"
#include "FileUtils.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
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

namespace specnets
{
  //-----------------------------------------------------------------------------
  bool performFilterPairs(ParameterList & ip,
                          SpecSet & inputSpectra,
                          SpecSet & inputSpectraMS2,
                          SpectrumPairSet & filteredPairs,
                          vector<TwoValues<float> > & ratios,
                          vector<TwoValues<float> > & means,
                          vector<float> & varTerms,
                          list<vector<float> > & alignStats,
                          vector<vector<float> > & specStats,
                          std::vector<unsigned int> & idxKept,
                          std::vector<TwoValues<float> > & pvalues,
                          bool gridExecutionFlag,
                          bool resume)
  {
    ParameterList filterpairsParams;
    filterpairsParams.addIfExists(ip, "AMINO_ACID_MASSES");
    filterpairsParams.addIfExists(ip, "TOLERANCE_PEAK");
    filterpairsParams.addIfExists(ip, "TOLERANCE_PM");
    filterpairsParams.addIfExists(ip, "PARTIAL_OVERLAPS");
    filterpairsParams.addIfExists(ip, "AA_DIFF_COUNT");
    filterpairsParams.addIfExists(ip, "MIN_OVERLAP_AREA");
    filterpairsParams.addIfExists(ip, "MIN_RATIO");
    filterpairsParams.addIfExists(ip, "MIN_SHIFT");
    filterpairsParams.addIfExists(ip, "MAX_SHIFT");
    filterpairsParams.addIfExists(ip, "PAIRS_MATCH_MODE");
    filterpairsParams.addIfExists(ip, "PAIRS_MIN_COSINE");
    filterpairsParams.addIfExists(ip, "PAIRS_PROJECTED_COSINE");

    filterpairsParams.addIfExists(ip, "GRID_TYPE");
    filterpairsParams.addIfExists(ip, "GRID_NUMNODES");
    filterpairsParams.addIfExists(ip, "GRID_NUMCPUS");
    filterpairsParams.addIfExists(ip, "GRID_EXE_DIR");
    filterpairsParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
    filterpairsParams.addIfExists(ip, "GRID_PARAMS");
    
    filterpairsParams.addIfExists(ip, "ALIGNGF");
    filterpairsParams.addIfExists(ip, "ALIGNGF_MAX_PVALUE");
    filterpairsParams.addIfExists(ip, "INSTRUMENT_TYPE");
    filterpairsParams.addIfExists(ip, "CHARGE_SEPARATE_PAIRS");

    filterpairsParams.addIfExists(ip, "TAG_FILTER");
    filterpairsParams.addIfExists(ip, "INPUT_TAGS");
    filterpairsParams.addIfExists(ip, "MAX_TAG_SIZE");
    
    filterpairsParams.addIfExists(ip, "MAX_PVALUE");
    filterpairsParams.addIfExists(ip, "RESOLUTION");
    filterpairsParams.addIfExists(ip, "TAGS_MATCH_FLANK");
    filterpairsParams.addIfExists(ip, "TAGS_MATCH_COUNT");

    filterpairsParams.addIfExists(ip, "PROJECT_DIR");
    filterpairsParams.addIfExists(ip, "DEBUG_FILTERPAIRS_SPLIT");

    filterpairsParams.addIfDoesntExist("USE_MIN_DIST_57", "1");
    filterpairsParams.addIfDoesntExist("SPEC_TYPE_MSMS", "0");

    filterpairsParams.setValue("MIN_MATCHED_PEAKS", ip.getValue("MIN_MATCHED_PEAKS"));
    filterpairsParams.setValue("MAX_SHIFT", ip.getValue("MAX_MOD_MASS"));

    filterpairsParams.setValue("GRID_DATA_DIR", getProjPath(ip, ALIGNS_DIR));
    filterpairsParams.setValue("GRID_DATA_DIR_PARAMS", getProjPath(ip, ALIGNS_DIR));
    filterpairsParams.setValue("GRID_DATA_DIR_IN", getProjPath(ip, ALIGNS_DIR));
    filterpairsParams.setValue("GRID_DATA_DIR_OUT", getProjPath(ip, ALIGNS_DIR));
    filterpairsParams.setValue("GRID_DATA_DIR_INTERMEDIATE", getProjPath(ip, ALIGNS_DIR));

    filterpairsParams.setValue("OUTPUT_SPECTRA_PATH",
                                getProjPath(ip, "aligns"));
    filterpairsParams.setValue("OUTPUT_ALIGNS", "pairs_raw.bin");
    filterpairsParams.setValue("OUTPUT_MEANS", "means.bin");
    filterpairsParams.setValue("OUTPUT_VARIANCE", "vars.bin");
    filterpairsParams.setValue("OUTPUT_RATIOS", "ratios.bin");

    filterpairsParams.setValue("DEBUG_PARAMS", "1");
    filterpairsParams.writeToFile("debug_filterpairs.params");

    DEBUG_TRACE;

    stringstream aux;
    filterpairsParams.print(aux);
    DEBUG_MSG(aux.str());

    DEBUG_TRACE;

    if (filterpairsParams.getValue("PAIRS_MATCH_MODE", "") == "cosine")
    {
      filterpairsParams.setValue("OUTPUT_ALIGNS",
                                 getProjPath(ip, "./aligns/pairs_cosine.bin"));
      /*
       for (unsigned int i = 0; i < inputSpectraMS2.size(); i++)
       {
       for (unsigned int j = 0; j < inputSpectraMS2[i].size(); j++)
       inputSpectraMS2[i][j][1] = sqrt(inputSpectraMS2[i][j][1]);
       inputSpectraMS2[i].normalize2();
       }
       */
    }

    ExecFilterPairs moduleFilterPairs(filterpairsParams,
                                      &inputSpectra,
                                      &inputSpectraMS2,
                                      &filteredPairs,
                                      &ratios,
                                      &means,
                                      &varTerms,
                                      &alignStats,
                                      &specStats,
                                      &idxKept,
                                      &pvalues);

    string errorString;
    bool isValid = moduleFilterPairs.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

    bool returnStatus;
    if (!ip.exists("GRID_NUMNODES") or ip.getValueInt("GRID_NUMNODES") <= 0)
    {
      returnStatus = moduleFilterPairs.invoke();
    }
    else
    {
      DEBUG_TRACE;
      int numNodes = ip.getValueInt("GRID_NUMNODES");

      string gridType = ip.getValue("GRID_TYPE");
      if (gridType == "pbs")
      {
        ParallelPbsExecution exec(&moduleFilterPairs,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "sge")
      {
        ParallelSgeExecution exec(&moduleFilterPairs,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "debug")
      {
        ParallelSgeExecution exec(&moduleFilterPairs,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume,
            true);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "threaded")
      {
        ParallelThreadedExecution exec(&moduleFilterPairs);
        returnStatus = exec.invoke(numNodes);
      }
    }

    // Test for return status
    TEST_RETURN_STATUS("moduleFilterPairs");

    DEBUG_TRACE;

    moduleFilterPairs.saveOutputData();

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performFilterAligns(ParameterList & ip,
                           SpectrumPairSet & filteredPairs,
                           SpecSet & prmSpectra,
                           PeptideSpectrumMatchSet & filterPsmSet,
                           vector<TwoValues<float> > & ratios,
                           vector<TwoValues<float> > & means,
                           vector<float> & varTerms,
                           std::vector<unsigned int> & idxKept,
                           std::vector<TwoValues<float> > & pvalues)
  {
    ParameterList filteralignsParams;
    filteralignsParams.addIfExists(ip, "TOLERANCE_PM");
    filteralignsParams.addIfExists(ip, "TOLERANCE_PEAK");
    filteralignsParams.addIfExists(ip, "MIN_RATIO");
    filteralignsParams.addIfExists(ip, "MAX_PVALUE");
    filteralignsParams.addIfExists(ip, "FILTER_TRIGS");
    filteralignsParams.addIfExists(ip, "ALIGNGF");
    filteralignsParams.addIfExists(ip, "MAX_COMPONENT_SIZE");
    filteralignsParams.addIfExists(ip, "MAX_UNIQUE_MASS_NUMBER");
    filteralignsParams.addIfExists(ip, "FILTERSTARPAIRS_EDGE_FDR");

    filteralignsParams.setValue("INPUT_SPECTRA_PATH",
                                getProjPath(ip, "aligns"));
    filteralignsParams.setValue("INPUT_ALIGNS",   "pairs_raw.bin");
    filteralignsParams.setValue("INPUT_MEANS",    "means.bin");
    filteralignsParams.setValue("INPUT_VARIANCE", "vars.bin");
    filteralignsParams.setValue("INPUT_RATIOS",   "ratios.bin");

    filteralignsParams.setValue("OUTPUT_SPECTRA_PATH",
                                getProjPath(ip, "aligns"));
    filteralignsParams.setValue("OUTPUT_ALIGNS", "pairs.bin");

    filteralignsParams.setValue("DEBUG_PARAMS", "1");
    filteralignsParams.writeToFile("debug_filteraligns.params");

    DEBUG_TRACE;

    stringstream aux;
    filteralignsParams.print(aux);
    DEBUG_MSG(aux.str());

    DEBUG_TRACE;

    ExecFilterAligns moduleFilterAligns(filteralignsParams,
                                        &filteredPairs,
                                        &prmSpectra,
                                        &filterPsmSet,
                                        &ratios,
                                        &means,
                                        &varTerms,
                                        &idxKept,
                                        &pvalues);

    string errorString;
    bool isValid = moduleFilterAligns.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

    bool returnStatus = moduleFilterAligns.invoke();

    // Test for return status
    TEST_RETURN_STATUS("moduleFilterAligns");

    DEBUG_TRACE;

    // Save output data
    moduleFilterAligns.saveOutputData();

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performAlignment(ParameterList & ip,
                        AAJumps & jumps,
                        SpecSet &inputSpectra,
                        SpectrumPairSet &inputPairs,
                        SpecSet &pairAlignments,
                        SpecSet &starSpectraOnly,
                        SpecSet &starSpectra,
                        vector<unsigned int> &alignedSpectra)
  {
    ParameterList alignParams;
    alignParams.addIfExists(ip, "TOLERANCE_PEAK");
    alignParams.addIfExists(ip, "TOLERANCE_PM");
    alignParams.addIfExists(ip, "RESOLUTION");
    alignParams.addIfExists(ip, "PARTIAL_OVERLAPS");
    alignParams.addIfExists(ip, "MAX_AA_JUMP");
    alignParams.addIfExists(ip, "PENALTY_PTM");
    alignParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

    alignParams.addIfExists(ip, "PROJECT_DIR");
    alignParams.setValue("OUTPUT_SPECTRA_PATH",
                         getProjPath(ip, SPECTRA_DIR));

    alignParams.setValue("PENALTY_PTM_PEAKS", "-1.0");
    alignParams.setValue("OUTPUT_SPECS", "aligns_specs.pklbin");
    alignParams.setValue("OUTPUT_STARS", "stars_only.pklbin");
    alignParams.setValue("OUTPUT_STARS_INDEX", "stars_indices.bin");
    alignParams.setValue("OUTPUT_STARS_ALL", "stars.pklbin");
    if (ip.exists("AMINO_ACID_MASSES"))
    {
      alignParams.setValue("AMINO_ACID_MASSES",
                           "../" + ip.getValue("AMINO_ACID_MASSES"));
    }

    alignParams.setValue("DEBUG_PARAMS", "1");
    alignParams.writeToFile("debug_align.params");

    DEBUG_TRACE;
    ExecAlignment moduleAlignment(alignParams,
                                  &inputSpectra,
                                  &inputPairs,
                                  &pairAlignments,
                                  &starSpectraOnly,
                                  &starSpectra,
                                  &alignedSpectra);

    string errorString;
    bool isValid = moduleAlignment.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

  #if 1
    bool returnStatus = moduleAlignment.invoke();
  #else
    DEBUG_TRACE;
    ParallelThreadedExecution exec(&moduleAlignment);
    bool returnStatus = exec.invoke(1, 1);
  #endif

    // Test for return status
    TEST_RETURN_STATUS("ExecAlignment");

    DEBUG_TRACE;
    // Save output data
    moduleAlignment.saveOutputData();

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performFilterStarPairs(ParameterList & ip,
                              SpectrumPairSet &inputPairs,
                              SpecSet & starSpectra,
                              vector<vector<float> > & ratios,
                              SpecSet & matchedPeaks)
  {
    ParameterList filterstarpairsParams;
    filterstarpairsParams.addIfExists(ip, "TOLERANCE_PEAK");
    filterstarpairsParams.addIfExists(ip, "TOLERANCE_PM");
    filterstarpairsParams.addIfExists(ip, "MAX_MOD_MASS");
    filterstarpairsParams.addIfExists(ip, "PARTIAL_OVERLAPS");
    filterstarpairsParams.addIfExists(ip, "MIN_RATIO");
    filterstarpairsParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
    filterstarpairsParams.addIfExists(ip, "MAX_AA_JUMP");
    filterstarpairsParams.addIfExists(ip, "PENALTY_PTM");
    filterstarpairsParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

    filterstarpairsParams.addIfExists(ip, "PROJECT_DIR");
    filterstarpairsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    filterstarpairsParams.addIfExists(ip, "INPUT_FILTER_PSMS");
    filterstarpairsParams.addIfExists(ip, "INPUT_FILTER_MSGFDB_PSMS");
    filterstarpairsParams.addIfExists(ip, "ALIGNGF");
    filterstarpairsParams.addIfExists(ip, "FILTERSTARPAIRS_BY_PEAKINTENSITY");

    filterstarpairsParams.setValue("SPEC_TYPE_MSMS", "0");
    filterstarpairsParams.setValue("OUTPUT_ALIGNS",
                                   getProjPath(ip, "aligns/pairs_stars.bin"));

    filterstarpairsParams.setValue("DEBUG_PARAMS", "1");
    filterstarpairsParams.writeToFile(getProjPath(ip,
                                                  "debug_filterstarpairs.params"));

    DEBUG_TRACE;
    ExecFilterStarPairs moduleFilterStarPairs(filterstarpairsParams,
                                              &inputPairs,
                                              &starSpectra,
                                              &ratios,
                                              &matchedPeaks);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleFilterStarPairs.validateParams(errorString);
    TEST_VALID;

  #if 1
    bool returnStatus = moduleFilterStarPairs.invoke();
  #else
    ParallelSgeExecution exec(&moduleFilterStarPairs);
    //ParallelThreadedExecution exec(&moduleFilterStarPairs);
    bool returnStatus = exec.invoke(1, 1);
  #endif

    // Test for return status
    TEST_RETURN_STATUS("ExecFilterStarPairs");

    // Save output data
    returnStatus = moduleFilterStarPairs.saveOutputData();
    //  TEST_SAVE_OUPUT_DATA("ExecFilterStarPairs");

    return true;
  }

} // end name specnets

