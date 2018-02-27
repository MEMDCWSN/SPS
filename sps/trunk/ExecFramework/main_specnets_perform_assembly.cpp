//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_perform_assembly.h"

// Module Includes
#include "Logger.h"
#include "ExecAssembly.h"
#include "ExecFilterContigPairs.h"
#include "ExecMetaAssembly.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"

// Specnets Includes
#include "clusters.h"
#include "ClusterData.h"
#include "abruijn.h"
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
  bool performAssembly(ParameterList & ip,
                       SpecSet & starSpectra,
                       SpectrumPairSet & starPairs,
                       Clusters & contigShifts,
                       abinfo_t & contigAbinfo)
  {
    ParameterList assemblyParams;
    assemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
    assemblyParams.addIfExists(ip, "TOLERANCE_PM");
    assemblyParams.addIfExists(ip, "MAX_MOD_MASS");
    assemblyParams.addIfExists(ip, "MAX_AA_JUMP");
    assemblyParams.addIfExists(ip, "PENALTY_PTM");
    assemblyParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
    assemblyParams.addIfExists(ip, "PENALTY_SAME_VERTEX");
    assemblyParams.addIfExists(ip, "PENALTY_PTM");
    assemblyParams.addIfExists(ip, "PARALLEL_PATHS");
    assemblyParams.addIfExists(ip, "ADD_ENDPOINTS");

    assemblyParams.addIfExists(ip, "PROJECT_DIR");

    assemblyParams.setValue("NO_SEQUENCING", "0");
    assemblyParams.setValue("EDGE_SCORE_TYPE", "1");
    assemblyParams.setValue("SPEC_TYPE_MSMS", "0");
    assemblyParams.setValue("GRAPH_TYPE", "2");
    assemblyParams.setValue("OUTPUT_COMPLETE_ABRUIJN", "");
    assemblyParams.setValue("PATH_MIN_PEAKS",
                            ip.getValue("SPSPATH_MIN_NUM_PEAKS"));
    assemblyParams.setValue("PATH_MIN_SPECS",
                            ip.getValue("SPSPATH_MIN_NUM_SPECS"));
    assemblyParams.setValue("MIN_EDGES_TO_COMPONENT",
                            ip.getValue("SPS_MIN_EDGES_TO_COMPONENT"));

    assemblyParams.setValue("OUTPUT_SPECTRA_PATH", getProjPath(ip, "assembly"));
    assemblyParams.setValue("OUTPUT_CLUSTERS", "path_spectra_as_cluster.txt");
    assemblyParams.setValue("OUTPUT_SPECS", "sps_seqs.pklbin");
    assemblyParams.setValue("OUTPUT_ABINFO", "component_info.bin");

    assemblyParams.setValue("DEBUG_PARAMS", "1");
    assemblyParams.writeToFile(getProjPath(ip, "debug_assembly.params"));

    DEBUG_TRACE;
    ExecAssembly moduleAssembly(assemblyParams,
                                &starSpectra,
                                &starPairs,
                                &contigShifts,
                                &contigAbinfo);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleAssembly.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleAssembly.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecAssembly");

    returnStatus = moduleAssembly.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecAssembly");

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performContigAlignment(ParameterList & ip,
                              Clusters & contigs,
                              SpectrumPairSet & contigPairs)
  {
    ParameterList alignParams;
    alignParams.addIfExists(ip, "TOLERANCE_PEAK");
    alignParams.addIfExists(ip, "TOLERANCE_PM");
    alignParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
    alignParams.addIfExists(ip, "MIN_METACONTIG_SCORE");

    alignParams.setValue("INPUT_SPECTRA",
                         getProjPath(ip, "assembly/sps_seqs.pklbin"));

    alignParams.setValue("OUTPUT_CONTIG_ALIGNS",
                         getProjPath(ip, "aligns/sps_seqs_aligns.bin"));

    alignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    DEBUG_TRACE;
    ExecFilterContigPairs moduleAlign(alignParams,
                                      &contigs.consensus,
                                      &contigPairs);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleAlign.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleAlign.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecFilterContigPairs");

    returnStatus = moduleAlign.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecFilterContigPairs");

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performMetaAssembly(ParameterList & ip,
                           Clusters & contigs,
                           SpectrumPairSet & contigAligns,
                           abinfo_t & contigAbruijn,
                           SpecSet & starSpectra,
                           Clusters & outputMetaContigs,
                           abinfo_t & metaContigAbruijn)
  {
    ParameterList assemblyParams;
    assemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
    assemblyParams.addIfExists(ip, "TOLERANCE_PM");
    assemblyParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
    assemblyParams.addIfExists(ip, "MIN_METACONTIG_SCORE");
    assemblyParams.addIfExists(ip, "MIN_METACONTIG_SIZE");

    assemblyParams.setValue("OUTPUT_SPECTRA_PATH", getProjPath(ip, "assembly"));

    assemblyParams.setValue("OUTPUT_CLUSTERS", "meta_path_spectra_as_cluster.txt");
    assemblyParams.setValue("OUTPUT_SPECTRA", "sps_seqs.pklbin");
    assemblyParams.setValue("OUTPUT_ABINFO", "component_info.bin");

    assemblyParams.writeToFile(getProjPath(ip, "debug_metaAssembly.params"));

    DEBUG_TRACE;
    ExecMetaAssembly moduleAssembly(assemblyParams,
                                    &contigs.consensus,
                                    &contigAligns,
                                    &contigAbruijn,
                                    &starSpectra,
                                    &outputMetaContigs,
                                    &metaContigAbruijn);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleAssembly.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleAssembly.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecMetaAssembly");

    returnStatus = moduleAssembly.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecMetaAssembly");

    return true;
  }


} // end name specnets

