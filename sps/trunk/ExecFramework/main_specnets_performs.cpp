//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_performs.h"

// Module Includes
#include "AlignmentPenaltyBased.h"
#include "Logger.h"
#include "ExecMsCluster.h"
#include "ExecGenoMS.h"
#include "ExecDeconvoluteMS2.h"
#include "ExecReportSpsplot.h"
#include "ExecReportSPSStats.h"
#include "ExecPrmScoring.h"
#include "ExecPrmClustering.h"
#include "ExecProtProtAlign.h"
#include "ExecHomologyAssembly.h"
#include "ExecMainSpecnets.h"
#include "FileUtils.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
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

namespace specnets
{
  /**
   * Performs the MS2 deconvolution step. Since this occurs before clustering, we need to keep
   *    the spectra in separate files.
   * @param ip list of user parameters
   * @param inputSpecs set of input SpecSets
   * @param outputSpecs set of output SpecSets
   * @return true if successful, false if not
   */
  bool performMS2Deconv(ParameterList & ip,
                        vector<string>& inputSpecs,
                        vector<string>& outputSpecs)
  {
    ParameterList deconvParams;
    deconvParams.addIfExists(ip, "MAX_KLDiv");
    deconvParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    // Need the input iso envelope. It should be in the execution directory
    string isoModelPath = getPath(ip.getValue("EXE_DIR"), "resources", false);
    isoModelPath = getPath(isoModelPath, DEFAULT_DECONV_MODEL, false);

    outputSpecs.resize(inputSpecs.size());
    for (unsigned int i = 0; i < inputSpecs.size(); i++)
    {
      FilenameManager mngr(inputSpecs[i]);
      string outFile = getPath(mngr.path, mngr.filename, false);
      outFile += "_z.pklbin";
      outputSpecs[i] = outFile;
      deconvParams.setValue("OUTPUT_SPECTRA", outputSpecs[i]);
      deconvParams.setValue("INPUT_SPECTRA", inputSpecs[i]);
      deconvParams.setValue("INPUT_ISO_ENV", isoModelPath);
      IsoEnvelope isoModel;

      SpecSet inputSpecs;
      SpecSet outputSpecs;
      ExecDeconvoluteMS2 module(deconvParams,
                                &inputSpecs,
                                &isoModel,
                                &outputSpecs);

      if (!module.loadInputData())
      {
        return false;
      }

      if (!module.invoke())
      {
        return false;
      }

      if (!module.saveOutputData())
      {
        return false;
      }

      DEBUG_VAR(outputSpecs.size());
    }

    string deconvFilesPath = 
        getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_DECONV_FILES_LIST;

    if (!writeFileIndex(deconvFilesPath, outputSpecs))
    {
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performMsCluster(ParameterList & ip,
                        vector<string>& inputFilesList,
                        SpecSet &clustSpectra)
  {
    DEBUG_TRACE;

    ParameterList msclusterParams;
    msclusterParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
    if (ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST) != CLUSTER_MSCLUST)
    {
      msclusterParams.setValue("CLUSTER_MIN_SIZE", "0");
    }
    msclusterParams.addIfExists(ip, "EXE_DIR");
    msclusterParams.addIfExists(ip, "TOLERANCE_PEAK");
    msclusterParams.addIfExists(ip, "TOLERANCE_PM");
    msclusterParams.addIfExists(ip, "TOLERANCE_PM_PPM");
    msclusterParams.addIfExists(ip, "CLUSTER_MODEL");
    msclusterParams.addIfExists(ip, "MIN_SPECTRUM_QUALITY");
    msclusterParams.addIfExists(ip, "RESET_MSCLUSTER_PMS");
    msclusterParams.addIfExists(ip, "MSCLUSTER_MIX_PROB");
    msclusterParams.addIfExists(ip, "CLUST_RANK_FILTER");
    msclusterParams.addIfExists(ip, "REMOVABLE_BAD_CLUSTER_SIZE");

    msclusterParams.setValue("PROJECT_DIR", getProjPath(ip, SPECTRA_DIR));

    string outputPaths = SPECS_MS_FILE + ";" + SPECS_MS_MGF_FILE;
    msclusterParams.setValue("OUTPUT_SPECTRA", outputPaths);
    msclusterParams.setValue("OUTPUT_CLUSTERS", SPECS_MS_CLUST_FILE);

    DEBUG_TRACE;
    stringstream aux;
    msclusterParams.print(aux);
    DEBUG_MSG(aux.str());
    ClusterSet clusters;

    ExecMsCluster moduleMsCluster(msclusterParams,
        &inputFilesList,
        &clusters,
        &clustSpectra);

    if (!moduleMsCluster.invoke())
    {
      return false;
    }

    if (!moduleMsCluster.saveOutputData())
    {
      return false;
    }

    DEBUG_VAR(clustSpectra.size());

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performScoring(ParameterList & ip,
                      SpecSet & inputSpectra,
                      SpecSet & outputSpectra,
                      bool gridExecutionFlag,
                      bool resume)
  {
    DEBUG_TRACE;

    ParameterList scoringParams;
    scoringParams.addIfExists(ip, "EXE_DIR");
    scoringParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");
    scoringParams.addIfExists(ip, "PEPNOVO_EXE_DIR");
    scoringParams.addIfExists(ip, "PEPNOVO_MODEL_DIR");
    scoringParams.addIfExists(ip, "TOLERANCE_PEAK");
    scoringParams.addIfExists(ip, "TOLERANCE_PM");
    scoringParams.addIfExists(ip, "TOLERANCE_PM_PPM");
    scoringParams.addIfExists(ip, "MIN_SPECTRUM_QUALITY");
    scoringParams.addIfExists(ip, "GUESS_CHARGE");
    scoringParams.addIfExists(ip, "CORRECT_PM");
    scoringParams.addIfExists(ip, "PEPNOVO_MODEL");
    scoringParams.addIfExists(ip, "INSTRUMENT_TYPE");

    scoringParams.addIfExists(ip, "PEPNOVO_MODEL");
    scoringParams.addIfExists(ip, "SKIP_PEPNOVO_INVOKE");
    scoringParams.addIfExists(ip, "PEPNOVO_PTMS");
    scoringParams.addIfExists(ip, "PEPNOVO_OUTDIR");

    scoringParams.addIfExists(ip, "GRID_TYPE");
    scoringParams.addIfExists(ip, "GRID_NUMNODES");
    scoringParams.addIfExists(ip, "GRID_NUMCPUS");
    scoringParams.addIfExists(ip, "GRID_EXE_DIR");
    scoringParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
    scoringParams.addIfExists(ip, "GRID_PARAMS");
    
    if (ip.getValue("PAIRS_MATCH_MODE", "") != "cosine") {
      scoringParams.addIfExists(ip, "SCORING_COMPRESS_RESULTS");
      scoringParams.addIfExists(ip, "TOLERANCE_PEAK");
      scoringParams.addIfExists(ip, "SPEC_TYPE_MSMS");
    }

    scoringParams.addIfExists(ip, "PROJECT_DIR");
    DEBUG_VAR(ip.getValue("PROJECT_DIR"));
    scoringParams.setValue("GRID_DATA_DIR_OUT", 
                           getProjPath(ip, SPECTRA_DIR));
    string gridDataDir = getProjPath(ip, SPECTRA_GRID_SCORING_DIR);
    DEBUG_VAR(gridDataDir);
    scoringParams.setValue("GRID_DATA_DIR", gridDataDir);
    scoringParams.setValue("GRID_DATA_DIR_IN", gridDataDir);
    scoringParams.setValue("GRID_DATA_DIR_PARAMS", gridDataDir);
    scoringParams.setValue("GRID_DATA_DIR_INTERMEDIATE", gridDataDir);

    /*scoringParams.setValue("INPUT_SPECTRA",
     getProjPath(ip, "spectra/specs_ms.pklbin"));
     */
    // Free memory for PepNovo
    inputSpectra.clear();

    string ms2SpecPath = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_MS_FILE;
    scoringParams.setValue("PEPNOVO_OUTDIR", getProjPath(ip, SPECTRA_DIR));
    scoringParams.setValue("INPUT_SPECTRA", ms2SpecPath);
    scoringParams.setValue("OUTPUT_SPECTRA", SPECS_SCORED_FILE);
    scoringParams.setValue("PEPNOVO_INPUT_MGF",
        getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_PEPNOVO_MGF_FILE);
    scoringParams.setValue("PEPNOVO_OUTPUT_PRMS", SPECS_PEPNOVO_PRMS_FILE);
    scoringParams.setValue("ENFORCE_DA_TOL", "1");

    scoringParams.setValue("DEBUG_PARAMS", "1");
    //  specProtAlignParams.writeToFile(getProjPath(ip, "debug_specprotalign.params"));
    scoringParams.writeToFile(getProjPath(ip, "debug_scoring.params"));

    ExecPrmScoring moduleScoring(scoringParams, &inputSpectra, &outputSpectra);

    string errorString;
    bool isValid = moduleScoring.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

    //moduleScoring.loadInputData();

    if (!moduleScoring.loadInputData())
    {
      return false;
    }

    DEBUG_TRACE;

    bool returnStatus;
    if (!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") <= 0 ||
        ip.exists("NO_GRID_SCORING"))
    {
      returnStatus = moduleScoring.invoke();
    }
    else
    {
      DEBUG_TRACE;
      bool res = mkdir_if_not_exist(gridDataDir.c_str());
      if (res)
      {
        DEBUG_MSG("Made directory \'"<< gridDataDir << "\'");
      }

      int numNodes = ip.getValueInt("GRID_NUMNODES");

      string gridType = ip.getValue("GRID_TYPE");
      if (gridType == "pbs")
      {
        ParallelPbsExecution exec(&moduleScoring,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "sge")
      {
        ParallelSgeExecution exec(&moduleScoring,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "debug")
      {
        ParallelSgeExecution exec(&moduleScoring,
            gridExecutionFlag,
            !gridExecutionFlag,
            resume,
            true);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "threaded")
      {
        ParallelThreadedExecution exec(&moduleScoring);
        returnStatus = exec.invoke(numNodes);
      }
    }

    // Test for return status
    TEST_RETURN_STATUS("moduleScoring");

    DEBUG_TRACE;

    if (!moduleScoring.saveOutputData())
    {
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performPrmClustering(ParameterList & ip,
                            SpecSet & inputSpectra,
                            SpecSet & outputSpectra)
  {
    DEBUG_TRACE;

    ParameterList scoringParams;
    scoringParams.addIfExists(ip, "EXE_DIR");
    scoringParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");
    scoringParams.addIfExists(ip, "TOLERANCE_PEAK");
    scoringParams.addIfExists(ip, "TOLERANCE_PM");
    scoringParams.addIfExists(ip, "TOLERANCE_PM_PPM");
    scoringParams.addIfExists(ip, "CLUSTER_MODEL");
    scoringParams.addIfExists(ip, "MIN_SPECTRUM_QUALITY");
    scoringParams.addIfExists(ip, "INSTRUMENT_TYPE");
    scoringParams.addIfExists(ip, "MIN_SILAC_COSINE");

    scoringParams.setValue("INPUT_CLUSTERS",
        getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_MS_CLUST_FILE);

    if (ip.getValueInt("MERGE_SAME_PREC", 0) > 0
        || (ip.getValueInt("CLUSTER_MIN_SIZE", 0) >= 1
            && ip.getValue("CLUSTER_TOOL", "") == CLUSTER_PRMS))
    {
      scoringParams.addIfExists(ip, "MERGE_SAME_PREC");
      scoringParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
      scoringParams.setValue("OUTPUT_CLUSTERS", SPECS_SCORED_CLUST_FILE);
    }

    scoringParams.addIfExists(ip, "PRM_RANK_FILTER");
    scoringParams.addIfExists(ip, "NUM_CONSECUTIVE");
    scoringParams.addIfExists(ip, "PRM_CLUSTER_RATIO");
    scoringParams.addIfExists(ip, "BOOST_SILAC_PRMS");
    scoringParams.addIfExists(ip, "SILAC_SCAN_RANGE");
    scoringParams.addIfExists(ip, "FILTER_NONSILAC_PRMS");

    if (scoringParams.getValueInt("BOOST_SILAC_PRMS", 0) > 0)
    {
      scoringParams.setValue("OUTPUT_SILAC_PAIRS", 
          getProjPath(ip, ALIGNS_DIR) + "/" + ALIGNS_SILAC_FILE);
    }

    const string unclusteredScoredPath = 
         getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_SCORED_UNCLUST_FILE;
    if (!ExecMergeConvert::saveSpecsetMultiple(unclusteredScoredPath,
            &inputSpectra))
    {
      ERROR_MSG("Failed to save to \'" << unclusteredScoredPath
          << "\'");
      return false;
    }

    scoringParams.addIfExists(ip, "PROJECT_DIR");
    /*scoringParams.setValue("INPUT_SPECTRA",
     getProjPath(ip, "spectra/specs_ms.pklbin"));
     */
    // Free memory for PepNovo
    inputSpectra.clear();

    string ms2SpecPath = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_MS_FILE;
    scoringParams.setValue("PEPNOVO_OUTDIR", getProjPath(ip, SPECTRA_DIR));
    scoringParams.setValue("INPUT_SPECTRA_MS", ms2SpecPath);
    scoringParams.setValue("INPUT_SPECTRA", unclusteredScoredPath);
    scoringParams.setValue("OUTPUT_SPECTRA_DIR", getProjPath(ip, SPECTRA_DIR));
    scoringParams.setValue("OUTPUT_SPECTRA", SPECS_SCORED_FILE);

    //  specProtAlignParams.writeToFile(getProjPath(ip, "debug_specprotalign.params"));
    scoringParams.writeToFile(getProjPath(ip, "debug_prmclustering.params"));

    ClusterSet outputClusters;

    ExecPrmClustering moduleClustering(scoringParams, &outputSpectra, &outputClusters);

    string errorString;
    bool isValid = moduleClustering.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

    //moduleScoring.loadInputData();

    if (!moduleClustering.loadInputData())
    {
      return false;
    }

    DEBUG_TRACE;

    bool returnStatus = moduleClustering.invoke();

    // Test for return status
    TEST_RETURN_STATUS("moduleClustering");

    DEBUG_TRACE;

    if (!moduleClustering.saveOutputData())
    {
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performStatProtSeqs(ParameterList & ip)
  {
    if (!ip.exists("FASTA_DATABASE") || !ip.exists("INPUT_SPEC_IDS"))
    {
      return false;
    }

    ParameterList statsParams;
    statsParams.addIfExists(ip, "TOLERANCE_PEAK");
    statsParams.setValue("INPUT_CONTIGS",
                         getProjPath(ip, "assembly/sps_seqs.pklbin"));
    statsParams.setValue("INPUT_STARS", getProjPath(ip, "spectra/stars.pklbin"));
    statsParams.setValue("INPUT_CONTIG_ABINFO",
                         getProjPath(ip, "assembly/component_info.bin"));
    statsParams.setValue("INPUT_MIDX",
                         getProjPath(ip, "homology/contigs_midx.pklbin"));
    statsParams.setValue("INPUT_MP", getProjPath(ip, "homology/contigs_mp.bin"));
    statsParams.setValue("INPUT_REF_INDICES",
                         getProjPath(ip, "spectra/contigs_indices.bin"));
    statsParams.setValue("INPUT_FASTA", ip.getValue("FASTA_DATABASE"));
    string aux = ip.getValue("EXE_DIR");
    aux += "/dancik_model.txt";
    statsParams.setValue("INPUT_ION_TYPES", aux);
    statsParams.addIfExists(ip, "INPUT_SPEC_IDS");
    statsParams.addIfExists(ip, "INPUT_CLUST_SPEC_IDS");
    statsParams.setValue("SPEC_ID_FORMAT",
                         ip.getValue("SPEC_ID_FORMAT", "MSGFDB"));

    statsParams.setValue("INPUT_FILE_INDEX",
        getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_USER_INPUT_FILES_LIST);
    statsParams.setValue("SPS_PROJECT_DIR", "");

    struct stat buf;
    if (stat(getProjPath(ip, "spectra/clusterData.bin").c_str(), &buf) != -1)
    {
      statsParams.setValue("INPUT_CLUSTERS_DIR", getProjPath(ip, ""));
    }
    statsParams.setValue("INPUT_FILE_MAPPING",
        getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_INPUT_MAPPING);

    statsParams.setValue("INPUT_SCAN_REF_FILES",
        getProjPath(ip, SPECTRA_DIR) + "/" + DEFAULT_BIN_FILES_FILENAME);

    if (ip.exists("STATS_MIN_CONTIG_AA_TAG"))
    {
      statsParams.setValue("MIN_CONTIG_AA_TAG",
                           ip.getValue("STATS_MIN_CONTIG_AA_TAG"));
    }

    if (ip.exists("STATS_MIN_CONTIG_DB_MP"))
    {
      statsParams.setValue("MIN_CONTIG_DB_MP",
                           ip.getValue("STATS_MIN_CONTIG_DB_MP"));
    }

    if (ip.exists("STATS_ENDS_CHOP"))
    {
      statsParams.setValue("ENDS_CHOP", ip.getValue("STATS_ENDS_CHOP"));
    }

    if (ip.exists("STATS_TARGET_PROTEINS"))
    {
      statsParams.setValue("TARGET_PROTEINS",
                           ip.getValue("STATS_TARGET_PROTEINS"));
    }
    else
    {
      statsParams.setValue("TARGET_PROTEINS", "");
    }

    statsParams.setValue("TABLE_DELIM", "\t");
    statsParams.setValue("OUTPUT_SPS_STATS_FILE",
                         getProjPath(ip, "ReportData/tableStatsCummulative.tsv"));
    statsParams.setValue("OUTPUT_CONTIG_STATS_FILE",
                         getProjPath(ip, "ReportData/tableStatsContig.tsv"));

    statsParams.writeToFile(getProjPath(ip, "debug_statProtSeqs.params"));

    ExecReportSPSStats statsModule(statsParams);
    if (!statsModule.loadInputData())
    {
      return false;
    }
    if (!statsModule.invoke())
    {
      return false;
    }
    if (!statsModule.saveOutputData())
    {
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performProtProtAlign(ParameterList & ip,
                            unsigned int refProtIdx,
                            set<unsigned int> dbIndexes,
                            DB_fasta & db)
  {
    DEBUG_TRACE;
    ParameterList protProtAlignParams;

    protProtAlignParams.addIfExists(ip, "PROJECT_DIR");
    protProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    //protProtAlignParams.addIfExists(ip,"CLUSTALW_EXE_DIR");
    if (ip.exists("EXE_DIR"))
    protProtAlignParams.setValue("CLUSTALW_EXE_DIR", ip.getValue("EXE_DIR"));

    ostringstream tmp;
    tmp << refProtIdx;
    protProtAlignParams.setValue("REFERENCE_PROTEIN_IDX", tmp.str());

    protProtAlignParams.addIfExists(ip, "CLUSTALW_MINSCORE");
    protProtAlignParams.addIfDoesntExist("CLUSTALW_MINSCORE", "250");

    protProtAlignParams.setValue("CLUSTALW_FILENAME_PREFIX",
        getProjPath(ip, "homology/cwseqs_"));
    protProtAlignParams.setValue("CLUSTALW_INDEX",
        getProjPath(ip, "homology/cwindex.txt"));
    protProtAlignParams.writeToFile("debug_protprotalign.params");

    DEBUG_TRACE;
    vector<string> alnFileNames;
    ExecProtProtAlign moduleProtProtAlign(protProtAlignParams,
        &db,
        &dbIndexes,
        &alnFileNames);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleProtProtAlign.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleProtProtAlign.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecProtProtAlign");
    // Save output data
    returnStatus = moduleProtProtAlign.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecProtProtAlign");

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performHomologyAssembly(ParameterList & ip,
                               SpecSet & spectra,
                               DB_fasta & db,
                               SpecSet & contigShifts,
                               SpecSet & contigMatchedIndices,
                               vector<vector<int> > & contigMatchedProts)
  {
    DEBUG_TRACE;
    ParameterList homologyAssemblyParams;
    homologyAssemblyParams.setValue("INPUT_HOMOLOGIES",
        getProjPath(ip, "homology/cwindex.txt"));
    //  homologyAssemblyParams.setValue("INPUT_SPECS_NAMES",      getProjPath(ip, "homology/ref_sps_names.txt") );
    homologyAssemblyParams.setValue("SPEC_TYPE_MSMS", "0");
    homologyAssemblyParams.setValue("GRAPH_TYPE", "2");
    homologyAssemblyParams.setValue("MIN_CONTIG_SET", "1");
    homologyAssemblyParams.setValue("EDGE_SCORE_TYPE", "1");
    homologyAssemblyParams.setValue("OUTPUT_SPECS",
        getProjPath(ip,
            "homology/homglue_matches.pklbin"));
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MIDX",
        getProjPath(ip,
            "homology/homglue_ref_midx.pklbin"));
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MP",
        getProjPath(ip, "homology/homglue_ref_mp.bin"));
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MIDX",
        getProjPath(ip,
            "homology/homglue_matches_midx.pklbin"));
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MP",
        getProjPath(ip,
            "homology/homglue_matches_mp.bin"));
    homologyAssemblyParams.setValue("OUTPUT_CSV",
        getProjPath(ip, "homglue_matches.txt"));
    homologyAssemblyParams.addIfExists(ip, "RESOLUTION");
    homologyAssemblyParams.addIfExists(ip, "MAX_MOD_MASS");
    homologyAssemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
    homologyAssemblyParams.addIfExists(ip, "TOLERANCE_PM");
    homologyAssemblyParams.addIfExists(ip, "MAX_AA_JUMP");
    homologyAssemblyParams.addIfExists(ip, "PENALTY_PTM");
    homologyAssemblyParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

    homologyAssemblyParams.addIfExists(ip, "PROJECT_DIR");
    homologyAssemblyParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    homologyAssemblyParams.writeToFile(getProjPath(ip, "debug_homology.params"));

    DEBUG_TRACE;
    ExecHomologyAssembly moduleHomologyAssembly(homologyAssemblyParams,
        &spectra,
        &db,
        &contigShifts,
        &contigMatchedIndices,
        &contigMatchedProts);
    DEBUG_TRACE;

    string errorString;
    bool isValid = moduleHomologyAssembly.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = moduleHomologyAssembly.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecHomologyAssembly");

    returnStatus = moduleHomologyAssembly.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecHomologyAssembly");

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performReport(ParameterList & ip)
  {
    DEBUG_MSG("Entering performReport()");

    string execMode = ip.getValue("EXECUTION_MODE", "sps");

    float resolution = ip.getValueDouble("RESOLUTION", 0.1);

    // Get current directory -- needed for absolute path composition
    string currentWorkindDir;
    currentWorkindDir = getCurrentDirectory();

    // auxiliary string used for file name composition using project name
    string aux;

    // Exec Spsplot section

    ParameterList reportSpsplotParams;

    reportSpsplotParams.addIfExists(ip, "PROJECT_DIR");
    reportSpsplotParams.addIfExists(ip, "EXE_DIR");
    reportSpsplotParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    reportSpsplotParams.addIfExists(ip, "REPORT_JOB");
    reportSpsplotParams.addIfExists(ip, "REPORT_USER");
    reportSpsplotParams.addIfExists(ip, "REPORT_PWD");
    // add dynamic report project directory (target upon relocation)
    reportSpsplotParams.addIfExists(ip, "REPORT_SERVER");
    reportSpsplotParams.addIfExists(ip, "REPORT_DIR_SERVER");
    reportSpsplotParams.addIfExists(ip, "REPORT_CELLS_PER_LINE");
    // specify if MS/MS images are shown
    reportSpsplotParams.addIfExists(ip, "REPORT_MSMS_IMAGES");
    // add PDF reports status
    reportSpsplotParams.addIfExists(ip, "REPORT_PDF");

    // add mass shift for report drawing
    reportSpsplotParams.addIfExists(ip, "REPORT_MASS_SHIFT");
    // add mass shift for report drawing
    reportSpsplotParams.addIfExists(ip, "MS2_SCORING_MODEL");
    // add amino acids file
    reportSpsplotParams.addIfExists(ip, "AMINO_ACID_MASSES");
    // add mass shift for report drawing
    //reportSpsplotParams.addIfExists(ip, "REPORT_MASS_SHIFT_PRM");
    // add mass shift for report drawing
    //reportSpsplotParams.addIfExists(ip, "REPORT_MS2_SCORING_MODEL_PRM");

    // add dynamic reports status
    //reportSpsplotParams.addIfExists(ip, "REPORT_DYNAMIC");
    if (ip.getValueInt("REPORT_DYNAMIC") == 1)
      reportSpsplotParams.setValue("REPORT_DYNAMIC", "1");

    if (ip.getValueInt("CLUSTER_MIN_SIZE") == 0
        && ip.getValueInt("MERGE_SAME_PREC", 0) > 1)
      reportSpsplotParams.setValue("NO_CLUSTERS", "1");

    // allow the realign in protein coverage pages
    if (ip.getValueInt("REPORT_REALIGN") == 1)
      reportSpsplotParams.setValue("ALLOW_REALIGN", "1");
    else
      reportSpsplotParams.setValue("ALLOW_REALIGN", "0");

    // allow interactivity with the server
    if (ip.getValueInt("REPORT_NOSERVER") == 1)
      reportSpsplotParams.setValue("DYNAMIC", "0");
    else
      reportSpsplotParams.setValue("DYNAMIC", "1");

    stringstream saux;
    saux << resolution;
    reportSpsplotParams.setValue("RESOLUTION", saux.str());

    int tool = 1;

    DEBUG_VAR(ip.getValueInt("GENOMS_FLAG"));
    DEBUG_VAR(ip.getValueInt("MERGE_FLAG"));

    if (ip.getValueInt("GENOMS_FLAG") == 1)
      tool = 2;

    if (ip.getValueInt("MERGE_FLAG") == 1)
      tool = 3;

    DEBUG_VAR(tool);

    reportSpsplotParams.setValue("TOOL", parseInt(tool));

    DEBUG_VAR(reportSpsplotParams.getValueInt("TOOL"));

    // add number of cpus
    if (ip.exists("CPUS"))
    {
      reportSpsplotParams.addIfExists(ip, "CPUS");
      //} else if(ip.exists("GRID_NUMCPUS")) {
      //  reportSpsplotParams.setValue("CPUS", ip.getValue("GRID_NUMCPUS") );
    }

    // Get report dir -- for the new reports
    string report_dir;
    string report_dir_aux = "report";
    if (ip.exists("REPORT_DIR"))
      report_dir_aux = ip.getValue("REPORT_DIR");

    // If report directory is a relative path, make it absolute
    if (report_dir_aux[0] != '/')
    {
      report_dir = currentWorkindDir;
      if (report_dir[report_dir.size() - 1] != '/')
        report_dir += '/';
    }
    report_dir += report_dir_aux;

    // set the report directory parameter
    reportSpsplotParams.setValue("OUTDIR", report_dir);

    aux = currentWorkindDir;
    reportSpsplotParams.setValue("FONT_PATH", aux);

    //aux = currentWorkindDir;
    //aux += "/spectra/stars_only.pklbin";
    aux = "spectra/stars_only.pklbin";
    reportSpsplotParams.setValue("OUTPUT_STARS", aux);
    //aux = currentWorkindDir;
    //aux += "/spectra/stars_indices.bin";
    aux = "spectra/stars_indices.bin";
    reportSpsplotParams.setValue("OUTPUT_STARS_INDEX", aux);
    //aux = currentWorkindDir ; aux += "/spectra/stars.pklbin";
    //reportSpsplotParams.setValue("OUTPUT_STARS_ALL",    aux);

    aux = ip.getValue("FASTA_DATABASE");
    string aux2;
    if (aux[0] != '/')
    {
      aux2 = currentWorkindDir;
      if (aux2[aux2.size() - 1] != '/')
        aux2 += '/';
    }
    aux2 += aux;
    reportSpsplotParams.setValue("FILE_FASTA", aux);

    if (execMode == string("sps"))
    {
      DEBUG_MSG("sps selected");

      //aux = currentWorkindDir;
      //aux += "/homology/ref_sps_names.txt";
      aux = "homology/ref_sps_names.txt";
      reportSpsplotParams.setValue("FILE_REFINDEX", aux);
      //aux = currentWorkindDir;
      //aux += "/assembly/component_info.bin";
      aux = "assembly/component_info.bin";
      reportSpsplotParams.setValue("FILE_COMP", aux);
      //aux = currentWorkindDir;
      //aux += "/assembly/sps_seqs.pklbin";
      aux = "assembly/sps_seqs.pklbin";
      reportSpsplotParams.setValue("FILE_SEQS", aux);

      // for the old reports
      //aux = currentWorkindDir;
      //aux += "/homology/contigs_mp_all.bin";
      aux = "homology/contigs_mp_all.bin";
      reportSpsplotParams.setValue("FILE_MP", aux);
      //aux = currentWorkindDir;
      //aux += "/homology/contigs_midx_all.pklbin";
      aux = "homology/contigs_midx_all.pklbin";
      reportSpsplotParams.setValue("FILE_MIDX", aux);

      // for the new reports
      //aux = currentWorkindDir;
      //aux += "/homology/contigs_mp.bin";
      aux = "homology/contigs_mp.bin";
      reportSpsplotParams.setValue("FILE_MP2", aux);
      //aux = currentWorkindDir;
      //aux += "/homology/contigs_midx.pklbin";
      aux = "homology/contigs_midx.pklbin";
      reportSpsplotParams.setValue("FILE_MIDX2", aux);

      //aux = currentWorkindDir;
      //aux += "/homology/homglue_ref_mp.bin";
      aux = "homology/homglue_ref_mp.bin";
      reportSpsplotParams.setValue("FILE_REFMP", aux);
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_ref_midx.pklbin";
      aux = "homology/homglue_ref_midx.pklbin";
      reportSpsplotParams.setValue("FILE_REFMIDX", aux);

      //aux = currentWorkindDir;
      //aux += "/spectra/contigs.pklbin";
      aux = "homology/contigs.pklbin";
      reportSpsplotParams.setValue("MATCHED_CONTIGS", aux);
      //aux = currentWorkindDir;
      //aux += "/spectra/contigs_indices.bin";
      aux = "homology/contigs_indices.bin";
      reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX", aux);

    }
    else
    {

      DEBUG_MSG("snets selected");

      //aux = currentWorkindDir;
      //aux += "/homology/ref_snets_names.txt";
      aux = "homology/ref_snets_names.txt";
      reportSpsplotParams.setValue("FILE_REFINDEX", aux);
      //aux = currentWorkindDir;
      //aux += "/homology/specnets_abruijn.bin";
      aux = "homology/specnets_abruijn.bin";
      reportSpsplotParams.setValue("FILE_COMP", aux);
      //    aux = currentWorkindDir ; aux += "/specnets/snets_specs.pklbin";
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches.pklbin";
      aux = "homology/homglue_matches.pklbin";
      reportSpsplotParams.setValue("FILE_SEQS", aux);

      // for the old reports

      //    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_mp.bin";
      aux = "homology/homglue_matches_mp.bin";
      reportSpsplotParams.setValue("FILE_MP", aux);
      //    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_midx.pklbin";
      aux = "homology/homglue_matches_midx.pklbin";
      reportSpsplotParams.setValue("FILE_MIDX", aux);

      // for the new reports

      //    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_mp.bin";
      aux = "homology/homglue_matches_mp.bin";
      reportSpsplotParams.setValue("FILE_MP2", aux);
      //    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_midx.pklbin";
      aux = "homology/homglue_matches_midx.pklbin";
      reportSpsplotParams.setValue("FILE_MIDX2", aux);

      // Same as FILE_MP / FILE_MIDX - no homology mapping for spectral networks
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_mp.bin";
      aux = "homology/homglue_matches_mp.bin";
      reportSpsplotParams.setValue("FILE_REFMP", aux);
      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches_midx.pklbin";
      aux = "homology/homglue_matches_midx.pklbin";
      reportSpsplotParams.setValue("FILE_REFMIDX", aux);

      //aux = currentWorkindDir;
      //aux += "/homology/homglue_matches.pklbin";
      aux = "homology/homglue_matches.pklbin";
      reportSpsplotParams.setValue("MATCHED_CONTIGS", aux);
      //    aux = currentWorkindDir ; aux += "/spectra/contigs_indices.bin";
      //    reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX",     aux);
    }

    //aux = currentWorkindDir;
    //aux += "/spectra/stars.pklbin";
    aux = "spectra/stars.pklbin";
    reportSpsplotParams.setValue("FILE_STARS", aux);
    //aux = currentWorkindDir;
    //aux += "/spectra/specs_ms.pklbin";
    if (ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST) == CLUSTER_PRMS)
      aux = SPECTRA_DIR + "/" + SPECS_SCORED_FILE;
      aux = SPECTRA_DIR + "/" + SPECS_MS_FILE;
    reportSpsplotParams.setValue("FILE_MS", aux);
    //aux = currentWorkindDir;
    //aux += "/spectra/input_index.txt";
    aux = "spectra/input_index.txt";
    reportSpsplotParams.setValue("FILE_INDEX", aux);

    //aux = currentWorkindDir;
    //aux += "/spectra/out/clust/clusters_0_";
    aux = "spectra/out/clust/clusters_0_";
    reportSpsplotParams.setValue("FILE_CLUSTER", aux);
    //aux = currentWorkindDir;
    //aux += '/';
    //aux += DEFAULT_INPUT_FILES_LIST;
    aux = DEFAULT_INPUT_FILES_LIST;
    reportSpsplotParams.setValue("FILE_CLUSTERMS", aux);
    //aux = currentWorkindDir;
    //aux += '/';
    //aux += DEFAULT_BIN_FILES_FILENAME;
    aux = DEFAULT_BIN_FILES_FILENAME;
    reportSpsplotParams.setValue("FILE_CLUSTERSCAN", aux);

    // write params to debug file
    reportSpsplotParams.writeToFile(getProjPath(ip, "debug_report.params"));

    // instatiate reports module
    DEBUG_TRACE;
    ExecReportSpsplot moduleReportSpsplot(reportSpsplotParams);
    DEBUG_TRACE;

    // test parameters
    string errorString;
    bool isValid = moduleReportSpsplot.validateParams(errorString);
    TEST_VALID;

    DEBUG_MSG("Invoking reports module");

    bool returnStatus = moduleReportSpsplot.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecReportSpsplot");

    string source = ip.getValue("EXE_DIR");
    source += "/resources/css";
    string dest = ip.getValue("REPORT_DIR");
    CopyDirectory(source, dest);

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performExecMainSpecnets(ParameterList & ip,
                               SpecSet * msSpectra,
                               SpecSet * scoredSpectra,
                               SpecSet * starSpectra,
                               SpectrumPairSet * pairs,
                               DB_fasta * db,
                               PenaltyMatrix * penaltyMatrixBlosum,
                               PenaltyMatrix * penaltyMatrixMods,
                               PeptideSpectrumMatchSet * psms,
                               PeptideSpectrumMatchSet * origPsms,
                               SpecSet * psms_spectra,
                               SpecSet * psms_midx,
                               vector<vector<int> > * psms_mp,
                               SpecSet * snets_contigs,
                               SpecSet * snets_midx,
                               vector<vector<int> > * snets_mp)
  {
    ParameterList specnetsParams;
    specnetsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");
    specnetsParams.setValue("RESOLUTION", ip.getValue("RESOLUTION", "0.1"));
    specnetsParams.setValue("MIN_MOD_MASS", ip.getValue("MIN_MOD_MASS", "-100"));
    specnetsParams.setValue("MAX_MOD_MASS", ip.getValue("MAX_MOD_MASS", "100"));
    specnetsParams.setValue("TOLERANCE_PEAK",
                            ip.getValue("TOLERANCE_PEAK", "0.45"));
    specnetsParams.setValue("TOLERANCE_PM", ip.getValue("TOLERANCE_PM", "1.5"));

    //alignment params
    specnetsParams.addIfExists(ip, "MAX_PARSIMONY");
    specnetsParams.addIfExists(ip, "MAX_NUM_MODS");
    specnetsParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
    specnetsParams.addIfExists(ip, "MATCHES_PER_MOD");
    specnetsParams.addIfExists(ip, "MIN_MATCHED_PEAKS_DB");
    specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT");
    specnetsParams.addIfExists(ip, "ENFORCE_ENDPEAKS");
    specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
    specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
    specnetsParams.addIfExists(ip, "MIN_PENALTY_FREQUENCY");
    specnetsParams.addIfExists(ip, "MIN_PEAK_EQUIVALENTS");
    specnetsParams.addIfExists(ip, "MAX_PEAK_EQUIVALENTS");
    specnetsParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
    specnetsParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");

    // Tagsearch params
    specnetsParams.addIfExists(ip, "INSPECT_PSMS");
    specnetsParams.addIfExists(ip, "INPUT_SPEC_IDS");
    specnetsParams.addIfExists(ip, "INPUT_CLUSTERS_DIR");
    specnetsParams.setValue("OUTPUT_PEPTIDES", "spectra/tagsearch_peptides.txt");
    specnetsParams.addIfExists(ip, "OUTPUT_RAW_RESULTS");
    specnetsParams.addIfExists(ip, "OUTPUT_SNETS_PTM_MATRIX");
    specnetsParams.addIfExists(ip, "OUTPUT_SPECS_PROJ");
    specnetsParams.addIfExists(ip, "OUTPUT_ALIGNS");
    specnetsParams.addIfExists(ip, "OUTPUT_ANNOTINFO");

    // SVM params
    specnetsParams.setValue("STARS_SCORED_PRM_OFFSET", "-1.0072763");
    specnetsParams.setValue("STARS_SCORED_SRM_OFFSET", "-1.0072763");
    specnetsParams.setValue("SCAN_ZERO_INDEX", "1");
    specnetsParams.addIfExists(ip, "SCAN_ZERO_INDEX");
    specnetsParams.addIfExists(ip, "MS2_SCORING_MODEL"); // To become "<trunk>/resources/dancik_model.txt
    specnetsParams.addIfExists(ip, "SPECS_MS_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_ms.txt
    specnetsParams.addIfExists(ip, "SPECS_SCORED_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_prm.txt
    specnetsParams.addIfExists(ip, "STARS_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_stars.txt
    specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE1"); // To become "<trunk>/resources/HEK293_charge1_model_svm_keys_range.txt
    specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE2"); // To become "<trunk>/resources/HEK293_charge2_model_svm_keys_range.txt
    specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE3"); // To become "<trunk>/resources/HEK293_charge3_model_svm_keys_range.txt
    specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE1"); // To become "<trunk>/resources/HEK293_charge1_model_SVM.model
    specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE2"); // To become "<trunk>/resources/HEK293_charge2_model_SVM.model
    specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE3"); // To become "<trunk>/resources/HEK293_charge3_model_SVM.model

    // Spectral networks params
    specnetsParams.setValue("SPECNETS_PROJ_TYPE",
                            ip.getValue("SPECNETS_PROJ_TYPE", "all")); // Possible values are "all" / "matched" to retain all/matched-only PRMs during propagation
    specnetsParams.setValue("SPECNETS_USE_SVM",
                            ip.getValue("SPECNETS_USE_SVM", "0"));

    // FDR params
    specnetsParams.setValue("INPUT_RESULTS_TYPE", "snets");
    specnetsParams.setValue("OUTPUT_FDR_RESULTS", "spectra/psms_fdr.txt");

    DEBUG_TRACE;
    SpectrumPairSet filteredPairs;

    specnetsParams.writeToFile("debug_specnets.params");

    MS2ScoringModel model;
    if (ip.exists("MS2_SCORING_MODEL"))
    {
      if (!model.LoadModel(ip.getValue("MS2_SCORING_MODEL").c_str()))
      {
        DEBUG_MSG("Could not load " << ip.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }
    DEBUG_MSG("SCORING MODEL LOADED");

    DEBUG_TRACE;
    ExecMainSpecnets mainSpecnets(specnetsParams,
                                  msSpectra,
                                  scoredSpectra,
                                  starSpectra,
                                  pairs,
                                  &model,
                                  db,
                                  penaltyMatrixBlosum,
                                  penaltyMatrixMods,
                                  psms,
                                  origPsms,
                                  psms_spectra,
                                  psms_midx,
                                  psms_mp,
                                  snets_contigs,
                                  snets_midx,
                                  snets_mp);
    DEBUG_TRACE;

    string errorString;
    bool isValid = mainSpecnets.validateParams(errorString);
    TEST_VALID;

    bool returnStatus = mainSpecnets.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecMainSpecnets");

    // Saving output data
    returnStatus = mainSpecnets.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecMainSpecnets");

    ofstream spsIndex(getProjPath(ip, "homology/ref_snets_names.txt").c_str(),
                      ios_base::out | ios_base::binary);
    for (unsigned int i = 0; i < snets_contigs->size(); i++)
      spsIndex << "snet:" << i + 1 << endl;
    spsIndex.close();

    //---------------------------------------------------------------------------
    // REPORT STAGE
    //---------------------------------------------------------------------------

    ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                     ios_base::out | ios_base::binary);
    spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
        << ip.getValue("TOLERANCE_PM") << "\n";
    spsProj.close();

    if (!performReport(ip))
    {
      ERROR_MSG("Problem encountered during Report stage");
      //exit(-100 - STAGE_REPORT);
      exit(-100);  // Lars changed so copy of define was not needed
    }

    return true;
  }

} // end name specnets

