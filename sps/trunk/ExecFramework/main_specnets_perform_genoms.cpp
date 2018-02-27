//
#include "main_specnets_defs.h"
#include "main_specnets_helpers.h"
#include "main_specnets_perform_genoms.h"

// Module Includes
#include "abruijn.h"
#include "ClusterData.h"
#include "ExecGenoMS.h"
#include "Logger.h"
#include "ParameterList.h"
#include "StatusFile.h"

// Specnets Includes
#include "utils.h"  // stringSplit
#include "SpecSet.h"

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
  bool generateRelaunchScript(ParameterList & ip)
  {
    string relauncherScriptName = "relauncher.sh";
    FILE* myfile;
    myfile = fopen(relauncherScriptName.c_str(), "rb");
    if (myfile)
    {
      fclose(myfile);
      return true;
    }

    string dbFileName = getProjPath(ip, "relaunch.protid.fasta");
    string paramsFile = "relaunch.params";
    string exeDir = ip.getValue("EXE_DIR", ".");

    string logLevel = ip.getValue("LOG_LEVEL", "0");
    if (strlen(logLevel.c_str()) == 0)
      logLevel = "0";
    //Create the command

    //Change to the project directory
    string relauncherScript = "cd " + getProjPath(ip, ".") + "\n";

    relauncherScript += "echo \"Running\" > status.txt\n";

    //If we were currently running in merge mode, then we need to move the old CSPS files back to their old places
    // if(ip.getValueInt("MERGE_FLAG",0) == 1)
    myfile = fopen("assembly/csps.sps_seqs.pklbin", "rb");
    if (myfile)
    {
      fclose(myfile);
      relauncherScript +=
          "mv assembly/csps.sps_seqs.pklbin assembly/sps_seqs.pklbin\n";
      relauncherScript +=
          "mv assembly/csps.component_info.bin assembly/component_info.bin\n";
      //relauncherScript +=
      //    "mv homology/csps.contigs_midx_all.pklbin homology/contigs_midx_all.pklbin\n";
      relauncherScript +=
          "mv homology/csps.contigs_midx.pklbin homology/contigs_midx.pklbin\n";
      relauncherScript +=
          "mv homology/csps.contigs_mp_all.bin homology/contigs_mp_all.bin\n";
      relauncherScript +=
          "mv homology/csps.contigs_mp.bin homology/contigs_mp.bin\n";
      relauncherScript +=
          "mv homology/csps.homglue_matches_midx.pklbin homology/homglue_matches_midx.pklbin\n";
      relauncherScript +=
          "mv homology/csps.homglue_matches_mp.bin homology/homglue_matches_mp.bin\n";
      relauncherScript +=
          "mv homology/csps.homglue_matches.pklbin homology/homglue_matches.pklbin\n";
      relauncherScript +=
          "mv homology/csps.homglue_ref_midx.pklbin homology/homglue_ref_midx.pklbin\n";
      relauncherScript +=
          "mv homology/csps.homglue_ref_mp.bin homology/homglue_ref_mp.bin\n";

      ifstream ifile("homology/csps.ref_sps_names.txt");
      if (ifile)
      {
        relauncherScript +=
            "mv homology/csps.ref_sps_names.txt homology/ref_sps_names.txt\n";
      }

      relauncherScript +=
          "mv spectra/csps.contigs_indices.bin spectra/contigs_indices.bin\n";
      relauncherScript +=
          "mv spectra/csps.contigs.pklbin spectra/contigs.pklbin\n";
      relauncherScript += "mv spectra/csps.stars.pklbin spectra/stars.pklbin\n";
    }

    //execute the command
    relauncherScript += exeDir + "/main_specnets " + paramsFile;
    relauncherScript += " -i tagsearch";
    relauncherScript += " -ll " + logLevel;
    relauncherScript += " -lf "
        + ip.getValue("LOG_FILE_NAME", "relauncher_log.txt");
    int val = ip.getValueInt("GRID_EXECUTION");
    if (val != 0)
      relauncherScript += " -g";
    //val = ip.getValueInt("GENOMS_FLAG");
    //if(val != 0)
    //  relauncherScript += " -q";
    //val = ip.getValueInt("MERGE_FLAG");
    //if(val != 0)
    //  relauncherScript += " -m";
    relauncherScript += "\n";
    //Write the relauncher.sh
    FILE * f = fopen(relauncherScriptName.c_str(), "wb");

    if (f == NULL)
    {
      ERROR_MSG("Unable to open relauncher file for writing");
      return false;
    }

    val = fwrite(relauncherScript.c_str(),
                 sizeof(char),
                 strlen(relauncherScript.c_str()),
                 f);
    if (val != strlen(relauncherScript.c_str()))
    {
      ERROR_MSG("Problem encountered writing relauncher file");
      ERROR_MSG(val);
      ERROR_MSG(strlen(relauncherScript.c_str()));
      return false;
    }
    fclose(f);

    //Write the updated params file
    ip.setValue("FASTA_DATABASE", dbFileName);
    ip.setValue("GENOMS_FLAG", "0");
    ip.setValue("MERGE_FLAG", "0");
    ip.writeToFile(getProjPath(ip, paramsFile));

    DEBUG_MSG("Wrote relauncher.sh");
    return true;
  }

  //-----------------------------------------------------------------------------
  bool performGenoMS(ParameterList & ip)
  {
    DEBUG_TRACE;

    ParameterList genoMSParams;

    genoMSParams.addIfExists(ip, "EXE_DIR");
    genoMSParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

    genoMSParams.addIfExists(ip, "TOLERANCE_PEAK");
    genoMSParams.addIfExists(ip, "TOLERANCE_PM");
    genoMSParams.addIfExists(ip, "PEAK_PENALTY");

    genoMSParams.addIfExists(ip, "ALL2ALL_SIMILARITY");
    genoMSParams.addIfExists(ip, "HMM_LATE_ADD");
    genoMSParams.addIfExists(ip, "FDR_CUTOFF");
    genoMSParams.addIfExists(ip, "MUTATION_MODE");
    genoMSParams.addIfExists(ip, "DBCOMBINED");
    genoMSParams.addIfExists(ip, "DBROOTNAME");
    genoMSParams.addIfExists(ip, "GENOMESEQ");
    genoMSParams.addIfExists(ip, "DIGEST");
    genoMSParams.addIfExists(ip, "TEMPLATECONSTRAINTFILE");
    genoMSParams.addIfExists(ip, "FIXEDMOD");
    genoMSParams.addIfExists(ip, "PROJECT_DIR");
    genoMSParams.addIfExists(ip, "RUN_DBSEARCH");

    //new params for CID, ETD, HCD, or pairs/triplets
    genoMSParams.addIfExists(ip, "NUM_CONSECUTIVE");
    genoMSParams.addIfExists(ip, "INSTRUMENT_TYPE");
    genoMSParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
    genoMSParams.addIfExists(ip, "CLUSTER_TOOL");
    genoMSParams.addIfExists(ip, "MERGE_SAME_PREC");

    //genoMSParams.setValue("PRMS",ip.getValue("PEPNOVO_OUTPUT_PRMS","spectra/specs_scored.prms"));
    //genoMSParams.setValue("PRMS", "spectra/pepnovo_out_CID.prms");
    string scoredMgfFile = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_SCORED_MGF_FILE;
    genoMSParams.setValue("PRMS", scoredMgfFile);
    //genoMSParams.setValue("SPECTRA",     "spectra/specs_ms_pepnovo.mgf");
    //genoMSParams.setValue("SPECTRA", "spectra/pepnovo_in_CID.mgf");
    string msMgfFile = getProjPath(ip, SPECTRA_DIR) + "/" + SPECS_MS_MGF_FILE;
    genoMSParams.setValue("SPECTRA", msMgfFile);
    genoMSParams.setValue("OUTPUT_FILE", "genoMS.out");
    genoMSParams.setValue("GENERATE_REPORTS", "1");
    genoMSParams.setValue("LOG_FILE", "genoMS.log");

    DEBUG_TRACE;
    stringstream aux;
    genoMSParams.print(aux);
    DEBUG_MSG(aux.str());
    DEBUG_TRACE;
    ExecGenoMS moduleGenoMS(genoMSParams);

    string errorString;
    bool isValid = moduleGenoMS.validateParams(errorString);
    TEST_VALID;

    DEBUG_TRACE;

    bool returnStatus = moduleGenoMS.invoke();
    // Test for return status
    TEST_RETURN_STATUS("moduleGenoMS");

    DEBUG_TRACE;

    return true;
  }

  //-----------------------------------------------------------------------------
  bool performMergeOfCSPSAndGenoMS(ParameterList & ip)
  {
    ///////////////////////////////////////////////
    //First we duplicate the input spectrum file
    ///////////////////////////////////////////////
    bool debug = false;
    DEBUG_TRACE;

    unsigned int cspsContigCount;
    unsigned int genomsContigCount;

    spsReports::ClusterData clusterData;

    SpecSet inputSpectra;
    SpecSet other;
    unsigned int specCount; // total (unique spectra)

    if (inputSpectra.loadPklBin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
      return false;
    }

    DEBUG_MSG("First we have " << inputSpectra.size() << " inputSpectra!");
    if (other.loadPklBin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);

    specCount = inputSpectra.size();

    DEBUG_MSG("Now we have " << specCount << " inputSpectra!");
    if (inputSpectra.savePklBin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving doubled spectra to ./spectra/specs_ms.pklbin");
      return false;
    }
    if (debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/specs_ms.mgf");
    }

    ///////////////////////////////////////////////
    //Also duplicate the clusterData.bin
    ///////////////////////////////////////////////
    if ((ip.exists("CLUSTER_TOOL")
        && ip.getValue("CLUSTER_TOOL", "") == CLUSTER_PRMS)
        || (ip.exists("CLUSTER_MIN_SIZE")
            && ip.getValueInt("CLUSTER_MIN_SIZE") <= 0))
    {

      DEBUG_MSG("Doubling input_mapping.bin");
      vector<vector<int> > map;
      string fName = getProjPath(ip, DEFAULT_INPUT_MAPPING);
      Load_binArray(fName.c_str(), map);
      DEBUG_MSG("Loaded double");
      Save_binArrayDouble(fName.c_str(), map);
      DEBUG_MSG("Saved double");

    }
    else
    {
      DEBUG_MSG("Doubling clusterData.bin");
      clusterData.loadDataDouble(".");
      DEBUG_MSG("Loaded double");
      clusterData.saveData("..");
      DEBUG_MSG("Saved double");
    }
    ///////////////////////////////////////////////
    //Merge the star spectra
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging star spectra, csps + genoms");
    if (inputSpectra.loadPklBin("./spectra/csps.stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/csps.stars.pklbin");
      return false;
    }

    if (other.loadPklBin("./spectra/genoms.stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/genoms.stars.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./spectra/stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving doubled star spectra to ./spectra/stars.pklbin");
      return false;
    }

    DEBUG_MSG("Made " << inputSpectra.size() << " star spectra!");
    if (debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/stars.mgf");
    }
    ///////////////////////////////////////////////
    //Merge the contig spectra
    //////////////////////////////////////////////
    DEBUG_MSG("Merging contig spectra");
    if (inputSpectra.loadPklBin("./assembly/csps.sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./assembly/csps.sps_seqs.pklbin");
      return false;
    }

    cspsContigCount = inputSpectra.size();

    if (other.loadPklBin("./assembly/genoms.sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./assembly/genoms.sps_seqs.pklbin");
      return false;
    }
    genomsContigCount = other.size();
    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./assembly/sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving merged contig spectra to ./assembly/sps_seqs.pklbin");
      return false;
    }
    DEBUG_MSG("Made " << inputSpectra.size() << " contig spectra!");
    DEBUG_MSG("CSPS has contigs 0-" << (cspsContigCount - 1));
    DEBUG_MSG("GenoMS has contigs " << cspsContigCount << "-" << (genomsContigCount + cspsContigCount - 1));

    if (debug)
    {
      inputSpectra.SaveSpecSet_mgf("./assembly/sps_seqs.mgf");
    }
    ///////////////////////////////////////////////
    //Merge abruin infos
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging abruijn info");
    abinfo_t info1;
    abinfo_t info2;
    string abinfoOutFile = "./assembly/component_info.bin";

    Load_abinfo("./assembly/csps.component_info.bin", info1);
    Load_abinfo("./assembly/genoms.component_info.bin", info2);

    if (!merge_abinfo(info1, info2, abinfoOutFile.c_str(), specCount / 2))
    {
      ERROR_MSG("Error merging abruijn info");
      return false;
    }
    if (debug)
    {
      abinfo_t info;
      Load_abinfo("./assembly/component_info.bin", info);
      dumpAbInfo("./assembly/component_info.txt", info);
    }
    ///////////////////////////////////////////////
    //merge contig_mp.bin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging contig_mp_all.bin");
    vector<vector<int> > *m_contigs_mp_sub1;
    vector<vector<int> > *m_contigs_mp_sub2;

    m_contigs_mp_sub1 = new vector<vector<int> >();
    m_contigs_mp_sub2 = new vector<vector<int> >();

    DEBUG_MSG("Initialized vectors..");
    // load the data
    if (Load_binArray<int, vector>("./homology/csps.contigs_mp.bin",
                                   *m_contigs_mp_sub1) < 0)
    {

      return false;
    }

    if (Load_binArray<int, vector>("./homology/genoms.contigs_mp.bin",
                                   *m_contigs_mp_sub2) < 0)
    {

      return false;
    }
    DEBUG_MSG("Loaded contig_mp_all vectors");
    if (Save_doubleBinArray<int>("./homology/contigs_mp.bin",
                                 *m_contigs_mp_sub1,
                                 *m_contigs_mp_sub2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/contigs_mp.bin");
      return false;

    }

    DEBUG_MSG("CSPS contig_mp.bin = " << (*m_contigs_mp_sub1).size() << "x" << (*m_contigs_mp_sub1)[0].size());
    DEBUG_MSG("GenoMS contig_mp.bin = " << (*m_contigs_mp_sub2).size() << "x" << (*m_contigs_mp_sub2)[0].size());
    if (debug)
    {

      cout << "DEBUG: homology/contigs_mp.bin" << endl;

      Load_binArray<int, vector>("./homology/contigs_mp.bin", *m_contigs_mp_sub1);
      for (int i = 0; i < (*m_contigs_mp_sub1).size(); ++i)
      {
        for (int j = 0; j < (*m_contigs_mp_sub1)[0].size(); ++j)
        {
          cout << " " << (*m_contigs_mp_sub1)[i][j];
        }
        cout << endl;
      }

    }
    ///////////////////////////////////////////////
    //Merge contigs_midx.pklbin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging contigs_midx.pklbin");
    if (inputSpectra.loadPklBin("./homology/csps.contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx.pklbin");
      return false;
    }

    if (other.loadPklBin("./homology/genoms.contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./homology/contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx.pklbin");
      return false;
    }
    if (debug)
    {
      inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx.mgf");
    }

    ///////////////////////////////////////////////
    //merge contig_mp_all.bin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging contig_mp_all.bin");
    vector<vector<int> > *m_contigs_mp1;
    vector<vector<int> > *m_contigs_mp2;

    m_contigs_mp1 = new vector<vector<int> >();
    m_contigs_mp2 = new vector<vector<int> >();

    DEBUG_MSG("Initialized vectors..");
    // load the data
    if (Load_binArray<int, vector>("./homology/csps.contigs_mp_all.bin",
                                   *m_contigs_mp1) < 0)
    {

      return false;
    }

    if (Load_binArray<int, vector>("./homology/genoms.contigs_mp_all.bin",
                                   *m_contigs_mp2) < 0)
    {

      return false;
    }
    DEBUG_MSG("Loaded contig_mp_all vectors");
    if (Save_doubleBinArray<int>("./homology/contigs_mp_all.bin",
                                 *m_contigs_mp1,
                                 *m_contigs_mp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/contigs_mp_all.bin");
      return false;

    }

    DEBUG_MSG("CSPS contig_mp_all.bin = " << (*m_contigs_mp1).size() << "x" << (*m_contigs_mp1)[0].size());
    DEBUG_MSG("GenoMS contig_mp_all.bin = " << (*m_contigs_mp2).size() << "x" << (*m_contigs_mp2)[0].size());
    if (debug)
    {

      cout << "DEBUG: homology/contigs_mp_all.bin" << endl;

      Load_binArray<int, vector>("./homology/contigs_mp_all.bin", *m_contigs_mp1);
      for (int i = 0; i < (*m_contigs_mp1).size(); ++i)
      {
        for (int j = 0; j < (*m_contigs_mp1)[0].size(); ++j)
        {
          cout << " " << (*m_contigs_mp1)[i][j];
        }
        cout << endl;
      }

    }
    ///////////////////////////////////////////////
    //Merge contigs_midx_all.pklbin
    ///////////////////////////////////////////////
    /*DEBUG_MSG("Merging contigs_midx_all.pklbin");
     if (inputSpectra.loadPklBin("./homology/csps.contigs_midx_all.pklbin") <= 0)
     {
     ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx_all.pklbin");
     return false;
     }

     if (other.loadPklBin("./homology/genoms.contigs_midx_all.pklbin") <= 0)
     {
     ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx_all.pklbin");
     return false;
     }

     inputSpectra.appendSpecSet(other);
     if (inputSpectra.savePklBin("./homology/contigs_midx_all.pklbin") <= 0)
     {
     ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx_all.pklbin");
     return false;
     }
     if (debug)
     {
     inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx_all.mgf");
     }
     */
    ///////////////////////////////////////////////
    //Merge contigs.pklbin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging contigs.pklbin");
    if (inputSpectra.loadPklBin("./spectra/csps.contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/csps.contigs.pklbin");
      return false;
    }
    unsigned int cspsMappedContigCount = inputSpectra.size();
    if (other.loadPklBin("./spectra/genoms.contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/genoms.contigs.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./spectra/contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig spectra to ./spectra/contigs.pklbin");
      return false;
    }
    DEBUG_MSG("Made " << inputSpectra.size() << " mapped contig spectra!");
    DEBUG_MSG("CSPS has contigs 0-" << (cspsMappedContigCount - 1));
    DEBUG_MSG("GenoMS has contigs " << cspsMappedContigCount << "-" << (inputSpectra.size() - 1));
    if (debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/contigs.mgf");
    }
    ///////////////////////////////////////////////
    //merge contig_indices.bin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging contig_indices.bin");

    vector<vector<int> > *m_contigIndex1;
    vector<vector<int> > *m_contigIndex2;

    m_contigIndex1 = new vector<vector<int> >();
    m_contigIndex2 = new vector<vector<int> >();

    // load the data
    if (Load_binArray<int, vector>("./spectra/csps.contigs_indices.bin",
                                   *m_contigIndex1) < 0)
    {

      ERROR_MSG("Error loading binary array ./spectra/csps.contigs_indices.bin");
      return false;
    }

    if (Load_binArray<int, vector>("./spectra/genoms.contigs_indices.bin",
                                   *m_contigIndex2) < 0)
    {

      ERROR_MSG("Error loading binary array ./spectra/genoms.contigs_indices.bin");
      return false;
    }
    DEBUG_MSG("CSPS contig_indices.bin = " << (*m_contigIndex1).size() << "x" << (*m_contigIndex1)[0].size());
    DEBUG_MSG("GenoMS contig_indices.bin = " << (*m_contigIndex2).size() << "x" << (*m_contigIndex2)[0].size());

    //Update component indices from genoms result
    for (int i = 0; i < m_contigIndex2->size(); ++i)
    {
      (*m_contigIndex2)[i][0] += cspsContigCount;

    }

    if (Save_doubleBinArray<int>("./spectra/contigs_indices.bin",
                                 *m_contigIndex1,
                                 *m_contigIndex2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./spectra/contigs_indices.bin");
      return false;

    }

    if (debug)
    {

      cout << "DEBUG: spectra/contigs_indices.bin" << endl;

      Load_binArray<int, vector>("./spectra/contigs_indices.bin",
                                 *m_contigIndex1);
      for (int i = 0; i < (*m_contigIndex1).size(); ++i)
      {
        for (int j = 0; j < (*m_contigIndex1)[0].size(); ++j)
        {
          cout << " " << (*m_contigIndex1)[i][j];
        }
        cout << endl;
      }

    }

    ///////////////////////////////////////////////
    //merge homglue_ref_mp.bin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging homglue_ref_mp.bin");
    vector<vector<int> > *m_homglueRefMp1;
    vector<vector<int> > *m_homglueRefMp2;

    m_homglueRefMp1 = new vector<vector<int> >();
    m_homglueRefMp2 = new vector<vector<int> >();

    // load the data
    if (Load_binArray<int, vector>("./homology/csps.homglue_ref_mp.bin",
                                   *m_homglueRefMp1) < 0)
    {

      return false;
    }

    if (Load_binArray<int, vector>("./homology/genoms.homglue_ref_mp.bin",
                                   *m_homglueRefMp2) < 0)
    {

      return false;
    }

    DEBUG_MSG("CSPS homglue_ref_mp.bin = " << (*m_homglueRefMp1).size() << "x" << (*m_homglueRefMp1)[0].size());
    DEBUG_MSG("GenoMS homglue_ref_mp.bin = " << (*m_homglueRefMp2).size() << "x" << (*m_homglueRefMp2)[0].size());

    if (Save_doubleBinArray("./homology/homglue_ref_mp.bin",
                            *m_homglueRefMp1,
                            *m_homglueRefMp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/homglue_ref_mp.bin");
      return false;

    }

    if (debug)
    {

      cout << "DEBUG: homology/homglue_ref_mp.bin" << endl;

      Load_binArray<int, vector>("./homology/homglue_ref_mp.bin",
                                 *m_homglueRefMp1);
      for (int i = 0; i < (*m_homglueRefMp1).size(); ++i)
      {
        for (int j = 0; j < (*m_homglueRefMp1)[0].size(); ++j)
        {
          cout << " " << (*m_homglueRefMp1)[i][j];
        }
        cout << endl;
      }

    }

    ///////////////////////////////////////////////
    //Merge homglue_ref_midx.pklbin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging homglue_ref_midx.pklbin");
    if (inputSpectra.loadPklBin("./homology/csps.homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_ref_midx.pklbin");
      return false;
    }

    if (other.loadPklBin("./homology/genoms.homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_ref_midx.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./homology/homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_ref_midx.pklbin");
      return false;
    }

    if (debug)
    {
      DEBUG_MSG("homglue_ref_midx.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_ref_midx.mgf");
    }

    ///////////////////////////////////////////////
    //Merge homglue_matches.pklbin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging homglue_matches.pklbin");
    if (inputSpectra.loadPklBin("./homology/csps.homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches.pklbin");
      return false;
    }

    if (other.loadPklBin("./homology/genoms.homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./homology/homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches.pklbin");
      return false;
    }
    if (debug)
    {
      DEBUG_MSG("homglue_matches.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches.mgf");
    }

    ///////////////////////////////////////////////
    //merge homglue_matches_mp.bin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging homglue_matches_mp.bin");
    vector<vector<int> > *m_homglueMatchMp1;
    vector<vector<int> > *m_homglueMatchMp2;

    m_homglueMatchMp1 = new vector<vector<int> >();
    m_homglueMatchMp2 = new vector<vector<int> >();

    // load the data
    if (Load_binArray<int, vector>("./homology/csps.homglue_matches_mp.bin",
                                   *m_homglueMatchMp1) < 0)
    {

      return false;
    }

    if (Load_binArray<int, vector>("./homology/genoms.homglue_matches_mp.bin",
                                   *m_homglueMatchMp2) < 0)
    {

      return false;
    }

    DEBUG_MSG("CSPS homglue_matches_mp.bin = " << (*m_homglueMatchMp1).size() << "x" << (*m_homglueMatchMp1)[0].size());
    DEBUG_MSG("GenoMS homglue_matches_mp.bin = " << (*m_homglueMatchMp2).size() << "x" << (*m_homglueMatchMp2)[0].size());

    if (Save_doubleBinArray("./homology/homglue_matches_mp.bin",
                            *m_homglueMatchMp1,
                            *m_homglueMatchMp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/homglue_matches_mp.bin");
      return false;

    }

    ///////////////////////////////////////////////
    //Merge homglue_matches_midx.pklbin
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging homglue_matches.pklbin");
    if (inputSpectra.loadPklBin("./homology/csps.homglue_matches_midx.pklbin")
        <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches_midx.pklbin");
      return false;
    }

    if (other.loadPklBin("./homology/genoms.homglue_matches_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches_midx.pklbin");
      return false;
    }

    inputSpectra.appendSpecSet(other);
    if (inputSpectra.savePklBin("./homology/homglue_matches_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches_midx.pklbin");
      return false;
    }

    if (debug)
    {
      DEBUG_MSG("homglue_matches_midx.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches_midx.mgf");
    }

    ///////////////////////////////////////////////
    //merge ref_sps_names.txt
    ///////////////////////////////////////////////
    DEBUG_MSG("Merging ref_sps_names.txt");
    vector<string> genoMScontigNames;
    int i, j;
    int contigNum = cspsContigCount + 1;
    //The GenoMS contig numbers need to be updated
    readFilesFromFile("./homology/genoms.ref_sps_names.txt", genoMScontigNames);

    FILE * fp = fopen("./homology/genoms.ref_sps_names.txt", "wb");
    //We assume the names are GenoMS:i, where i is the contig number
    for (i = 0; i < genoMScontigNames.size(); ++i)
    {
      DEBUG_MSG("Orig genoMS contig Name: " << genoMScontigNames[i]);
      string newName = "";
      for (j = 0; j < genoMScontigNames[i].size(); ++j)
      {
        if (genoMScontigNames[i][j] == ':')
          break;
        else
          newName += genoMScontigNames[i][j];
      }
      newName += ":";
      newName += intToString(contigNum);
      newName += "\n";
      genoMScontigNames[i] = newName;
      DEBUG_MSG("New genoMS contig Name: " << genoMScontigNames[i]);

      fwrite((void *)(genoMScontigNames[i].c_str()),
             genoMScontigNames[i].size(),
             1,
             fp);
      contigNum += 1;

    }

    fclose(fp);

    //Only concatenate if the csps version exists.
    ifstream tmpFile("./homology/csps.ref_sps_names.txt");
    if (tmpFile)
    {
      if (!concatenateFiles("./homology/csps.ref_sps_names.txt",
                            "./homology/genoms.ref_sps_names.txt",
                            "./homology/ref_sps_names.txt"))
        return false;
    }

    return true;

  }

  void performGenoMSRename(string & statusFileName)
  {

    //We need to rename all of our results files
    if (rename("./spectra/stars.pklbin", "./spectra/genoms.stars.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./assembly/sps_seqs.pklbin",
               "./assembly/genoms.sps_seqs.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./assembly/component_info.bin",
               "./assembly/genoms.component_info.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_mp.bin", "./homology/genoms.contigs_mp.bin")
        < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_midx.pklbin",
               "./homology/genoms.contigs_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_mp_all.bin",
               "./homology/genoms.contigs_mp_all.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_midx_all.pklbin",
               "./homology/genoms.contigs_midx_all.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./spectra/contigs.pklbin", "./spectra/genoms.contigs.pklbin")
        < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./spectra/contigs_indices.bin",
               "./spectra/genoms.contigs_indices.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_ref_mp.bin",
               "./homology/genoms.homglue_ref_mp.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_ref_midx.pklbin",
               "./homology/genoms.homglue_ref_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches.pklbin",
               "./homology/genoms.homglue_matches.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches_mp.bin",
               "./homology/genoms.homglue_matches_mp.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches_midx.pklbin",
               "./homology/genoms.homglue_matches_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/ref_sps_names.txt",
               "./homology/genoms.ref_sps_names.txt") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }

    if (rename("./protid.fasta", "./genoms.protid.fasta") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
  }

} // end name specnets

