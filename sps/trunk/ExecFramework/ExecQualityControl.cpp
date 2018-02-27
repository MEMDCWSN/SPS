/*
 * ExecQualityControl.cpp
 *
 *  Created on: Jul 8, 2013
 *      Author: aguthals
 */

#include "ExecQualityControl.h"
#include "FileUtils.h"

void ExecQualityControl::loadMS2ScoringModels(const string &modelDir,
                                              MS2ScoringModelSet &modelSet)
{
  string pathCID = getPath(modelDir, "model_cid.txt", false);
  string pathETD = getPath(modelDir, "model_etd.txt", false);
  string pathPRM = getPath(modelDir, "model_prm.txt", false);

  MS2ScoringModel scoringModel;
  Spectrum::FragType typeCID = Spectrum::FragType_CID;
  Spectrum::FragType typeHCD = Spectrum::FragType_HCD;
  Spectrum::FragType typeETD = Spectrum::FragType_ETD;
  Spectrum::FragType typePRM = Spectrum::FragType_PRM;
  if (scoringModel.LoadModel(pathCID.c_str()))
  {
    modelSet.addModel(Spectrum::activationToString(typeCID), scoringModel);
    modelSet.addModel(Spectrum::activationToString(typeHCD), scoringModel);
  }
  else
  {
    WARN_MSG("Failed to locate CID scoring model from " << pathCID);
  }

  if (scoringModel.LoadModel(pathETD.c_str()))
  {
    modelSet.addModel(Spectrum::activationToString(typeETD), scoringModel);
  }
  else
  {
    WARN_MSG("Failed to locate ETD scoring model from " << pathETD);
  }

  if (scoringModel.LoadModel(pathPRM.c_str()))
  {
    modelSet.addModel(Spectrum::activationToString(typePRM), scoringModel);
  }
  else
  {
    WARN_MSG("Failed to locate PRM scoring model from " << pathPRM);
  }
}

ExecQualityControl::ExecQualityControl(void) :
  ExecReportSPSStats(), m_paramsSPS(0x0), m_spectraRawMS2(0x0),
      m_spectraClustMS2(0x0), m_spectraScored(0x0), m_ms2ScoringModels(0x0),
      m_statsSpectraRawMS2(0x0), m_statsSpectraClustMS2(0x0),
      m_statsSpectraScored(0x0), m_statsSpectraStar(0x0), m_statsContigs(0x0),
      m_statsSummaryRawMS2(0x0), m_statsSummaryClustMS2(0x0),
      m_statsSummaryScored(0x0), m_statsSummaryStar(0x0),
      m_statsSummaryContigs(0x0)
{
  m_name = "ExecQualityControl";
  m_type = "ExecQualityControl";
}

ExecQualityControl::ExecQualityControl(const ParameterList & inputParams) :
  ExecReportSPSStats(inputParams), m_paramsSPS(0x0), m_spectraRawMS2(0x0),
      m_spectraClustMS2(0x0), m_spectraScored(0x0), m_ms2ScoringModels(0x0),
      m_statsSpectraRawMS2(0x0), m_statsSpectraClustMS2(0x0),
      m_statsSpectraScored(0x0), m_statsSpectraStar(0x0), m_statsContigs(0x0),
      m_statsSummaryRawMS2(0x0), m_statsSummaryClustMS2(0x0),
      m_statsSummaryScored(0x0), m_statsSummaryStar(0x0),
      m_statsSummaryContigs(0x0)
{
  m_name = "ExecQualityControl";
  m_type = "ExecQualityControl";
}

ExecQualityControl::~ExecQualityControl(void)
{
  if (ownInput)
  {
    delete m_paramsSPS;
    delete m_spectraRawMS2;
    delete m_spectraClustMS2;
    delete m_spectraScored;
    delete m_ms2ScoringModels;
  }
  if (ownOutput)
  {
    delete m_statsSpectraRawMS2;
    delete m_statsSpectraClustMS2;
    delete m_statsSpectraScored;
    delete m_statsSpectraStar;
    delete m_statsContigs;

    delete m_statsSummaryRawMS2;
    delete m_statsSummaryClustMS2;
    delete m_statsSummaryScored;
    delete m_statsSummaryStar;
    delete m_statsSummaryContigs;
  }
}

ExecBase * ExecQualityControl::clone(const ParameterList & inputParams) const
{
  return new ExecQualityControl(inputParams);
}

bool ExecQualityControl::invoke(void)
{
  if (ownOutput)
  {
    if (m_statsSpectraRawMS2 == 0x0)
    {
      m_statsSpectraRawMS2 = new OutputTable;
    }
    if (m_statsSpectraClustMS2 == 0x0)
    {
      m_statsSpectraClustMS2 = new OutputTable;
    }
    if (m_statsSpectraScored == 0x0)
    {
      m_statsSpectraScored = new OutputTable;
    }
    if (m_statsSpectraStar == 0x0)
    {
      m_statsSpectraStar = new OutputTable;
    }
    if (m_statsContigs == 0x0)
    {
      m_statsContigs = new MappedContigSetTable;
    }

    if (m_statsSummaryRawMS2 == 0x0)
    {
      m_statsSummaryRawMS2 = new OutputTable;
    }
    if (m_statsSummaryClustMS2 == 0x0)
    {
      m_statsSummaryClustMS2 = new OutputTable;
    }
    if (m_statsSummaryScored == 0x0)
    {
      m_statsSummaryScored = new OutputTable;
    }
    if (m_statsSummaryStar == 0x0)
    {
      m_statsSummaryStar = new OutputTable;
    }
    if (m_statsSummaryContigs == 0x0)
    {
      m_statsSummaryContigs = new MappedSPSStatTable;
    }
  }

  DEBUG_VAR(m_spectraRawIDs->size());
  DEBUG_VAR(m_spectraClustIDs->size());
  // Call parent invoke to annotate contigs
  if (!ExecReportSPSStats::invoke())
  {
    return false;
  }

  DEBUG_TRACE;

  // PSMs to annotate spectra
  PeptideSpectrumMatchSet PSMsRawMS2; // raw MS/MS
  PeptideSpectrumMatchSet PSMsClustMS2; // clustered MS/MS (specs_ms.pklbin)
  PeptideSpectrumMatchSet PSMsScored; // scored PRM (specs_scored.pklbin)
  PeptideSpectrumMatchSet PSMsStars; // star spectra (stars.pklbin)

  // Scan # and filename mappings to origin MS/MS spectra
  map<int, list<pair<int, string> > >* rawScanInfo;
  map<int, list<pair<int, string> > >* clustScanInfo;
  map<int, list<pair<int, string> > >* scoredScanInfo;
  map<int, list<pair<int, string> > >* starScanInfo;

  // Ions to annotate for each fragmentation mode
  vector<string> ionTypes(4);
  ionTypes[(unsigned int)Spectrum::FragType_CID] = "b,y,b++,y++";
  ionTypes[(unsigned int)Spectrum::FragType_HCD] = "b,y,b++,y++";
  ionTypes[(unsigned int)Spectrum::FragType_ETD] = "c,z,z.,z.++";
  ionTypes[(unsigned int)Spectrum::FragType_PRM] = "prm,srm+H2O";

  AAJumps jumps(1);

  // Annotate raw MS/MS spectra
  DEBUG_MSG("Annotating raw MS/MS spectra ...");
  PSMsRawMS2 = *m_spectraRawIDs;
  rawScanInfo = &m_rawScanInfo;
  PSMsRawMS2.addSpectra(m_spectraRawMS2);
  for (unsigned int i = 0; i < PSMsRawMS2.size(); i++)
  {
    if (PSMsRawMS2[i]->m_spectrum)
    {
      PSMsRawMS2[i]->annotate(PSMsRawMS2[i]->m_annotation,
                              ionTypes[(unsigned int)PSMsRawMS2[i]->m_spectrum->msFragType],
                              *m_ms2ScoringModels,
                              0,
                              0,
                              jumps);
    }
  }

  DEBUG_TRACE;

  // Annotate clustered MS/MS spectra
  DEBUG_MSG("Annotating clustered MS/MS spectra ...");
  if (m_paramsSPS->getValueInt("CLUSTER_MIN_SIZE", 0) > 0
      && m_paramsSPS->getValue("CLUSTER_TOOL", "MSCluster") == "MSCluster")
  {
    PSMsClustMS2 = *m_spectraClustIDs;
    clustScanInfo = &m_clustScanInfo;
  }
  else
  {
    PSMsClustMS2 = *m_spectraRawIDs;
    clustScanInfo = &m_rawScanInfo;
  }
  PSMsClustMS2.addSpectra(m_spectraClustMS2);
  for (unsigned int i = 0; i < PSMsClustMS2.size(); i++)
  {
    if (PSMsClustMS2[i]->m_spectrum)
    {
      PSMsClustMS2[i]->annotate(PSMsClustMS2[i]->m_annotation,
                                ionTypes[(unsigned int)PSMsClustMS2[i]->m_spectrum->msFragType],
                                *m_ms2ScoringModels,
                                0,
                                0,
                                jumps);
    }
  }

  DEBUG_TRACE;

  // Annotate scored spectra
  DEBUG_MSG("Annotating scored spectra ...");
  if (m_paramsSPS->getValueInt("CLUSTER_MIN_SIZE", 0) > 0)
  {
    PSMsScored = *m_spectraClustIDs;
    scoredScanInfo = &m_clustScanInfo;
  }
  else
  {
    PSMsScored = *m_spectraRawIDs;
    scoredScanInfo = &m_rawScanInfo;
  }
  PSMsScored.addSpectra(m_spectraScored);
  Spectrum::FragType prmType = Spectrum::FragType_PRM;
  for (unsigned int i = 0; i < PSMsScored.size(); i++)
  {
    if (PSMsScored[i]->m_spectrum)
    {
      PSMsScored[i]->annotate(PSMsScored[i]->m_annotation,
                              ionTypes[(unsigned int)Spectrum::FragType_PRM],
                              m_ms2ScoringModels->getModel(Spectrum::activationToString(prmType)),
                              0,
                              0,
                              jumps);
    }
  }

  DEBUG_TRACE;

  // Annotate star spectra
  DEBUG_MSG("Annotating star spectra ...");
  if (m_paramsSPS->getValueInt("CLUSTER_MIN_SIZE", 0) > 0)
  {
    PSMsStars = *m_spectraClustIDs;
    starScanInfo = &m_clustScanInfo;
  }
  else
  {
    PSMsStars = *m_spectraRawIDs;
    starScanInfo = &m_rawScanInfo;
  }
  PSMsStars.addSpectra(m_spectraStars);
  for (unsigned int i = 0; i < PSMsStars.size(); i++)
  {
    if (PSMsStars[i]->m_spectrum)
    {
      PSMsStars[i]->annotate(PSMsStars[i]->m_annotation,
                             ionTypes[(unsigned int)Spectrum::FragType_PRM],
                             m_ms2ScoringModels->getModel(Spectrum::activationToString(prmType)),
                             0,
                             0,
                             jumps);
    }
  }

  DEBUG_TRACE;

  // Fill in per-spectrum result tables
  prepareStatsSpectra(PSMsRawMS2, *rawScanInfo, *m_statsSpectraRawMS2);
  prepareStatsSpectra(PSMsClustMS2, *clustScanInfo, *m_statsSpectraClustMS2);
  prepareStatsSpectra(PSMsScored, *scoredScanInfo, *m_statsSpectraScored);
  prepareStatsSpectra(PSMsStars, *starScanInfo, *m_statsSpectraStar);

  // Fill in spectral summary result tables
  prepareStatsSummary(PSMsRawMS2, *rawScanInfo, *m_statsSummaryRawMS2);
  prepareStatsSummary(PSMsClustMS2, *clustScanInfo, *m_statsSummaryClustMS2);
  prepareStatsSummary(PSMsScored, *scoredScanInfo, *m_statsSummaryScored);
  prepareStatsSummary(PSMsStars, *starScanInfo, *m_statsSummaryStar);

  // Prepare per-contig result table
  m_statsContigs->mapped_sps_proj = m_mappedProj;

  // Prepare contig summary result table
  m_statsSummaryContigs->mapped_sps_proj = m_mappedProj;

  return true;
}

bool ExecQualityControl::loadInputData(void)
{

  if (ownInput)
  {
    if (!m_paramsSPS)
    {
      m_paramsSPS = new ParameterList;
    }
    if (!m_spectraRawMS2)
    {
      m_spectraRawMS2 = new SpecSet;
    }
    if (!m_spectraClustMS2)
    {
      m_spectraClustMS2 = new SpecSet;
    }
    if (!m_spectraScored)
    {
      m_spectraScored = new SpecSet;
    }
    if (!m_ms2ScoringModels)
    {
      m_ms2ScoringModels = new MS2ScoringModelSet;
    }
  }

  string scoringModelDir = m_params.getValue("MS2_SCORING_MODEL_DIR");
  ExecQualityControl::loadMS2ScoringModels(scoringModelDir, *m_ms2ScoringModels);

  string specIDFormat = m_params.getValue("SPEC_ID_FORMAT", "msgfdb");
  std::transform(specIDFormat.begin(),
                 specIDFormat.end(),
                 specIDFormat.begin(),
                 ::tolower);

  string spsDir = m_params.getValue("SPS_PROJECT_DIR");
  string spectraDir = getPath(spsDir, "spectra", false);

  string
      spsParamsFile = m_params.getValue("SPS_PARAMS_FILE",
                                        getPath(spsDir, "sps.params", false));
  if (!m_paramsSPS->readFromFile(spsParamsFile))
  {
    ERROR_MSG("Failed to read from \'" << spsParamsFile << "\'");
    return false;
  }

  string rawSpecFiles = getPath(spectraDir, "pklbin_files.txt", false);
  vector<string> rawSpecFileList;
  if (!readFilesFromFile(rawSpecFiles, rawSpecFileList))
  {
    ERROR_MSG("Failed to read from \'" << rawSpecFiles << "\'");
    return false;
  }

  m_spectraRawMS2->resize(0);
  for (unsigned int fIdx = 0; fIdx < rawSpecFileList.size(); fIdx++)
  {
    string locFile = getPath(spsDir, rawSpecFileList[fIdx], false);
    SpecSet locSpecs;
    DEBUG_MSG("Loading \'" << locFile << "\'");
    if (!locSpecs.Load(locFile.c_str()))
    {
      ERROR_MSG("Failed to read from \'" << locFile << "\'");
      return false;
    }
    DEBUG_VAR(locSpecs.size());
    m_spectraRawMS2->swapAppendSpecSet(locSpecs, false);
  }
  DEBUG_VAR(m_spectraRawMS2->size());

  // Scans need to be set to 1-based indices so the PSMs can be matched (their scans are also corrected)
  for (unsigned int i = 0; i < m_spectraRawMS2->size(); i++)
  {
    (*m_spectraRawMS2)[i].scan = i + 1;
  }

  string clustSpecFile = getPath(spectraDir, "specs_ms.pklbin", false);
  DEBUG_MSG("Loading \'" << clustSpecFile << "\'");
  if (!m_spectraClustMS2->Load(clustSpecFile.c_str()))
  {
    ERROR_MSG("Failed to read from \'" << clustSpecFile << "\'");
    return false;
  }
  DEBUG_VAR(m_spectraClustMS2->size());

  // Scans need to be set to 1-based indices so the PSMs can be matched (their scans are also corrected)
  for (unsigned int i = 0; i < m_spectraClustMS2->size(); i++)
  {
    (*m_spectraClustMS2)[i].scan = i + 1;
  }

  string scoredSpecFile = getPath(spectraDir, "specs_scored.pklbin", false);
  DEBUG_MSG("Loading \'" << scoredSpecFile << "\'");
  if (!m_spectraScored->Load(scoredSpecFile.c_str()))
  {
    ERROR_MSG("Failed to read from \'" << scoredSpecFile << "\'");
    return false;
  }
  DEBUG_VAR(m_spectraScored->size());

  // Scans need to be set to 1-based indices so the PSMs can be matched (their scans are also corrected)
  for (unsigned int i = 0; i < m_spectraScored->size(); i++)
  {
    (*m_spectraScored)[i].scan = i + 1;
  }

  string assemblyDir = getPath(spsDir, "assembly", false);
  string homologyDir = getPath(spsDir, "homology", false);
  string alignsDir = getPath(spsDir, "aligns", false);

  // ADD ADDITIONAL INPUT FILE LOADS HERE


  // fill in the appropriate parameters so we can call ExecReportSPSStats::loadInputData()
  ParameterList spsStatsParams;
  spsStatsParams.setValue("SPS_PROJECT_DIR", spsDir);
  spsStatsParams.setValue("INPUT_CONTIGS", getPath(assemblyDir,
                                                   "sps_seqs.pklbin",
                                                   false));
  spsStatsParams.setValue("INPUT_STARS", getPath(spectraDir,
                                                 "stars.pklbin",
                                                 false));
  spsStatsParams.setValue("INPUT_CONTIG_ABINFO", getPath(assemblyDir,
                                                         "component_info.bin",
                                                         false));
  spsStatsParams.setValue("INPUT_MIDX", getPath(homologyDir,
                                                "contigs_midx.pklbin",
                                                false));
  spsStatsParams.setValue("INPUT_MP", getPath(homologyDir,
                                              "contigs_mp.bin",
                                              false));
  spsStatsParams.setValue("INPUT_REF_INDICES", getPath(spectraDir,
                                                       "contigs_indices.bin",
                                                       false));
  spsStatsParams.setValue("INPUT_FASTA",
                          m_paramsSPS->getValue("FASTA_DATABASE"));
  spsStatsParams.setValue("INPUT_ION_TYPES", getPath(scoringModelDir,
                                                     "model_prm.txt",
                                                     false));
  spsStatsParams.setValue("SPEC_ID_FORMAT", specIDFormat);
  spsStatsParams.setValue("INPUT_SPEC_IDS", m_params.getValue("INPUT_SPEC_IDS"));
  spsStatsParams.addIfExists(m_params, "INPUT_CLUST_SPEC_IDS");

  struct stat buf;
  string clustFilePath = getPath(spectraDir, "clusterData.bin", false);
  if (stat(clustFilePath.c_str(), &buf) != -1)
  {
    spsStatsParams.setValue("INPUT_CLUSTERS_DIR", spsDir);
  }
  spsStatsParams.setValue("INPUT_FILE_INDEX", getPath(spectraDir,
                                                      "input_index.txt",
                                                      false));
  spsStatsParams.setValue("INPUT_FILE_MAPPING", getPath(spectraDir,
                                                        "input_mapping.bin",
                                                        false));
  spsStatsParams.setValue("INPUT_SCAN_REF_FILES", getPath(spectraDir,
                                                          "bin_files.txt",
                                                          false));
  spsStatsParams.setValue("TARGET_PROTEINS",
                          m_params.getValue("TARGET_PROTEINS", ""));

  spsStatsParams.addIfExists(m_params, "MIN_CONTIG_AA_TAG");
  spsStatsParams.addIfExists(m_params, "MIN_CONTIG_DB_MP");
  spsStatsParams.addIfExists(m_params, "ENDS_CHOP");
  spsStatsParams.addIfExists(*m_paramsSPS, "TOLERANCE_PEAK");

  ParameterList myParamsCopy(m_params);
  m_params = spsStatsParams;
  bool res = ExecReportSPSStats::loadInputData();
  if (!res)
  {
    return false;
  }
  m_params = myParamsCopy;

  float peakTol = m_paramsSPS->getValueFloat("TOLERANCE_PEAK");
  m_spectraRawMS2->setPeakTolerance(peakTol, false);
  m_spectraClustMS2->setPeakTolerance(peakTol, false);
  m_spectraScored->setPeakTolerance(peakTol, false);

  return true;
}

bool ExecQualityControl::saveInputData(std::vector<std::string> & filenames)
{
  return false;
}

bool ExecQualityControl::saveOutputData(void)
{
  string outDir = m_params.getValue("OUTPUT_STATS_DIR", "");
  const char* delim = "\t";

  string outFile = "per-spectrum_raw_MS2.tsv";
  outFile = m_params.getValue("OUTPUT_PER_SPECTRUM_STATS_RAW", getPath(outDir,
                                                                       outFile,
                                                                       false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSpectraRawMS2->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "per-spectrum_clust_MS2.tsv";
  outFile = m_params.getValue("OUTPUT_PER_SPECTRUM_STATS_CLUST",
                              getPath(outDir, outFile, false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSpectraClustMS2->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "per-spectrum_scored_PRM.tsv";
  outFile = m_params.getValue("OUTPUT_PER_SPECTRUM_STATS_SCORED",
                              getPath(outDir, outFile, false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSpectraScored->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "per-spectrum_star_PRM.tsv";
  outFile = m_params.getValue("OUTPUT_PER_SPECTRUM_STATS_STAR",
                              getPath(outDir, outFile, false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSpectraStar->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "per-contig.tsv";
  outFile = m_params.getValue("OUTPUT_SPECTRA_STATS_CONTIG", getPath(outDir,
                                                                     outFile,
                                                                     false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsContigs->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "summary_raw_MS2.tsv";
  outFile = m_params.getValue("OUTPUT_SUMMARY_STATS_RAW", getPath(outDir,
                                                                  outFile,
                                                                  false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSummaryRawMS2->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "summary_clust_MS2.tsv";
  outFile = m_params.getValue("OUTPUT_SUMMARY_STATS_CLUST", getPath(outDir,
                                                                    outFile,
                                                                    false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSummaryClustMS2->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "summary_scored_PRM.tsv";
  outFile = m_params.getValue("OUTPUT_SUMMARY_STATS_SCORED", getPath(outDir,
                                                                     outFile,
                                                                     false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSummaryScored->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "summary_star_PRM.tsv";
  outFile = m_params.getValue("OUTPUT_SUMMARY_STATS_STAR", getPath(outDir,
                                                                   outFile,
                                                                   false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSummaryStar->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  outFile = "summary_contig.tsv";
  outFile = m_params.getValue("OUTPUT_SUMMARY_STATS_CONTUG", getPath(outDir,
                                                                     outFile,
                                                                     false));
  DEBUG_MSG("Saving statistics to " << outFile << " ...");
  if (!m_statsSummaryContigs->printToCSV(outFile.c_str(), delim))
  {
    ERROR_MSG("Failed to save to " << outFile);
    return false;
  }

  return true;
}

bool ExecQualityControl::loadOutputData(void)
{
  return false;
}

vector<ExecBase *>
const & ExecQualityControl::split(int numSplit)
{
  m_subModules.resize(0);
  return m_subModules;
}

// -------------------------------------------------------------------------
bool ExecQualityControl::merge(void)
{
  return false;
}

bool ExecQualityControl::validateParams(std::string & error)
{
  m_isValid = false;
  VALIDATE_PARAM_EXIST("SPS_PROJECT_DIR");
  VALIDATE_PARAM_EXIST("MS2_SCORING_MODEL_DIR");
  //VALIDATE_PARAM_EXIST("SPS_PARAMS_FILE");
  VALIDATE_PARAM_EXIST("INPUT_SPEC_IDS");
  //VALIDATE_PARAM_EXIST("INPUT_CLUST_SPEC_IDS");
  m_isValid = true;
  return true;
}

void ExecQualityControl::prepareStatsSpectra(const PeptideSpectrumMatchSet &psms,
                                             const map<int, list<pair<int,
                                                 string> > > &scanInfo,
                                             OutputTable &statsTable)
{
  const unsigned int numCols = 11;
  unsigned int row = 0;
  unsigned int col = 0;
  statsTable.values.resize(psms.size() + 1);
  statsTable.values[row].resize(numCols);
  statsTable.setValue(row, col++, "Index", true);
  statsTable.setValue(row, col++, "Fragmentation", true);
  statsTable.setValue(row, col++, "Charge", true);
  statsTable.setValue(row, col++, "Peptide", true);
  statsTable.setValue(row, col++, "Score", true);
  statsTable.setValue(row, col++, "P-value", true);
  statsTable.setValue(row, col++, "%ExplainedIntensity", true);
  statsTable.setValue(row, col++, "%ObservedBreaks", true);
  statsTable.setValue(row, col++, "#Peaks", true);
  statsTable.setValue(row, col++, "Scan", true);
  statsTable.setValue(row, col++, "Filename", true);
  row++;
  col = 0;

  SpectrumAnnotStatistics stats;
  string allIons = "all";

  for (unsigned int i = 0; i < psms.size(); i++)
  {
    statsTable.values[row].resize(numCols);
    const psmPtr &psm = psms[i];
    statsTable.setValue(row, col++, parseInt(psm->m_scanNum - 1));
    string scans = "";
    string fileNames = "";
    const list<pair<int, string> > &childScans = scanInfo.at(psm->m_scanNum);
    for (list<pair<int, string> >::const_iterator cIt = childScans.begin(); cIt
        != childScans.end(); cIt++)
    {
      scans += parseInt(cIt->first);
      fileNames += cIt->second;
      scans += "/";
      fileNames += "/";
    }
    scans.erase(scans.length() - 1, 1);
    fileNames.erase(fileNames.length() - 1, 1);
    statsTable.setValue(row,
                        col++,
                        Spectrum::activationToString(psm->m_spectrum->msFragType));
    statsTable.setValue(row, col++, parseInt(psm->m_spectrum->parentCharge));
    statsTable.setValue(row, col++, psm->m_annotation);
    statsTable.setValue(row, col++, parseFloat(psm->m_score, 5));
    statsTable.setValue(row, col++, parseDouble(psm->m_pValue, 5));
    float explainedIntensity = 0;//stats.percentExplainedIntensity(*psm, allIons);
    statsTable.setValue(row, col++, parseFloat(explainedIntensity, 2));
    float observedBreaks = 0;//stats.observedBreaks(*psm, allIons);
    statsTable.setValue(row, col++, parseFloat(observedBreaks, 2));
    statsTable.setValue(row, col++, parseInt(psm->m_spectrum->size()));
    statsTable.setValue(row, col++, scans);
    statsTable.setValue(row, col++, fileNames);
    row++;
    col = 0;
  }
}

void ExecQualityControl::prepareStatsSummary(const PeptideSpectrumMatchSet &psms,
                                             const map<int, list<pair<int,
                                                 string> > > &scanInfo,
                                             OutputTable &statsTable)
{
  const unsigned int numCols = 2;
  const unsigned int numRows = 4;

  statsTable.values.resize(numRows);
  unsigned int row = 0;
  unsigned int col = 0;
  string allIons = "all";
  SpectrumAnnotStatistics stats;
  statsTable.values[row].resize(numCols);

  statsTable.setValue(row, col++, "#PSMs", true);
  statsTable.setValue(row, col++, parseInt(psms.size()));
  row++;
  col = 0;

  statsTable.values[row].resize(numCols);
  statsTable.setValue(row, col++, "#Non-empty PSMs", true);
  unsigned int numSpecs = stats.numNonEmptySpectra(psms);
  statsTable.setValue(row, col++, parseInt(numSpecs));
  row++;
  col = 0;

  statsTable.values[row].resize(numCols);
  statsTable.setValue(row, col++, "%ExplainedIntensity", true);
  float explainedIntensity = 0;//stats.percentExplainedIntensity(psms, allIons);
  statsTable.setValue(row, col++, parseFloat(explainedIntensity, 2));
  row++;
  col = 0;

  statsTable.values[row].resize(numCols);
  statsTable.setValue(row, col++, "%ObservedBreaks", true);
  float observedBreaks = 0;//stats.observedBreaks(psms, allIons);
  statsTable.setValue(row, col++, parseFloat(observedBreaks, 2));
  row++;
  col = 0;

}
