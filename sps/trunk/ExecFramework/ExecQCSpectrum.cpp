/*
 * ExecQCSpectrum.cpp
 *
 *  Created on: Aug 29, 2013
 *      Author: aguthals
 */

#include "ExecQCSpectrum.h"
#include "ExecMergeConvert.h"

using namespace std;

namespace specnets
{
  void ExecQCSpectrum::loadMS2ScoringModels(const string &modelDir,
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

  ExecQCSpectrum::ExecQCSpectrum() :
    ExecBase(), m_ms2ScoringModels(0x0), m_inputSpectra(0x0),
        m_inputClusters(0x0), m_inputPSMs(0x0), m_statsParams(0x0),
        m_outputPSMs(0x0), m_spectraStats(0x0), m_spectraStatsHeader(0x0),
        ownInput(true), ownOutput(true)
  {
    m_ms2ScoringModels = new MS2ScoringModelSet;
    m_inputSpectra = new SpecSet;
    m_inputClusters = new ClusterSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_outputPSMs = new PeptideSpectrumMatchSet;
    m_statsParams = new SpectrumAnnotParameterList;
    m_spectraStats = new std::vector<vector<float> >;
    m_spectraStatsHeader = new std::vector<string>;
    m_name = "ExecQCSpectrum";
    m_type = "ExecQCSpectrum";
  }

  ExecQCSpectrum::ExecQCSpectrum(const ParameterList & inputParams) :
    ExecBase(inputParams), m_ms2ScoringModels(0x0), m_inputSpectra(0x0),
        m_inputClusters(0x0), m_inputPSMs(0x0), m_statsParams(0x0),
        m_outputPSMs(0x0), m_spectraStats(0x0), m_spectraStatsHeader(0x0),
        ownInput(true), ownOutput(true)
  {
    m_ms2ScoringModels = new MS2ScoringModelSet;
    m_inputSpectra = new SpecSet;
    m_inputClusters = new ClusterSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_outputPSMs = new PeptideSpectrumMatchSet;
    m_statsParams = new SpectrumAnnotParameterList;
    m_spectraStats = new std::vector<vector<float> >;
    m_spectraStatsHeader = new std::vector<string>;
    m_name = "ExecQCSpectrum";
    m_type = "ExecQCSpectrum";
  }

  ExecQCSpectrum::ExecQCSpectrum(const ParameterList & inputParams,
                                 MS2ScoringModelSet *ms2ScoringModels,
                                 SpecSet *inputSpectra,
                                 ClusterSet *inputClusters,
                                 PeptideSpectrumMatchSet *inputPSMs,
                                 SpectrumAnnotParameterList * statsParams,
                                 PeptideSpectrumMatchSet *outputPSMs,
                                 std::vector<vector<float> > * spectraStats,
                                 std::vector<string> * spectraStatsHeader) :
    ExecBase(inputParams), m_ms2ScoringModels(ms2ScoringModels),
        m_inputSpectra(inputSpectra), m_inputClusters(inputClusters),
        m_inputPSMs(inputPSMs), m_statsParams(statsParams),
        m_outputPSMs(outputPSMs), m_spectraStats(spectraStats),
        m_spectraStatsHeader(spectraStatsHeader), ownInput(false),
        ownOutput(false)
  {
    m_name = "ExecQCSpectrum";
    m_type = "ExecQCSpectrum";
  }

  ExecQCSpectrum::~ExecQCSpectrum(void)
  {
    if (ownInput)
    {
      delete m_ms2ScoringModels;
      delete m_inputSpectra;
      delete m_inputClusters;
      delete m_inputPSMs;
      delete m_statsParams;
    }

    if (ownOutput)
    {
      delete m_outputPSMs;
      delete m_spectraStats;
      delete m_spectraStatsHeader;
    }
  }

  ExecBase * ExecQCSpectrum::clone(const ParameterList & inputParams) const
  {
    return new ExecQCSpectrum(inputParams);
  }

  bool ExecQCSpectrum::invoke(void)
  {
    m_outputPSMs->operator =(*m_inputPSMs);
    if (m_inputClusters->size() > 0)
    {
      map<int, list<pair<int, string> > > clusterInfo;
      m_inputClusters->getScanToFileMapping(clusterInfo);
      m_outputPSMs->cluster(clusterInfo);
      m_outputPSMs->addSpectra(m_inputSpectra, true);
    }
    else
    {
      m_outputPSMs->addSpectraByFilename(m_inputSpectra, true);
    }

    float prmOffset = m_params.getValueFloat("PRM_OFFSET", 0);
    float srmOffset = m_params.getValueFloat("SRM_OFFSET", 0);

    vector<psmPtr>::iterator resultIterator;

    SpectrumAnnotStatistics stats;

    AAJumps jumps(1);
    string ionTypes = "all";

    for (resultIterator = m_outputPSMs->m_psmSet.begin(); resultIterator
        != m_outputPSMs->m_psmSet.end(); resultIterator++)
    {
      PeptideSpectrumMatch* psm = resultIterator->get();
      /*
       DEBUG_VAR(psm->m_annotation);
       DEBUG_TRACE;
       DEBUG_VAR(psm->m_spectrum->parentMass);
       */
      //annotate
      psm->annotate(psm->m_annotation,
                    ionTypes,
                    *m_ms2ScoringModels,
                    prmOffset,
                    srmOffset,
                    jumps);
    }

    bool addedHeaders = false;
    if (m_params.getValueInt("REPORT_CUMMULATIVE", 0) > 0)
    {
      vector<float> currStatistics(m_statsParams->m_params.size() + 1);
      currStatistics[0] = -1.0;
      addedHeaders = true;
      m_spectraStatsHeader->push_back("Scan#");
      for (int i = 0; i < m_statsParams->m_params.size(); i++)
      {
        SpectrumAnnotParameter * currParam = &(m_statsParams->m_params[i]);

        //add to header
        m_spectraStatsHeader->push_back(currParam->statisticName);

        if (currParam->statistic.compare("%explained intensity") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // we're ignoring this parameter
            //do nothing
          }
          else
          {
            float explainedIntensity =
                stats.percentExplainedIntensity(*m_outputPSMs, *currParam);
            currStatistics[i + 1] = explainedIntensity;
          }
        }
        else if (currParam->statistic.compare("%explained peaks") == 0)
        {
          float explainedPeaks = stats.percentExplainedPeaks(*m_outputPSMs,
                                                             *currParam);
          currStatistics[i + 1] = explainedPeaks;
        }
        else if (currParam->statistic.compare("%observed ions") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // check if we're ignoring this parameter
            // do nothing
          }
          else
          {
            float observedIons = stats.percentObservedIons(*m_outputPSMs,
                                                           *currParam);
            currStatistics[i + 1] = observedIons;
          }
        }
        else if (currParam->statistic.compare("%observed breaks") == 0)
        {
          float observedBreaks =
              stats.observedBreaks(*m_outputPSMs, *currParam);
          currStatistics[i + 1] = observedBreaks;
        }
        else if (currParam->statistic.compare("") != 0)
        {
          currStatistics[i + 1] = -1.0;
        }
      }
      m_spectraStats->push_back(currStatistics);
    }
    //annotate spectra and output stats
    for (resultIterator = m_outputPSMs->m_psmSet.begin(); resultIterator
        != m_outputPSMs->m_psmSet.end(); resultIterator++)
    {
      PeptideSpectrumMatch& psm = *(resultIterator->get());

      int scanNum = psm.m_scanNum;

      Spectrum * currSpec = psm.m_spectrum;

      if (currSpec == NULL)
      {
        WARN_MSG("Spectrum for scan " << psm.m_scanNum << " not defined!");
        continue;
      }

      vector<float> currStatistics(m_statsParams->m_params.size() + 1); //vector to hold stats for current spectra

      currStatistics[0] = (float)scanNum;

      //add to header
      if (!addedHeaders && resultIterator == m_outputPSMs->m_psmSet.begin())
      {
        m_spectraStatsHeader->push_back("Scan#");
      }

      for (int i = 0; i < m_statsParams->m_params.size(); i++)
      {
        SpectrumAnnotParameter * currParam = &(m_statsParams->m_params[i]);

        //add to header
        if (!addedHeaders && resultIterator == m_outputPSMs->m_psmSet.begin())
        {
          m_spectraStatsHeader->push_back(currParam->statisticName);
        }

        if (currParam->statistic.compare("%explained intensity") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // we're ignoring this parameter
            //do nothing
          }
          else
          {
            float explainedIntensity =
                stats.percentExplainedIntensity(psm, *currParam);
            currStatistics[i + 1] = explainedIntensity;
          }
        }
        else if (currParam->statistic.compare("number of gaps") == 0)
        {
          int numGaps = psm.countGaps();
          currStatistics[i + 1] = (float)numGaps;
        }
        else if (currParam->statistic.compare("%explained peaks") == 0)
        {
          float explainedPeaks = stats.percentExplainedPeaks(psm, *currParam);
          currStatistics[i + 1] = explainedPeaks;
        }
        else if (currParam->statistic.compare("%observed ions") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // check if we're ignoring this parameter
            // do nothing
          }
          else
          {
            float observedIons = stats.percentObservedIons(psm, *currParam);
            currStatistics[i + 1] = observedIons;
          }
        }
        else if (currParam->statistic.compare("total peaks") == 0)
        {
          int totalPeaks = stats.totalPeaks(psm);
          currStatistics[i + 1] = (float)totalPeaks;
        }
        else if (currParam->statistic.compare("parent mass error ppm") == 0)
        {
          float parentMassErrorPpm =
              stats.parentMassErrorPPM(psm, psm.m_charge);
          //parentMassErrorPpm = abs(parentMassErrorPpm / (1 - peakTol)); //output error adjusted by tolerance
          currStatistics[i + 1] = parentMassErrorPpm;
        }
        else if (currParam->statistic.compare("parent mass error da") == 0)
        {
          float parentMassErrorDa = stats.parentMassErrorDa(psm, psm.m_charge);
          //parentMassErrorDa = abs(parentMassErrorDa / (1 - peakTol)); //output error adjusted by tolerance
          currStatistics[i + 1] = parentMassErrorDa;
        }
        else if (currParam->statistic.compare("%observed breaks") == 0)
        {
          float observedBreaks = stats.observedBreaks(psm, *currParam);
          currStatistics[i + 1] = observedBreaks;
        }
        else if (currParam->statistic.compare("%observed difference") == 0)
        {
          vector<string> ionList;
          string
              ionNames =
                  currParam->getFragSpecificIonNames(Spectrum::activationToString(currSpec->msFragType));
          stringSplit(ionNames, ionList, ",");
          SpectrumAnnotParameter p1;
          SpectrumAnnotParameter p2;
          p1.ionNames = ionList[0];
          p2.ionNames = ionList[1];
          float firstValue = stats.percentObservedIons(psm, p1);
          float secondValue = stats.percentObservedIons(psm, p2);
          float observedDifference = abs(firstValue - secondValue);
          currStatistics[i + 1] = observedDifference;
        }
        else if (currParam->statistic.compare("%explained difference") == 0)
        {
          vector<string> ionList;
          string
              ionNames =
                  currParam->getFragSpecificIonNames(Spectrum::activationToString(currSpec->msFragType));
          stringSplit(ionNames, ionList, ",");
          SpectrumAnnotParameter p1;
          SpectrumAnnotParameter p2;
          p1.ionNames = ionList[0];
          p2.ionNames = ionList[1];
          float firstValue = stats.percentExplainedIntensity(psm, p1);
          float secondValue = stats.percentExplainedIntensity(psm, p2);
          float explainedDifference = abs(firstValue - secondValue);
          currStatistics[i + 1] = explainedDifference;
        }
        else if (currParam->statistic.compare("") != 0)
        {
          DEBUG_MSG("Unknown parameter type: <" << currParam->statistic << ">");
        }
      }
      addedHeaders = true;
      m_spectraStats->push_back(currStatistics);
    }
    return true;
  }

  bool ExecQCSpectrum::loadInputData(void)
  {
    string psmType = m_params.getValue("INPUT_PSMS_TYPE", "");
    std::transform(psmType.begin(), psmType.end(), psmType.begin(), ::tolower);

    if (psmType.length() == 0)
    {
      if (!m_inputPSMs->loadFromFile(m_params.getValue("INPUT_PSMS").c_str()))
      {
        ERROR_MSG("Failed to load PSMs in native format from \'"
            << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }
    }
    else if (psmType == "msgfplus")
    {
      if (!m_inputPSMs->loadMSGFPlusResultsFile(m_params.getValue("INPUT_PSMS").c_str()))
      {
        ERROR_MSG("Failed to load PSMs in MSGFPlus format from \'"
            << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }
    }
    else if (psmType == "msgfdb")
    {
      if (!m_inputPSMs->loadMSGFDBResultsFile(m_params.getValue("INPUT_PSMS").c_str()))
      {
        ERROR_MSG("Failed to load PSMs in MSGFDB format from \'"
            << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }
    }
    else if (psmType == "moda")
    {
      if (!m_inputPSMs->loadModaResultsFile(m_params.getValue("INPUT_PSMS").c_str()))
      {
        ERROR_MSG("Failed to load PSMs in ModA format from \'"
            << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }
    }
    else if (psmType == "InspecT")
    {
      if (!m_inputPSMs->loadInspectResultsFile(m_params.getValue("INPUT_PSMS").c_str()))
      {
        ERROR_MSG("Failed to load PSMs in InspecT format from \'"
            << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }
    }
    else
    {
      ERROR_MSG("Unrecognized PSM format \'" << psmType << "\'");
      return false;
    }

    string exeDir = m_params.getValue("EXE_DIR");

    if (!ExecMergeConvert::loadSpecsetMultiple(exeDir,
                                               m_params.getValue("INPUT_SPECTRA"),
                                               m_inputSpectra))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_SPECTRA"));
      return false;
    }

    if (m_params.getValueBool("SET_SCAN_NUMS", false))
    {
      DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
      for (int i = 0; i < m_inputSpectra->size(); i++)
      {
        (*m_inputSpectra)[i].scan = (unsigned int)(i + 1);
      }
    }

    if (m_params.exists("INPUT_CLUSTERS"))
    {
      if (!m_inputClusters->loadBinaryFile(m_params.getValue("INPUT_CLUSTERS").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CLUSTERS"));
        return false;
      }
    }

    /*
     for (unsigned int i = 0; i < m_inputSpectra->size(); i++)
     {
     DEBUG_MSG("Spectrum " << i << ": " << (*m_inputSpectra)[i].getUniqueID());
     }
     abort();
     */

    string rsrcDir = getPath(exeDir, "resources", false);

    loadMS2ScoringModels(rsrcDir, *m_ms2ScoringModels);

    DEBUG_VAR(m_inputSpectra->size());
    if (m_params.exists("TOLERANCE_PEAK_PPM"))
    {
      float ppmTol = m_params.getValueFloat("TOLERANCE_PEAK_PPM");
      m_inputSpectra->setPeakTolerance(ppmTol, true);
    }
    else if (m_params.exists("TOLERANCE_PEAK"))
    {
      float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
      m_inputSpectra->setPeakTolerance(peakTol, false);
    }

    DEBUG_TRACE;

    //Load in statistics generation parameters
    if (m_params.exists("STATISTICS_CONFIG"))
    {
      if (!m_statsParams->loadSpectrumAnnotFile(m_params.getValue("STATISTICS_CONFIG").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("STATISTICS_CONFIG"));
        return false;
      }
    }

    DEBUG_TRACE;

    return true;
  }

  bool ExecQCSpectrum::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  bool ExecQCSpectrum::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_PSMS"))
    {
      if (!m_outputPSMs->saveToFile(m_params.getValue("OUTPUT_PSMS").c_str(),
                                    true))
      {
        ERROR_MSG("Failed to save output PSMs to \'"
            << m_params.getValue("OUTPUT_PSMS") << "\'!");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_SPECTRA_STATS"))
    {
      string delim = "\t";

      ofstream outputFile;

      DEBUG_MSG("Opening output stats file...");
      outputFile.open(m_params.getValue("OUTPUT_SPECTRA_STATS").c_str(),
                      ios::out | ios::trunc | ios::binary);
      if (outputFile.fail())
      {
        ERROR_MSG("Unable to open stats file! "
            << m_params.getValue("OUTPUT_SPECTRA_STATS"));
        return false;
      }

      if (outputFile.is_open() && outputFile.good())
      {
        //output header

        outputFile << (*m_spectraStatsHeader)[0];
        for (int i = 1; i < m_spectraStatsHeader->size(); i++)
        {
          outputFile << delim << (*m_spectraStatsHeader)[i];
        }
        outputFile << endl;

        for (int i = 0; i < m_spectraStats->size(); i++)
        {
          outputFile << (*m_spectraStats)[i][0]; //scan number

          for (int j = 1; j < (*m_spectraStats)[i].size(); j++)
          {
            outputFile << delim << (*m_spectraStats)[i][j];
          }
          outputFile << endl;
        }
      }
      else
      {
        ERROR_MSG("Unable to open file!");
        return false;
      }
    }
    return true;
  }

  bool ExecQCSpectrum::loadOutputData(void)
  {
    return false;
  }

  vector<ExecBase *>
  const & ExecQCSpectrum::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecQCSpectrum::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecQCSpectrum::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("INPUT_PSMS");
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
    VALIDATE_PARAM_EXIST("EXE_DIR");
    m_isValid = true;
    return true;
  }
}
