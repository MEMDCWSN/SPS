// Header Include
#include "ExecSvmStatistics.h"

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecSvmStatistics::ExecSvmStatistics(void) :
    ownInput(true), ownOutput(true), ownStats(true), ownModels(true),
        m_starSpectraStatsHeader(0x0), m_starSpectraStats(0x0),
        m_specsMsSpectraStatsHeader(0x0), m_specsMsSpectraStats(0x0),
        m_specsScoredSpectraStatsHeader(0x0), m_specsScoredSpectraStats(0x0),
        m_peptideResults(0x0), m_peptideOutputResults(0x0),
        m_svmScaleParamsCharge1(0x0), m_svmScaleParamsCharge2(0x0),
        m_svmScaleParamsCharge3(0x0), m_svmModelCharge1(0x0),
        m_svmModelCharge2(0x0), m_svmModelCharge3(0x0), m_msAnnotParams(0x0),
        m_starAnnotParams(0x0), m_specsScoredAnnotParams(0x0)
  {
    m_name = "ExecSvmStatistics";
    m_type = "ExecSvmStatistics";
  }

  // -------------------------------------------------------------------------
  ExecSvmStatistics::ExecSvmStatistics(const ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(true), ownOutput(true), ownStats(true),
        ownModels(true), m_starSpectraStatsHeader(0x0),
        m_starSpectraStats(0x0), m_specsMsSpectraStatsHeader(0x0),
        m_specsMsSpectraStats(0x0), m_specsScoredSpectraStatsHeader(0x0),
        m_specsScoredSpectraStats(0x0), m_peptideResults(0x0),
        m_peptideOutputResults(0x0), m_svmScaleParamsCharge1(0x0),
        m_svmScaleParamsCharge2(0x0), m_svmScaleParamsCharge3(0x0),
        m_svmModelCharge1(0x0), m_svmModelCharge2(0x0), m_svmModelCharge3(0x0),
        m_msAnnotParams(0x0), m_starAnnotParams(0x0),
        m_specsScoredAnnotParams(0x0)
  {
    m_name = "ExecSvmStatistics";
    m_type = "ExecSvmStatistics";
  }

  // -------------------------------------------------------------------------
  ExecSvmStatistics::ExecSvmStatistics(const ParameterList & inputParams,
                                       std::vector<string> * starSpectraStatsHeader,
                                       std::vector<vector<float> > * starSpectraStats,
                                       std::vector<string> * specsMsSpectraStatsHeader,
                                       std::vector<vector<float> > * specsMsSpectraStats,
                                       std::vector<string> * specsScoredSpectraStatsHeader,
                                       std::vector<vector<float> > * specsScoredSpectraStats,
                                       PeptideSpectrumMatchSet * peptideResults,
                                       PeptideSpectrumMatchSet * peptideOutputResults) :
    ExecBase(inputParams), ownInput(false), ownStats(false), ownOutput(false),
        ownModels(true), m_starSpectraStatsHeader(starSpectraStatsHeader),
        m_starSpectraStats(starSpectraStats),
        m_specsMsSpectraStatsHeader(specsMsSpectraStatsHeader),
        m_specsMsSpectraStats(specsMsSpectraStats),
        m_specsScoredSpectraStatsHeader(specsScoredSpectraStatsHeader),
        m_specsScoredSpectraStats(specsScoredSpectraStats),
        m_peptideResults(peptideResults),
        m_peptideOutputResults(peptideOutputResults),
        m_svmScaleParamsCharge1(0x0), m_svmScaleParamsCharge2(0x0),
        m_svmScaleParamsCharge3(0x0), m_svmModelCharge1(0x0),
        m_svmModelCharge2(0x0), m_svmModelCharge3(0x0), m_msAnnotParams(0x0),
        m_starAnnotParams(0x0), m_specsScoredAnnotParams(0x0)
  {
    m_name = "ExecSvmStatistics";
    m_type = "ExecSvmStatistics";
  }
  // -------------------------------------------------------------------------
  ExecSvmStatistics::ExecSvmStatistics(const ParameterList & inputParams,
                                       PeptideSpectrumMatchSet * peptideResults,
                                       PeptideSpectrumMatchSet * peptideOutputResults) :
    ExecBase(inputParams), ownInput(false), ownOutput(false), ownStats(true),
        ownModels(true), m_peptideResults(peptideResults),
        m_peptideOutputResults(peptideOutputResults),
        m_starSpectraStatsHeader(0x0), m_starSpectraStats(0x0),
        m_specsMsSpectraStatsHeader(0x0), m_specsMsSpectraStats(0x0),
        m_specsScoredSpectraStatsHeader(0x0), m_specsScoredSpectraStats(0x0),
        m_svmScaleParamsCharge1(0x0), m_svmScaleParamsCharge2(0x0),
        m_svmScaleParamsCharge3(0x0), m_svmModelCharge1(0x0),
        m_svmModelCharge2(0x0), m_svmModelCharge3(0x0), m_msAnnotParams(0x0),
        m_starAnnotParams(0x0), m_specsScoredAnnotParams(0x0)
  {
    m_name = "ExecSvmStatistics";
    m_type = "ExecSvmStatistics";
  }

  // -------------------------------------------------------------------------
  ExecSvmStatistics::~ExecSvmStatistics(void)
  {
    if (ownStats)
    {
      if (m_starSpectraStatsHeader)
      {
        delete m_starSpectraStatsHeader;
        m_starSpectraStatsHeader = 0x0;
      }
      if (m_starSpectraStats)
      {
        delete m_starSpectraStats;
        m_starSpectraStats = 0x0;
      }
      if (m_specsMsSpectraStatsHeader)
      {
        delete m_specsMsSpectraStatsHeader;
        m_specsMsSpectraStatsHeader = 0x0;
      }
      if (m_specsMsSpectraStats)
      {
        delete m_specsMsSpectraStats;
        m_specsMsSpectraStats = 0x0;
      }
      if (m_specsScoredSpectraStatsHeader)
      {
        delete m_specsScoredSpectraStatsHeader;
        m_specsScoredSpectraStatsHeader = 0x0;
      }
      if (m_specsScoredSpectraStats)
      {
        delete m_specsScoredSpectraStats;
        m_specsScoredSpectraStats = 0x0;
      }
    }
    if (ownInput)
    {
      if (m_peptideResults)
      {
        delete m_peptideResults;
        m_peptideResults = 0x0;
      }
    }
    if (ownModels)
    {
      if (m_svmScaleParamsCharge1)
      {
        delete m_svmScaleParamsCharge1;
        m_svmScaleParamsCharge1 = 0x0;
      }
      if (m_svmScaleParamsCharge2)
      {
        delete m_svmScaleParamsCharge2;
        m_svmScaleParamsCharge2 = 0x0;
      }
      if (m_svmScaleParamsCharge3)
      {
        delete m_svmScaleParamsCharge3;
        m_svmScaleParamsCharge3 = 0x0;
      }
      if (m_svmModelCharge1)
      {
        delete m_svmModelCharge1;
        m_svmModelCharge1 = 0x0;
      }
      if (m_svmModelCharge2)
      {
        delete m_svmModelCharge2;
        m_svmModelCharge2 = 0x0;
      }
      if (m_svmModelCharge3)
      {
        delete m_svmModelCharge3;
        m_svmModelCharge3 = 0x0;
      }
      if (m_msAnnotParams)
      {
        delete m_msAnnotParams;
        m_msAnnotParams = 0x0;
      }
      if (m_starAnnotParams)
      {
        delete m_starAnnotParams;
        m_starAnnotParams = 0x0;
      }
      if (m_specsScoredAnnotParams)
      {
        delete m_specsScoredAnnotParams;
        m_specsScoredAnnotParams = 0x0;
      }
    }
    if (ownOutput)
    {
      if (m_peptideOutputResults)
      {
        delete m_peptideOutputResults;
        m_peptideOutputResults = 0x0;
      }
    }
  }
  // -------------------------------------------------------------------------
  ExecBase * ExecSvmStatistics::clone(const ParameterList & inputParams) const
  {
    return new ExecSvmStatistics(inputParams);
  }
  // -------------------------------------------------------------------------

  bool compareKeys(vector<unsigned int> i, vector<unsigned int> j)
  {
    return i[2] < j[2];
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::buildStatsVector(vector<vector<unsigned int> > &keyMapping,
                                           vector<float> &outputVector,
                                           int lineIndex)
  {
    //index 0 - 1 for ms stats, 2 for stars stats, 3 for scored
    //index 1 - index within spectrum statistics
    //index 2 - index within key mapping
    for (int i = 0; i < keyMapping.size(); i++)
    {
      unsigned int statsIndex = keyMapping[i][1];
      unsigned int keyIndex = keyMapping[i][2];
      //DEBUG_VAR(lineIndex);
      //DEBUG_VAR(statsIndex);
      //DEBUG_VAR(keyIndex);

      switch (keyMapping[i][0])
      {
      //MS stats
      case 1:
        //DEBUG_MSG("ms spectra" << (*m_specsMsSpectraStats)[lineIndex][statsIndex]);
        outputVector[keyIndex - 1]
            = (*m_specsMsSpectraStats)[lineIndex][statsIndex];
        break;
        //Star stats
      case 2:
        //DEBUG_MSG("star spectra " << (*m_starSpectraStats)[lineIndex][statsIndex]);
        outputVector[keyIndex - 1]
            = (*m_starSpectraStats)[lineIndex][statsIndex];
        break;
        //Specs Scored stats
      case 3:
        //DEBUG_MSG("scored spectra " << (*m_specsScoredSpectraStats)[lineIndex][statsIndex])
        outputVector[keyIndex - 1]
            = (*m_specsScoredSpectraStats)[lineIndex][statsIndex];
        break;
      default:
        WARN_MSG("Unknown statistics file mapping! File key: " << keyMapping[i][0])
        ;
        return false;
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  void ExecSvmStatistics::buildKeyMap(vector<vector<unsigned int> > &keyMapping,
                                      SvmModel * currModel)
  {
    //DEBUG_VAR(currModel->getKeyIndexSize());
    addKeys(keyMapping, m_specsMsSpectraStatsHeader, currModel, 1);
    addKeys(keyMapping, m_starSpectraStatsHeader, currModel, 2);
    addKeys(keyMapping, m_specsScoredSpectraStatsHeader, currModel, 3);
    sort(keyMapping.begin(), keyMapping.end(), compareKeys);
  }

  // -------------------------------------------------------------------------
  void ExecSvmStatistics::addKeys(vector<vector<unsigned int> > &keyMapping,
                                  vector<string> * statisticHeader,
                                  SvmModel * currModel,
                                  int fileEnum)
  {
    for (int i = 0; i < statisticHeader->size(); i++)
    {
      unsigned int keyIndex;
      vector<unsigned int> currKey(3);

      string keyName = (*statisticHeader)[i];
      //DEBUG_VAR(keyName);

      if (currModel->getKeyIndexByName(keyName, keyIndex))
      {
        currKey[0] = fileEnum;
        currKey[1] = i;
        currKey[2] = keyIndex;
        keyMapping.push_back(currKey);
      }
      else
      {
        //        WARN_MSG("Statistic generated which is not used in model: " << (*statisticHeader)[i]);
      }
    }
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::getSvmScore(Spectrum * msSpectra,
                                      Spectrum * scoredSpectra,
                                      Spectrum * starSpectra,
                                      PeptideSpectrumMatch &psm,
                                      MS2ScoringModel * model)
  {
    //make sure stats parameters are loaded
    if (m_msAnnotParams == NULL)
    {
      ownModels = true;
      //statistics params
      m_msAnnotParams = new SpectrumAnnotParameterList();
      m_starAnnotParams = new SpectrumAnnotParameterList();
      m_specsScoredAnnotParams = new SpectrumAnnotParameterList();

      if (!loadAnnotStatisticsParams())
      {
        ERROR_MSG("Unable to load annotation statistics");
        return false;
      }
    }

    if (m_starSpectraStatsHeader == 0x0) // if stats vectors not yet defined.
    {
      m_starSpectraStats = new vector<vector<float> > ();
      m_starSpectraStatsHeader = new vector<string> ();
      m_specsMsSpectraStats = new vector<vector<float> > ();
      m_specsMsSpectraStatsHeader = new vector<string> ();
      m_specsScoredSpectraStats = new vector<vector<float> > ();
      m_specsScoredSpectraStatsHeader = new vector<string> ();
      ownStats = true;
    }

    DEBUG_VAR(m_starSpectraStats->size());
    DEBUG_VAR(m_starSpectraStatsHeader->size());
    DEBUG_VAR(m_specsMsSpectraStats->size());
    DEBUG_VAR(m_specsMsSpectraStatsHeader->size());
    DEBUG_VAR(m_specsScoredSpectraStats->size());
    DEBUG_VAR(m_specsScoredSpectraStatsHeader->size());

    //clear statistics
    m_starSpectraStats->resize(0);
    m_starSpectraStatsHeader->resize(0);
    m_specsMsSpectraStats->resize(0);
    m_specsMsSpectraStatsHeader->resize(0);
    m_specsScoredSpectraStats->resize(0);
    m_specsScoredSpectraStatsHeader->resize(0);

    ParameterList statsParams;
    //set generic values for all statistics generation
    statsParams.addIfExists(m_params, "INPUT_RESULTS_TYPE");
    statsParams.addIfExists(m_params, "TOLERANCE_PEAK");

    SpecSet currSpecSet;
    currSpecSet.resize(1);

    PeptideSpectrumMatchSet currPSMSet;
    currPSMSet.resize(1);

    psmPtr currPsm(new PeptideSpectrumMatch);
    *currPsm = psm;
    currPSMSet[0] = currPsm;
    currSpecSet[0] = *msSpectra;
    currPSMSet[0]->m_spectrum = &(currSpecSet[0]);

    //generate statistics
    ExecStatistics msStatsModule(statsParams,
                                 &currSpecSet,
                                 model,
                                 m_msAnnotParams,
                                 &currPSMSet,
                                 m_specsMsSpectraStats,
                                 m_specsMsSpectraStatsHeader,
                                 0x0,
                                 0x0);

    DEBUG_MSG("Executing statistics: ");

    if (!msStatsModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    currSpecSet[0] = *scoredSpectra;
    currPSMSet[0]->m_spectrum = &(currSpecSet[0]);

    if (m_params.exists("STARS_SCORED_PRM_OFFSET"))
    {
      statsParams.setValue("PRM_OFFSET",
                           m_params.getValue("STARS_SCORED_PRM_OFFSET"));
    }

    if (m_params.exists("STARS_SCORED_SRM_OFFSET"))
    {
      statsParams.setValue("SRM_OFFSET",
                           m_params.getValue("STARS_SCORED_SRM_OFFSET"));
    }

    ExecStatistics scoredStatsModule(statsParams,
                                     &currSpecSet,
                                     model,
                                     m_specsScoredAnnotParams,
                                     &currPSMSet,
                                     m_specsScoredSpectraStats,
                                     m_specsScoredSpectraStatsHeader,
                                     0x0,
                                     0x0);

    if (!scoredStatsModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    currSpecSet[0] = *starSpectra;
    currPSMSet[0]->m_spectrum = &(currSpecSet[0]);

    //set parameters for stars stats generation
    statsParams.setValue("STATISTICS_CONFIG",
                         m_params.getValue("STARS_STATISTICS_CONFIG"));

    if (m_params.exists("STARS_SCORED_SRM_OFFSET"))
    {
      statsParams.setValue("SRM_OFFSET",
                           m_params.getValue("STARS_SCORED_SRM_OFFSET"));
    }

    ExecStatistics starModule(statsParams,
                              &currSpecSet,
                              model,
                              m_starAnnotParams,
                              &currPSMSet,
                              m_starSpectraStats,
                              m_starSpectraStatsHeader,
                              0x0,
                              0x0);

    if (!starModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    PeptideSpectrumMatchSet * resultPtr = m_peptideResults;
    PeptideSpectrumMatchSet * outputPtr = m_peptideOutputResults;
    //temporarily set currPSMSet to current psm
    m_peptideResults = &currPSMSet;
    PeptideSpectrumMatchSet outputTemp;
    m_peptideOutputResults = &outputTemp;

    if (!invoke())
    {
      DEBUG_MSG("Unable to generate SVM results!");
      return false;
    }

    m_peptideResults = resultPtr;
    m_peptideOutputResults = outputPtr;

    psm.m_score = outputTemp[0]->m_score;
    if (psm.m_score != psm.m_score)
    {
      for (int i = 0; i < m_starSpectraStats->size(); i++)
      {
        for (int j = 0; j < (*m_starSpectraStats)[i].size(); j++)
        {
          DEBUG_VAR((*m_starSpectraStats)[i][j]);
        }
      }
      for (int i = 0; i < m_starSpectraStatsHeader->size(); i++)
      {
        DEBUG_VAR((*m_starSpectraStatsHeader)[i]);
      }
      DEBUG_VAR(m_specsMsSpectraStats->size());
      for (int i = 0; i < m_specsMsSpectraStats->size(); i++)
      {
        DEBUG_VAR((*m_specsMsSpectraStats)[i].size());
        for (int j = 0; j < (*m_specsMsSpectraStats)[i].size(); j++)
        {
          DEBUG_VAR((*m_specsMsSpectraStats)[i][j]);
        }
      }
      DEBUG_VAR(m_specsMsSpectraStatsHeader->size());
      for (int i = 0; i < m_specsMsSpectraStatsHeader->size(); i++)
      {
        DEBUG_VAR((*m_specsMsSpectraStatsHeader)[i]);
      }
      for (int i = 0; i < m_specsScoredSpectraStats->size(); i++)
      {
        for (int j = 0; j < (*m_specsScoredSpectraStats)[i].size(); j++)
        {
          DEBUG_VAR((*m_specsScoredSpectraStats)[i][j]);
        }
      }
      for (int i = 0; i < m_specsScoredSpectraStatsHeader->size(); i++)
      {
        DEBUG_VAR((*m_specsScoredSpectraStatsHeader)[i]);
      }
      DEBUG_VAR(msSpectra->parentMZ);
      DEBUG_VAR(msSpectra->parentCharge);
      DEBUG_VAR(msSpectra->parentMass);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::invokeWithStatistics(SpecSet * specsMsSpectra,
                                               SpecSet * specsScoredSpectra,
                                               SpecSet * starSpectra,
                                               MS2ScoringModel * model)
  {
    if (m_starSpectraStatsHeader == 0x0) // if stats vectors not yet defined.
    {
      m_starSpectraStats = new vector<vector<float> > ();
      m_starSpectraStatsHeader = new vector<string> ();
      m_specsMsSpectraStats = new vector<vector<float> > ();
      m_specsMsSpectraStatsHeader = new vector<string> ();
      m_specsScoredSpectraStats = new vector<vector<float> > ();
      m_specsScoredSpectraStatsHeader = new vector<string> ();
      ownStats = true;
    }

    if (m_msAnnotParams == NULL)
    {
      ownModels = true;
      //statistics params
      m_msAnnotParams = new SpectrumAnnotParameterList();
      m_starAnnotParams = new SpectrumAnnotParameterList();
      m_specsScoredAnnotParams = new SpectrumAnnotParameterList();

      if (!loadAnnotStatisticsParams())
      {
        ERROR_MSG("Unable to load annotation statistics");
        return false;
      }
    }

    DEBUG_TRACE;
    DEBUG_VAR(m_starAnnotParams->size());
    DEBUG_VAR(m_msAnnotParams->size());
    DEBUG_VAR(m_specsScoredAnnotParams->size());

    ParameterList statsParams;
    //set generic values for all statistics generation
    statsParams.addIfExists(m_params, "INPUT_RESULTS_TYPE");
    statsParams.addIfExists(m_params, "TOLERANCE_PEAK");

    PeptideSpectrumMatchSet msPeptideResults = *m_peptideResults;
    msPeptideResults.addSpectra(specsMsSpectra, false);
    DEBUG_VAR(msPeptideResults.size());

    //generate statistics
    ExecStatistics msStatsModule(statsParams,
                                 specsMsSpectra,
                                 model,
                                 m_msAnnotParams,
                                 &msPeptideResults,
                                 m_specsMsSpectraStats,
                                 m_specsMsSpectraStatsHeader,
                                 0x0,
                                 0x0);

    DEBUG_MSG("Executing statistics: ");

    if (!msStatsModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    if (m_params.exists("STARS_SCORED_PRM_OFFSET"))
    {
      statsParams.setValue("PRM_OFFSET",
                           m_params.getValue("STARS_SCORED_PRM_OFFSET"));
    }

    if (m_params.exists("STARS_SCORED_SRM_OFFSET"))
    {
      statsParams.setValue("SRM_OFFSET",
                           m_params.getValue("STARS_SCORED_SRM_OFFSET"));
    }

    PeptideSpectrumMatchSet scoredPeptideResults = *m_peptideResults;

    scoredPeptideResults.addSpectra(specsScoredSpectra, false);

    ExecStatistics scoredStatsModule(statsParams,
                                     specsScoredSpectra,
                                     model,
                                     m_specsScoredAnnotParams,
                                     &scoredPeptideResults,
                                     m_specsScoredSpectraStats,
                                     m_specsScoredSpectraStatsHeader,
                                     0x0,
                                     0x0);

    if (!scoredStatsModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    //set parameters for stars stats generation
    statsParams.setValue("STATISTICS_CONFIG",
                         m_params.getValue("STARS_STATISTICS_CONFIG"));

    if (m_params.exists("STARS_SCORED_SRM_OFFSET"))
    {
      statsParams.setValue("SRM_OFFSET",
                           m_params.getValue("STARS_SCORED_SRM_OFFSET"));
    }

    PeptideSpectrumMatchSet starPeptideResults = *m_peptideResults;
    starPeptideResults.addSpectra(starSpectra, false);

    ExecStatistics starModule(statsParams,
                              starSpectra,
                              model,
                              m_starAnnotParams,
                              &starPeptideResults,
                              m_starSpectraStats,
                              m_starSpectraStatsHeader,
                              0x0,
                              0x0);

    if (!starModule.invoke())
    {
      DEBUG_MSG("Unable to generate statistics!");
      return false;
    }

    DEBUG_TRACE;

    if (!invoke())
    {
      DEBUG_MSG("Unable to generate SVM results!");
      return false;
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::invoke(void)
  {
    DEBUG_VAR(m_starSpectraStats->size());
    DEBUG_VAR(m_specsMsSpectraStats->size());
    DEBUG_VAR(m_specsScoredSpectraStats->size());

    //check to make sure input stats are the same length
    if (m_starSpectraStats->size() != m_specsMsSpectraStats->size()
        || m_specsScoredSpectraStats->size() != m_specsMsSpectraStats->size())
    {
      ERROR_MSG("Input star,ms or specs_scored statistic vectors not the same length!");
      return false;
    }

    DEBUG_TRACE;
    if (m_peptideOutputResults == NULL)
    {
      ownOutput = true;
      m_peptideOutputResults = new PeptideSpectrumMatchSet;
    }

    //Load in model files
    if (m_svmScaleParamsCharge1 == NULL)
    {
      ownModels = true;
      m_svmScaleParamsCharge1 = new SvmScaleParameterList();
      m_svmScaleParamsCharge2 = new SvmScaleParameterList();
      m_svmScaleParamsCharge3 = new SvmScaleParameterList();

      m_svmModelCharge1 = new SvmModel();
      m_svmModelCharge2 = new SvmModel();
      m_svmModelCharge3 = new SvmModel();

      if (!loadModels())
      {
        ERROR_MSG("Unable to load model files");
        return false;
      }
    }
    DEBUG_VAR(m_svmScaleParamsCharge1->size());
    DEBUG_VAR(m_svmModelCharge1->size());
    DEBUG_VAR(m_svmModelCharge1->getKeyIndexSize());
    DEBUG_VAR(m_svmScaleParamsCharge2->size());
    DEBUG_VAR(m_svmModelCharge2->size());
    DEBUG_VAR(m_svmModelCharge2->getKeyIndexSize());
    DEBUG_VAR(m_svmScaleParamsCharge3->size());
    DEBUG_VAR(m_svmModelCharge3->getKeyIndexSize());
    DEBUG_VAR(m_svmModelCharge3->size());

    //maps key names to which file they belong to and
    //which index in the file they are contained in.
    //index 0 - 1 for ms stats, 2 for stars stats, 3 for scored
    //index 1 - index within spectrum statistics
    //index 2 - index within key mapping
    vector<vector<unsigned int> > keyMapping;

    m_peptideOutputResults->resize(m_starSpectraStats->size());

    buildKeyMap(keyMapping, m_svmModelCharge1);

    for (int i = 0; i < m_specsMsSpectraStats->size(); i++)
    {
      int lastKey = keyMapping.back()[2];
      //DEBUG_VAR(lastKey);

      vector<float> statsVector(lastKey); //Last key index should be largest.
      vector<float> scaledVector(lastKey);

      psmPtr resultLine = (*m_peptideResults)[i];
      psmPtr outputLine(new PeptideSpectrumMatch);
      *outputLine = *resultLine;
      DEBUG_VAR(resultLine->m_annotation);

      double svmScore; //this is only a double because it is easier to compare to libsvm stats.

      switch (resultLine->m_charge)
      {
      case 1:
        buildStatsVector(keyMapping, statsVector, i);
        m_svmModelCharge1->scaleVector(statsVector, scaledVector);
        svmScore = m_svmModelCharge1->getSvmScore(scaledVector);
        for (int j = 0; j < statsVector.size(); j++)
        {
          DEBUG_VAR(statsVector[j]);
        }
        break;
      case 2:
        buildStatsVector(keyMapping, statsVector, i);
        m_svmModelCharge2->scaleVector(statsVector, scaledVector);
        svmScore = m_svmModelCharge2->getSvmScore(scaledVector);
        for (int j = 0; j < scaledVector.size(); j++)
        {
          DEBUG_VAR(scaledVector[j]);
        }
        break;
      default:
        buildStatsVector(keyMapping, statsVector, i);
        m_svmModelCharge3->scaleVector(statsVector, scaledVector);
        svmScore = m_svmModelCharge3->getSvmScore(scaledVector);
        for (int j = 0; j < statsVector.size(); j++)
        {
          DEBUG_VAR(statsVector[j]);
        }
        break;
      }

      DEBUG_VAR(svmScore);
      outputLine->m_score = svmScore;
      m_peptideOutputResults->m_psmSet[i] = outputLine;
    }
    DEBUG_TRACE;
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::loadModels(void)
  {
    if (m_params.exists("SVM_SCALE_PARAMS_CHARGE1"))
    {
      if (!m_svmScaleParamsCharge1->loadSvmKeyRange(m_params.getValue("SVM_SCALE_PARAMS_CHARGE1").c_str()))
      {
        DEBUG_MSG("Unable to load key file! " << m_params.getValue("SVM_SCALE_PARAMS_CHARGE1"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_SCALE_PARAMS_CHARGE1 parameter!");
      return false;
    }

    m_svmModelCharge1->setKeyRanges(m_svmScaleParamsCharge1);

    if (m_params.exists("SVM_SCALE_PARAMS_CHARGE2"))
    {
      if (!m_svmScaleParamsCharge2->loadSvmKeyRange(m_params.getValue("SVM_SCALE_PARAMS_CHARGE2").c_str()))
      {
        DEBUG_MSG("Unable to load key file! " << m_params.getValue("SVM_SCALE_PARAMS_CHARGE2"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_SCALE_PARAMS_CHARGE2 parameter!");
      return false;
    }

    m_svmModelCharge2->setKeyRanges(m_svmScaleParamsCharge2);

    if (m_params.exists("SVM_SCALE_PARAMS_CHARGE3"))
    {
      if (!m_svmScaleParamsCharge3->loadSvmKeyRange(m_params.getValue("SVM_SCALE_PARAMS_CHARGE3").c_str()))
      {
        DEBUG_MSG("Unable to load key file! " << m_params.getValue("SVM_SCALE_PARAMS_CHARGE3"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_SCALE_PARAMS_CHARGE3 parameter!");
      return false;
    }

    m_svmModelCharge3->setKeyRanges(m_svmScaleParamsCharge3);

    if (m_params.exists("SVM_MODEL_CHARGE1"))
    {
      if (!m_svmModelCharge1->loadSvmModel(m_params.getValue("SVM_MODEL_CHARGE1").c_str()))
      {
        DEBUG_MSG("Unable to load SVM Model: " << m_params.getValue("SVM_MODEL_CHARGE1"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_MODEL_CHARGE1 parameter!");
      return false;
    }

    if (m_params.exists("SVM_MODEL_CHARGE2"))
    {
      if (!m_svmModelCharge2->loadSvmModel(m_params.getValue("SVM_MODEL_CHARGE2").c_str()))
      {
        DEBUG_MSG("Unable to load SVM Model: " << m_params.getValue("SVM_MODEL_CHARGE2"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_MODEL_CHARGE2 parameter!");
      return false;
    }

    if (m_params.exists("SVM_MODEL_CHARGE3"))
    {
      if (!m_svmModelCharge3->loadSvmModel(m_params.getValue("SVM_MODEL_CHARGE3").c_str()))
      {
        DEBUG_MSG("Unable to load SVM Model: " << m_params.getValue("SVM_MODEL_CHARGE3"));
        return false;
      }
    }
    else
    {
      ERROR_MSG("Missing SVM_MODEL_CHARGE3 parameter!");
      return false;
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::loadAnnotStatisticsParams(void)
  {
    if (!m_msAnnotParams->loadSpectrumAnnotFile(m_params.getValue("SPECS_MS_STATISTICS_CONFIG").c_str()))
    {
      ERROR_MSG("Unable to read " << m_params.getValue("SPECS_MS_STATISTICS_CONFIG"));
      return false;
    }

    if (!m_specsScoredAnnotParams->loadSpectrumAnnotFile(m_params.getValue("SPECS_SCORED_STATISTICS_CONFIG").c_str()))
    {
      ERROR_MSG("Unable to read " << m_params.getValue("SPECS_SCORED_STATISTICS_CONFIG"));
      return false;
    }

    if (!m_starAnnotParams->loadSpectrumAnnotFile(m_params.getValue("STARS_STATISTICS_CONFIG").c_str()))
    {
      ERROR_MSG("Unable to read " << m_params.getValue("STARS_STATISTICS_CONFIG"));
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::loadInputData(void)
  {
    //Load in statistics if they haven't been passed in by
    //another function
    if (m_starSpectraStats == 0x0)
    {
      ownInput = true;
      m_starSpectraStats = new vector<vector<float> > ();
      m_starSpectraStatsHeader = new vector<string> ();
      m_specsMsSpectraStats = new vector<vector<float> > ();
      m_specsMsSpectraStatsHeader = new vector<string> ();
      m_specsScoredSpectraStats = new vector<vector<float> > ();
      m_specsScoredSpectraStatsHeader = new vector<string> ();
      m_svmScaleParamsCharge1 = new SvmScaleParameterList();
      m_svmScaleParamsCharge2 = new SvmScaleParameterList();
      m_svmScaleParamsCharge3 = new SvmScaleParameterList();

      m_svmModelCharge1 = new SvmModel();
      m_svmModelCharge2 = new SvmModel();
      m_svmModelCharge3 = new SvmModel();

      m_msAnnotParams = new SpectrumAnnotParameterList();
      m_starAnnotParams = new SpectrumAnnotParameterList();
      m_specsScoredAnnotParams = new SpectrumAnnotParameterList();

      ownStats = true;

      m_peptideResults = new PeptideSpectrumMatchSet;
    }

    //Load models
    DEBUG_MSG("Load models");

    if (!loadModels())
    {
      ERROR_MSG("Unable to read model files!");
      return false;
    }

    //Load in inspect / specnets results
    if (m_params.exists("INPUT_RESULTS"))
    {
      DEBUG_MSG("Input results: " << m_params.getValue("INPUT_RESULTS"));
      if (m_params.getValue("INPUT_RESULTS_TYPE").compare("inspect") == 0)
      {
        if (!m_peptideResults->loadInspectResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                      m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
      else
      {
        if (!m_peptideResults->loadSpecnetsResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                       m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
    }

    DEBUG_MSG("InputLoaded")

    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    vector<int> requiredHeaderIndex(requiredHeader.size());

    if (m_params.exists("INPUT_SPECTRA_STATS_STARS"))
    {
      if (!DelimitedTextReader::loadDelimitedFile(m_params.getValue("INPUT_SPECTRA_STATS_STARS").c_str(),
                                                  ",",
                                                  "#",
                                                  *m_starSpectraStatsHeader,
                                                  *m_starSpectraStats,
                                                  requiredHeader,
                                                  requiredHeaderIndex))
      {
        DEBUG_MSG("Unable to load stars stats! " << m_params.getValue("INPUT_SPECTRA_STATS_STARS"));
        return false;
      }
    }

    DEBUG_MSG("LOADED STARS");

    if (m_params.exists("INPUT_SPECTRA_STATS_SPECS_MS"))
    {
      if (!DelimitedTextReader::loadDelimitedFile(m_params.getValue("INPUT_SPECTRA_STATS_SPECS_MS").c_str(),
                                                  ",",
                                                  "#",
                                                  *m_specsMsSpectraStatsHeader,
                                                  *m_specsMsSpectraStats,
                                                  requiredHeader,
                                                  requiredHeaderIndex))
      {
        DEBUG_MSG("Unable to load ms stats! " << m_params.getValue("INPUT_SPECTRA_STATS_SPECS_MS"));
        return false;
      }
    }

    DEBUG_MSG("LOADED MS");

    if (m_params.exists("INPUT_SPECTRA_STATS_SPECS_SCORED"))
    {
      if (!DelimitedTextReader::loadDelimitedFile(m_params.getValue("INPUT_SPECTRA_STATS_SPECS_SCORED").c_str(),
                                                  ",",
                                                  "#",
                                                  *m_specsScoredSpectraStatsHeader,
                                                  *m_specsScoredSpectraStats,
                                                  requiredHeader,
                                                  requiredHeaderIndex))
      {
        DEBUG_MSG("Unable to load scored stats! " << m_params.getValue("INPUT_SPECTRA_STATS_SPECS_SCORED"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPECTRA_SVM_SCORES"))
    {
      if (!m_peptideOutputResults->saveToFile(m_params.getValue("OUTPUT_SPECTRA_SVM_SCORES").c_str(),
                                              true))
      {
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::loadOutputData(void)
  {
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecSvmStatistics::split(int numSplit)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::merge(void)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::validateStatisticsParams(std::string & error)
  {
    m_isValid = false;
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE1");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE2");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE3");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE1");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE2");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE3");
    VALIDATE_PARAM_EXIST("SPECS_MS_STATISTICS_CONFIG");
    VALIDATE_PARAM_EXIST("SPECS_SCORED_STATISTICS_CONFIG");
    VALIDATE_PARAM_EXIST("STARS_STATISTICS_CONFIG");
    m_isValid = true;
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecSvmStatistics::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE1");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE2");
    VALIDATE_PARAM_EXIST("SVM_SCALE_PARAMS_CHARGE3");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE1");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE2");
    VALIDATE_PARAM_EXIST("SVM_MODEL_CHARGE3");
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA_STATS_SPECS_SCORED");
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA_STATS_SPECS_MS");
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA_STATS_STARS");

    m_isValid = true;
    return true;
  }

}
