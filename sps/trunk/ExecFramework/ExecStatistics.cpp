// Header Include
#include "ExecStatistics.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "utils.h"
#include "ExecMergeConvert.h"

// System Includes
#include <stdio.h>

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(void) :
    m_spectra(0x0), m_model(0x0), m_statsParams(0x0),
        m_peptideResults(0x0), m_spectraStats(0x0), m_spectraStatsHeader(0x0),
        m_DaHistogram(0x0), m_PPMHistogram(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(const ParameterList & inputParams) :
    ExecBase(inputParams), m_spectra(0x0), m_model(0x0), m_statsParams(0x0),
        m_peptideResults(0x0), m_spectraStats(0x0), m_spectraStatsHeader(0x0),
        m_DaHistogram(0x0), m_PPMHistogram(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(const ParameterList & inputParams,
                                 SpecSet * spectra,
                                 MS2ScoringModel * model,
                                 SpectrumAnnotParameterList * statsParams,
                                 PeptideSpectrumMatchSet * peptideResults,
                                 vector<vector<float> > * spectraStats,
                                 vector<string> * spectraHeader,
                                 OutputTable *DaHistogram,
                                 OutputTable *PPMHistogram) :
    ExecBase(inputParams), m_spectra(spectra), m_model(model),
        m_statsParams(statsParams), m_peptideResults(peptideResults),
        m_spectraStats(spectraStats), m_spectraStatsHeader(spectraHeader),
        m_DaHistogram(DaHistogram), m_PPMHistogram(PPMHistogram),
        ownInput(false), ownOutput(false)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::~ExecStatistics(void)
  {
    if (ownInput)
    {
      if (m_spectra)
      {
        delete m_spectra;
        m_spectra = 0x0;
      }
      if (m_model)
      {
        delete m_model;
        m_model = 0x0;
      }
      if (m_statsParams)
      {
        delete m_statsParams;
        m_statsParams = 0x0;
      }
      if (m_peptideResults)
      {
        delete m_peptideResults;
        m_peptideResults = 0x0;
      }
    }

    if (ownOutput)
    {
      if (m_spectraStats)
      {
        delete m_spectraStats;
        m_spectraStats = 0x0;
      }
      if (m_spectraStatsHeader)
      {
        delete m_spectraStatsHeader;
        m_spectraStatsHeader = 0x0;
      }
      if (m_DaHistogram)
      {
        delete m_DaHistogram;
      }
      if (m_PPMHistogram)
      {
        delete m_PPMHistogram;
      }
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecStatistics::clone(const ParameterList & inputParams) const
  {
    return new ExecStatistics(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::invoke(void)
  {
    if (m_spectraStats == 0x0)
    {
      ownOutput = true;
      m_spectraStats = new vector<vector<float> > ();
      m_spectraStatsHeader = new vector<string> ();
    }
    if (m_DaHistogram == 0x0 and ownOutput)
    {
      m_DaHistogram = new OutputTable;
    }
    if (m_PPMHistogram == 0x0 and ownOutput)
    {
      m_PPMHistogram = new OutputTable;
    }

    //Make sure required parameters are set:
    if (!m_spectra or m_spectra->size() == 0)
    {
      DEBUG_MSG("Missing input spectra!");
      return false;
    }

    if (!m_peptideResults or !m_statsParams)
    {
      DEBUG_MSG("Missing stats configuration or annotations!")
      return false;
    }

    if (m_statsParams->m_params.size() == 0)
    {
      DEBUG_MSG("No defined stats parameters!");
      return false;
    }

    float prmOffset = m_params.getValueFloat("PRM_OFFSET");
    float srmOffset = m_params.getValueFloat("SRM_OFFSET");
    float peakTol = 0.05;
    if (m_params.exists("TOLERANCE_PEAK"))
    {
      peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
    }

    SpecSet processedSpecs;
    vector<pair<int, int> > loadedIndices;
    ExecMergeConvert loader(m_params,
                            &loadedIndices,
                            m_spectra,
                            &processedSpecs);

    if (!loader.invoke())
    { // Let MergeConvert do all the pre-processing
      return false;
    }

    m_spectra->operator =(processedSpecs);

    vector<psmPtr>::iterator resultIterator;

    SpectrumAnnotStatistics stats;

    AAJumps jumps(1);
    string ionTypes = "all";

    for (resultIterator = m_peptideResults->m_psmSet.begin(); resultIterator
        != m_peptideResults->m_psmSet.end(); resultIterator++)
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
                    *m_model,
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
                stats.percentExplainedIntensity(*m_peptideResults, *currParam);
            currStatistics[i + 1] = explainedIntensity;
          }
        }
        else if (currParam->statistic.compare("fragment mass error da") == 0
            and m_DaHistogram)
        {
          float maxDaError;
          if (m_params.exists("MAX_DA_ERROR"))
          {
            maxDaError = m_params.getValueFloat("MAX_DA_ERROR");
          }
          else
          {
            maxDaError = 0.06;
            for (unsigned int specIdx = 0; specIdx < m_spectra->size(); specIdx++)
            {
              for (unsigned int peakIdx = 0; peakIdx
                  < (*m_spectra)[specIdx].size(); peakIdx++)
              {
                maxDaError = max(maxDaError,
                                 (*m_spectra)[specIdx].getTolerance(peakIdx));
              }
            }
          }
          vector<TwoValues<float> > bins;
          stats.fragmentMassErrorDaHistogram(*m_peptideResults,
                                             *currParam,
                                             maxDaError,
                                             bins);
          m_DaHistogram->loadHistogram(bins);
        }
        else if (currParam->statistic.compare("fragment mass error ppm") == 0
            and m_PPMHistogram)
        {
          float maxPPMError;
          if (m_params.exists("MAX_DA_ERROR"))
          {
            maxPPMError = m_params.getValueFloat("MAX_DA_ERROR");
          }
          else
          {
            maxPPMError = 0.06;
            for (unsigned int specIdx = 0; specIdx < m_spectra->size(); specIdx++)
            {
              for (unsigned int peakIdx = 0; peakIdx
                  < (*m_spectra)[specIdx].size(); peakIdx++)
              {
                maxPPMError = max(maxPPMError, (float)(1000000.0
                    * (*m_spectra)[specIdx].getTolerance(peakIdx)
                    / (*m_spectra)[specIdx][peakIdx][0]));
              }
            }
          }
          vector<TwoValues<float> > bins;
          stats.fragmentMassErrorPPMHistogram(*m_peptideResults,
                                              *currParam,
                                              maxPPMError,
                                              bins);
          m_PPMHistogram->loadHistogram(bins);
        }
        else if (currParam->statistic.compare("%explained peaks") == 0)
        {
          float explainedPeaks = stats.percentExplainedPeaks(*m_peptideResults,
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
            float observedIons = stats.percentObservedIons(*m_peptideResults,
                                                           *currParam);
            currStatistics[i + 1] = observedIons;
          }
        }
        else if (currParam->statistic.compare("%observed breaks") == 0)
        {
          float observedBreaks = stats.observedBreaks(*m_peptideResults,
                                                      *currParam);
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
    for (resultIterator = m_peptideResults->m_psmSet.begin(); resultIterator
        != m_peptideResults->m_psmSet.end(); resultIterator++)
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
      if (!addedHeaders && resultIterator == m_peptideResults->m_psmSet.begin())
      {
        m_spectraStatsHeader->push_back("Scan#");
      }

      for (int i = 0; i < m_statsParams->m_params.size(); i++)
      {
        SpectrumAnnotParameter * currParam = &(m_statsParams->m_params[i]);

        //add to header
        if (!addedHeaders && resultIterator
            == m_peptideResults->m_psmSet.begin())
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
          parentMassErrorPpm = abs(parentMassErrorPpm / (1 - peakTol)); //output error adjusted by tolerance
          currStatistics[i + 1] = parentMassErrorPpm;
        }
        else if (currParam->statistic.compare("parent mass error da") == 0)
        {
          float parentMassErrorDa = stats.parentMassErrorDa(psm, psm.m_charge);
          parentMassErrorDa = abs(parentMassErrorDa / (1 - peakTol)); //output error adjusted by tolerance
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
  // -------------------------------------------------------------------------
  bool loadSpecSet(string & path, string & fileExt, SpecSet & currSet)
  {
    //load in file based on file ext.
    if (fileExt.compare("pklbin") == 0)
    {
      std::stringstream ss;
      ss << path.substr(0, path.find_last_of('.')) << ".bin";
      string binFile = ss.str();

      DEBUG_VAR(binFile);

      if (!currSet.loadPklBin(path.c_str()))
      {
        return false;
      }

      DEBUG_VAR(currSet.size());

      vector<vector<int> > data;
      if (fileExists(binFile))
      {
        if (!Load_binArray(binFile.c_str(), data))
        {
          return false;
        }

        DEBUG_VAR(data.size());

        for (int i = 0; i < data.size(); i++)
        {
          unsigned int scanNum = (unsigned int)data[i][0];
          short msLevel = (short)data[i][1];

          if (msLevel != 1 && msLevel != 2 && msLevel != 3)
          {
            DEBUG_VAR(scanNum);
            DEBUG_VAR(currSet[i].scan);
            DEBUG_VAR(msLevel);
          }

          currSet[i].msLevel = msLevel;
          currSet[i].scan = scanNum;
        }
      }
    }
    else if (fileExt.compare("mgf") == 0)
    {
      if (!currSet.LoadSpecSet_mgf(path.c_str()))
      {
        return false;
      }
      //reset scan numbers since we ignore values from mgf.
      //DELETE ME IN THE FUTURE! I'M HERE TO DEAL WITH SCAN
      //NUMBERS THAT DON'T MATCH INSPECT!
      for (int i = 0; i < currSet.size(); i++)
      {
        currSet[i].scan = i + 1;
        currSet[i].msLevel = 2;
      }
      DEBUG_VAR(currSet.size());
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool loadInspectSpectrum(const char * specSetListFile,
                           PeptideSpectrumMatchSet &psms,
                           SpecSet &outputResults,
                           bool setScanNums)
  {
    vector<vector<string> > spectraFileList;
    map<string, unsigned int> spectraFileListHeader;
    vector<string> requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(specSetListFile,
                                                "\t",
                                                "#",
                                                spectraFileListHeader,
                                                spectraFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load spectrum files!" << specSetListFile);
      return false;
    }

    vector<int> outputResultIndex(psms.m_psmSet.size(), -1);
    outputResults.resize(psms.m_psmSet.size());

    for (int i = 0; i < spectraFileList.size(); i++)
    {
      //get file path
      string path = spectraFileList[i][requiredHeaderIndex[0]];
      DEBUG_VAR(path);

      string filename;
      string fileExt;
      extractFileName(path, filename, false); //pull out stripped file name
      extractFileExt(path, fileExt); // pull out file extension
      DEBUG_VAR(fileExt);
      DEBUG_VAR(filename);

      SpecSet currSet;

      if (!loadSpecSet(path, fileExt, currSet))
      {
        ERROR_MSG("Unable to load file! " << path);
        return false;
      }

      if (setScanNums)
      {
        DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
        for (int i = 0; i < currSet.size(); i++)
        {
          currSet[i].scan = (unsigned int)(i + 1);
        }
      }

      for (int j = 0; j < psms.m_psmSet.size(); j++)
      {
        psmPtr currPsm = psms.m_psmSet[j];
        string psmFilename;
        extractFileName(currPsm->m_spectrumFile, psmFilename, false);

        if (psmFilename.compare(filename) == 0)
        {
          Spectrum * currScan = currSet.getScan(currPsm->m_scanNum);
          if (currScan != NULL)
          {
            outputResults[j] = *currScan;
            outputResultIndex[j] = outputResults.size() - 1;
          }
          else
          {
            DEBUG_VAR(currPsm->m_scanNum);
            DEBUG_TRACE;
          }
        }
      }
    }

    for (int i = 0; i < outputResultIndex.size(); i++)
    {
      if (outputResultIndex[i] > -1)
      {
        psms.m_psmSet[i]->m_spectrum = &(outputResults[i]);
      }
    }

    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecStatistics::loadInputData(void)
  {
    //Load in spectra if they haven't been passed in by
    //another function
    if (m_spectra == 0x0)
    {
      ownInput = true;
      m_spectra = new SpecSet;
      m_statsParams = new SpectrumAnnotParameterList;
      m_peptideResults = new PeptideSpectrumMatchSet;
      m_model = new MS2ScoringModel;
    }

    //Load in inspect / specnets results
    if (m_params.exists("INPUT_RESULTS"))
    {
      DEBUG_MSG("Input results: " << m_params.getValue("INPUT_RESULTS"));
      string resultsType(m_params.getValue("INPUT_RESULTS_TYPE"));
      std::transform(resultsType.begin(),
                     resultsType.end(),
                     resultsType.begin(),
                     ::tolower);
      if (resultsType.compare("inspect") == 0)
      {
        if (!m_peptideResults->loadInspectResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                      m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
      else if (resultsType.compare("msgfdb") == 0)
      {
        if (!m_peptideResults->loadMSGFDBResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
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

    if (m_peptideResults->size() == 0)
    {
      DEBUG_MSG("No peptide annotation results!");
      return false;
    }

    DEBUG_VAR(m_peptideResults->size());

    if (m_params.exists("INPUT_SPECS"))
    {
      if (!m_spectra->LoadSpecSet_pkl(m_params.getValue("INPUT_SPECS").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS"));
        return false;
      }

    }
    else if (m_params.exists("INPUT_SPECS_PKLBIN"))
    {
      if (m_spectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
          <= 0)
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS_PKLBIN"));
        return false;
      }
      if (m_params.getValueBool("SET_SCAN_NUMS"))
      {
        DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
        for (int i = 0; i < m_spectra->size(); i++)
        {
          (*m_spectra)[i].scan = (unsigned int)(i + 1);
        }
      }
      m_peptideResults->addSpectra(m_spectra, false);
    }
    else if (m_params.exists("INPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::loadSpecset(m_params.getValue("EXE_DIR"),
                                         m_params.getValue("INPUT_SPECTRA"),
                                         m_spectra))
      {
        return false;
      }
    }
    else if (m_params.exists("INPUT_SPECTRA_FILE_LIST"))
    {
      if (!loadInspectSpectrum(m_params.getValue("INPUT_SPECTRA_FILE_LIST").c_str(),
                               *m_peptideResults,
                               *m_spectra,
                               m_params.getValueBool("SET_SCAN_NUMS")))
      {
        ERROR_MSG("Unable to load spectrum files!" << m_params.getValue("INPUT_SPECTRA_FILE_LIST"));
        return false;
      }
    }

    if (m_spectra->size() == 0)
    {
      DEBUG_MSG("Input spectra size is 0!");
      return false;
    }

    if (m_params.getValueBool("SET_SCAN_NUMS"))
    {
      DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
      for (int i = 0; i < m_spectra->size(); i++)
      {
        (*m_spectra)[i].scan = (unsigned int)(i + 1);
      }
    }
    m_peptideResults->addSpectra(m_spectra, false);

    DEBUG_VAR(m_spectra->size());
    if (m_params.exists("TOLERANCE_PEAK_PPM"))
    {
      float ppmTol = m_params.getValueFloat("TOLERANCE_PEAK_PPM");
      m_spectra->setPeakTolerance(ppmTol, true);
    }
    else if (m_params.exists("TOLERANCE_PEAK"))
    {
      float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
      m_spectra->setPeakTolerance(peakTol, false);
    }

    //Load in statistics generation parameters
    if (m_params.exists("STATISTICS_CONFIG"))
    {
      if (!m_statsParams->loadSpectrumAnnotFile(m_params.getValue("STATISTICS_CONFIG").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("STATISTICS_CONFIG"));
        return false;
      }
    }

    if (m_statsParams->size() == 0)
    {
      DEBUG_MSG("No statistics asked for!");
      return false;
    }

    if (m_params.exists("MS2_SCORING_MODEL"))
    {
      if (!m_model->LoadModel(m_params.getValue("MS2_SCORING_MODEL").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }

    if (m_model->probs.size() == 0)
    {
      DEBUG_MSG("No model parameters!");
      return false;
    }

    DEBUG_VAR(m_statsParams->size());

    //associate results with spectra
    // m_peptideResults->addSpectra(m_spectra);

    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecStatistics::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPECTRA_STATS_BIN"))
    {
      DEBUG_MSG("Outputting statistics...");
      if (!Save_binArray(m_params.getValue("OUTPUT_SPECTRA_STATS_BIN").c_str(),
                         *m_spectraStats))
      {
        ERROR_MSG("Unable to write bin file!");
        return false;
      }
      DEBUG_MSG("done...");
    }

    if (m_params.exists("OUTPUT_DA_ERROR_HISTOGRAM") && m_DaHistogram)
    {
      m_DaHistogram->printToCSV(m_params.getValue("OUTPUT_DA_ERROR_HISTOGRAM").c_str(),
                                ",");
    }

    if (m_params.exists("OUTPUT_PPM_ERROR_HISTOGRAM") && m_PPMHistogram)
    {
      m_PPMHistogram->printToCSV(m_params.getValue("OUTPUT_PPM_ERROR_HISTOGRAM").c_str(),
                                 ",");
    }

    if (m_params.exists("OUTPUT_SPECTRA_STATS"))
    {
      ofstream outputFile;

      DEBUG_MSG("Opening output stats file...");
      outputFile.open(m_params.getValue("OUTPUT_SPECTRA_STATS").c_str(),
                      ios::out | ios::trunc | ios::binary);
      if (outputFile.fail())
      {
        ERROR_MSG("Unable to open stats file! " << m_params.getValue("OUTPUT_SPECTRA_STATS"));
        return false;
      }

      if (outputFile.is_open() && outputFile.good())
      {
        //output header

        outputFile << (*m_spectraStatsHeader)[0];
        for (int i = 1; i < m_spectraStatsHeader->size(); i++)
        {
          outputFile << "," << (*m_spectraStatsHeader)[i];
        }
        outputFile << endl;

        for (int i = 0; i < m_spectraStats->size(); i++)
        {
          outputFile << (*m_spectraStats)[i][0]; //scan number

          for (int j = 1; j < (*m_spectraStats)[i].size(); j++)
          {
            outputFile << "," << (*m_spectraStats)[i][j];
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

  // -------------------------------------------------------------------------
  bool ExecStatistics::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::loadOutputData(void)
  {
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecStatistics::split(int numSplit)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::merge(void)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::validateParams(std::string & error)
  {
    m_isValid = false;

    //VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");

    m_isValid = true;
    return true;
  }
}
