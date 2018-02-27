// Header Includes
#include "ExecPrmScoring.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "ExecMergeConvert.h"
#include "ExecFilterPairs.h"
#include "DeconvSpectrum.h"
#include "PairedSpecSet.h"
#include "UnionFind.h"
#include "prm_alignment.h"
#include "mzrange.h"

// External Includes
#include "Specific.h"

// System Includes
#include <stdlib.h>
#include <fstream>

using namespace std;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(void) :
      ownInput(true), m_inputSpectra(0x0), ownOutput(true),
          m_outputSpectra(0x0), enforceDaTolerance(false)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_inputSpectra = new SpecSet;
    m_outputSpectra = new SpecSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
          ownOutput(true), m_outputSpectra(0x0), enforceDaTolerance(false)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_inputSpectra = new SpecSet;
    m_outputSpectra = new SpecSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams,
                                 SpecSet * inputSpectra,
                                 SpecSet * outputSpectra) :
      ExecBase(inputParams), ownInput(false), m_inputSpectra(inputSpectra),
          ownOutput(false), m_outputSpectra(outputSpectra),
          enforceDaTolerance(false)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams,
                                 SpecSet * inputSpectra) :
      ExecBase(inputParams), ownInput(false), m_inputSpectra(inputSpectra),
          ownOutput(true), m_outputSpectra(0x0), enforceDaTolerance(false)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_outputSpectra = new SpecSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::~ExecPrmScoring(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
    }

    if (ownOutput)
    {
      delete m_outputSpectra;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecPrmScoring::clone(const ParameterList & inputParams) const
  {
    return new ExecPrmScoring(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::invoke(void)
  {
    if (m_inputSpectra == 0)
    {
      ERROR_MSG("No input spectra!");
      return false;
    }

    DEBUG_TRACE;
    if (m_params.getValueBool("DEBUG_PARAMS")) {
      stringstream aux;
      m_params.print(aux);
      DEBUG_MSG(aux.str());
    }
    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)m_inputSpectra->size()) - 1);

    SpecSet inputSpecsSubset;
    SpecSet *inputMS2Use = m_inputSpectra;
    DEBUG_VAR(m_inputSpectra->size());
    if (idxStart > 0 || idxEnd < ((int)m_inputSpectra->size()) - 1)
    {
      inputSpecsSubset.resize(idxEnd - idxStart + 1);
      for (int i = idxStart; i <= idxEnd; i++)
      {
        inputSpecsSubset[i - idxStart] = (*m_inputSpectra)[i];
      }
      inputMS2Use = &inputSpecsSubset;
    }
    DEBUG_VAR(inputMS2Use->size());

    if (inputMS2Use->size() == 0)
    {
      WARN_MSG("No input spectra, terminating invoke without error");
      return true;
    }

    int c13Offset = m_params.getValueInt("PM_OFFSET_C13", 0);

    if (c13Offset != 0)
    {
      double pmOffset = AAJumps::massC_Iso * ((double)c13Offset);
      DEBUG_MSG("Offsetting parent masses by " << pmOffset);
      for (unsigned int i = 0; i < inputMS2Use->size(); i++)
      {
        (*inputMS2Use)[i].setParentMass((*inputMS2Use)[i].parentMass
            + pmOffset);
      }
    }

    const bool guessedPM = m_params.getValueBool("CORRECT_PM", false);
    DEBUG_VAR(guessedPM);
    const bool guessedCharge = m_params.getValueBool("GUESS_CHARGE", false);
    DEBUG_VAR(guessedCharge);

    SpecSet headerInfo;
    headerInfo.resize(inputMS2Use->size());
    map<int, int> scanToIdx;
    for (int i = 0; i < inputMS2Use->size(); i++)
    {
      headerInfo[i].copyNP((*inputMS2Use)[i]);
      headerInfo[i].resize(0);

      if (scanToIdx.count(headerInfo[i].scan) > 0)
      {
        ERROR_MSG("Found duplicate scan " << headerInfo[i].scan <<", all scans must be unique");
        return false;
      }
      scanToIdx[headerInfo[i].scan] = i;
    }

    // if true, just enforce Da tolerance rather than capturing PPM tolerance from input MS/MS spectra
    enforceDaTolerance = m_params.getValueInt("ENFORCE_DA_TOL", 1) > 0;

    SpecSet prmSpecs;

    if (m_params.exists("PEPNOVO_MODEL") && m_params.exists("PEPNOVO_INPUT_MGF")
        && m_params.exists("PEPNOVO_OUTPUT_PRMS"))
    {
      try
      {
        bool res = invokePepNovo(m_params, *inputMS2Use, prmSpecs);
        if (!res)
        {
          throw 20;
        }
      }
      catch (exception& e)
      {
        ERROR_MSG("Exception occurred: " << e.what());
        return false;
      }
      DEBUG_TRACE;
    }
    else
    {
      // extract spectra for each fragmentation mode
      SpecSet CIDms2specs;
      SpecSet HCDms2specs;
      SpecSet ETDms2specs;
      inputMS2Use->swapExtractSpectra(CIDms2specs, Spectrum::FragType_CID);
      inputMS2Use->swapExtractSpectra(HCDms2specs, Spectrum::FragType_HCD);
      inputMS2Use->swapExtractSpectra(ETDms2specs, Spectrum::FragType_ETD);

      DEBUG_VAR(CIDms2specs.size());
      DEBUG_VAR(HCDms2specs.size());
      DEBUG_VAR(ETDms2specs.size());

      int idxUse = 0;
      bool execSuccess = true;

      if (CIDms2specs.size() > 0)
      {

        SpecSet CIDprms(CIDms2specs.size());

        // generate PRM spectra w/ CID model
        bool res = invokePepNovoCID(m_params, CIDms2specs, CIDprms);
        if (!res)
        {
          return false;
        }

        CIDms2specs.clear();

        // append CID PRMs
        prmSpecs.appendSpecSet(CIDprms, false);

      }
      if (ETDms2specs.size() > 0)
      {
        SpecSet ETDprms(ETDms2specs.size());

        // generate PRM spectra with ETD model
        bool res = invokePepNovoETD(m_params, ETDms2specs, ETDprms);
        if (!res)
        {
          return false;
        }

        ETDms2specs.clear();

        prmSpecs.appendSpecSet(ETDprms, false);

      }
      if (HCDms2specs.size() > 0)
      {
        SpecSet HCDprms(HCDms2specs.size());

        // generate PRM spectra with HCD model
        bool res = invokePepNovoHCD(m_params, HCDms2specs, HCDprms);
        DEBUG_TRACE;
        if (!res)
        {
          return false;
        }

        HCDms2specs.clear();

        prmSpecs.appendSpecSet(HCDprms, false);
      }
    }

    DEBUG_VAR(prmSpecs.size());

    m_outputSpectra->resize(prmSpecs.size());

    for (int i = 0; i < prmSpecs.size(); i++)
    {
      int originSpecIdx = scanToIdx[prmSpecs[i].scan];
      (*m_outputSpectra)[originSpecIdx] = prmSpecs[i];
    }

    DEBUG_TRACE;

    int numEmpty = 0;

    for (int i = 0; i < m_outputSpectra->size(); i++)
    {
      const int curCharge =
          (guessedCharge) ? (*m_outputSpectra)[i].parentCharge :
              headerInfo[i].parentCharge;
      const float curMass =
          (guessedPM) ? (*m_outputSpectra)[i].parentMass :
              headerInfo[i].parentMass;

      (*m_outputSpectra)[i].copyNP(headerInfo[i]);

      if (guessedPM || guessedCharge)
      {
        (*m_outputSpectra)[i].parentMass = curMass;
        (*m_outputSpectra)[i].setCharge(curCharge);
      }

      if ((*m_outputSpectra)[i].size() == 0)
      {
        numEmpty++;
      }
    }

    DEBUG_VAR(m_outputSpectra->size());
    DEBUG_MSG("Have " << numEmpty << " empty spectra");

    bool compressResults = m_params.getValueBool("SCORING_COMPRESS_RESULTS");
    if (compressResults) {
      // Compress the PRM spectra vector so that less memory is required
      DEBUG_TRACE;
      float peakTol = m_params.getValueDouble("TOLERANCE_PEAK", 0.5);
      bool specTypeMSMS = m_params.getValueBool("SPEC_TYPE_MSMS", false);
      float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
      unsigned int lastUsed = 0;
      for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
      {
          size_t peakSize = (*m_outputSpectra)[i].size();
          if (peakSize != 0) {
              if (lastUsed != i) {
                  (*m_outputSpectra)[lastUsed] = (*m_outputSpectra)[i];
              }
              lastUsed++;
          }
      }
      m_outputSpectra->resize(lastUsed);
      DEBUG_VAR(m_outputSpectra->size());
    }


    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::loadInputData(void)
  {

    if (m_params.exists("INPUT_SPECTRA"))
    {
      if (!m_inputSpectra->Load(m_params.getValue("INPUT_SPECTRA").c_str()))
      {
        ERROR_MSG("Failed to load from \'" << m_params.getValue("INPUT_SPECTRA") << "\'");
        return false;
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_SPECTRA?");
      return false;
    }

    // If re-loading the subset of spectra on a child node, must start processing at zero
    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)m_inputSpectra->size()) - 1);

    if (idxEnd - idxStart + 1 == (int)m_inputSpectra->size())
    {
      m_params.setValue("IDX_END", parseInt(m_inputSpectra->size() - 1));
      m_params.setValue("IDX_START", parseInt(0));
    }

    /*
     m_scanToIdx.resize(m_inputSpectra->size());
     m_scanToInputScan.resize(m_inputSpectra->size());
     unsigned int total = 0;
     for (unsigned int i = 0; i < m_loader->m_recordedInput.size(); i++)
     {
     for (unsigned int j = 0; j < m_loader->m_recordedInput[i].second; j++)
     {
     m_scanToInputScan[total] = (*m_inputSpectra)[total].scan;
     (*m_inputSpectra)[total].scan = total;
     m_scanToIdx[total].first = i;
     m_scanToIdx[total].second = j;
     }
     }
     */

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR_IN");
    if (dataDir.empty()) {
      dataDir = ".";
    }
    DEBUG_VAR(dataDir);

    string spectraFilename = m_params.getValue("INPUT_SPECTRA");

    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)m_inputSpectra->size()) - 1);

    SpecSet spectraSubSet(idxEnd - idxStart + 1);
    for (int i = idxStart; i <= idxEnd; i++)
    {
      spectraSubSet[i - idxStart] = (*m_inputSpectra)[i];
    }

    if (!ExecMergeConvert::saveSpecsetMultiple(spectraFilename, &spectraSubSet))
    {
      ERROR_MSG("Failed to save to \'" << spectraFilename << "\'");
      return false;
    }

    string paramDir = m_params.getValue("GRID_DATA_DIR_PARAMS");
    if (paramDir.empty()) {
      paramDir = ".";
    }
    string baseDirectory = paramDir + "/";
    string baseFilename = baseDirectory + getName();
    string paramFilename = baseFilename + ".params";
    DEBUG_VAR(paramFilename);

    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(spectraFilename);

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::saveOutputData(void)
  {
    string outDir = m_params.getValue("GRID_DATA_DIR_OUT");
    if (outDir.empty()) {
      outDir = ".";
    }

    if (m_params.exists("OUTPUT_SPECTRA")) {
      string saveName = outDir + "/" + m_params.getValue("OUTPUT_SPECTRA");
      if (!ExecMergeConvert::saveSpecsetMultiple(saveName,
                                                 m_outputSpectra)) {
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::loadOutputData(void)
  {
    string outDir = m_params.getValue("GRID_DATA_DIR_INTERMEDIATE");
    if (outDir.empty()) {
      outDir = ".";
    }

    if (m_params.exists("OUTPUT_SPECTRA")) {
      string fileName = outDir + "/" + m_params.getValue("OUTPUT_SPECTRA");
      if (!m_outputSpectra->Load(fileName.c_str())) {
        ERROR_MSG("Failed to load from [" << fileName << "]");
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPrmScoring::split(int numSplit)
  {
    m_subModules.resize(0);
    if (m_inputSpectra->size() == 0)
    {
      if (!loadInputData())
      {
        return m_subModules;
      }
    }

    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)m_inputSpectra->size()) - 1);
    DEBUG_VAR(idxStart);
    DEBUG_VAR(idxEnd);

    DEBUG_VAR(numSplit);

    if (numSplit == 0)
    {
      return m_subModules;
    }

    double totalComputations = 0;
    vector<double> computationsPerSpec(m_inputSpectra->size());
    for (int i = idxStart; i <= idxEnd; i++)
    {
      computationsPerSpec[i] = (*m_inputSpectra)[i].parentMass;
      totalComputations += computationsPerSpec[i];
    }
    double computationsPerModule = totalComputations / ((double)numSplit);

    DEBUG_VAR(totalComputations);
    DEBUG_VAR(computationsPerModule);

    int lastIdx1 = -1;
    int curIdxStart = idxStart;
    int curSplit = 0;
    m_subModules.resize(numSplit, (ExecBase*)0);
    double runningComputations = 0;
    double runningComputationsPerModule = 0;

    string dataDirIn = m_params.getValue("GRID_DATA_DIR_IN");
    if (dataDirIn.empty())
    {
      dataDirIn = ".";
    }
    string dataDirIntermediate = m_params.getValue("GRID_DATA_DIR_INTERMEDIATE");
    if (dataDirIntermediate.empty())
    {
      dataDirIntermediate = ".";
    }

    for (int i = idxStart; i <= idxEnd; i++)
    {
      runningComputations += computationsPerSpec[i];
      runningComputationsPerModule += computationsPerSpec[i];

      if ((runningComputationsPerModule >= computationsPerModule
          && curSplit < numSplit + 1) || i == idxEnd)
      {
        ParameterList childParams(m_params);
        childParams.removeParam("GRID_EXECUTION"); // necessary for Proteosafe

        childParams.setValue("GRID_DATA_DIR_OUT", dataDirIntermediate);
        
        childParams.setValue("IDX_END", parseInt(i));
        childParams.setValue("IDX_START", parseInt(curIdxStart));

        DEBUG_MSG("Job " << curSplit << " has range [" << curIdxStart << "," << i << "]");

        ExecBase * theClone = new ExecPrmScoring(childParams, m_inputSpectra);

        // Give it a new name based on the split
        theClone->setName(makeName(m_name, curSplit));

        // Have to set up the output files also so the params will be correct on reload
        string baseName = theClone->getName();
        //DEBUG_VAR(baseName);

        theClone->m_params.setValue("INPUT_SPECTRA", 
                dataDirIn + "/" + baseName + "_input_spectra.pklbin");
        theClone->m_params.setValue("OUTPUT_SPECTRA", 
                baseName + "_output_spectra.pklbin");

        m_subModules[curSplit++] = theClone;
        curIdxStart = i + 1;

        runningComputationsPerModule = 0;

        if (curSplit < numSplit)
        {
          computationsPerModule = (totalComputations - runningComputations)
              / (((double)numSplit) - ((double)curSplit));
        }
      }
    }
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::merge(void)
  {
    m_outputSpectra->resize(0);

    int numSpecsLoaded = 0;

    for (int nSplit = 0; nSplit < m_subModules.size(); nSplit++)
    {
      if (m_subModules[nSplit] == 0)
      {
        ERROR_MSG("ExecBase module at index " << nSplit << " is zero");
        return false;
      }
      ExecPrmScoring* module =
          dynamic_cast<ExecPrmScoring*>(m_subModules[nSplit]);

      if (module == 0)
      {
        ERROR_MSG("Failed to cast module \'" << m_subModules[nSplit]->getName() << "\' to ExecPrmScoring");
        return false;
      }

      for (int i = 0; i < module->m_outputSpectra->size(); i++)
      {
        m_outputSpectra->push_back((*(module->m_outputSpectra))[i]);
      }

      numSpecsLoaded += module->m_outputSpectra->size();
    }
    DEBUG_VAR(m_outputSpectra->size());
      
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("CORRECT_PM");
    VALIDATE_PARAM_EXIST("GUESS_CHARGE");
    if (m_params.exists("PEPNOVO_EXE_DIR"))
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_EXE_DIR");
    }
    else
    {
      VALIDATE_PARAM_EXIST("EXE_DIR");
    }
    if (m_params.exists("PEPNOVO_MODEL"))
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_INPUT_MGF");
      VALIDATE_PARAM_EXIST("PEPNOVO_OUTPUT_PRMS");
    }
    else
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_OUTDIR");
    }

    m_isValid = true;
    return true;
  }

  bool ExecPrmScoring::invokePepNovoCID(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra)
  {

    DEBUG_MSG("Invoking PepNovo for input CID spectra...");

    string ftType = "FT";
    bool highAcc = ftType == params.getValue("INSTRUMENT_TYPE", "IT");

    ParameterList pepnovoParams(params);

    if (highAcc)
    {
      pepnovoParams.setValue("PEPNOVO_MODEL", "CID_FT_Tryp_z");
    }
    else
    {
      pepnovoParams.setValue("PEPNOVO_MODEL", "CID_IT_TRYP");
    }

    if ((!enforceDaTolerance) && highAcc)
    {
      pepnovoParams.setValue("TOLERANCE_PEAK", "0.04");
    }
    else if ((!enforceDaTolerance) && (!highAcc))
    {
      pepnovoParams.setValue("TOLERANCE_PEAK", "0.5");

    }

    string inFile(pepnovoParams.getValue("PEPNOVO_OUTDIR"));
    inFile = this->getTemporaryFilename(inFile, "pepnovo_in_CID.mgf");
    pepnovoParams.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(pepnovoParams.getValue("PEPNOVO_OUTDIR"));
    outFile = this->getTemporaryFilename(outFile, "pepnovo_out_CID.prms");
    pepnovoParams.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(pepnovoParams, inputSpectra, outputSpectra);

    if (!res)
    {
      return false;
    }

    return true;
  }

  bool ExecPrmScoring::invokePepNovoHCD(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra)
  {

    DEBUG_MSG("Invoking PepNovo for input HCD spectra...");

    ParameterList pepnovoParams(params);

    string ftType = "FT";

    pepnovoParams.setValue("TOLERANCE_PEAK", "0.04");

    SpecSet HCDms2charge2(inputSpectra.size());
    unsigned int charge2Idx = 0;
    SpecSet HCDms2charge3(inputSpectra.size());
    unsigned int charge3Idx = 0;

    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      // extract charge 2 and charge >= 3 spectra so they can be scored with different models
      Spectrum* hcdSpec = &inputSpectra[i];
      if (hcdSpec->parentCharge <= 2)
      {
        HCDms2charge2[charge2Idx] = *hcdSpec;
        ++charge2Idx;
      }
      else
      {
        HCDms2charge3[charge3Idx] = *hcdSpec;
        HCDms2charge3[charge3Idx].setCharge(3);
        ++charge3Idx;
      }
      // Only keep one copy of each MS/MS spectrum in memory
      hcdSpec->resize(0);
    }
    inputSpectra.resize(0);
    HCDms2charge2.resize(charge2Idx);
    HCDms2charge3.resize(charge3Idx);

    SpecSet HCDprms2(HCDms2charge2.size());
    SpecSet HCDprms3(HCDms2charge3.size());

    DEBUG_VAR(HCDms2charge2.size());
    DEBUG_VAR(HCDms2charge3.size());

    bool returnVal = true;

    ParameterList charge2Params(pepnovoParams);
    ParameterList charge3Params(pepnovoParams);

    charge2Params.setValue("PEPNOVO_MODEL", "HCD_Tryp_z");
    charge3Params.setValue("PEPNOVO_MODEL", "HCD_LysC_z");

    if (m_params.getValueBool("USE_EXACT_PM", false))
    {
      charge2Params.setValue("PEPNOVO_MODEL", "HCD_z_corrPm");
      charge3Params.setValue("PEPNOVO_MODEL", "HCD_z_corrPm");
    }
    else if (m_params.getValue("INSTRUMENT_TYPE", "FT") == "qtof")
    {
      charge2Params.setValue("PEPNOVO_MODEL", "HCD_aLP_z");
      charge3Params.setValue("PEPNOVO_MODEL", "HCD_aLP_z");
    }

    string inFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    inFile = this->getTemporaryFilename(inFile, "pepnovo_in_HCD_2.mgf");
    charge2Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    outFile = this->getTemporaryFilename(outFile, "pepnovo_out_HCD_2.prms");
    charge2Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(charge2Params, HCDms2charge2, HCDprms2);

    if (!res)
    {
      return false;
    }

    inFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    inFile = this->getTemporaryFilename(inFile, "pepnovo_in_HCD_3.mgf");
    charge3Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    outFile = charge2Params.getValue("PEPNOVO_OUTDIR");
    outFile = this->getTemporaryFilename(outFile, "pepnovo_out_HCD_3.prms");
    charge3Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    res = invokePepNovo(charge3Params, HCDms2charge3, HCDprms3);

    if (!res)
    {
      return false;
    }

    outputSpectra = HCDprms2;
    outputSpectra.appendSpecSet(HCDprms3, false);

    // Restore input MS/MS spectra
    inputSpectra.resize(HCDms2charge2.size() + HCDms2charge3.size());
    for (unsigned int i = 0; i < HCDms2charge2.size(); i++)
    {
      inputSpectra[i] = HCDms2charge2[i];
      HCDms2charge2[i].resize(0);
    }
    HCDms2charge2.resize(0);

    for (unsigned int i = 0; i < HCDms2charge3.size(); i++)
    {
      inputSpectra[HCDms2charge2.size() + i] = HCDms2charge3[i];
      HCDms2charge3[i].resize(0);
    }
    return true;
  }

  bool ExecPrmScoring::invokePepNovoETD(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra)
  {

    DEBUG_MSG("Invoking PepNovo for input ETD spectra...");

    ParameterList pepnovoParams(params);

    string ftType = "FT";
    bool highAcc = ftType == pepnovoParams.getValue("INSTRUMENT_TYPE", "IT");

    SpecSet ETDms2charge2(inputSpectra.size());
    unsigned int charge2Idx = 0;
    SpecSet ETDms2charge3(inputSpectra.size());
    unsigned int charge3Idx = 0;

    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      // extract charge 2 and charge >= 3 spectra so they can be scored with different models
      Spectrum* etdSpec = &inputSpectra[i];
      if (etdSpec->parentCharge <= 2)
      {
        ETDms2charge2[charge2Idx] = *etdSpec;
        ++charge2Idx;
      }
      else
      {
        ETDms2charge3[charge3Idx] = *etdSpec;
        ++charge3Idx;
      }
      // Only keep one copy of each MS/MS spectrum in memory
      etdSpec->resize(0);
    }
    inputSpectra.resize(0);

    ETDms2charge2.resize(charge2Idx);
    ETDms2charge3.resize(charge3Idx);

    SpecSet ETDprms2(ETDms2charge2.size());
    SpecSet ETDprms3(ETDms2charge3.size());

    DEBUG_VAR(ETDms2charge2.size());
    DEBUG_VAR(ETDms2charge3.size());

    bool returnVal = true;

    ParameterList charge2Params(pepnovoParams);
    ParameterList charge3Params(pepnovoParams);

    if (highAcc)
    {
      charge2Params.setValue("PEPNOVO_MODEL", "ETD_FT_Tryp_z");
      charge3Params.setValue("PEPNOVO_MODEL", "ETD_FT_LysC_z");
    }
    else
    {
      charge2Params.setValue("PEPNOVO_MODEL", "ETD_IT_Tryp");
      charge3Params.setValue("PEPNOVO_MODEL", "ETD_IT_LysC");
    }

    string inFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    inFile = this->getTemporaryFilename(inFile, "pepnovo_in_ETD_2.mgf");
    charge2Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    outFile = this->getTemporaryFilename(outFile, "pepnovo_out_ETD_2.prms");
    charge2Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(charge2Params, ETDms2charge2, ETDprms2);

    if (!res)
    {
      return false;
    }

    inFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    inFile = this->getTemporaryFilename(inFile, "pepnovo_in_ETD_3.mgf");
    charge3Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    outFile = charge2Params.getValue("PEPNOVO_OUTDIR");
    outFile = this->getTemporaryFilename(outFile, "pepnovo_out_ETD_3.prms");
    charge3Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    res = invokePepNovo(charge3Params, ETDms2charge3, ETDprms3);

    if (!res)
    {
      return false;
    }

    outputSpectra = ETDprms2;
    outputSpectra.appendSpecSet(ETDprms3, false);

    // Restore input MS/MS spectra
    inputSpectra.resize(ETDms2charge2.size() + ETDms2charge3.size());
    for (unsigned int i = 0; i < ETDms2charge2.size(); i++)
    {
      inputSpectra[i] = ETDms2charge2[i];
      ETDms2charge2[i].resize(0);
    }
    ETDms2charge2.resize(0);

    for (unsigned int i = 0; i < ETDms2charge3.size(); i++)
    {
      inputSpectra[ETDms2charge2.size() + i] = ETDms2charge3[i];
      ETDms2charge3[i].resize(0);
    }
    return true;
  }

  bool ExecPrmScoring::invokePepNovo(ParameterList& params,
                                     SpecSet& inputSpectra,
                                     SpecSet& outputSpectra)
  {

    if (inputSpectra.size() == 0)
    {
      outputSpectra.resize(0);
      return true;
    }

    DEBUG_MSG("Invoking PepNovo ...");
    //DEBUG_VAR(inputSpectra.size());

    //DEBUG_VAR(params.getValue("TOLERANCE_PEAK"));

    string pepnovoMode;
    if (params.exists("MIN_SPECTRUM_QUALITY"))
    {
      pepnovoMode = "-min_filter_prob ";
      pepnovoMode += params.getValue("MIN_SPECTRUM_QUALITY");
    }
    else
    {
      pepnovoMode = "-no_quality_filter";
    }

    // TODO: Lower case the value?
    bool guessedPM = false;
    if (params.getValue("CORRECT_PM") == "yes")
    {
      pepnovoMode += " -correct_pm";
      guessedPM = true;
    }
    else
    {
      pepnovoMode += " -use_spectrum_mz";
    }

    // TODO: Lower case the value?
    bool guessedCharge = true;
    if (params.getValue("GUESS_CHARGE") == "no")
    {
      pepnovoMode += " -use_spectrum_charge";
      guessedCharge = false;

      for (unsigned int i = 0; i < inputSpectra.size(); i++)
      {
        if (inputSpectra[i].parentCharge == 0)
        {
          WARN_MSG("Found charge 0+ for spectrum " << i << " (scan " << inputSpectra[i].scan << "), resizing to peak list zero");
          inputSpectra[i].parentCharge = 2;
          inputSpectra[i].resize(0);
        }
      }
    }

    if (params.exists("PEPNOVO_PTMS"))
    {
      pepnovoMode += " -PTMs M+16:C+57:";
      pepnovoMode += params.getValue("PEPNOVO_PTMS");
    }
    else
    {
      pepnovoMode += " -PTMs M+16:C+57";
    }

    //DEBUG_VAR(pepnovoMode);
    string modelName = params.getValue("PEPNOVO_MODEL");

    string exeDir = params.getValue("EXE_DIR");
    string pepNovoModelDir(exeDir);

    // pepnovo models directory
    pepNovoModelDir = exeDir;
    if (params.exists("PEPNOVO_MODEL_DIR"))
    {
      pepNovoModelDir = params.getValue("PEPNOVO_MODEL_DIR");
    }
    else
    {
      pepNovoModelDir += "/resources/Models_pepnovo";
    }

    // pepnovo executable directory
    if (params.exists("PEPNOVO_EXE_DIR"))
    {
      exeDir = params.getValue("PEPNOVO_EXE_DIR");
    }

    string pepnovoCmd(exeDir);
    pepnovoCmd += "/PepNovo_bin -prm_norm -model_dir ";
    pepnovoCmd += pepNovoModelDir;
    pepnovoCmd += " -model ";
    pepnovoCmd += modelName;
    pepnovoCmd += " -file ";
    pepnovoCmd += params.getValue("PEPNOVO_INPUT_MGF");
    pepnovoCmd += " -fragment_tolerance ";
    pepnovoCmd += params.getValue("TOLERANCE_PEAK");
    if (params.exists("TOLERANCE_PM"))
    {
      pepnovoCmd += " -pm_tolerance ";
      pepnovoCmd += params.getValue("TOLERANCE_PM");
    }
    pepnovoCmd += " -digest NON_SPECIFIC ";
    pepnovoCmd += pepnovoMode;
    pepnovoCmd += " > ";
    pepnovoCmd += params.getValue("PEPNOVO_OUTPUT_PRMS");

    // Remove spectra with 5 or less peaks.
    checkSpecset(inputSpectra);

    map<int, int> scanRef;
    for (int i = 0; i < inputSpectra.size(); i++)
    {
      scanRef[inputSpectra[i].scan] = i;
    }

    // Prepare a copy of MS/MS spectra for PepNovo input

    vector<short> oldCharges(inputSpectra.size());

    for (int i = 0; i < inputSpectra.size(); i++)
    {
      Spectrum* spec = &inputSpectra[i];
      oldCharges[i] = spec->parentCharge;
      // ETD spectra have a charge 4 model while HCD spectra only have a charge 3 model
      if (modelName.find("ETD") != string::npos && spec->parentCharge > 4)
      {
        spec->setCharge(4);
      }
      else if (modelName.find("ETD") == string::npos && spec->parentCharge > 3)
      {
        spec->setCharge(3);
      }
    }

    DEBUG_VAR(enforceDaTolerance);

    ifstream infile;
    infile.open(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str(), ios::binary);

    if (params.getValueInt("SKIP_PEPNOVO_INVOKE", 0) == 0 || !infile.is_open())
    {
      if (infile.is_open())
      {
        infile.close();
      }

      // invoke PepNovo via command line
      DEBUG_MSG("Saving temporary file \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
      if (!inputSpectra.SaveSpecSet_mgf(params.getValue("PEPNOVO_INPUT_MGF").c_str()), 0, false)
      {
        return false;
      }

      DEBUG_VAR(pepnovoCmd);
      int status = spsSystem(pepnovoCmd.c_str());
      if (status != 0)
      {
        string errorString = "Executing ";
        errorString += exeDir;
        errorString += "/Pepnovo_bin!";
        ERROR_MSG(errorString);
        return false;
      }
    }
    if (infile.is_open())
    {
      infile.close();
    }

    SpecSet * pepnovoSpectra = new SpecSet(inputSpectra.size());

    DEBUG_MSG("Loading temporary file \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
    if (!pepnovoSpectra->LoadSpecSet_prmsv3(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str()))
    {
      return false;
    }

    if (params.getValueInt("SKIP_PEPNOVO_INVOKE", 0) == 0)
    {
      // remove temporary input/output files if they were generated
      if (remove(params.getValue("PEPNOVO_INPUT_MGF").c_str()) != 0)
      {
        int aa = errno;
        string aux = strerror(aa);
        ERROR_MSG("Failed to remove \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
        ERROR_MSG("error = " << aux);
        //        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
      }

      if (remove(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str()) != 0)
      {
        int aa = errno;
        string aux = strerror(aa);
        ERROR_MSG("Failed to remove \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
        ERROR_MSG("error = " << aux);
        //        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
      }
    }

    //DEBUG_VAR(pepnovoSpectra->size());

    outputSpectra.resize(inputSpectra.size());

    for (int i = 0; i < inputSpectra.size(); i++)
    {
      inputSpectra[i].setCharge(oldCharges[i]);
      outputSpectra[i].resize(0);
      outputSpectra[i].copyNP(inputSpectra[i]);
    }
    for (int i = 0; i < pepnovoSpectra->size(); i++)
    {
      int specIdx = scanRef[(*pepnovoSpectra)[i].scan];

      // Copy PepNovo output spectra

      outputSpectra[specIdx].mergeClosestPeaks((*pepnovoSpectra)[i], 0);

      if (guessedPM)
      {
        // Let PepNovo set parent mass if we asked it to guess it
        outputSpectra[specIdx].setParentMass((*pepnovoSpectra)[i].parentMass);
      }

      if (guessedCharge || outputSpectra[specIdx].parentCharge == 0)
      {
        // Let PepNovo set the parent charge if we asked it to guess it
        outputSpectra[specIdx].setCharge((*pepnovoSpectra)[i].parentCharge);
      }

      Spectrum* ms2Spec = &inputSpectra[specIdx];
      Spectrum* prmSpec = &outputSpectra[specIdx];

      prmSpec->setPeakTolerance(params.getValueFloat("TOLERANCE_PEAK"));
      prmSpec->filterLowIntensity(-1.0);
    }

    delete pepnovoSpectra;

    //DEBUG_VAR(outputSpectra.size());

    return true;
  }

  /*
   void ExecPrmScoring::setOriginMasses(Spectrum* prmSpec,
   vector<bool>* reversedPeaks,
   Spectrum* ms2Spec,
   vector<string>* prmOriginPeaks)
   {

   reversedPeaks->resize(prmSpec->size());

   if (prmSpec->size() == 0)
   {
   return;
   }

   Spectrum expectedPrmOffsets;
   expectedPrmOffsets.resize(8);
   expectedPrmOffsets[0][0] = 1.00728;
   expectedPrmOffsets[1][0] = -17.00328;
   expectedPrmOffsets[2][0] = -26.98710;
   expectedPrmOffsets[3][0] = -16.01927;
   expectedPrmOffsets[4][0] = 18.03430;
   expectedPrmOffsets[5][0] = 17.02232;
   expectedPrmOffsets[6][0] = -25.98199;
   expectedPrmOffsets[7][0] = 9.52079;
   expectedPrmOffsets.sortPeaks();

   Spectrum expectedSrmOffsets;
   expectedSrmOffsets.resize(5);
   expectedSrmOffsets[0][0] = 19.01840;
   expectedSrmOffsets[1][0] = 1.00784;
   expectedSrmOffsets[2][0] = 10.01284;
   expectedSrmOffsets[3][0] = 2.99000;
   expectedSrmOffsets[4][0] = 1.99864;
   expectedSrmOffsets.sortPeaks();

   prmSpec->setTolerance(0, 0.001);
   prmSpec->setTolerance(prmSpec->size() - 1, ms2Spec->parentMassTol);

   map<float, bool> reversedPeakRef;
   reversedPeakRef[prmSpec->front()->operator [](0)] = false;
   reversedPeakRef[prmSpec->back()->operator [](0)] = false;

   for (unsigned int i = 1; i < prmSpec->size() - 1; i++)
   {

   vector<string> originInfo(4);
   splitText((*prmOriginPeaks)[i].c_str(), originInfo, ",");
   float ms2mass = getFloat(originInfo[0].c_str());
   int orient = getInt(originInfo[1].c_str());
   int charge = getInt(originInfo[2].c_str());
   float massOffset = getFloat(originInfo[3].c_str());

   float srmOffset = 0;

   bool isSuf = false;

   if (orient == 0)
   { // from prefix ion
   int closestIdx = expectedPrmOffsets.findClosest(massOffset);
   if (MZRange::EqualWithinRange(massOffset,
   expectedPrmOffsets[closestIdx][0],
   0.001))
   {
   massOffset = expectedPrmOffsets[closestIdx][0];
   }
   }
   else
   { // suffix ion
   isSuf = true;
   int closestIdx = expectedSrmOffsets.findClosest(massOffset);
   if (MZRange::EqualWithinRange(massOffset,
   expectedSrmOffsets[closestIdx][0],
   0.001))
   {
   massOffset = expectedSrmOffsets[closestIdx][0];
   }
   srmOffset =
   (ms2Spec->msFragType == Spectrum::FragType_ETD) ?
   0.0 - AAJumps::massNH : AAJumps::massH2O;

   }

   int m2Idx = ms2Spec->findClosest(ms2mass);
   float massUse = (*ms2Spec)[m2Idx][0];

   if (charge > 1)
   {
   massUse = DeconvSpectrum::GetMonoisotopicMass(massUse, charge);
   }
   (*prmSpec)[i][0] = massUse - massOffset + srmOffset;
   prmSpec->setTolerance(i, ms2Spec->getTolerance(m2Idx) * ((float)charge));

   reversedPeakRef[(*prmSpec)[i][0]] = isSuf;

   }
   prmSpec->sortPeaks();

   for (unsigned int i = 0; i < prmSpec->size(); i++)
   {
   (*reversedPeaks)[i] = reversedPeakRef[(*prmSpec)[i][0]];
   }
   }

   void ExecPrmScoring::setOriginTolerances(Spectrum* prmSpec,
   Spectrum* ms2Spec,
   vector<string>* prmOriginPeaks)
   {
   if (prmSpec->size() == 0)
   {
   return;
   }

   prmSpec->setTolerance(0, 0.001);
   prmSpec->setTolerance(prmSpec->size() - 1, ms2Spec->parentMassTol);

   for (unsigned int i = 1; i < prmSpec->size() - 1; i++)
   {

   vector<string> originInfo(4);
   splitText((*prmOriginPeaks)[i].c_str(), originInfo, ",");
   float ms2mass = getFloat(originInfo[0].c_str());
   int orient = getInt(originInfo[1].c_str());
   int charge = getInt(originInfo[2].c_str());

   int m2Idx = ms2Spec->findClosest(ms2mass);

   float locms2Tol = (ms2Spec->getTolerance(m2Idx) + ((float)charge))
   + 0.001;

   if (orient == 0)
   { // from prefix ion
   prmSpec->setTolerance(i, locms2Tol);
   }
   else
   { // suffix ion
   prmSpec->setTolerance(i, locms2Tol + ms2Spec->parentMassTol);
   }
   }
   }

   */

  void ExecPrmScoring::checkSpecset(SpecSet &inOutSpecs)
  {
    DEBUG_MSG("Checking sepecset. Size: " << inOutSpecs.size());

    int numResized = 0;

    for (int i = 0; i < inOutSpecs.size(); i++)
    {
      //DEBUG_VAR(inOutSpecs[i].size());
      if (inOutSpecs[i].size() <= 8)
      {
        //DEBUG_MSG("Spectrum #" << i << " has " << inOutSpecs[i].size() << " peaks. Resizing to 0.");
        inOutSpecs[i].resize(0);
        numResized++;
      }
    }

    DEBUG_MSG("Checking sepecset done, emptied " << numResized << " spectra with less than 9 peaks");
  }

} // namespace specnets
