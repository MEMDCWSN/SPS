/*
 * ExecTrainMassError.cpp
 *
 *  Created on: Nov 15, 2013
 *      Author: aguthals
 */

#include "ExecTrainMassError.h"
#include "utils.h"
#include "ExecMergeConvert.h"
#include "OutputTable.h"
#include "FdrPeptide.h"
//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#  include <tr1/unordered_set>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#  include <unordered_set>
#endif

using namespace std;

namespace specnets
{

  const int ExecTrainMassError::MIN_ALLOWED_PSMS = 1000;

  ExecTrainMassError::ExecTrainMassError(void) :
      ExecBase(), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputSpectraPRM(0x0), m_inputPSMs(0x0), m_outputModelTarget(0x0), m_outputModelDecoy(0x0)
  {
    m_name = "ExecTrainMassError";
    m_type = "ExecTrainMassError";
    m_inputSpectra = new SpecSet;
    m_inputSpectraPRM = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_outputModelTarget = new MassErrorModel;
    m_outputModelDecoy = new MassErrorModel;
  }

  ExecTrainMassError::ExecTrainMassError(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputSpectraPRM(0x0), m_inputPSMs(0x0), m_outputModelTarget(0x0), m_outputModelDecoy(0x0)
  {
    m_name = "ExecTrainMassError";
    m_type = "ExecTrainMassError";
    m_inputSpectra = new SpecSet;
    m_inputSpectraPRM = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_outputModelTarget = new MassErrorModel;
    m_outputModelDecoy = new MassErrorModel;
  }

  ExecTrainMassError::ExecTrainMassError(const ParameterList & inputParams,
                                         SpecSet * inputSpectra,
                                         SpecSet * inputSpectraPRM,
                                         PeptideSpectrumMatchSet * inputPSMs) :
      ExecBase(inputParams), ownInput(false), ownOutput(true), m_inputSpectra(inputSpectra), m_inputSpectraPRM(inputSpectraPRM), m_inputPSMs(inputPSMs), m_outputModelTarget(0x0), m_outputModelDecoy(0x0)
  {
    m_name = "ExecTrainMassError";
    m_type = "ExecTrainMassError";
    m_outputModelTarget = new MassErrorModel;
    m_outputModelDecoy = new MassErrorModel;
  }

  ExecTrainMassError::ExecTrainMassError(const ParameterList & inputParams,
                                         SpecSet * inputSpectra,
                                         SpecSet * inputSpectraPRM,
                                         PeptideSpectrumMatchSet * inputPSMs,
                                         MassErrorModel * outputModelTarget,
                                         MassErrorModel * outputModelDecoy) :
      ExecBase(inputParams), ownInput(false), ownOutput(false), m_inputSpectra(inputSpectra), m_inputSpectraPRM(inputSpectraPRM), m_inputPSMs(inputPSMs), m_outputModelTarget(outputModelTarget), m_outputModelDecoy(outputModelDecoy)
  {
    m_name = "ExecTrainMassError";
    m_type = "ExecTrainMassError";
  }

  ExecTrainMassError::~ExecTrainMassError(void)
  {
    if (ownInput)
    {
      delete m_inputPSMs;
      delete m_inputSpectra;
      delete m_inputSpectraPRM;
    }
    if (ownOutput)
    {
      delete m_outputModelTarget;
      delete m_outputModelDecoy;
    }
  }

  ExecBase * ExecTrainMassError::clone(const ParameterList & inputParams) const
  {
    return new ExecTrainMassError(inputParams);
  }

  bool ExecTrainMassError::invoke(void)
  {
    const string decoyID = m_params.getValue("DECOY_PROTEIN_ID");

    PeptideSpectrumMatchSet topPsms;

    int numTotalPSMs = 0;
    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      numTotalPSMs += ((*m_inputSpectra)[i].psmList.size() > 0) ? 1 : 0;
    }

    if (numTotalPSMs < MIN_ALLOWED_PSMS)
    {
      WARN_MSG("Number of PSMs (" << numTotalPSMs << ") is insufficient to train mass error model");
      return false;
    }

    topPsms.resize(numTotalPSMs);
    int idxUse = 0;

    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      double minPval = 2.0;
      const psmPtr *bestPsm = 0x0;

      for (list<psmPtr>::const_iterator pIt =
          (*m_inputSpectra)[i].psmList.begin();
          pIt != (*m_inputSpectra)[i].psmList.end(); pIt++)
      {
        if ((*pIt)->m_pValue < minPval)
        {
          bestPsm = &(*pIt);
          minPval = (*pIt)->m_pValue;
        }
      }

      (*m_inputSpectra)[i].psmList.clear();
      if (bestPsm != 0x0)
      {
        topPsms[idxUse++] = *bestPsm;
      }
    }

    topPsms.resize(idxUse);

    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      (*m_inputSpectra)[i].psmList.clear();
    }

    for (int i = 0; i < topPsms.size(); i++)
    {
      topPsms[i]->m_isDecoy = topPsms[i]->m_protein.find(decoyID)
          != string::npos;
    }

    PeptideSpectrumMatchSet topPsmsTDA;
    //PeptideSpectrumMatchSet nonNetPSMsTDA;

    DEBUG_TRACE;

    FdrPeptide::calculatePValues(topPsms,
                                 topPsmsTDA,
                                 1,
                                 FDR_SORT_PVALUE,
                                 true);

    DEBUG_TRACE;

    int numFDR = 0;
    double maxCutoff = 0;
    tr1::unordered_set<string> foundPeptides;

    for (int i = 0; i < topPsmsTDA.size(); i++)
    {
      if (topPsmsTDA[i]->m_pepFdr <= 0.01 && !topPsmsTDA[i]->m_isDecoy)
      {
        foundPeptides.insert(topPsmsTDA[i]->m_annotation);
      }

      if (topPsmsTDA[i]->m_fdr <= 0.01 && !topPsmsTDA[i]->m_isDecoy)
      {
        numFDR++;
        maxCutoff = max(maxCutoff, topPsmsTDA[i]->m_pValue);
      }
    }

    DEBUG_MSG("Have " << numFDR << " PSMs at 1% spectrum-level FDR");
    DEBUG_MSG("Have " << foundPeptides.size() << " peptides at 1% peptide-level FDR");
    DEBUG_MSG("Maximum probability cutoff is " << maxCutoff << " at 1% spectrum-level FDR");

    if (numFDR < MIN_ALLOWED_PSMS)
    {
      WARN_MSG("Number of PSMs (" << numFDR << ") is insufficient to train mass error model");
      return false;
    }

    topPsmsTDA.addMostSpectra(m_inputSpectra, true);

    AAJumps jumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);

    DEBUG_TRACE;
    m_outputModelTarget->initialize(*m_inputSpectra,
                                    *m_inputSpectraPRM,
                                    jumps,
                                    false);
    DEBUG_TRACE;
    m_outputModelDecoy->initialize(*m_inputSpectra,
                                   *m_inputSpectraPRM,
                                   jumps,
                                   true);
    DEBUG_TRACE;

    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      (*m_inputSpectra)[i].psmList.clear();
    }

    return true;
  }

  bool ExecTrainMassError::loadInputData(void)
  {
    if (m_params.exists("INPUT_SPECTRA_MS"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA_MS from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\' ...");
      if (!m_inputSpectra->Load(m_params.getValue("INPUT_SPECTRA_MS").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputSpectra->size());
      for (int i = 0; i < m_inputSpectra->size(); i++)
      {
        (*m_inputSpectra)[i].psmList.clear();
      }
    }

    if (m_params.exists("INPUT_SPECTRA"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA from \'" << m_params.getValue("INPUT_SPECTRA") << "\' ...");
      if (!m_inputSpectraPRM->Load(m_params.getValue("INPUT_SPECTRA").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputSpectraPRM->size());
      for (int i = 0; i < m_inputSpectraPRM->size(); i++)
      {
        (*m_inputSpectraPRM)[i].psmList.clear();
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!!!");
      return false;
    }

    if (m_params.exists("INPUT_PSMS"))
    {
      const string psmType = m_params.getValue("INPUT_PSMS_TYPE", "");
      const string psmFile = m_params.getValue("INPUT_PSMS");

      if (!m_inputPSMs->Load(psmFile, psmType))
      {
        ERROR_MSG("Failed to load PSMs from \'" << m_params.getValue("INPUT_PSMS") << "\'");
        return false;
      }

      int numAdded = m_inputPSMs->addMostSpectra(m_inputSpectra);

      if (numAdded == 0)
      {
        ERROR_MSG("No PSMs could be matched to input spectra!!!");
        return false;
      }
    }

    if (m_inputPSMs->size() == 0)
    {
      ERROR_MSG("Input PSMs size is 0!!!");
      return false;
    }

    return true;
  }

  bool ExecTrainMassError::saveOutputData(void)
  {
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_MODEL_TARGET"))
    {
      const string filename = m_params.getValue("OUTPUT_MODEL_TARGET");
      DEBUG_MSG("Saving target mass error model to \'" << filename << "\'");
      if (!m_outputModelTarget->saveBinaryFile(filename))
      {
        ERROR_MSG("Failed to save target mass error model to \'"
            << filename << "\'");
        return false;
      }
    }
    if (m_params.exists("OUTPUT_MODEL_DECOY"))
    {
      const string filename = m_params.getValue("OUTPUT_MODEL_DECOY");
      DEBUG_MSG("Saving decoy mass error model to \'" << filename << "\'");
      if (!m_outputModelDecoy->saveBinaryFile(filename))
      {
        ERROR_MSG("Failed to save decoy mass error model to \'"
            << filename << "\'");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_ERROR_HISTOGRAM"))
    {
      OutputTable histogram;

      const string filename = m_params.getValue("OUTPUT_ERROR_HISTOGRAM");
      DEBUG_MSG("Saving TDA error histogram to \'" << filename << "\'");
      MassErrorModel::getTDAHistogram(*m_outputModelTarget, *m_outputModelDecoy, histogram);
      if (!histogram.printToCSV(filename.c_str(), "\t"))
      {
        ERROR_MSG("Failed to save TDA error histogram to \'"
            << filename << "\'");
        return false;
      }
    }
    return true;
  }

  bool ExecTrainMassError::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  bool ExecTrainMassError::loadOutputData(void)
  {
    return false;
  }

  vector<ExecBase*> const & ExecTrainMassError::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  bool ExecTrainMassError::merge(void)
  {
    return false;
  }

  bool ExecTrainMassError::validateParams(std::string & error)
  {
    VALIDATE_PARAM_EXIST("DECOY_PROTEIN_ID");
    return true;
  }
}

