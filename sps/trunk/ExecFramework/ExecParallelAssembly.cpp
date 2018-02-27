/*
 * ExecParallelAssembly.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: aguthals
 */

// Header Includes
#include "ExecParallelAssembly.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes

using namespace std;
using namespace specnets;

namespace specnets
{

  ExecParallelAssembly::ExecParallelAssembly(void) :
    ExecBase(), m_inputAbruijn(0x0), m_stars(0x0), ownInput(true),
        m_contigShifts(0x0), m_outputAbruijn(0x0), ownOutput(true)
  {
    m_name = "ExecParallelAssembly";
    m_type = "ExecParallelAssembly";
  }

  // -------------------------------------------------------------------------

  ExecParallelAssembly::ExecParallelAssembly(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputAbruijn(0x0), m_stars(0x0), ownInput(true),
        m_contigShifts(0x0), m_outputAbruijn(0x0), ownOutput(true)
  {
    m_name = "ExecParallelAssembly";
    m_type = "ExecParallelAssembly";
  }

  // -------------------------------------------------------------------------

  ExecParallelAssembly::ExecParallelAssembly(const ParameterList & inputParams,
                                             abinfo_t * inputAbruijn,
                                             SpecSet * inputStars,
                                             Clusters * outputContigShifts,
                                             abinfo_t * outputAbruijn) :
    ExecBase(inputParams), m_inputAbruijn(inputAbruijn), m_stars(inputStars),
        ownInput(false), m_contigShifts(outputContigShifts),
        m_outputAbruijn(outputAbruijn), ownOutput(false)
  {
    m_name = "ExecParallelAssembly";
    m_type = "ExecParallelAssembly";
  }

  ExecParallelAssembly::ExecParallelAssembly(const ParameterList & inputParams,
                                             abinfo_t * inputAbruijn,
                                             SpecSet * inputStars) :
    ExecBase(inputParams), m_inputAbruijn(inputAbruijn), m_stars(inputStars),
        ownInput(false), m_contigShifts(0x0), m_outputAbruijn(0x0),
        ownOutput(true)
  {
    m_name = "ExecParallelAssembly";
    m_type = "ExecParallelAssembly";
  }

  // -------------------------------------------------------------------------

  ExecParallelAssembly::~ExecParallelAssembly(void)
  {
    if (ownInput)
    {

      if (m_inputAbruijn)
        delete m_inputAbruijn;
      if (m_stars)
        delete m_stars;
    }
    if (ownOutput)
    {
      if (m_contigShifts)
        delete m_contigShifts;
      if (m_outputAbruijn)
        delete m_outputAbruijn;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecParallelAssembly::clone(const ParameterList & inputParams) const
  {
    return new ExecParallelAssembly(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::invoke(void)
  {
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM");

    abruijn::AbruijnGraph2::REVERSE_STARS = m_params.getValueInt("REVERSE_STARS", 1) > 0;

    abruijn::AbruijnGraph2::ENFORCE_B_ENDPTS = m_params.getValueInt("ENFORCE_B_ENDPTS", 0) > 0;

    abruijn::AbruijnGraph2::PEAK_TOLERANCE = (float)m_params.getValueDouble("TOLERANCE_PEAK");

    DEBUG_VAR(abruijn::AbruijnGraph2::MAX_NUM_JUMPS);

    if (!AAJumps::globalJumpsInitialized())
    {
      AAJumps::initializeGlobalJumps(abruijn::AbruijnGraph2::MAX_NUM_JUMPS,
                                     0.01,
                                     0);
    }

    const int debugIdx = -1;

    for (abinfo_t::iterator abIt = m_inputAbruijn->begin(); abIt
        != m_inputAbruijn->end(); abIt++)
    {
      unsigned int i = abIt->first;

      if (i >= m_contigShifts->consensus.size())
      {
        m_contigShifts->consensus.resize(i + 1);
      }

      if (m_inputAbruijn->count(i) == 0)
      {
        DEBUG_MSG("Skipping component " << i);
        continue;
      }

      if (debugIdx >= 0 && debugIdx != i)
      {
        continue;
      }
      pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
          double> > > > outAbinfo;

      DEBUG_MSG("Recomputing consensus for component " << i << " ...");

      abruijn::AbruijnGraph2::DEBUG_EXPANSION = false;//= i == debugIdx;

      abruijn::AbruijnGraph2 abGraph(m_inputAbruijn->at(i),
                                    *m_stars,
                                    &outAbinfo);

      if (outAbinfo.second.size() == 0)
      {
        m_contigShifts->consensus[i].resize(0);
        continue;
      }

      (*m_outputAbruijn)[i] = outAbinfo;

      abGraph.outputConsensusSpectrum(m_contigShifts->consensus[i]);

      m_contigShifts->consensus[i].setPeakTolerance(peakTol);

      if (i == debugIdx)
      {
        abGraph.saveGraphviz("postExpand.dot");
        DEBUG_MSG("Consensus spectrum for component " << i << ":");
        m_contigShifts->consensus[i].output(cerr);
        DEBUG_MSG("Abinfo for component " << i << ":");
        stringstream logMsg;
        logMsg << "Spectra: ";
        for (unsigned int j = 0; j < outAbinfo.first.first.size(); j++)
        {
          logMsg << outAbinfo.first.first[j] << "("
              << outAbinfo.first.second[j] << "), ";
        }
        DEBUG_VAR(logMsg.str());
        for (unsigned int j = 0; j < outAbinfo.second.size(); j++)
        {
          stringstream logMsg2;
          logMsg2 << "Vertex " << j << ": ";
          for (unsigned int k = 0; k < outAbinfo.second[j].first.size(); k++)
          {
            logMsg2 << outAbinfo.second[j].first[k] << "("
                << outAbinfo.second[j].second[k] << "), ";
          }
          DEBUG_VAR(logMsg2.str());
        }
      }
    }

    m_contigShifts->initializeFromAbinfo(m_contigShifts->consensus,
                                         *m_outputAbruijn,
                                         *m_stars,
                                         peakTol,
                                         peakTol);
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::loadInputData(void)
  {
    if (ownInput)
    {

      if (!m_inputAbruijn)
      {
        m_inputAbruijn = new abinfo_t;
      }
      if (!m_stars)
      {
        m_stars = new SpecSet();
      }
    }

    m_inputAbruijn->clear();
    m_stars->resize(0);
    if (ownOutput)
    {
      if (!m_contigShifts)
      {
        m_contigShifts = new Clusters;
      }
      if (!m_outputAbruijn)
      {
        m_outputAbruijn = new abinfo_t;
      }
    }
    m_contigShifts->resize(0);
    m_outputAbruijn->clear();

    if (m_params.exists("AMINO_ACID_MASSES"))
    {
      AAJumps tmpJumps(-1);
      tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
    }

    if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue("EXE_DIR"),
                                               m_params.getValue("INPUT_STARS"),
                                               m_stars) or m_stars->size() == 0)
    {
      ERROR_MSG("Error reading input star spectra from "<< m_params.getValue("INPUT_STARS"));
      return false;
    }

    if (!Load_abinfo(m_params.getValue("INPUT_ABINFO").c_str(), *m_inputAbruijn))
    {
      ERROR_MSG("Error reading abinfo from \'"<<m_params.getValue("INPUT_ABINFO")<<"\'");
      return false;
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_CLUSTERS"))
    {
      if (m_contigShifts->Save(m_params.getValue("OUTPUT_CLUSTERS").c_str())
          != 0)
      {
        ERROR_MSG("Failed to save clusters to \'" << m_params.getValue("OUTPUT_CLUSTERS") << "\'");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::saveSpecsetMultiple(m_params.getValue("OUTPUT_SPECTRA"),
                                                 &m_contigShifts->consensus))
      {
        ERROR_MSG("Failed to save meta-contig spectra to \'" << m_params.getValue("OUTPUT_SPECTRA") << "\'");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_ABINFO"))
    {
      if (!Save_abinfo_v1_0(m_params.getValue("OUTPUT_ABINFO").c_str(),
                            *m_outputAbruijn))
      {
        ERROR_MSG("Failed to save abinfo to \'" << m_params.getValue("OUTPUT_ABINFO") << "\'");
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecParallelAssembly::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecParallelAssembly::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}
