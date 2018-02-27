/*
 * ExecMetaAssembly.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: aguthals
 */

// Header Includes
#include "ExecMetaAssembly.h"
#include "ExecParallelAssembly.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes

using namespace std;
using namespace specnets;

namespace specnets
{

  ExecMetaAssembly::ExecMetaAssembly(void) :
    ExecBase(), m_spectra(0x0), m_spectrumPairs(0x0), m_inputAbruijn(0x0),
        m_stars(0x0), ownInput(true), m_contigShifts(0x0),
        m_outputAbruijn(0x0), ownOutput(true), m_overlaps(0x0),
        m_prot_match(0x0), m_fasta(0x0)
  {
    m_name = "ExecMetaAssembly";
    m_type = "ExecMetaAssembly";
  }

  // -------------------------------------------------------------------------

  ExecMetaAssembly::ExecMetaAssembly(const ParameterList & inputParams) :
    ExecBase(inputParams), m_spectra(0x0), m_spectrumPairs(0x0),
        m_inputAbruijn(0x0), m_stars(0x0), ownInput(true), m_contigShifts(0x0),
        m_outputAbruijn(0x0), ownOutput(true), m_overlaps(0x0),
        m_prot_match(0x0), m_fasta(0x0)
  {
    m_name = "ExecMetaAssembly";
    m_type = "ExecMetaAssembly";
  }

  // -------------------------------------------------------------------------

  ExecMetaAssembly::ExecMetaAssembly(const ParameterList & inputParams,
                                     SpecSet * spectra,
                                     SpectrumPairSet * spectrumPairs,
                                     abinfo_t * inputAbruijn,
                                     SpecSet * inputStars,
                                     Clusters * outputContigShifts,
                                     abinfo_t * outputAbruijn) :
    ExecBase(inputParams), m_spectra(spectra), m_spectrumPairs(spectrumPairs),
        m_inputAbruijn(inputAbruijn), m_stars(inputStars), ownInput(false),
        m_contigShifts(outputContigShifts), m_outputAbruijn(outputAbruijn),
        ownOutput(false), m_overlaps(0x0), m_prot_match(0x0), m_fasta(0x0)
  {
    m_name = "ExecMetaAssembly";
    m_type = "ExecMetaAssembly";
  }

  ExecMetaAssembly::ExecMetaAssembly(const ParameterList & inputParams,
                                     SpecSet * spectra,
                                     SpectrumPairSet * spectrumPairs,
                                     abinfo_t * inputAbruijn,
                                     SpecSet * inputStars) :
    ExecBase(inputParams), m_spectra(spectra), m_spectrumPairs(spectrumPairs),
        m_inputAbruijn(inputAbruijn), m_stars(inputStars), ownInput(false),
        m_contigShifts(0x0), m_outputAbruijn(0x0), ownOutput(true),
        m_overlaps(0x0), m_prot_match(0x0), m_fasta(0x0)
  {
    m_name = "ExecMetaAssembly";
    m_type = "ExecMetaAssembly";
  }

  // -------------------------------------------------------------------------

  ExecMetaAssembly::~ExecMetaAssembly(void)
  {
    if (ownInput)
    {
      if (m_spectra)
        delete m_spectra;
      if (m_spectrumPairs)
        delete m_spectrumPairs;
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

    if (m_overlaps)
    {
      delete m_overlaps;
    }

    if (m_prot_match)
    {
      delete m_prot_match;
    }

    if (m_fasta)
    {
      delete m_fasta;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecMetaAssembly::clone(const ParameterList & inputParams) const
  {
    return new ExecMetaAssembly(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::invoke(void)
  {
  
    bool noOp = m_params.getValueBool("META_ASSEMBLY_NO_OP", false);
    if (noOp) {
      DEBUG_MSG("META_ASSEMBLY_NO_OP");
      // Everything gets handled in the saveOutputData()
      return true;
    } else {
      // Clear these in case they were set by ProteoSAFe to the original
      m_contigShifts->resize(0);
    }
  
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM");

    float minScore = (float)m_params.getValueDouble("MIN_METACONTIG_SCORE");
    int minCompSize = m_params.getValueInt("MIN_METACONTIG_SIZE");
    unsigned int minMatchedPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS");
    bool parallelPaths = m_params.getValueInt("PARALLEL_PATHS", 0) > 0;

    m_stars->addZPMpeaks(peakTol, 0, true);
    m_stars->setPeakTolerance(pmTol, false);
    m_spectra->setPeakTolerance(pmTol, false);

    DEBUG_VAR(m_spectra->size());
    DEBUG_VAR(m_spectrumPairs->size());

    list<SpectrumPair> filteredPairs;
    unsigned int idxUse = 0;

    for (unsigned int i = 0; i < m_spectrumPairs->size(); i++)
    {
      if (min((*m_spectrumPairs)[i].score1, (*m_spectrumPairs)[i].score2)
          < minScore)
      {
        continue;
      }
      PRMAlignment alignmentObj;

      if ((*m_spectrumPairs)[i].spec1 >= m_spectra->size()
          || (*m_spectrumPairs)[i].spec2 >= m_spectra->size())
      {
        ERROR_MSG("Alignment index is out of range!!");
        DEBUG_VAR((*m_spectrumPairs)[i].spec1);
        DEBUG_VAR((*m_spectrumPairs)[i].spec2);
        DEBUG_VAR(m_spectra->size());
        return false;
      }
      alignmentObj.setSpec1(&(*m_spectra)[(*m_spectrumPairs)[i].spec1]);
      alignmentObj.setSpec2(&(*m_spectra)[(*m_spectrumPairs)[i].spec2]);
      pair<int, pair<float, float> > shiftScore =
          alignmentObj.getShiftScore((*m_spectrumPairs)[i].shift1, peakTol, 0);
      if (shiftScore.first < minMatchedPeaks)
      {
        continue;
      }
      filteredPairs.push_back((*m_spectrumPairs)[i]);
    }

    DEBUG_MSG("Retained " << filteredPairs.size() << " of " << m_spectrumPairs->size() << " input pairs");

    /*
     // Initialize overlap graph
     ContigNetwork overlapGraph(*m_spectra, filteredPairs, pmTol, pmTol);

     DEBUG_TRACE;

     // assemble meta-contigs
     int numMerged = overlapGraph.assembleIteratively(minScore, minMatchedPeaks);
     */

    CombineContigsParams params;
    params.contigs = m_spectra;
    params.star_spectra = m_stars;
    params.in_abinfo = m_inputAbruijn;
    params.contig_alignments = &filteredPairs;

    params.contig_peak_idx = m_overlaps;
    params.contig_prot_idx = m_prot_match;
    params.proteins = m_fasta;

    params.parent_mass_tol = m_params.getValueFloat("TOLERANCE_PM", 0.5);
    params.peak_tol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);
    params.resolution = m_params.getValueFloat("RESOLUTION", 0.01);
    params.min_component_size = minCompSize;
    params.min_matched_peaks_contig_align = minMatchedPeaks;
    params.contig_overlap_score = minScore;

    set<int> protein_idx_ignore;
    list<string> strIdxs;

    if (!splitText(m_params.getValue("CONTAMINANT_PROTEINS", "").c_str(),
                   strIdxs,
                   ":"))
      return false;

    for (list<string>::iterator strIt = strIdxs.begin(); strIt != strIdxs.end(); strIt++)
      protein_idx_ignore.insert(atoi((*strIt).c_str()));

    params.protein_idx_ignore = &protein_idx_ignore;

    CombineContigs comp;
    comp.construct(&params);
    comp.combineEdges(0);

    *m_outputAbruijn = params.out_abinfo;

    m_contigShifts->initializeFromAbinfo(params.consensus,
                                         *m_outputAbruijn,
                                         *m_stars,
                                         params.parent_mass_tol,
                                         params.parent_mass_tol);

    if (parallelPaths)
    {
      ParameterList pParams(m_params);
      pParams.setValue("ENFORCE_B_ENDPTS", "1");
      pParams.setValue("REVERSE_STARS", "1");
      pParams.setValue("PATH_EXPAND_LIMIT", "2");

      Clusters newShifts;
      newShifts = *m_contigShifts;
      abinfo_t newAbruijn;

      DEBUG_MSG("Invoking ParallelAssembly ...");
      ExecParallelAssembly pAssembly(pParams,
                                     m_outputAbruijn,
                                     m_stars,
                                     &newShifts,
                                     &newAbruijn);

      if (!pAssembly.invoke())
      {
        ERROR_MSG("Failed to invoke ParallelAssembly!!!");
        abort();
      }

      m_outputAbruijn->operator =(newAbruijn);
      m_contigShifts->operator =(newShifts);
    }

    for (unsigned int i = 0; i < m_contigShifts->consensus.size(); i++)
    {
      m_contigShifts->consensus[i].scan = i + 1;
    }

    DEBUG_VAR(m_outputAbruijn->size());

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::loadInputData(void)
  {
    if (ownInput)
    {
      if (!m_spectra)
      {
        m_spectra = new SpecSet();
      }
      if (!m_spectrumPairs)
      {
        m_spectrumPairs = new SpectrumPairSet();
      }
      if (!m_inputAbruijn)
      {
        m_inputAbruijn = new abinfo_t;
      }
      if (!m_stars)
      {
        m_stars = new SpecSet();
      }
    }
    m_spectra->resize(0);
    m_spectrumPairs->resize(0);
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
                                               m_params.getValue("INPUT_SPECTRA"),
                                               m_spectra) or m_spectra->size()
        == 0)
    {
      ERROR_MSG("Error reading input contig spectra from "<< m_params.getValue("INPUT_SPECTRA"));
      return false;
    }

    if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue("EXE_DIR"),
                                               m_params.getValue("INPUT_STARS"),
                                               m_stars) or m_stars->size() == 0)
    {
      ERROR_MSG("Error reading input star spectra from "<< m_params.getValue("INPUT_STARS"));
      return false;
    }

    if (!m_spectrumPairs->loadFromBinaryFile(m_params.getValue("INPUT_CONTIG_ALIGNS")))
    {
      ERROR_MSG("Error reading contig pairs from "<<m_params.getValue("INPUT_CONTIG_ALIGNS"));
      return false;
    }

    if (!Load_abinfo(m_params.getValue("INPUT_ABINFO").c_str(), *m_inputAbruijn))
    {
      ERROR_MSG("Error reading abinfo from \'"<<m_params.getValue("INPUT_ABINFO")<<"\'");
      return false;
    }

    if (m_params.exists("INPUT_MIDX") && m_params.exists("INPUT_MP")
        && m_params.exists("INPUT_FASTA"))
    {
      if (!m_overlaps)
      {
        m_overlaps = new SpecSet;
      }

      if (!m_prot_match)
      {
        m_prot_match = new vector<vector<int> > ;
      }

      if (!m_fasta)
      {
        m_fasta = new DB_fasta;
      }

      if (!ExecMergeConvert::loadSpecset(m_params.getValue("EXE_DIR"),
                                         m_params.getValue("INPUT_MIDX"),
                                         m_overlaps))
      {
        return false;
      }

      DEBUG_MSG("Loading matched protein idxs from " << m_params.getValue("INPUT_MP") << " ...");
      if (!Load_binArray(m_params.getValue("INPUT_MP").c_str(), *m_prot_match))
      {
        ERROR_MSG("Failed to load " << m_params.getValue("INPUT_MP"));
        return false;
      }

      DEBUG_MSG("Loading fasta proteins from " << m_params.getValue("INPUT_FASTA") << " ...");
      if (!m_fasta->Load(m_params.getValue("INPUT_FASTA").c_str()))
      {
        ERROR_MSG("Failed to load " << m_params.getValue("INPUT_FASTA"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::saveOutputData(void)
  {
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }

    bool noOp = m_params.getValueBool("META_ASSEMBLY_NO_OP", false);
  
    if (m_params.exists("OUTPUT_CLUSTERS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_CLUSTERS");
      DEBUG_MSG("Saving clusters to: [" << fileName << "]");
      // In the case of NO-OP the m_contigShifts should be the original
      if (m_contigShifts->Save(fileName.c_str()) != 0) {
        ERROR_MSG("Failed to save clusters to [" << fileName << "[");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_SPECTRA")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_SPECTRA");
      DEBUG_MSG("Saving spectra to: [" << fileName << "]");
      if (noOp) {
        if (!ExecMergeConvert::saveSpecsetMultiple(fileName,
                                                   m_spectra)) {
          ERROR_MSG("Failed to save meta-contig spectra to [" << fileName << "]");
          return false;
        }
      } else {
        if (!ExecMergeConvert::saveSpecsetMultiple(fileName,
                                                   &m_contigShifts->consensus)) {
          ERROR_MSG("Failed to save meta-contig spectra to [" << fileName << "]");
          return false;
        }
      }
    }

    if (m_params.exists("OUTPUT_ABINFO")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_ABINFO");
      DEBUG_MSG("Saving abinfo to: [" << fileName << "]");
      if (noOp) {
        if (!Save_abinfo_v1_0(fileName.c_str(), *m_inputAbruijn)) {
          ERROR_MSG("Failed to save abinfo to [" << fileName << "]");
          return false;
        }
      } else {
        if (!Save_abinfo_v1_0(fileName.c_str(), *m_outputAbruijn)) {
          ERROR_MSG("Failed to save abinfo to [" << fileName << "]");
          return false;
        }
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecMetaAssembly::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecMetaAssembly::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("MIN_METACONTIG_SIZE");
    VALIDATE_PARAM_EXIST("MIN_METACONTIG_SCORE");
    VALIDATE_PARAM_EXIST("MIN_MATCHED_PEAKS");
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}
