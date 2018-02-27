// Headed Include
#include "ExecFilterStarPairs.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "SetMerger.h"
#include "SpectralPairs.h"

// External Includes
//#include "abruijn.h"
#include "alignment_scoring.h"
//#include "filters.h"
//#include "graph.h"

// System Includes
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
using namespace specnets;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecFilterStarPairs::ExecFilterStarPairs(void) :
    m_inputPairs(0x0), m_starSpectra(0x0), ownInput(true), m_ratios(0x0),
        m_matchedPeaks(0x0), ownOutput(true)
  {
    m_name = "ExecFilterStarPairs";
    m_type = "ExecFilterStarPairs";
  }

  // -------------------------------------------------------------------------
  ExecFilterStarPairs::ExecFilterStarPairs(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputPairs(0x0), m_starSpectra(0x0),
        ownInput(true), m_ratios(0x0), m_matchedPeaks(0x0), ownOutput(true)
  {
    m_name = "ExecFilterStarPairs";
    m_type = "ExecFilterStarPairs";
  }

  SpectrumPairSet * m_inputPairs; //! Input partial alignments (or almost same peptides)
  SpecSet * m_starSpectra; //! Input spectra

  // -------------------------------------------------------------------------
  ExecFilterStarPairs::ExecFilterStarPairs(const ParameterList & inputParams,
                                           SpectrumPairSet * inputPairs,
                                           SpecSet * starSpectra,
                                           vector<vector<float> > * ratios,
                                           SpecSet * matchedPeaks) :
    ExecBase(inputParams), m_inputPairs(inputPairs), m_starSpectra(starSpectra), 
        ownInput(false), 
        m_ratios(ratios), m_matchedPeaks(matchedPeaks), ownOutput(false)
  {
    m_name = "ExecFilterStarPairs";
    m_type = "ExecFilterStarPairs";
  }

  // -------------------------------------------------------------------------
  ExecFilterStarPairs::~ExecFilterStarPairs(void)
  {
    if (ownInput)
    {
      delete m_inputPairs;
      delete m_starSpectra;
    }

    if (ownOutput)
    {
      delete m_ratios;
      delete m_matchedPeaks;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFilterStarPairs::clone(const ParameterList & inputParams) const
  {
    return new ExecFilterStarPairs(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::invoke(void)
  {
    DEBUG_TRACE;
    if (m_params.getValueBool("DEBUG_PARAMS")) {
      stringstream aux;
      m_params.print(aux);
      DEBUG_MSG(aux.str());
    }

    DEBUG_VAR(m_starSpectra->size());
    DEBUG_VAR(m_inputPairs->size());

    float pmTol   = m_params.getValueFloat("TOLERANCE_PM", 1.0);
    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);
    bool filterByPeakIntensity = m_params.getValueBool("FILTERSTARPAIRS_BY_PEAKINTENSITY", true);
   
    if (filterByPeakIntensity) {
      m_matchedPeaks->resize(2 * m_inputPairs->size());
      m_ratios->resize(m_inputPairs->size());
      for (unsigned int i = 0; i < m_ratios->size(); i++)
      {
        (*m_ratios)[i].resize(3);
        (*m_ratios)[i][0] = 0;
        (*m_ratios)[i][1] = 0;
        (*m_ratios)[i][2] = 0;
      }

      DEBUG_TRACE;

      // Separate m_inputPairs into isolated networks (connected components)
      SetMerger components(m_starSpectra->size());
      SpectrumPairSet foo;
      components.createSets(m_starSpectra->size(),
                            2,
                            foo,
                            *m_inputPairs);

      float penalty_ptm             = m_params.getValueFloat("PENALTY_PTM");
      float penalty_sameVert        = m_params.getValueFloat("PENALTY_SAME_VERTEX");
      int maxAAjump                 = m_params.getValueInt("MAX_AA_JUMP", 0);
      float maxModMass              = m_params.getValueFloat("MAX_MOD_MASS", 100.0);
      float minRatio                = m_params.getValueFloat("MIN_RATIO", -1);
      unsigned int minMatchedPeaks  = m_params.getValueInt("MIN_MATCHED_PEAKS", 0);
      int specType                  = m_params.getValueInt("SPEC_TYPE_MSMS", 0);
      float ionOffset               = specType ? AAJumps::massHion : 0;
      
      DEBUG_VAR(penalty_ptm);
      DEBUG_VAR(penalty_sameVert);
      DEBUG_VAR(maxAAjump);
      DEBUG_VAR(maxModMass);
      DEBUG_VAR(peakTol);
      DEBUG_VAR(minRatio);
      DEBUG_VAR(minMatchedPeaks);
      DEBUG_VAR(specType);
      DEBUG_VAR(ionOffset);

      // Allow overide of the other parameters with FilterStarPairs specific params
      if (m_params.exists("FILTERSTARPAIRS_MIN_ALIGN_RATIO")) {
        minRatio = m_params.getValueFloat("FILTERSTARPAIRS_MIN_ALIGN_RATIO", -1);
        DEBUG_VAR(minRatio);
      }
      if (m_params.exists("FILTERSTARPAIRS_MIN_MATCHED_PEAKS")) {
        minMatchedPeaks = m_params.getValueFloat("FILTERSTARPAIRS_MIN_MATCHED_PEAKS", 0);
        DEBUG_VAR(minMatchedPeaks);
      }
  
      vector<SpectrumPairSet> cAligns;
      vector<vector<int> > cAligns_idx; // To allow going back from components to original order
      components.splitAligns(*m_inputPairs, cAligns, cAligns_idx);

      DEBUG_TRACE;

      vector<float> maxSpecScores(m_starSpectra->size());
      vector<bool> specFlipped(m_starSpectra->size());
      for (unsigned int i = 0; i < m_starSpectra->size(); i++)
      {
        maxSpecScores[i] = 0;
        for (unsigned int j = 0; j < (*m_starSpectra)[i].size(); j++)
        {
          maxSpecScores[i] += (*m_starSpectra)[i][j][1];
        }
        (*m_starSpectra)[i].addZPMpeaks(peakTol, ionOffset, true);
        specFlipped[i] = false;
      }

      DEBUG_TRACE;

      vector<vector<TwoValues<int> > > matches;
      vector<vector<float> > cRatios;
      vector<float> modPos;
      SpecSet matchedPeaks(2 * m_inputPairs->size());
      
      DEBUG_VAR(cAligns.size());
      
      for (unsigned int cIdx = 0; cIdx < cAligns.size(); cIdx++)
      {
        //DEBUG_VAR(cIdx);
        
        SplitPairs(*m_starSpectra,
                   cAligns[cIdx],
                   peakTol,
                   pmTol,
                   maxAAjump,
                   maxModMass,
                   penalty_sameVert,
                   penalty_ptm,
                   matches,
                   specFlipped,
                   modPos,
                   minMatchedPeaks,
                   0,
                   false,
                   false,
                   &cRatios);

        // Copy alignment statistics to the corresponding global pair indices (as opposed to per-connected-component pair indices)
        for (unsigned int pairIdx = 0; pairIdx < cAligns[cIdx].size(); pairIdx++)
        {
          //DEBUG_VAR(pairIdx);
          
          for (unsigned int i = 0; i < 3; i++)
            (*m_ratios)[cAligns_idx[cIdx][pairIdx]][i] = cRatios[pairIdx][i];

          unsigned int mpIdx = 2 * cAligns_idx[cIdx][pairIdx];
          matchedPeaks[mpIdx].resize(matches[pairIdx].size());
          matchedPeaks[mpIdx + 1].resize(matches[pairIdx].size());
          unsigned int s1 = cAligns[cIdx][pairIdx].spec1, s2 =
              cAligns[cIdx][pairIdx].spec2;

  #if 0
          if(s1==0 and s2==1)
          {
            ERROR_MSG(" * pairIdx = "<<cAligns_idx[cIdx][pairIdx]<<", mpIdx = "<<mpIdx);
          }
  #endif

          for (unsigned int matchIdx = 0; matchIdx < matches[pairIdx].size(); matchIdx++)
          {
            matchedPeaks[mpIdx][matchIdx]
                = (*m_starSpectra)[s1][matches[pairIdx][matchIdx][0]];
            matchedPeaks[mpIdx + 1][matchIdx]
                = (*m_starSpectra)[s2][matches[pairIdx][matchIdx][1]];
  #if 0

            if(s1==0 and s2==1)
            {
              ERROR_MSG(" * "<<matches[pairIdx][matchIdx][0]<<"/("<<matchedPeaks[mpIdx][matchIdx][0]<<","<<matchedPeaks[mpIdx][matchIdx][1]<<") - "
              <<matches[pairIdx][matchIdx][1]<<"/("<<matchedPeaks[mpIdx+1][matchIdx][0]<<","<<matchedPeaks[mpIdx+1][matchIdx][1]<<")");
            }
  #endif

          }
        }

      //DEBUG_TRACE;

  #if 0
        for(unsigned int pairIdx=0; pairIdx<cAligns[cIdx].size(); pairIdx++)
        {
          float score1=0, score2=0; int s1=cAligns[cIdx][pairIdx].spec1, s2=cAligns[cIdx][pairIdx].spec2;
          for(unsigned int i=0; i<matches[pairIdx].size(); i++)
          {
            if(matches[pairIdx][i][0]>=0) score1+=(*m_starSpectra)[s1][matches[pairIdx][i][0]][1];
            if(matches[pairIdx][i][1]>=0) score2+=(*m_starSpectra)[s2][matches[pairIdx][i][1]][1];
          }
          ratios[cAligns_idx[cIdx][pairIdx]][0] = score1/maxSpecScores[s1];
          ratios[cAligns_idx[cIdx][pairIdx]][1] = score2/maxSpecScores[s2];
          ratios[cAligns_idx[cIdx][pairIdx]][2] = matches[pairIdx].size();
        }
  #endif

      }

      DEBUG_TRACE;

      // Keep only the spectral pairs matching >= MIN_RATIO percent intensity in both aligned spectra
      unsigned int keptPairs = 0;
      SpectrumPair curPair;
      for (unsigned int idxRatio = 0; idxRatio < m_ratios->size(); idxRatio++)
      {
        //DEBUG_VAR(idxRatio);
        
        curPair = (*m_inputPairs)[idxRatio];
        unsigned int extraMatchPeaks = 0; // ASP pairs always match either/both PRMs at 0/18 or/and PM-19/PM-1 so these have higher numbers of peaks that need to match
        if (abs(curPair.shift1) <= peakTol)
          extraMatchPeaks++; // matching 0/18 does not count towards achieving minMatchedPeaks
        if (abs(curPair.shift2) <= peakTol)
          extraMatchPeaks++; // matching PM-19/PM-1 does not count towards achieving minMatchedPeaks

  #if 0
        if((*m_inputPairs)[idxRatio].spec1==9630 and (*m_inputPairs)[idxRatio].spec2==9642)
        {
          ERROR_MSG(" ----- (*m_inputPairs)[idxRatio].shift1 = "<<(*m_inputPairs)[idxRatio].shift1<<", peakTol = "<<peakTol<<", (abs((*m_inputPairs)[idxRatio].shift1)<=peakTol) = "<<(abs((*m_inputPairs)[idxRatio].shift1)<=peakTol));
          ERROR_MSG(" ----- idxRatio = "<<idxRatio<<", keptPairs = "<<keptPairs<<", ratios[idxRatio][2] = "<<ratios[idxRatio][2]<<", minMatchedPeaks = "<<minMatchedPeaks<<", extraMatchPeaks = "<<extraMatchPeaks<<", pass = "<<(ratios[idxRatio][2]>=minMatchedPeaks+extraMatchPeaks));
        }
  #endif

        if ((*m_ratios)[idxRatio][0] >= minRatio && (*m_ratios)[idxRatio][1]
            >= minRatio && (*m_ratios)[idxRatio][2] >= minMatchedPeaks
            + extraMatchPeaks)
        {
          (*m_inputPairs)[keptPairs] = (*m_inputPairs)[idxRatio];
          for (unsigned int pivot = 0; pivot < 3; pivot++)
            (*m_ratios)[keptPairs][pivot] = (*m_ratios)[idxRatio][pivot];
          keptPairs++;
        }
      }

      DEBUG_TRACE;
    
      m_inputPairs->resize(keptPairs);
      m_ratios->resize(keptPairs);

      DEBUG_MSG("Retained " << keptPairs << " pairs");
      
    } // if (filterByPeakIntensity) {

    DEBUG_VAR(m_inputPairs->size());

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::loadInputData(void)
  {
    DEBUG_TRACE;

    ownInput = true;
    m_inputPairs = new SpectrumPairSet();
    m_starSpectra = new SpecSet();
    ownOutput = true;
    m_ratios = new vector<vector<float> > ();
    m_matchedPeaks = new SpecSet();

    string starsFilename = m_params.getValue("INPUT_SPECS_PKLBIN");
    DEBUG_MSG("Loading stars file [" << starsFilename << "]...");
    if (!m_starSpectra->loadPklBin(starsFilename.c_str()))
    {
      ERROR_MSG("Could not load " << starsFilename);
      return false;
    }

    DEBUG_MSG("Loading stars complete. Num stars [" << m_starSpectra->size() << "]");

    string alignFilename = m_params.getValue("INPUT_ALIGNS");
    DEBUG_MSG("Loading aligns file [" << alignFilename << "]...");
    if (!m_inputPairs->loadFromBinaryFile(alignFilename))
    {
      ERROR_MSG("Could not load " << alignFilename);
      return false;
    }
    
    DEBUG_MSG("Loading alignments complete. Num pairs [" << m_inputPairs->size() << "]");

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::saveInputData(std::vector<std::string> & filenames)
  {
    DEBUG_TRACE;

    std::string pairsFilename = getName() + "_pairs.bin";
    m_params.setValue("INPUT_ALIGNS", pairsFilename);
    if (!fileExists(pairsFilename))
    {
      m_inputPairs->saveToBinaryFile(pairsFilename);
    }

    DEBUG_VAR(m_inputPairs->size());

    DEBUG_TRACE;

    string starsFilename = getName() + "_stars.pklbin";
    m_params.setValue("INPUT_SPECS_PKLBIN", starsFilename);
    if (!fileExists(starsFilename))
    {
      m_starSpectra->savePklBin(starsFilename.c_str());
    }

    DEBUG_VAR(m_starSpectra->size());

    // Have to set up the output files also so the params will be correct on reload
    m_params.setValue("OUTPUT_ALIGNS", getName() + "_pairs_stars.bin");

    //SpecSet m_inputSpectra; // the input spectra
    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(pairsFilename);
    filenames.push_back(starsFilename);

    DEBUG_TRACE;

    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::saveOutputData(void)
  {
    m_inputPairs->saveToBinaryFile(m_params.getValue("OUTPUT_ALIGNS"));

    if (m_params.exists("OUTPUT_RATIOS"))
    {
      Save_binArray(m_params.getValue("OUTPUT_RATIOS").c_str(), *m_ratios);
    }
    if (m_params.exists("OUTPUT_SPECS"))
    {
      m_matchedPeaks->savePklBin(m_params.getValue("OUTPUT_SPECS").c_str());
    }
    else if (m_params.exists("OUTPUT_SPECS_PKLBIN"))
    {
      m_matchedPeaks->savePklBin(m_params.getValue("OUTPUT_SPECS_PKLBIN").c_str());
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::loadOutputData(void)
  {
    m_inputPairs->loadFromBinaryFile(m_params.getValue("OUTPUT_ALIGNS"));
    DEBUG_VAR(m_inputPairs->size());

    if (m_params.exists("OUTPUT_RATIOS"))
    {
      if (Load_binArray(m_params.getValue("OUTPUT_RATIOS").c_str(), *m_ratios) <= 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_RATIOS"));
        return false;
      }
      DEBUG_VAR(m_ratios->size());
    }
    if (m_params.exists("OUTPUT_SPECS"))
    {
      if (!m_matchedPeaks->loadPklBin(m_params.getValue("OUTPUT_SPECS").c_str()))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_SPECS"));
        return false;
      }
      DEBUG_VAR(m_matchedPeaks->size());
    }

    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFilterStarPairs::split(int numSplit)
	{
		m_subModules.resize(0);
		return m_subModules;
	}

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterStarPairs::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("PENALTY_PTM");
    VALIDATE_PARAM_EXIST("PENALTY_SAME_VERTEX");
//    VALIDATE_PARAM_EXIST("MAX_AA_JUMP");
//    VALIDATE_PARAM_EXIST("MAX_MOD_MASS");
//    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
//    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
//    VALIDATE_PARAM_EXIST("MIN_RATIO");
//    VALIDATE_PARAM_EXIST("MIN_MATCHED_PEAKS");
//    VALIDATE_PARAM_EXIST("SPEC_TYPE_MSMS");

    m_isValid = true;
    return true;
  }

} // namespace specnets
