// Header Includes
#include "ExecAlignment.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes
#include "alignment_scoring.h"
#include "SpectralPairs.h"
using namespace specnets;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecAlignment::ExecAlignment(void) :
    m_inputSpectra(0x0), m_inputPairs(0x0), ownInput(true),
        m_pairAlignments(0x0), m_starSpectraOnly(0x0), m_starSpectra(0x0),
        m_alignedSpectra(0x0), ownOutput(true)
  {
    m_name = "ExecAlignment";
    m_type = "ExecAlignment";
  }

  // -------------------------------------------------------------------------
  ExecAlignment::ExecAlignment(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_inputPairs(0x0),
        ownInput(true), m_pairAlignments(0x0), m_starSpectraOnly(0x0),
        m_starSpectra(0x0), m_alignedSpectra(0x0), ownOutput(true)
  {
    m_name = "ExecAlignment";
    m_type = "ExecAlignment";
  }

  // -------------------------------------------------------------------------
  ExecAlignment::ExecAlignment(const ParameterList & inputParams,
                               SpecSet * inputSpectra,
                               SpectrumPairSet * inputPairs,
                               SpecSet * outputPairAlignments,
                               SpecSet * outputStarSpectraOnly,
                               SpecSet * outputStarSpectra,
                               vector<unsigned int> * outputAlignedSpectra) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra),
        m_inputPairs(inputPairs), ownInput(false),
        m_pairAlignments(outputPairAlignments),
        m_starSpectraOnly(outputStarSpectraOnly),
        m_starSpectra(outputStarSpectra),
        m_alignedSpectra(outputAlignedSpectra), ownOutput(false)
  {
    m_name = "ExecAlignment";
    m_type = "ExecAlignment";
  }

  // -------------------------------------------------------------------------
  ExecAlignment::~ExecAlignment(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
      delete m_inputPairs;
    }
    if (ownOutput)
    {
      delete m_pairAlignments;
      delete m_starSpectraOnly;
      delete m_starSpectra;
      delete m_alignedSpectra;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecAlignment::clone(const ParameterList & inputParams) const
  {
    return new ExecAlignment(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecAlignment::invoke(void)
  {
    DEBUG_TRACE;
    if (m_params.getValueBool("DEBUG_PARAMS")) {
      stringstream aux;
      m_params.print(aux);
      DEBUG_MSG(aux.str());
    }

    float penalty_sameVert = m_params.getValueFloat("PENALTY_SAME_VERTEX");
    DEBUG_VAR(penalty_sameVert);

    //	int aaDiff = m_params.getValueInt(paramStrs[4]);
    //	double minShift = m_params.getValueDouble(paramStrs[5]);
    //	double maxShift = m_params.getValueDouble(paramStrs[6]);

    float penalty_ptm = m_params.getValueFloat("PENALTY_PTM", 0.0);
    DEBUG_VAR(penalty_ptm);

    // If specified, indicates the minimum number peaks on each side of a modification to _not_ make it a terminal modification
    //   This is encoded by computing the average prm score per spectrum (at least 57 Da apart) and multiplying it by PENALTY_PTM_PEAKS to obtain the penalty (per spectrum)
    float penalty_ptm_peaks = m_params.getValueFloat("PENALTY_PTM_PEAKS", -1.0);
    DEBUG_VAR(penalty_ptm_peaks);

    int startBaseIdx = m_params.getValueInt("IDX_START", 0);
    int endBaseIdx = m_params.getValueInt("IDX_END", -1);
    int startStarIdx = m_params.getValueInt("STARS_START", 0);
    int endStarIdx = m_params.getValueInt("STARS_END", 0);
    int maxAAjump = m_params.getValueInt("MAX_AA_JUMP", 0);
    DEBUG_VAR(startBaseIdx);
    DEBUG_VAR(endBaseIdx);
    DEBUG_VAR(startStarIdx);
    DEBUG_VAR(endStarIdx);
    DEBUG_VAR(maxAAjump);

    //    vector<int> baseSpecIdx(endBaseIdx-startBaseIdx+1);  for (int i=0; i<endBaseIdx-startBaseIdx+1; i++) baseSpecIdx[i]=startBaseIdx+i;
    //	float minRatio = m_params.exists("MIN_RATIO") ? (float) m_params.getValueDouble("MIN_RATIO") : 0;

    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);
    float pmTol = m_params.getValueFloat("TOLERANCE_PM", 1);
    float resolution = m_params.getValueFloat("RESOLUTION", 0.1);
    int specType = m_params.getValueInt("SPEC_TYPE_MSMS", 0);
    float ionOffset = specType ? AAJumps::massHion : 0;
    DEBUG_VAR(peakTol);
    DEBUG_VAR(pmTol);
    DEBUG_VAR(resolution);
    DEBUG_VAR(specType);
    DEBUG_VAR(ionOffset);

    bool alignPA = m_params.getValueBool("PARTIAL_OVERLAPS", false);
    DEBUG_VAR(alignPA);

    DEBUG_TRACE;
    vector<float> penalties(m_inputSpectra->size());
    // Indices of matched peaks for partial overlaps (default sizes should avoid most resize() operations)
    vector<int> idx1dp(2048), idx2dp(2048);

    DEBUG_TRACE;
    float score;
    for (unsigned int specIdx = 0; specIdx < m_inputSpectra->size(); specIdx++) {
      (*m_inputSpectra)[specIdx].parentMass += 2 * ionOffset; // Quick fix to support alignment of MS/MS spectra
      if (penalty_ptm_peaks > 0) {
        score = ScoreOverlap6((*m_inputSpectra)[specIdx],
                              (*m_inputSpectra)[specIdx],
                              0,
                              peakTol,
                              idx1dp,
                              idx2dp,
                              AAJumps::minAAmass);
        DEBUG_VAR(score);
        if (idx1dp.size() == 0) {
          penalties[specIdx] = 0;
        } else {
          penalties[specIdx] = -penalty_ptm_peaks * score
              / ((float)idx1dp.size());
        }
      }
      (*m_inputSpectra)[specIdx].addZPMpeaks(peakTol, ionOffset, true);
    }

    DEBUG_TRACE;
    if (m_inputPairs->size() == 0) {
      DEBUG_MSG("No work to do... exiting.");
      return true;
    }

    if (endBaseIdx <= 0) {
      endBaseIdx = m_inputPairs->size() - 1;
    }

    int idxResultSpecs = 0;
    int spec1;
    int spec2;
    float shift;
    float curPenalty;

    //---------------------------------------------------------------------
    // Only perform this first section if we don't already have the results
    //---------------------------------------------------------------------
    if (!m_params.exists("INPUT_SPECS_PAIRS")) {
      DEBUG_TRACE;
      m_pairAlignments->resize(2 * m_inputPairs->size()); // 2 results per pair of spectra
      vector<float> modPositions;

      DEBUG_TRACE;
      modPositions.resize(m_inputPairs->size()); // Keeps track of the mass where the mod was placed

      DEBUG_MSG("Computing pairs from " << startBaseIdx << " to " << endBaseIdx << "...");

      for (unsigned int idxPair = startBaseIdx; 
           idxPair <= endBaseIdx; 
           idxPair++, idxResultSpecs += 2) {
           
        if (idxPair % 5000 == 0) {
          DEBUG_MSG("Processing pairs " << idxPair << " -> " << endBaseIdx << 
                    ", " << 100.0 * ((float)idxPair - (float)startBaseIdx) / 
                    ((float)endBaseIdx - (float)startBaseIdx) << "% completed");
        }

        int spec1 = (*m_inputPairs)[idxPair].spec1;
        int spec2 = (*m_inputPairs)[idxPair].spec2;
        float shift = (*m_inputPairs)[idxPair].shift1;

        (*m_pairAlignments)[idxResultSpecs].copyNP((*m_inputSpectra)[spec1]);
        (*m_pairAlignments)[idxResultSpecs + 1].copyNP((*m_inputSpectra)[spec2]);
        if (penalty_ptm_peaks > 0) {
          curPenalty = penalties[spec1] + penalties[spec2];
        } else {
          curPenalty = penalty_ptm;
        }

        // Same-peptide pairs can't be used to compute spectral stars
        // 		if (abs(shift) <= 2 * pmTol &&
        //        abs((*m_inputSpectra)[spec2].parentMass + shift - (*m_inputSpectra)[spec1].parentMass) <= 2 * pmTol)
        //      continue;
        if (!alignPA || (abs(shift) <= 2 * pmTol)
            || (abs((*m_inputSpectra)[spec2].parentMass + shift
                - (*m_inputSpectra)[spec1].parentMass) <= 2 * pmTol)) {
          // Common start or common end
          //DEBUG_MSG("Current pair is (" << spec1 + 1 << ", " << spec2 + 1 << ")");
          modPositions[idxPair]
              = SpectrumAlignment(&(*m_inputSpectra)[spec1],
                            &(*m_inputSpectra)[spec2],
                            peakTol,
                            &(*m_pairAlignments)[idxResultSpecs],
                            &(*m_pairAlignments)[idxResultSpecs + 1],
                            maxAAjump,
                            penalty_sameVert,
                            curPenalty);
          //for (unsigned int matchIdx=0; matchIdx<(*m_pairAlignments)[idxResultSpecs].size(); matchIdx++)
          //	cerr<<" * ("<<(*m_pairAlignments)[idxResultSpecs][matchIdx][0]<<","<<(*m_pairAlignments)[idxResultSpecs][matchIdx][1]<<") - ("
          //		<<(*m_pairAlignments)[idxResultSpecs+1][matchIdx][0]<<","<<(*m_pairAlignments)[idxResultSpecs+1][matchIdx][1]<<")" << endl;
        } else {
          // Partial overlap
          ScoreOverlap6((*m_inputSpectra)[spec1],
                        (*m_inputSpectra)[spec2],
                        shift,
                        peakTol,
                        idx1dp,
                        idx2dp,
                        AAJumps::minAAmass);

          (*m_pairAlignments)[idxResultSpecs].resize(idx1dp.size());
          for (int idxMatch = 0; idxMatch < idx1dp.size(); idxMatch++) {
            (*m_pairAlignments)[idxResultSpecs][idxMatch]
                = (*m_inputSpectra)[spec1][idx1dp[idxMatch]];
          }
          (*m_pairAlignments)[idxResultSpecs + 1].resize(idx2dp.size());
          for (int idxMatch = 0; idxMatch < idx2dp.size(); idxMatch++) {
            (*m_pairAlignments)[idxResultSpecs + 1][idxMatch]
                = (*m_inputSpectra)[spec2][idx2dp[idxMatch]];
          }
          modPositions[idxPair] = 0;
        }
      }

      // Quick fix to support alignment of MS/MS spectra
      for (unsigned int specIdx = 0; specIdx < m_pairAlignments->size(); specIdx++) {
        (*m_pairAlignments)[specIdx].parentMass -= 2 * ionOffset;
      }

    } // if ( !m_params.exists("INPUT_SPECS_PAIRS") )

    DEBUG_TRACE;

    //---------------------------------------------------------------------
    // Second step: Can be run using saved intermediate results
    //---------------------------------------------------------------------
    if (m_params.exists("OUTPUT_STARS") || m_params.exists("OUTPUT_STARS_ALL")) {
      unsigned int numStars = 0;
      SpecSet curSpecs;
      vector<unsigned short> numPairs(m_inputSpectra->size());

      for (unsigned int specIdx = 0; specIdx < numPairs.size(); specIdx++) {
        numPairs[specIdx] = 0;
      }

      if (endStarIdx <= 0 || endStarIdx >= m_inputSpectra->size()) {
        endStarIdx = m_inputSpectra->size() - 1;
      }

      DEBUG_MSG("Computing stars with indices from " << startStarIdx << 
                " to " << endStarIdx << "...");

      for (unsigned int idxPair = startBaseIdx; idxPair <= endBaseIdx; idxPair++) {
        spec1 = (*m_inputPairs)[idxPair].spec1;
        spec2 = (*m_inputPairs)[idxPair].spec2;
        shift = (*m_inputPairs)[idxPair].shift1;

        if (abs(shift) <= 2 * pmTol && abs((*m_inputSpectra)[spec2].parentMass
            + shift - (*m_inputSpectra)[spec1].parentMass) <= 2 * pmTol)  {
          continue; // Same-peptide pairs can't be used to compute spectral stars
        }
        if ((spec1 < startStarIdx || spec1 > endStarIdx) && (spec2
            < startStarIdx || spec2 > endStarIdx))  {
          continue; // Star-spectrum will be computed by another instance
        }

        if (spec1 >= startStarIdx && spec1 <= endStarIdx) {
          if (numPairs[spec1] == 0) {
            numStars++;
          }
          numPairs[spec1]++;
        }
        if (spec2 >= startStarIdx && spec2 <= endStarIdx) {
          if (numPairs[spec2] == 0) {
            numStars++;
          }
          numPairs[spec2]++;
        }
      }

      m_starSpectraOnly->resize(numStars);
      m_alignedSpectra->resize(numStars);

      DEBUG_MSG(numStars << " stars in range...");

      unsigned int starIdx = 0;
      for (unsigned int specIdx = startStarIdx; specIdx <= endStarIdx; specIdx++) {
        if (numPairs[specIdx] == 0) {
          continue;
        }
        (*m_alignedSpectra)[starIdx] = specIdx;
        if (starIdx % 100 == 0) {
          DEBUG_MSG("Computing star spectrum " << starIdx + 1 << "/" << numStars);
        }

        // Find all spectral pairs for this spectrum
        curSpecs.resize(numPairs[specIdx]);
        numPairs[specIdx] = 0;
        idxResultSpecs = 0;
        for (unsigned int idxPair = startBaseIdx; idxPair <= endBaseIdx; idxPair++) {
          spec1 = (*m_inputPairs)[idxPair].spec1;
          spec2 = (*m_inputPairs)[idxPair].spec2;
          shift = (*m_inputPairs)[idxPair].shift1;

          // Same-peptide pairs can't be used to compute spectral stars
          //      if ( abs(shift) <= 2 * pmTol &&
          //           abs((*m_inputSpectra)[spec2].parentMass + shift - (*m_inputSpectra)[spec1].parentMass) <= 2 * pmTol)
          //        continue;

          if (spec1 == (int)specIdx) {
            curSpecs[numPairs[specIdx]++] = (*m_pairAlignments)[2 * idxPair];
          } else if (spec2 == (int)specIdx) {
            curSpecs[numPairs[specIdx]++]
                = (*m_pairAlignments)[2 * idxPair + 1];
          }
          if (numPairs[specIdx] == curSpecs.size()) {
            break;
          }
        }
        if (numPairs[specIdx] == 1) {
          (*m_starSpectraOnly)[starIdx] = curSpecs[0];
        } else {
          ComputeSpectralStars(curSpecs,
                               (*m_starSpectraOnly)[starIdx],
                               peakTol,
                               resolution);
        }
        starIdx++;
      }

      DEBUG_MSG("Done");
    }

    // Fill-in spectra that were not converted to spectral stars
    DEBUG_VAR(m_starSpectraOnly->size());
    m_starSpectra->resize(m_inputSpectra->size());
    DEBUG_VAR(m_starSpectra->size());

    // Copy the star spectra
    for (unsigned int starIdx = 0; starIdx < m_alignedSpectra->size(); starIdx++)
    {
      (*m_starSpectra)[(*m_alignedSpectra)[starIdx]]
          = (*m_starSpectraOnly)[starIdx];
    }
    // Copy all other (unchanged) spectra
    for (unsigned int n = 0; n < m_starSpectra->size(); n++) {
      if ((*m_starSpectra)[n].parentMass == 0.0) {
        (*m_starSpectra)[n] = (*m_inputSpectra)[n];
      }
      (*m_starSpectra)[n].copyNP((*m_inputSpectra)[n]);
    }

    DEBUG_VAR(m_pairAlignments->size());
    DEBUG_VAR(m_starSpectraOnly->size());
    DEBUG_VAR(m_starSpectra->size());
    DEBUG_VAR(m_alignedSpectra->size());

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecAlignment::loadInputData(void)
  {
    ownInput = true;
    m_inputSpectra = new SpecSet();
    m_inputPairs = new SpectrumPairSet();
    ownOutput = true;
    m_pairAlignments = new SpecSet();
    m_starSpectraOnly = new SpecSet();
    m_starSpectra = new SpecSet();
    m_alignedSpectra = new vector<unsigned int> ();

    //---------------------------------
    // Load spectrum data
    //---------------------------------
    if (m_params.exists("INPUT_SPECS"))
    {
      if (!m_inputSpectra->LoadSpecSet_pkl(m_params.getValue("INPUT_SPECS").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_SPECS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_SPECS_PKLBIN"))
    {
      if (!m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_SPECS_PKLBIN"));
        return false;
      }
    }
    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!");
      return false;
    }

    DEBUG_MSG("Loading specs complete. Num specs [" << m_inputSpectra->size() << "]");

    string alignFilename = m_params.getValue("INPUT_ALIGNS");
    if (!m_inputPairs->loadFromBinaryFile(alignFilename))
    {
      ERROR_MSG("Could not load " << alignFilename);
      return false;
    }
    
    DEBUG_MSG("Loading alignments complete. Num pairs [" << m_inputPairs->size() << "]");

    //--------------------------------------------
    // Load the intermediate results (if desired)
    //--------------------------------------------
    if (m_params.exists("INPUT_SPECS_PAIRS"))
    {
      string specPairsFilename = m_params.getValue("INPUT_SPECS_PAIRS");
      if (!m_pairAlignments->loadPklBin(specPairsFilename.c_str())
          || m_pairAlignments->size() != 2 * m_inputSpectra->size())
      {
        ERROR_MSG("Could not load " << specPairsFilename);
        return false;
      }
      
      DEBUG_MSG("Loading matched peaks complete. Num pairs : " << m_pairAlignments->size() / 2);
    }

    //---------------------------------
    // Load amino acid masses
    //---------------------------------
    AAJumps jumps(-1);
    if (m_params.exists("AMINO_ACID_MASSES"))
    {
      jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecAlignment::saveInputData(std::vector<std::string> & filenames)
  {
    std::string pairsFilename = getName() + "_pairs.bin";
    m_params.setValue("INPUT_ALIGNS", pairsFilename);
    if (!fileExists(pairsFilename))
    {
      if (!m_inputPairs->saveToBinaryFile(pairsFilename))
      {
        ERROR_MSG("Could not load " << pairsFilename);
        return false;
      }
    }

    //SpecSet m_inputSpectra; // the input spectra
    std::string spectraFilename = getName() + "_spectra.pklbin";
    if (!fileExists(spectraFilename))
    {
      if (!m_inputSpectra->savePklBin(spectraFilename.c_str()))
      {
        ERROR_MSG("Could not load " << spectraFilename);
        return false;
      }
    }
    m_params.setValue("INPUT_SPECS_PKLBIN", spectraFilename);

    // Have to set up the output files also so the params will be correct on reload
    m_params.setValue("OUTPUT_SPECS", getName() + "_pair_aligns.bin");
    m_params.setValue("OUTPUT_STARS", getName() + "_stars.bin");
    m_params.setValue("OUTPUT_STARS_ALL", getName() + "_stars_only.bin");
    m_params.setValue("OUTPUT_STARS_INDEX", getName() + "_stars_index.bin");

    //SpecSet m_inputSpectra; // the input spectra
    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(pairsFilename);
    filenames.push_back(spectraFilename);

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecAlignment::saveOutputData(void)
  {
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    DEBUG_VAR(dataDir);
    if (m_params.exists("OUTPUT_SPECS"))
    {
      DEBUG_VAR(m_pairAlignments->size());
      string fileName = dataDir + m_params.getValue("OUTPUT_SPECS");
      DEBUG_MSG("Saving spectra to: [" << fileName << "]");
      m_pairAlignments->savePklBin(fileName.c_str());
    }

    if (m_params.exists("OUTPUT_STARS"))
    {
      DEBUG_VAR(m_starSpectraOnly->size());
      string fileName = dataDir + m_params.getValue("OUTPUT_STARS");
      DEBUG_MSG("Saving stars to: [" << fileName << "]");
      m_starSpectraOnly->savePklBin(fileName.c_str());
    }

    if (m_params.exists("OUTPUT_STARS_ALL"))
    {
      DEBUG_VAR(m_starSpectra->size());
      string fileName = dataDir + m_params.getValue("OUTPUT_STARS_ALL");
      DEBUG_MSG("Saving stars all to: [" << fileName << "]");
      m_starSpectra->savePklBin(fileName.c_str());
    }

    if (m_params.exists("OUTPUT_STARS_INDEX"))
    {
      DEBUG_VAR(m_alignedSpectra->size());
      string fileName = dataDir + m_params.getValue("OUTPUT_STARS_INDEX");
      DEBUG_MSG("Saving stars index all to: [" << fileName << "]");
      Save_binArray(fileName.c_str(), *m_alignedSpectra);
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecAlignment::loadOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPECS"))
    {
      if (!m_pairAlignments->loadPklBin(m_params.getValue("OUTPUT_SPECS").c_str()))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_SPECS"));
        return false;
      }
      DEBUG_VAR(m_pairAlignments->size());
    }

    if (m_params.exists("OUTPUT_STARS"))
    {
      if (!m_starSpectraOnly->loadPklBin(m_params.getValue("OUTPUT_STARS").c_str()))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_STARS"));
        return false;
      }
      DEBUG_VAR(m_starSpectraOnly->size());
    }

    if (m_params.exists("OUTPUT_STARS_ALL"))
    {
      if (!m_starSpectra->loadPklBin(m_params.getValue("OUTPUT_STARS_ALL").c_str()))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_STARS_ALL"));
        return false;
      }
      DEBUG_VAR(m_starSpectra->size());
    }

    if (m_params.exists("OUTPUT_STARS_INDEX"))
    {
      if (Load_binArray(m_params.getValue("OUTPUT_STARS_INDEX").c_str(),
                        *m_alignedSpectra) <= 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_STARS_INDEX"));
        return false;
      }
      DEBUG_VAR(m_alignedSpectra->size());
    }

    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase *> const & ExecAlignment::split(int numSplit)
  {
	  m_subModules.resize(0);
	  return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecAlignment::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecAlignment::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("PENALTY_SAME_VERTEX");
    //VALIDATE_PARAM_EXIST("PENALTY_PTM");
    //VALIDATE_PARAM_EXIST("PENALTY_PTM_PEAKS");
    //VALIDATE_PARAM_EXIST("IDX_START");
    //VALIDATE_PARAM_EXIST("IDX_END");
    //VALIDATE_PARAM_EXIST("STARS_START");
    //VALIDATE_PARAM_EXIST("STARS_END");
    //VALIDATE_PARAM_EXIST("MAX_AA_JUMP");
    //VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    //VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    //VALIDATE_PARAM_EXIST("RESOLUTION");
    //VALIDATE_PARAM_EXIST("SPEC_TYPE_MSMS");
    //VALIDATE_PARAM_EXIST("PARTIAL_OVERLAPS");

    m_isValid = true;
    return true;
  }

} // namespace specnets
