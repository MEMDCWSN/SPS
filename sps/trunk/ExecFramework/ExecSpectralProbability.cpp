/*
 * ExecSpectralProbability.cpp
 *
 *  Created on: Nov 5, 2013
 *      Author: aguthals
 */

#include "ExecSpectralProbability.h"
#include "GFTable.h"
#include "alignment_scoring.h"
#include "ExecMergeConvert.h"
#include <time.h>

using namespace std;

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/memory>
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#  include <unordered_map>
#endif

namespace specnets
{

  ExecSpectralProbability::ExecSpectralProbability(void) :
      ExecBase(), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputPSMs(0x0), m_inputErrorModelTarget(0x0), m_inputErrorModelDecoy(0x0), m_inputMsSpectra(0x0), m_outputSpectra(0x0), m_outputPSMs(0x0), m_outputClusters(0x0), decoysRemoved(false), m_jumps(0x0)
  {
    m_name = "ExecSpectralProbability";
    m_type = "ExecSpectralProbability";
    m_inputSpectra = new vector<SpecSet*>;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputErrorModelTarget = new MassErrorModel;
    m_inputErrorModelDecoy = new MassErrorModel;
    m_inputMsSpectra = new SpecSet;
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new PeptideSpectrumMatchSet;
    m_outputClusters = new ClusterSet;
    m_jumps = new AAJumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);
    m_jumps->multiplyMasses(AA_ROUNDING);
  }

  ExecSpectralProbability::ExecSpectralProbability(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), ownOutput(true), m_inputSpectra(0x0), m_inputPSMs(0x0), m_inputErrorModelTarget(0x0), m_inputErrorModelDecoy(0x0), m_inputMsSpectra(0x0), m_outputSpectra(0x0), m_outputPSMs(0x0), m_outputClusters(0x0), decoysRemoved(false), m_jumps(0x0)
  {
    m_name = "ExecSpectralProbability";
    m_type = "ExecSpectralProbability";
    m_inputSpectra = new vector<SpecSet*>;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputMsSpectra = new SpecSet;
    m_inputErrorModelTarget = new MassErrorModel;
    m_inputErrorModelDecoy = new MassErrorModel;
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new PeptideSpectrumMatchSet;
    m_outputClusters = new ClusterSet;
    m_jumps = new AAJumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);
    m_jumps->multiplyMasses(AA_ROUNDING);
  }

  ExecSpectralProbability::ExecSpectralProbability(const ParameterList &inputParams,
                                                   vector<SpecSet*> *inputSpectra,
                                                   PeptideSpectrumMatchSet *inputPSMs,
                                                   MassErrorModel *inputErrorModelTarget,
                                                   MassErrorModel *inputErrorModelDecoy,
                                                   SpecSet *inputMsSpectra,
                                                   AAJumps *jumps) :
      ExecBase(inputParams), ownInput(false), ownOutput(true), m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs), m_inputErrorModelTarget(inputErrorModelTarget), m_inputErrorModelDecoy(inputErrorModelDecoy), m_inputMsSpectra(inputMsSpectra), m_outputSpectra(0x0), m_outputPSMs(0x0), m_outputClusters(0x0), decoysRemoved(false), m_jumps(jumps)
  {
    m_name = "ExecSpectralProbability";
    m_type = "ExecSpectralProbability";
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new PeptideSpectrumMatchSet;
    m_outputClusters = new ClusterSet;
  }

  ExecSpectralProbability::ExecSpectralProbability(const ParameterList &inputParams,
                                                   vector<SpecSet*> *inputSpectra,
                                                   PeptideSpectrumMatchSet *inputPSMs,
                                                   MassErrorModel *inputErrorModelTarget,
                                                   MassErrorModel *inputErrorModelDecoy,
                                                   SpecSet *inputMsSpectra,
                                                   AAJumps *jumps,
                                                   SpecSet *outputSpectra,
                                                   ClusterSet *outputClusters,
                                                   PeptideSpectrumMatchSet *outputPSMs) :
      ExecBase(inputParams), ownInput(false), ownOutput(false), m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs), m_inputErrorModelTarget(inputErrorModelTarget), m_inputErrorModelDecoy(inputErrorModelDecoy), m_inputMsSpectra(inputMsSpectra), m_outputSpectra(outputSpectra), m_outputPSMs(outputPSMs), m_outputClusters(outputClusters), decoysRemoved(false), m_jumps(jumps)
  {
    m_name = "ExecSpectralProbability";
    m_type = "ExecSpectralProbability";
  }

  ExecSpectralProbability::~ExecSpectralProbability(void)
  {
    if (ownInput)
    {
      for (int i = 0; i < m_inputSpectra->size(); i++)
      {
        delete (*m_inputSpectra)[i];
      }
      delete m_inputSpectra;
      delete m_inputPSMs;
      delete m_inputErrorModelTarget;
      delete m_inputErrorModelDecoy;
      delete m_inputMsSpectra;
      delete m_jumps;
    }
    if (ownOutput)
    {
      delete m_outputSpectra;
      delete m_outputPSMs;
      delete m_outputClusters;
    }
  }

  ExecBase * ExecSpectralProbability::clone(const ParameterList & inputParams) const
  {
    return new ExecSpectralProbability(inputParams);
  }

  bool ExecSpectralProbability::invoke(void)
  {
    if ((*m_inputSpectra)[0]->size() == 0)
    {
      WARN_MSG("Input spectra size is 0, exiting invoke()");
      return true;
    }

    AAJumps jumpsNoRound(1, 0.01, -1, AAJumps::NO_MODS, false, true);

    const float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.05);
    const float pmTol = m_params.getValueFloat("TOLERANCE_PM");

    const int minMatchedPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS", 0);

    const bool computeOptimalProb = m_params.getValueInt("COMPUTE_OPTIMAL", 0)
        > 0;

    const int maxPsmRank = m_params.getValueInt("MAX_RANK", 1);

    const float normalizedSpecScore =
        m_params.getValueFloat("NORMALIZED_SPECTRUM_SCORE", -1.0);

    const float ppmEdgeError = m_params.getValueFloat("TOLERANCE_PEAK_PPM")
        * 2.0;

    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)(*m_inputSpectra)[0]->size())
                                                - 1);

    const int debugScan = m_params.getValueInt("DEBUG_SCAN", -1);

    const bool parseSilac = m_params.getValueBool("PARSE_SILAC_PSMS", false);

    const bool ignoreMS2ErrorModel = (m_inputErrorModelTarget->size() == 0
        || m_inputErrorModelDecoy->size() == 0 || m_inputMsSpectra->size() == 0);

    removeDecoyPSMsMatchingTarget();

    if (m_inputSpectra->size() == 0 || (*m_inputSpectra)[0]->size() == 0)
    {
      WARN_MSG("Found no input spectra, exiting invoke()");
      m_outputSpectra->resize(0);
      m_outputClusters->resize(0);
      return true;
    }

    int maxPSMSz = 0;
    for (int i = idxStart; i <= idxEnd; i++)
    {
      maxPSMSz = max(maxPSMSz, (int)(*(*m_inputSpectra)[0])[i].psmList.size());
    }

#pragma omp critical
    {
      DEBUG_VAR(maxPSMSz);
      DEBUG_VAR(idxStart);
      DEBUG_VAR(idxEnd);
      DEBUG_VAR(m_inputSpectra->size());
      DEBUG_VAR((*m_inputSpectra)[0]->size());
      DEBUG_VAR(normalizedSpecScore);
      DEBUG_VAR(maxPsmRank);
      DEBUG_VAR(parseSilac);
    }

    m_outputSpectra->resize(0);
    m_outputClusters->resize((*m_inputSpectra)[0]->size());

    clock_t init, final;
    int numComputed = 0;
    init = clock();
    GFTable gfTable;

    Spectrum emptySpectrum;

    int lastProg = -1;
    double secondsPerIntersection = 0;

    for (int i = idxStart; i <= idxEnd; i++)
    {
      (*m_outputClusters)[i].resize(0);

      if ((*(*m_inputSpectra)[0])[i].psmList.size() == 0)
      {
        continue;
      }

      if (debugScan >= 0 && (*(*m_inputSpectra)[0])[i].scan != debugScan)
      {
        continue;
      }

      bool debug = (debugScan >= 0);

      int curPsmIdx = 0;

      vector<PeptideSpectrumMatch> bestPsms(maxPSMSz);
      SpecSet psmSpecs(maxPSMSz);
      vector<TwoValues<int> > rankedPSMs(maxPSMSz);

      for (list<psmPtr>::reverse_iterator pIt =
          (*(*m_inputSpectra)[0])[i].psmList.rbegin();
          pIt != (*(*m_inputSpectra)[0])[i].psmList.rend(); pIt++)
      {
        string befPeptide = (*pIt)->m_annotation;

        const float silacMod = (parseSilac) ? (*pIt)->removeSILACMod(false) : 0;

        const string &peptide = (*pIt)->m_annotation;

        PeptideSpectrumMatch bestPsm;
        int maxMatchScore = -1;
        Spectrum bestSpec;

        Spectrum prmMasses0;
        m_jumps->getPRMMasses(peptide, prmMasses0);

        Spectrum prmMassesNoRound;
        jumpsNoRound.getPRMMasses(peptide, prmMassesNoRound);

        float observedPM = (*(*m_inputSpectra)[0])[i].parentMass;

        for (int c13Idx = 0; c13Idx < m_inputSpectra->size(); c13Idx++)
        {
          Spectrum spec = (*(*m_inputSpectra)[c13Idx])[i];
          if (debug)
          {
            DEBUG_VAR(spec.parentMass);
            DEBUG_VAR(prmMassesNoRound.parentMass);
          }

          if (parseSilac)
          {
            spec.setParentMass(spec.parentMass - silacMod);
          }
          else if (m_inputSpectra->size() > 1
              && !isEqual(spec.parentMass, prmMassesNoRound.parentMass, pmTol))
          {
            continue;
          }

          spec.removeZPMpeaks();

          if (spec.size() == 0)
          {
            continue;
          }

          Spectrum specY = spec;
          specY.reverse(0);

          spec.roundPeaks(AA_ROUNDING, true, false);
          specY.roundPeaks(AA_ROUNDING, true, false);

          if (normalizedSpecScore > 0)
          {
            spec.normalize(normalizedSpecScore, 0.0);
            specY.normalize(normalizedSpecScore, 0.0);
          }

          float score0B = 0, score0Y = 0;

          /*
           if (debug)
           {
           DEBUG_TRACE;
           spec.output(cerr);
           DEBUG_TRACE;
           prmMasses0.output(cerr);
           }
           */

          vector<int> idxMatched1;
          vector<int> idxMatched2;
          FindMatchPeaksAll2(spec,
                             prmMasses0,
                             0,
                             peakTol,
                             idxMatched1,
                             idxMatched2);

          for (unsigned int p = 0; p < idxMatched1.size(); p++)
          {
            if (debug)
              DEBUG_MSG("Matched mass " << spec[idxMatched1[p]][0] << " (prefix idx = " << idxMatched2[p] << ")");
            score0B += spec[idxMatched1[p]][1];
          }

          if (score0B > maxMatchScore && idxMatched1.size() >= minMatchedPeaks)
          {
            maxMatchScore = score0B;
            bestPsm = **pIt;
            bestSpec = spec;
            bestPsm.m_useYendPts = 0;
          }

          FindMatchPeaksAll2(specY,
                             prmMasses0,
                             0,
                             peakTol,
                             idxMatched1,
                             idxMatched2);

          for (unsigned int p = 0; p < idxMatched1.size(); p++)
          {
            score0Y += specY[idxMatched1[p]][1];
          }

          if (score0Y > maxMatchScore && idxMatched1.size() >= minMatchedPeaks)
          {
            maxMatchScore = score0Y;
            bestPsm = **pIt;
            bestSpec = specY;
            bestPsm.m_useYendPts = 1;
          }
        }

        if (maxMatchScore >= 0)
        {
          bestPsms[curPsmIdx] = bestPsm;
          psmSpecs[curPsmIdx] = bestSpec;
          rankedPSMs[curPsmIdx][0] = maxMatchScore;
          rankedPSMs[curPsmIdx][1] = curPsmIdx;
          curPsmIdx++;
        }
      }

      for (int psmIdx = curPsmIdx; psmIdx < rankedPSMs.size(); psmIdx++)
      {
        rankedPSMs[psmIdx][0] = -1;
        rankedPSMs[psmIdx][1] = -1;
      }

      sort(rankedPSMs.begin(), rankedPSMs.end());

      bool addedCluster = false;

      int specCount = 0;
      for (int rankIdx = rankedPSMs.size() - 1;
          rankIdx >= 0 && (rankedPSMs.size() - rankIdx) <= maxPsmRank
              && rankedPSMs[rankIdx][1] >= 0; rankIdx--)
      {
        if (!addedCluster)
        {
          (*m_outputClusters)[i].m_index = i;
          (*m_outputClusters)[i].m_scan = (*(*m_inputSpectra)[0])[i].scan;
          (*m_outputClusters)[i].m_fileIndex =
              (*(*m_inputSpectra)[0])[i].fileIndex;
          (*m_outputClusters)[i].m_filename =
              (*(*m_inputSpectra)[0])[i].fileName;

          addedCluster = true;
        }

        int psmIdx = rankedPSMs[rankIdx][1];
        psmSpecs[psmIdx].scan = ((i + 1) * maxPSMSz) + specCount;
        specCount++;

        vector<float> prmMassesRounded;

        if (!m_jumps->getRoundedPRMMasses(bestPsms[psmIdx].m_annotation,
                                          prmMassesRounded))
        {
          ERROR_MSG("Cannot parse peptide sequence \'" << bestPsms[psmIdx].m_annotation << "\'");
          abort();
        }

        float newPm = prmMassesRounded.back() + round(AAJumps::massMH);
        psmSpecs[psmIdx].setParentMass(newPm);

        psmSpecs[psmIdx].setPeakTolerance(0);
        psmSpecs[psmIdx].setParentMassTol(0);
        psmSpecs[psmIdx].roundPeaks(1, true, false);

        Spectrum &msSpectrum =
            (ignoreMS2ErrorModel) ? emptySpectrum : (*m_inputMsSpectra)[i];

        gfTable.initialize(psmSpecs[psmIdx],
                           *m_jumps,
                           msSpectrum,
                           *m_inputErrorModelTarget,
                           *m_inputErrorModelDecoy,
                           ppmEdgeError); //&bestPsm.m_annotation);

        numComputed++;

        int score = gfTable.getMatchScore(bestPsms[psmIdx].m_annotation);

        gfTable.computeGeneratingFunction(score);

        double pValue = gfTable.getPValue(score);

        if (debug)
        {
          DEBUG_VAR(psmIdx);
          DEBUG_VAR(bestPsms[psmIdx].m_annotation);
          DEBUG_VAR(score);
          DEBUG_VAR(pValue);
        }

        bestPsms[psmIdx].m_strict_envelope_score = bestPsms[psmIdx].m_score;
        bestPsms[psmIdx].m_pValue = pValue;
        bestPsms[psmIdx].m_score = score;
        bestPsms[psmIdx].m_charge = (*(*m_inputSpectra)[0])[i].parentCharge;

        psmSpecs[psmIdx].psmList.clear();

        int specIdx = m_outputSpectra->size();
        bestPsms[psmIdx].m_scanNum = psmSpecs[psmIdx].scan;
        bestPsms[psmIdx].m_spectrumFile = psmSpecs[psmIdx].fileName;
        m_outputSpectra->push_back(psmSpecs[psmIdx]);

        (*m_outputClusters)[i].add(specIdx,
                                   psmSpecs[psmIdx].scan,
                                   psmSpecs[psmIdx].fileIndex,
                                   psmSpecs[psmIdx].fileName);

        psmPtr newPsm(new PeptideSpectrumMatch);
        *newPsm = bestPsms[psmIdx];

        (*m_outputSpectra)[specIdx].psmList.push_back(newPsm);
        m_outputPSMs->push_back(newPsm);

        //DEBUG_MSG("Spectrum " << i << " has ID " << psm.m_annotation << " with probability " << parseDoubleSci(psm.m_pValue, 2) << "; sec/spec = " << secondsPerIntersection);

        double curProg = ((((double)i) - ((double)idxStart)) * 100.0)
            / (((double)idxEnd) - ((double)idxStart));
        if ((int)curProg > lastProg)
        {
          lastProg = (int)curProg;

          final = clock();
          double seconds = (((double)final) - ((double)init))
              / ((double)CLOCKS_PER_SEC);
          secondsPerIntersection = seconds / ((double)numComputed);

#pragma omp critical
          {
            DEBUG_MSG(m_name << " ... " << lastProg << "%, sec/spec = " << secondsPerIntersection);
          }
        }
      }
    }

    DEBUG_TRACE;

    return true;
  }

  bool ExecSpectralProbability::loadInputData(void)
  {
    int debugScan = m_params.getValueInt("DEBUG_SCAN", -1);

    m_inputSpectra->resize(0);
    for (int c13 = 0; c13 <= 5; c13++)
    {
      string param = ParameterList::getInputSpectraC13Param(c13);
      if (m_params.exists(param))
      {
        DEBUG_MSG("Loading INPUT_SPECTRA from \'" << m_params.getValue(param) << "\' ...");
        SpecSet* tempSpecs = new SpecSet;
        m_inputSpectra->push_back(tempSpecs);
        if (!m_inputSpectra->back()->Load(m_params.getValue(param).c_str()))
        {
          ERROR_MSG("Failed to load spectra from \'" << m_params.getValue(param) << "\'");
          return false;
        }
        DEBUG_VAR(m_inputSpectra->back()->size());
      }
    }

    for (int c13 = -1; c13 >= -5; c13--)
    {
      string param = ParameterList::getInputSpectraC13Param(c13);
      if (m_params.exists(param))
      {
        DEBUG_MSG("Loading INPUT_SPECTRA from \'" << m_params.getValue(param) << "\' ...");
        SpecSet* tempSpecs = new SpecSet;
        m_inputSpectra->push_back(tempSpecs);
        if (!m_inputSpectra->back()->Load(m_params.getValue(param).c_str()))
        {
          ERROR_MSG("Failed to load spectra from \'" << m_params.getValue(param) << "\'");
          return false;
        }
        DEBUG_VAR(m_inputSpectra->back()->size());
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("No input spectra!!!");
      return false;
    }

    if (m_params.exists("INPUT_SPECTRA_MS"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA_MS from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\' ...");
      if (!m_inputMsSpectra->Load(m_params.getValue("INPUT_SPECTRA_MS").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputMsSpectra->size());

      DEBUG_MSG("Pre-processing MS/MS spectra ...");
      for (int i = 0; i < m_inputMsSpectra->size(); i++)
      {
        if ((*m_inputMsSpectra)[i].scan != (*(*m_inputSpectra)[0])[i].scan)
        {
          ERROR_MSG("Mismatched scans for index " << i << ": " << (*m_inputMsSpectra)[i].scan << " != " << (*(*m_inputSpectra)[0])[i].scan);
          abort();
        }
        if (debugScan >= 0 && debugScan != (*m_inputMsSpectra)[i].scan)
        {
          continue;
        }

        (*m_inputMsSpectra)[i].rankFilterPeaks(7);
        (*m_inputMsSpectra)[i].setPeakTolerance(0);
        (*m_inputMsSpectra)[i].roundPeaks(AA_ROUNDING, true, false);

      }
      DEBUG_MSG("finished");
    }

    if (m_params.exists("INPUT_ERROR_MODEL_TARGET"))
    {
      string param = m_params.getValue("INPUT_ERROR_MODEL_TARGET");
      DEBUG_MSG("Loading INPUT_ERROR_MODEL_TARGET from \'" << param << "\' ...");
      if (!m_inputErrorModelTarget->loadBinaryFile(param))
      {
        ERROR_MSG("Failed to load target error model from \'" << param << "\'");
        return false;
      }
      DEBUG_VAR(m_inputErrorModelTarget->size());
    }

    if (m_params.exists("INPUT_ERROR_MODEL_DECOY"))
    {
      string param = m_params.getValue("INPUT_ERROR_MODEL_DECOY");
      DEBUG_MSG("Loading INPUT_ERROR_MODEL_DECOY from \'" << param << "\' ...");
      if (!m_inputErrorModelDecoy->loadBinaryFile(param))
      {
        ERROR_MSG("Failed to load target error model from \'" << param << "\'");
        return false;
      }
      DEBUG_VAR(m_inputErrorModelDecoy->size());
    }

    if (m_params.exists("INPUT_PSMS"))
    {
      const string psmType = m_params.getValue("INPUT_PSMS_TYPE", "");
      const string psmFiles = m_params.getValue("INPUT_PSMS");

      vector<string> fileList;
      splitText(psmFiles.c_str(), fileList, ";");

      m_inputPSMs->resize(0);

      for (int fIdx = 0; fIdx < fileList.size(); fIdx++)
      {
        PeptideSpectrumMatchSet tempPsms;
        DEBUG_MSG("Loading INPUT_PSMS from \'" << fileList[fIdx] << "\' ...");
        if (!tempPsms.Load(fileList[fIdx], psmType))
        {
          ERROR_MSG("Failed to load PSMs from \'" << fileList[fIdx] << "\'");
          return false;
        }
        DEBUG_VAR(tempPsms.size());

        int curIdx = m_inputPSMs->size();
        m_inputPSMs->resize(curIdx + tempPsms.size());
        for (int i = 0; i < tempPsms.size(); i++)
        {
          (*m_inputPSMs)[curIdx++] = tempPsms[i];
        }
      }

      DEBUG_VAR(m_inputPSMs->size());
      m_inputPSMs->addMostSpectra((*m_inputSpectra)[0]);
    }
    DEBUG_TRACE;

    return true;
  }

  bool ExecSpectralProbability::saveOutputData(void)
  {
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::saveSpecsetMultiple(m_params.getValue("OUTPUT_SPECTRA"),
              m_outputSpectra))
      {
        ERROR_MSG("Failed to save to \'" << m_params.getValue("OUTPUT_SPECTRA")
            << "\'");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_PSMS"))
    {
      const string psmFile = m_params.getValue("OUTPUT_PSMS");

      if (!m_outputPSMs->saveToFile(psmFile.c_str(), true))
      {
        ERROR_MSG("Failed to save to \'" << m_params.getValue("OUTPUT_PSMS")
            << "\'");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_CLUSTERS"))
    {
      if (!m_outputClusters->saveBinaryFile(m_params.getValue("OUTPUT_CLUSTERS")))
      {
        ERROR_MSG("Failed to save clusters to \'"
            << m_params.getValue("OUTPUT_CLUSTERS") << "\'");
        return false;
      }
    }
    return true;
  }

  bool ExecSpectralProbability::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
    /*
     string dataDir = m_params.getValue("GRID_DATA_DIR", ".");
     string baseDirectory = dataDir + "/";
     string baseFilename = baseDirectory + getName();
     string paramFilename = baseFilename + ".params";

     DEBUG_VAR(dataDir);
     DEBUG_VAR(paramFilename);
     filenames.push_back(paramFilename); // Parameter file MUST be first in vector

     string spectraFilename = baseDirectory + "input_specs_scored_";
     spectraFilename += parseInt(m_childID);
     spectraFilename += ".pklbin";

     if (!ExecMergeConvert::saveSpecsetMultiple(spectraFilename, m_inputSpectra))
     {
     ERROR_MSG("Failed to save to \'" << spectraFilename << "\'");
     return false;
     }

     m_params.setValue("INPUT_SPECTRA", spectraFilename);

     string psmsFilename = baseDirectory + "input_psms_";
     psmsFilename += parseInt(m_childID);
     psmsFilename += ".txt";

     if (!m_inputPSMs->saveToFile(psmsFilename.c_str(), true))
     {
     ERROR_MSG("Failed to save to \'" << psmsFilename << "\'");
     return false;
     }

     m_params.setValue("INPUT_PSMS", psmsFilename);

     m_params.writeToFile(paramFilename);

     return true;
     */
  }

  bool ExecSpectralProbability::loadOutputData(void)
  {
    return false;
    /*
     DEBUG_TRACE;
     if (m_params.exists("OUTPUT_SPECTRA"))
     {
     if (!m_outputSpectra->Load(m_params.getValue("OUTPUT_SPECTRA").c_str()))
     {
     ERROR_MSG("Failed to load from \'" << m_params.getValue("OUTPUT_SPECTRA") << "\'");
     return false;
     }
     }

     if (m_params.exists("OUTPUT_PSMS"))
     {
     const string psmFile = m_params.getValue("OUTPUT_PSMS");

     if (!m_outputPSMs->loadFromFile(psmFile.c_str()))
     {
     ERROR_MSG("Failed to load from \'" << m_params.getValue("OUTPUT_PSMS") << "\'");
     return false;
     }
     }
     return true;
     */
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*>
  const & ExecSpectralProbability::split(int numSplit)
  {
    m_subModules.resize(0);
    if (m_inputSpectra->size() == 0 || m_inputPSMs->size() == 0)
    {
      if (!loadInputData())
      {
        return m_subModules;
      }
    }

    removeDecoyPSMsMatchingTarget();

    const int idxStart = m_params.getValueInt("IDX_START", 0);
    const int idxEnd = m_params.getValueInt("IDX_END",
                                            ((int)(*m_inputSpectra)[0]->size())
                                                - 1);
    DEBUG_VAR(idxStart);
    DEBUG_VAR(idxEnd);

    DEBUG_VAR(numSplit);

    if (numSplit == 0)
    {
      return m_subModules;
    }

    double totalComputations = 0;
    vector<double> computationsPerSpec((*m_inputSpectra)[0]->size());
    for (int i = idxStart; i <= idxEnd; i++)
    {
      computationsPerSpec[i] = (*(*m_inputSpectra)[0])[i].parentMass
          * (*(*m_inputSpectra)[0])[i].getTotalIonCurrent() * 20.0;
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

    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }

    for (int i = idxStart; i <= idxEnd; i++)
    {
      runningComputations += computationsPerSpec[i];
      runningComputationsPerModule += computationsPerSpec[i];

      if ((runningComputationsPerModule >= computationsPerModule
          && curSplit < numSplit + 1) || i == idxEnd)
      {
        ParameterList childParams(m_params);
        childParams.setValue("IDX_END", parseInt(i));
        childParams.setValue("IDX_START", parseInt(curIdxStart));
        childParams.setValue("DECOY_PROTEIN_ID", "");

        DEBUG_MSG("Job " << curSplit << " has range [" << curIdxStart << "," << i << "]");

        ExecBase * theClone =
            new ExecSpectralProbability(childParams,
                                        m_inputSpectra,
                                        m_inputPSMs,
                                        m_inputErrorModelTarget,
                                        m_inputErrorModelDecoy,
                                        m_inputMsSpectra,
                                        m_jumps);

        theClone->setName(makeName(m_name, curSplit));

        ExecSpectralProbability* castedClone =
            dynamic_cast<ExecSpectralProbability*>(theClone);
        castedClone->decoysRemoved = true;

        string baseName = dataDir + "/" + theClone->getName();
        theClone->m_params.setValue("OUTPUT_SPECTRA",
                                    baseName + "_output_spectra.pklbin");
        theClone->m_params.setValue("OUTPUT_PSMS",
                                    baseName + "_output_psms.txt");
        theClone->m_params.setValue("OUTPUT_CLUSTERS",
                                    baseName + "_output_spectra.clust");

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
  bool ExecSpectralProbability::merge(void)
  {
    m_outputSpectra->resize(0);
    m_outputPSMs->resize(0);
    m_outputClusters->resize(0);

    int numSpecsLoaded = 0;

    for (int nSplit = 0; nSplit < m_subModules.size(); nSplit++)
    {
      if (m_subModules[nSplit] == 0)
      {
        ERROR_MSG("ExecBase module at index " << nSplit << " is zero");
        return false;
      }
      ExecSpectralProbability* module =
          dynamic_cast<ExecSpectralProbability*>(m_subModules[nSplit]);

      if (module == 0)
      {
        ERROR_MSG("Failed to cast module \'" << m_subModules[nSplit]->getName() << "\' to ExecSpectralProbability");
        return false;
      }

      for (int i = 0; i < module->m_outputClusters->size(); i++)
      {
        int clustIdx = m_outputClusters->size();
        m_outputClusters->push_back((*(module->m_outputClusters))[i]);
        (*m_outputClusters)[clustIdx].m_index = clustIdx;

        for (int j = 0; j < (*m_outputClusters)[clustIdx].size(); j++)
        {
          (*m_outputClusters)[clustIdx][j].m_index += numSpecsLoaded;
        }
      }

      for (int i = 0; i < module->m_outputSpectra->size(); i++)
      {
        m_outputSpectra->push_back((*(module->m_outputSpectra))[i]);
      }

      for (int i = 0; i < module->m_outputPSMs->size(); i++)
      {
        m_outputPSMs->push_back((*(module->m_outputPSMs))[i]);
      }

      numSpecsLoaded += module->m_outputSpectra->size();
    }
    return true;
  }

  bool ExecSpectralProbability::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK_PPM");

    m_isValid = true;
    return true;
  }

  void ExecSpectralProbability::removeDecoyPSMsMatchingTarget()
  {
    const string decoyID = m_params.getValue("DECOY_PROTEIN_ID", "");

    if (decoyID.length() > 0 && !decoysRemoved)
    {
      tr1::unordered_set<string> targetPeptides((*m_inputSpectra)[0]->size());

      int numTargets = 0;
      for (int i = 0; i < (*m_inputSpectra)[0]->size(); i++)
      {
        for (list<psmPtr>::const_iterator pIt =
            (*(*m_inputSpectra)[0])[i].psmList.begin();
            pIt != (*(*m_inputSpectra)[0])[i].psmList.end(); pIt++)
        {
          if ((*pIt)->m_protein.find(decoyID) == string::npos)
          {
            numTargets++;
            const string &peptide = (*pIt)->m_annotation;

            string peptideUse = AAJumps::stripMods(peptide);

            replaceAll(peptideUse, "I", "L");
            replaceAll(peptideUse, "K", "Q");

            targetPeptides.insert(peptideUse);
          }
        }
      }

      int numDecoyRemoved = 0;
      int numDecoys = 0;

      for (int i = 0; i < (*m_inputSpectra)[0]->size(); i++)
      {
        Spectrum tempSpec;
        for (list<psmPtr>::const_iterator pIt =
            (*(*m_inputSpectra)[0])[i].psmList.begin();
            pIt != (*(*m_inputSpectra)[0])[i].psmList.end(); pIt++)
        {
          const string &peptide = (*pIt)->m_annotation;
          if (decoyID.length() > 0
              && (*pIt)->m_protein.find(decoyID) != string::npos)
          {
            numDecoys++;
            string peptideUse = AAJumps::stripMods(peptide);
            replaceAll(peptideUse, "I", "L");
            replaceAll(peptideUse, "K", "Q");

            if (targetPeptides.count(peptideUse) == 0)
            {
              tempSpec.psmList.push_back(*pIt);
            }
            else
            {
              numDecoyRemoved++;
              continue;
            }
          }
          else
          {
            tempSpec.psmList.push_back(*pIt);
          }
        }
        (*(*m_inputSpectra)[0])[i].psmList = tempSpec.psmList;
      }
      decoysRemoved = true;
      DEBUG_MSG("Removed " << numDecoyRemoved << " decoy PSMs matching to target");
      DEBUG_VAR(numTargets);
      DEBUG_VAR(numDecoys);
    }
  }
}
