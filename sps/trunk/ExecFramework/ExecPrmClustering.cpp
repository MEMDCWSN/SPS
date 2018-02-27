// Header Includes
#include "ExecPrmClustering.h"

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
  ExecPrmClustering::ExecPrmClustering(void) :
      ownInput(true), m_inputSpectra(0x0), m_inputMS2Spectra(0x0),
          m_inputClusters(0x0), ownOutput(true), m_outputSpectra(0x0),
          m_outputClusters(0x0), enforceDaTolerance(false), m_scanToIdx(0),
          m_scanToInputScan(0)
  {
    m_name = "ExecPrmClustering";
    m_type = "ExecPrmClustering";
    m_inputSpectra = new SpecSet;
    m_inputClusters = new ClusterSet;
    m_outputSpectra = new SpecSet;
    m_outputClusters = new ClusterSet;
    m_inputMS2Spectra = new SpecSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmClustering::ExecPrmClustering(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
          m_inputMS2Spectra(0x0), m_inputClusters(0x0), ownOutput(true),
          m_outputSpectra(0x0), m_outputClusters(0x0),
          enforceDaTolerance(false), m_scanToIdx(0), m_scanToInputScan(0)
  {
    m_name = "ExecPrmClustering";
    m_type = "ExecPrmClustering";
    m_inputSpectra = new SpecSet;
    m_inputClusters = new ClusterSet;
    m_outputSpectra = new SpecSet;
    m_outputClusters = new ClusterSet;
    m_inputMS2Spectra = new SpecSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmClustering::ExecPrmClustering(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * inputMS2Spectra,
                                       ClusterSet * inputClusters,
                                       SpecSet * outputSpectra,
                                       ClusterSet * outputClusters) :
      ExecBase(inputParams), ownInput(false), m_inputSpectra(inputSpectra),
          m_inputMS2Spectra(inputMS2Spectra), m_inputClusters(inputClusters),
          ownOutput(false), m_outputSpectra(outputSpectra),
          m_outputClusters(outputClusters), enforceDaTolerance(false),
          m_scanToIdx(0), m_scanToInputScan(0)
  {
    m_name = "ExecPrmClustering";
    m_type = "ExecPrmClustering";
  }

  // -------------------------------------------------------------------------
  ExecPrmClustering::ExecPrmClustering(const ParameterList & inputParams,
                                       SpecSet * outputSpectra,
                                       ClusterSet * outputClusters) :
      ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
          m_inputMS2Spectra(0x0), m_inputClusters(0x0), ownOutput(false),
          m_outputSpectra(outputSpectra), m_outputClusters(outputClusters),
          enforceDaTolerance(false), m_scanToIdx(0), m_scanToInputScan(0)
  {
    m_name = "ExecPrmClustering";
    m_type = "ExecPrmClustering";
    m_inputSpectra = new SpecSet;
    m_inputMS2Spectra = new SpecSet;
    m_inputClusters = new ClusterSet;
  }

  // -------------------------------------------------------------------------
  ExecPrmClustering::~ExecPrmClustering(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
      delete m_inputMS2Spectra;
      delete m_inputClusters;
    }

    if (ownOutput)
    {
      delete m_outputSpectra;
      delete m_outputClusters;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecPrmClustering::clone(const ParameterList & inputParams) const
  {
    return new ExecPrmClustering(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecPrmClustering::invoke(void)
  {
    if (m_inputSpectra == 0)
    {
      ERROR_MSG("No input spectra!");
      return false;
    }
    
    // initialize converter/loader
    vector<pair<int, int> > loadedIndices;
    ExecMergeConvert loader(m_params,
                            &loadedIndices,
                            m_inputSpectra,
                            m_outputSpectra);

    // Do any pre-processing here
    if (!loader.invoke())
    {
      return false;
    }

    bool noOp = m_params.getValueBool("PRM_CLUSTERING_NO_OP", false);
    if (noOp) {
      DEBUG_MSG("PRM_CLUSTERING_NO_OP");
      *m_outputClusters = *m_inputClusters;
      return true;
    }

    SpecSet headerInfo;
    headerInfo.resize(m_inputSpectra->size());
    for (unsigned int i = 0; i < m_inputSpectra->size(); i++)
    {
      headerInfo[i].copyNP((*m_inputSpectra)[i]);
      headerInfo[i].resize(0);
    }

    m_scanToIdx.resize(m_outputSpectra->size());
    m_scanToInputScan.resize(m_outputSpectra->size());
    for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
    {
      m_scanToInputScan[i] = (*m_outputSpectra)[i].scan;
      (*m_outputSpectra)[i].scan = i;
      m_scanToIdx[i].first = 0;
      m_scanToIdx[i].second = i;
    }

    // if true, just enforce Da tolerance rather than capturing PPM tolerance from input MS/MS spectra
    enforceDaTolerance = m_params.getValueInt("ENFORCE_DA_TOL", 1) > 0;

    // if true, do clustering of PRM spectra
    bool mergeSamePrec = m_params.getValueInt("MERGE_SAME_PREC", 0) > 0;

    int clusterMinSz = m_params.getValueInt("CLUSTER_MIN_SIZE", 0);

    int numConsec = (mergeSamePrec) ? -1 : 1;
    DEBUG_VAR(numConsec);

    vector<vector<bool> >* reversedPeaks = (vector<vector<bool> >*)0;

    vector<vector<unsigned int> > pairedScans;
    m_outputSpectra->extractPairedScans(pairedScans, numConsec);

    DEBUG_VAR(pairedScans.size());

    // Keep spectra in the same order as input MS/MS
    SpecSet orderedSpecs(m_outputSpectra->size());
    SpecSet orderedMs2Specs(m_inputMS2Spectra->size());
    DEBUG_VAR(m_outputSpectra->size());
    DEBUG_VAR(m_inputMS2Spectra->size())

    if (m_outputSpectra->size() != m_inputMS2Spectra->size())
    {
      ERROR_MSG("Number of PRM spectra does not equal # of MS/MS spectra!!");
      abort();
    }

    DEBUG_TRACE;

    if (m_params.getValueInt("BOOST_SILAC_PRMS", 0) > 0)
    {
      vector<vector<unsigned int> >* newClusters = new vector<
          vector<unsigned int> >(pairedScans.size());
      boostSilacPairs(m_outputSpectra,
                      m_inputMS2Spectra,
                      pairedScans,
                      *newClusters);
      pairedScans = *newClusters;
      delete newClusters;
    }

    pair<unsigned, unsigned> refIdx;

    string outFilename("");
    if (m_params.exists("OUTPUT_SPECTRA"))
    {
      vector<string> files;
      splitText(m_params.getValue("OUTPUT_SPECTRA").c_str(), files, ";");
      FilenameManager mngr(files[0]);
      outFilename = mngr.getFilenameWithExtension();
    }

    if (mergeSamePrec || clusterMinSz > 0)
    {
      // cluster spectra from the same precursor
      DEBUG_TRACE;
      vector<vector<unsigned int> > newClusters(pairedScans.size());

      DEBUG_VAR(m_outputSpectra->size());

      mergeSamePrecursorStaged(m_outputSpectra,
          pairedScans,
          newClusters,
          mergeSamePrec,
          clusterMinSz,
          reversedPeaks);

      DEBUG_VAR(m_outputSpectra->size());
      m_outputClusters->resize(m_outputSpectra->size());

      for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
      {
        (*m_outputClusters)[i].initialize(i,
            i + 1,
            0,
            outFilename,
            0);

        //cout << "Cluster " << i << ": ";
        for (unsigned int j = 0; j < newClusters[i].size(); j++)
        {
          unsigned int specIdx = newClusters[i][j];

          if (m_inputClusters->size() > 0)
          {
            for (int clustChild = 0; clustChild < (*m_inputClusters)[specIdx].size(); clustChild++)
            {
              (*m_outputClusters)[i].push_back((*m_inputClusters)[specIdx][clustChild]);
            }
          }
          else
          {
            (*m_outputClusters)[i].add(specIdx,
                headerInfo[specIdx].scan,
                headerInfo[specIdx].fileIndex,
                headerInfo[specIdx].fileName);
          }

          //cout << "(" << refIdx.first << "," << refIdx.second << "); ";
        }
        //cout << endl;
      }
      DEBUG_TRACE;

    }
    else
    {
      DEBUG_TRACE;
      if (m_inputClusters->size() > 0)
      {
        m_outputClusters->operator =(*m_inputClusters);
      }
      else
      {
        // otherwise, just generate single-ton clusters
        m_outputClusters->resize(m_outputSpectra->size());
        for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
        {
          (*m_outputClusters)[i].initialize(i,
              headerInfo[i].scan,
              headerInfo[i].fileIndex,
              outFilename,
              1);
          (*m_outputClusters)[i][0].initialize(i,
              headerInfo[i].scan,
              headerInfo[i].fileIndex,
              headerInfo[i].fileName);
        }
      }
    }

    if (reversedPeaks)
    {
      delete reversedPeaks;
    }

    DEBUG_TRACE;

    int rankFilt = m_params.getValueInt("PRM_RANK_FILTER", -1);

    if (rankFilt > 0)
    {
      DEBUG_MSG("Rank-filtering PRM spectra with K = " << rankFilt);
    }

    for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
    {

      const int curCharge = headerInfo[i].parentCharge;
      const float curMass = headerInfo[i].parentMass;

      if (rankFilt > 0)
      {
        (*m_outputSpectra)[i].rankFilterPeaks(rankFilt);
      }

      if (mergeSamePrec || clusterMinSz > 0)
      {
        (*m_outputSpectra)[i].scan = (*m_outputClusters)[i].m_scan;
        (*m_outputSpectra)[i].fileName = (*m_outputClusters)[i].m_filename;
        (*m_outputSpectra)[i].fileIndex = (*m_outputClusters)[i].m_fileIndex;
      }
      else
      {
        (*m_outputSpectra)[i].copyNP(headerInfo[i]);
      }

      (*m_outputSpectra)[i].msFragType = Spectrum::FragType_PRM;
    }

    DEBUG_TRACE;

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmClustering::loadInputData(void)
  {

    vector<pair<int, int> > loadedIndices;

    ExecMergeConvert loader(m_params,
                            &loadedIndices,
                            m_inputSpectra,
                            m_outputSpectra);

    if (!loader.loadInputData())
    {
      return false;
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_SPECTRA?");
      return false;
    }

    if (m_params.exists("INPUT_SPECTRA_MS"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA_MS from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\' ...");
      if (!m_inputMS2Spectra->Load(m_params.getValue("INPUT_SPECTRA_MS").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA_MS") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputMS2Spectra->size());

    }

    if (m_params.exists("INPUT_CLUSTERS"))
    {
      if (!m_inputClusters->loadBinaryFile(m_params.getValue("INPUT_CLUSTERS")))
      {
        ERROR_MSG("Failed to load clusters from \'" << m_params.getValue("INPUT_CLUSTERS") << "\'");
        return false;
      }
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
  bool ExecPrmClustering::saveInputData(std::vector<std::string> & filenames)
  {
    //SpecSet m_inputSpectra; // the input spectra
    /*
     if (m_params.exists("INPUT_SPECTRA")) {
     if (!ExecMergeConvert::saveSpecset(m_params.getValue("INPUT_SPECTRA"),
     m_inputSpectra)) {
     return false;
     }
     }
     */
    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmClustering::saveOutputData(void)
  {
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    DEBUG_MSG("Output spectra path: [" << dataDir << "]");

    if (m_params.exists("OUTPUT_SPECTRA")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_SPECTRA");
      DEBUG_MSG("Saving spectra to: [" << fileName << "]");
      if (!ExecMergeConvert::saveSpecsetMultiple(fileName,
                                                 m_outputSpectra))
      {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_CLUSTERS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_CLUSTERS");
      DEBUG_MSG("Saving clusters to: [" << fileName << "]");
      if (!m_outputClusters->saveBinaryFile(fileName)) {
        ERROR_MSG("Failed to save clusters to \'" << m_params.getValue("OUTPUT_CLUSTERS") << "\'");
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmClustering::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPrmClustering::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmClustering::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmClustering::validateParams(std::string & error)
  {
    m_isValid = false;

    m_isValid = true;
    return true;
  }

  void ExecPrmClustering::boostSilacPairs(SpecSet* inOutSpecs,
                                          SpecSet* ms2Spectra,
                                          vector<vector<unsigned int> >& inClusteredScans,
                                          vector<vector<unsigned int> >& outClusteredScans)
  {
    float modMassR = 10.008269;
    float modMassK = 8.014199;
    unsigned int minNumMP = 6;

    unsigned int numSilacSpecs = 0;

    float minRatioClust = m_params.getValueFloat("PRM_CLUSTER_RATIO", 0.72);

    int scanRange = m_params.getValueInt("SILAC_SCAN_RANGE", 100);

    if (scanRange <= 0 || scanRange > inOutSpecs->size())
    {
      scanRange = inOutSpecs->size();
    }

    bool filterNonBoost = m_params.getValueInt("FILTER_NONSILAC_PRMS", 0) > 0;

    DEBUG_VAR(filterNonBoost);

    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

    DEBUG_VAR(inOutSpecs->size());

    vector<bool> usedSpecs(inOutSpecs->size(), false);
    vector<bool> rootSpecs(inOutSpecs->size(), false);

    PairedSpecSet cluster;

    SpecSet outSpecs(inOutSpecs->size());
    unsigned int idxUse = 0;
    unsigned int numPairs = 0;

    outClusteredScans = inClusteredScans;

    unsigned int specIdx = 0;
    unsigned int pairedIdx;

    vector<list<unsigned int> > lightIdxs;
    vector<list<unsigned int> > heavyIdxs;
    vector<list<unsigned int> > startIdxs;
    list<unsigned int> emptyList;

    SpecSet locMS2Specs(scanRange);
    vector<unsigned int> locIdxMap(scanRange);

    SpecSet emptySpecs(scanRange);
    vector<TwoValues<float> > ratios;
    vector<TwoValues<float> > means;
    vector<float> varTerms;
    list<vector<float> > alignStats;
    vector<vector<float> > specStats;
    std::vector<unsigned int> idxKept;
    std::vector<TwoValues<float> > pvalues;

    SpectrumPairSet allUsedPairs;

    ParameterList filterParams;
    filterParams.setValue("MAX_SHIFT", "12.0");
    filterParams.setValue("TOLERANCE_PEAK",
                          m_params.getValue("TOLERANCE_PEAK"));
    filterParams.setValue("TOLERANCE_PM", m_params.getValue("TOLERANCE_PM"));
    filterParams.setValue("PAIRS_MIN_COSINE",
                          m_params.getValue("MIN_SILAC_COSINE", "0.4"));
    filterParams.setValue("MIN_NUM_MATCHED_PEAKS", "4");
    filterParams.setValue("SPEC_TYPE_MSMS", "1");
    filterParams.setValue("PAIRS_MATCH_MODE", "cosine");

    float curPMass, nextPMass, curPMtol, nextPMtol;

    while (specIdx < inOutSpecs->size())
    {
      if (usedSpecs[specIdx] || (*inOutSpecs)[specIdx].size() == 0)
      {
        specIdx++;
        continue;
      }

      unsigned int specIdxLoc = specIdx;
      for (unsigned int i = 0; i < scanRange && specIdxLoc < inOutSpecs->size();
          i++)
      {
        locIdxMap[i] = specIdxLoc;
        if (usedSpecs[specIdxLoc] || (*inOutSpecs)[specIdxLoc].size() == 0)
        {
          locMS2Specs[i].resize(0);
          specIdxLoc++;
          continue;
        }
        //DEBUG_VAR((*ms2Spectra)[specIdxLoc].parentMass);
        //DEBUG_VAR((*inOutSpecs)[specIdxLoc].parentMass);

        locMS2Specs[i] = (*ms2Spectra)[specIdxLoc];

        specIdxLoc++;
      }
      if (specIdxLoc == inOutSpecs->size())
      {
        locMS2Specs.resize(specIdxLoc - specIdx + 1);
        locIdxMap.resize(specIdxLoc - specIdx + 1);
      }

      DEBUG_VAR(specIdx);

      // set of matched heavy/light pairs
      SpectrumPairSet filteredPairs;

      // Let ExecFilterPairs do cosine matching to find heavy/light pairs
      ExecFilterPairs moduleAlign(filterParams,
                                  &emptySpecs,
                                  &locMS2Specs,
                                  &filteredPairs,
                                  &ratios,
                                  &means,
                                  &varTerms,
                                  &alignStats,
                                  &specStats,
                                  &idxKept,
                                  &pvalues);

      if (!moduleAlign.invoke())
      {
        ERROR_MSG("Failed to invoke ExecFilterPairs!!!");
        abort();
      }

      DEBUG_VAR(filteredPairs.size());

      curPMass = (*inOutSpecs)[specIdx].parentMass;
      curPMtol = (*inOutSpecs)[specIdx].parentMassTol;
      vector<float> srmOffset;
      vector<float> lightPMass;
      vector<float> heavyPMass;

      if (scanRange == inOutSpecs->size())
      {
        startIdxs.assign(inOutSpecs->size(), emptyList);
        for (unsigned int i = 0; i < inOutSpecs->size(); i++)
        {
          startIdxs[i].push_back(i);
        }
        lightIdxs.assign(inOutSpecs->size(), emptyList);
        heavyIdxs.assign(inOutSpecs->size(), emptyList);
        srmOffset.assign(inOutSpecs->size(), -1.0);
        lightPMass.assign(inOutSpecs->size(), -1.0);
        heavyPMass.assign(inOutSpecs->size(), -1.0);
      }
      else
      {
        startIdxs.assign(1, emptyList);
        startIdxs[0].push_back(specIdx);
        lightIdxs.assign(1, emptyList);
        heavyIdxs.assign(1, emptyList);
        srmOffset.assign(1, -1.0);
        lightPMass.assign(1, -1.0);
        heavyPMass.assign(1, -1.0);
        rootSpecs[specIdx] = true;
        usedSpecs[specIdx] = true;
      }

      DEBUG_VAR(specIdx);
      DEBUG_VAR(curPMass);
      DEBUG_VAR(curPMtol);
      filteredPairs.sort_pairs();

      DEBUG_TRACE;
      for (unsigned int i = 0; i < filteredPairs.size(); i++)
      {
        unsigned int specIdx1 = locIdxMap[filteredPairs[i].spec1];
        unsigned int specIdx2 = locIdxMap[filteredPairs[i].spec2];
        unsigned int refIdx = 0;

        if (scanRange == inOutSpecs->size())
        {
          if (specIdx1 < specIdx2)
          {
            specIdx = specIdx1;
            pairedIdx = specIdx2;
          }
          else
          {
            specIdx = specIdx2;
            pairedIdx = specIdx1;
          }
          refIdx = specIdx;
          curPMass = (*inOutSpecs)[specIdx].parentMass;
          curPMtol = (*inOutSpecs)[specIdx].parentMassTol;
          nextPMass = (*inOutSpecs)[pairedIdx].parentMass;
          nextPMtol = (*inOutSpecs)[pairedIdx].parentMassTol;
        }
        else
        {
          if (specIdx1 == specIdx)
          {
            nextPMass = (*inOutSpecs)[specIdx2].parentMass;
            nextPMtol = (*inOutSpecs)[specIdx2].parentMassTol;
            pairedIdx = specIdx2;
          }
          else if (specIdx2 == specIdx)
          {
            nextPMass = (*inOutSpecs)[specIdx1].parentMass;
            nextPMtol = (*inOutSpecs)[specIdx1].parentMassTol;
            pairedIdx = specIdx1;
          }
          else
          {
            continue;
          }
        }

        if (usedSpecs[pairedIdx] || (usedSpecs[specIdx] && !rootSpecs[specIdx]))
        {
          continue;
        }

        bool samePM = MZRange::EqualWithinRange(nextPMass,
                                                curPMass,
                                                curPMtol + nextPMtol);
        bool nextHR = MZRange::EqualWithinRange(nextPMass,
                                                curPMass + modMassR,
                                                curPMtol + nextPMtol);
        bool nextHK = MZRange::EqualWithinRange(nextPMass,
                                                curPMass + modMassK,
                                                curPMtol + nextPMtol);
        bool nextLR = MZRange::EqualWithinRange(nextPMass + modMassR,
                                                curPMass,
                                                curPMtol + nextPMtol);
        bool nextLK = MZRange::EqualWithinRange(nextPMass + modMassK,
                                                curPMass,
                                                curPMtol + nextPMtol);

        bool foundSim = samePM || nextHR || nextHK || nextLR || nextLK;

        if (!foundSim)
        {
          continue;
        }

        /*
         DEBUG_VAR(specIdx);
         DEBUG_VAR(pairedIdx);
         DEBUG_VAR(filteredPairs[i].score1);
         DEBUG_VAR((*inOutSpecs)[specIdx].parentMass);
         DEBUG_VAR((*inOutSpecs)[pairedIdx].parentMass);
         */

        //DEBUG_VAR(filteredPairs[i].score1);
        if (srmOffset[refIdx] < 0)
        {
          if (samePM)
          {
            startIdxs[refIdx].push_back(pairedIdx);
          }
          else if (nextHR)
          {
            srmOffset[refIdx] = modMassR;
            lightIdxs[refIdx] = startIdxs[refIdx];
            heavyIdxs[refIdx].push_back(pairedIdx);
            lightPMass[refIdx] = curPMass;
            heavyPMass[refIdx] = curPMass + modMassR;
          }
          else if (nextHK)
          {
            srmOffset[refIdx] = modMassK;
            lightIdxs[refIdx] = startIdxs[refIdx];
            heavyIdxs[refIdx].push_back(pairedIdx);
            lightPMass[refIdx] = curPMass;
            heavyPMass[refIdx] = curPMass + modMassK;
          }
          else if (nextLR)
          {
            srmOffset[refIdx] = modMassR;
            heavyIdxs[refIdx] = startIdxs[refIdx];
            lightIdxs[refIdx].push_back(pairedIdx);
            heavyPMass[refIdx] = curPMass;
            lightPMass[refIdx] = curPMass - modMassR;
          }
          else
          {
            srmOffset[refIdx] = modMassK;
            heavyIdxs[refIdx] = startIdxs[refIdx];
            lightIdxs[refIdx].push_back(pairedIdx);
            heavyPMass[refIdx] = curPMass;
            lightPMass[refIdx] = curPMass - modMassK;
          }
          if (srmOffset[refIdx] >= 0)
          {
            rootSpecs[specIdx] = true;
            for (list<unsigned int>::const_iterator idxIt =
                lightIdxs[refIdx].begin(); idxIt != lightIdxs[refIdx].end();
                idxIt++)
            {
              usedSpecs[*idxIt] = true;
            }
            for (list<unsigned int>::const_iterator idxIt =
                heavyIdxs[refIdx].begin(); idxIt != heavyIdxs[refIdx].end();
                idxIt++)
            {
              usedSpecs[*idxIt] = true;
            }
          }
        }
        else if (MZRange::EqualWithinRange(nextPMass,
                                           heavyPMass[refIdx],
                                           curPMtol + nextPMtol))
        {
          heavyIdxs[refIdx].push_back(pairedIdx);
          usedSpecs[pairedIdx] = true;
        }
        else if (MZRange::EqualWithinRange(nextPMass,
                                           lightPMass[refIdx],
                                           curPMtol + nextPMtol))
        {
          lightIdxs[refIdx].push_back(pairedIdx);
          usedSpecs[pairedIdx] = true;
        }
        else
        {
          continue;
        }

        SpectrumPair newPair;
        newPair = filteredPairs[i];
        newPair.spec1 = specIdx1;
        newPair.spec2 = specIdx2;
        allUsedPairs.push_back(newPair);

      }

      if (scanRange != inOutSpecs->size()
          && (lightIdxs[0].size() == 0 || heavyIdxs[0].size() == 0))
      {
        if (!filterNonBoost)
        {
          outSpecs[specIdx] = (*inOutSpecs)[specIdx];
          outSpecs[specIdx].msFragType = Spectrum::FragType_CID;
        }

        specIdx++;
        continue;
      }

      for (unsigned int i = 0; i < srmOffset.size(); i++)
      {
        if (scanRange == inOutSpecs->size())
        {
          specIdx = i;
        }

        if (lightIdxs[i].size() == 0 || heavyIdxs[i].size() == 0)
        {
          if (!filterNonBoost)
          {
            outSpecs[specIdx] = (*inOutSpecs)[specIdx];
            outSpecs[specIdx].msFragType = Spectrum::FragType_CID;
          }

          continue;
        }

        numPairs++;

        SpecSet pairedSpecs(lightIdxs[i].size() + heavyIdxs[i].size());
        unsigned int locIdx = 0;
        vector<unsigned int> boostedIdxs(lightIdxs[i].size()
            + heavyIdxs[i].size());
        for (list<unsigned int>::const_iterator idxIt = lightIdxs[i].begin();
            idxIt != lightIdxs[i].end(); idxIt++)
        {
          pairedSpecs[locIdx] = (*inOutSpecs)[*idxIt];
          pairedSpecs[locIdx].msFragType = Spectrum::FragType_CID;
          boostedIdxs[locIdx] = *idxIt;
          locIdx++;
        }
        for (list<unsigned int>::const_iterator idxIt = heavyIdxs[i].begin();
            idxIt != heavyIdxs[i].end(); idxIt++)
        {
          pairedSpecs[locIdx] = (*inOutSpecs)[*idxIt];
          pairedSpecs[locIdx].msFragType = Spectrum::FragType_ETD;
          boostedIdxs[locIdx] = *idxIt;
          locIdx++;
        }

        numSilacSpecs += lightIdxs[i].size() + heavyIdxs[i].size();

        /*
         DEBUG_VAR(specIdx);
         DEBUG_VAR(lightIdxs.size());
         DEBUG_VAR(heavyIdxs.size());
         */

        cluster.initialize(&pairedSpecs, true, 0, peakTol, srmOffset[i]);

        SpecSet boostedSpecs;

        cluster.boostPRMs(boostedSpecs);

        //DEBUG_VAR(boostedSpecs.size());

        for (locIdx = 0; locIdx < boostedSpecs.size(); locIdx++)
        {
          if (filterNonBoost)
          {
            outSpecs[idxUse] = boostedSpecs[locIdx];
            outSpecs[idxUse].msFragType = Spectrum::FragType_CID;
            outClusteredScans[idxUse] = inClusteredScans[boostedIdxs[locIdx]];
            idxUse++;
          }
          else
          {
            outSpecs[boostedIdxs[locIdx]] = boostedSpecs[locIdx];
            outSpecs[boostedIdxs[locIdx]].msFragType = Spectrum::FragType_CID;
          }
        }
      }

      if (scanRange == inOutSpecs->size())
      {
        break;
      }

      specIdx++;
    }

    if (filterNonBoost)
    {
      DEBUG_VAR(idxUse);
      outSpecs.resize(idxUse);
      outClusteredScans.resize(idxUse);
    }
    inOutSpecs->operator =(outSpecs);

    DEBUG_VAR(inOutSpecs->size());

    DEBUG_MSG("Found " << numSilacSpecs << " silac spectra and " << numPairs << " components");

    if (m_params.exists("OUTPUT_SILAC_PAIRS"))
    {
      if (!allUsedPairs.saveToBinaryFile(m_params.getValue("OUTPUT_SILAC_PAIRS")))
      {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_SILAC_PAIRS"));
        return;
      }
    }
  }

  struct SortPairs : public std::binary_function<
      pair<float, TwoValues<unsigned int> >,
      pair<float, TwoValues<unsigned int> >, bool>
  {
    bool operator()(pair<float, TwoValues<unsigned int> > left,
                    pair<float, TwoValues<unsigned int> > right) const
    {
      return left.first > right.first;
    }
    ;
  };

  void ExecPrmClustering::mergeSamePrecursorStaged(SpecSet* inOutSpecs,
                                                   vector<vector<unsigned int> >& inPairedScans,
                                                   vector<vector<unsigned int> >& outClusteredScans,
                                                   bool mergeSamePrec,
                                                   int clusterMinSz,
                                                   vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Merging spectra from the same precursor with hierarchical method");

    // Minimum allowable overlapping score between clustered PRM spectra
    float minRatioClust = m_params.getValueFloat("PRM_CLUSTER_RATIO", 0.72);

    // How many peaks to keep in a +/- 56 Da radius in clustered PRM spectra
    unsigned int rankFiltK = m_params.getValueInt("PRM_RANK_FILTER", 3);

    // Minimum allowable number of PRM spectra within a cluster
    int minClustSize = clusterMinSz;

    DEBUG_VAR(rankFiltK);

    DEBUG_VAR(minRatioClust);

    // Keep a reference of scan # to lookup index
    map<unsigned int, unsigned int> scanToIdx;
    for (unsigned int i = 0; i < inOutSpecs->size(); i++)
    {
      scanToIdx[(*inOutSpecs)[i].scan] = i;
    }

    // Compute initial set of clusters (pairs or triples)

    // initial clustered spectra
    SpecSet* initialClusters = new SpecSet(inPairedScans.size());

    // keep track of clustered MS/MS indices
    vector<list<unsigned int> > clusterIdxs(initialClusters->size());

    // used to merge PRM spectra
    PairedSpecSet nextCluster;

    for (unsigned int i = 0; i < inPairedScans.size(); i++)
    {
      SpecSet pairedSpecs;
      vector<vector<bool> > pairedRev(0);
      //DEBUG_VAR(i);

      // Put triplet or paired spectra in a container
      for (unsigned int j = 0; j < inPairedScans[i].size(); j++)
      {
        //DEBUG_VAR(j);
        unsigned int specIdx = scanToIdx[inPairedScans[i][j]];
        pairedSpecs.push_back((*inOutSpecs)[specIdx]);
        //DEBUG_VAR((*inOutSpecs)[specIdx].parentMass);
        //DEBUG_VAR((*inOutSpecs)[specIdx].msFragType);
        //DEBUG_VAR((*inOutSpecs)[specIdx].size());

        if (reversedPeaks)
        {
          pairedRev.push_back((*reversedPeaks)[specIdx]);
        }
      }

      // also get the indices of MS/MS spectra
      unsigned int firstScan = inPairedScans[i][0];
      clusterIdxs[i].clear();
      clusterIdxs[i].push_back(scanToIdx[firstScan]);
      for (unsigned int j = 1; j < inPairedScans[i].size(); j++)
      {
        clusterIdxs[i].push_back(scanToIdx[inPairedScans[i][j]]);
      }

      if (reversedPeaks)
      {
        nextCluster.initialize(&pairedSpecs,
                               enforceDaTolerance,
                               &pairedRev,
                               m_params.getValueFloat("TOLERANCE_PEAK"));
      }
      else
      {
        nextCluster.initialize(&pairedSpecs,
                               enforceDaTolerance,
                               0,
                               m_params.getValueFloat("TOLERANCE_PEAK"));
      }

      // merge PRM spectra
      nextCluster.mergePRMs();
      nextCluster.getMergedSpectrum(&(*initialClusters)[i]);

      // rank filter peaks
      (*initialClusters)[i].rankFilterPeaks(rankFiltK);
    }

    DEBUG_VAR(initialClusters->size());
    DEBUG_VAR(inOutSpecs->size());

    // The final set of clusters after hierarchical clustering
    SpecSet finalClusters(initialClusters->size());

    finalClusters = *initialClusters;

    bool foundPair = true;
    PRMAlignment nextPair;
    vector<bool> mergedSpecs(finalClusters.size(), false);
    list<pair<float, TwoValues<unsigned int> > > candidatePairs;
    pair<float, TwoValues<unsigned int> > tempPair;
    unsigned int stage = 1;
    unsigned int numClusters = finalClusters.size();
    Spectrum sortedPrecursors;
    list<int> idxCheck;
    //map<MZRange, list<unsigned int> > precursorMap;

    if (minClustSize > 0)
    {

      while (foundPair)
      {
        DEBUG_MSG("Clustering stage " << stage << " ...");
        stage++;

        foundPair = false;
        candidatePairs.clear();

        DEBUG_MSG("Indexing precursors ...");

        // Use a sorted list (ie a spectrum) to index precursor masses
        sortedPrecursors.resize(finalClusters.size());
        unsigned int tempIdxUse = 0;

        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          Spectrum* spec1 = &finalClusters[i];

          if (spec1->size() == 0)
          {
            continue;
          }

          sortedPrecursors[tempIdxUse][0] = spec1->parentMass;
          sortedPrecursors[tempIdxUse][1] = (float)i;
          sortedPrecursors.setTolerance(tempIdxUse, spec1->parentMassTol);
          tempIdxUse++;
        }
        sortedPrecursors.resize(tempIdxUse);

        DEBUG_MSG("Sorting ...");
        sortedPrecursors.sortPeaks();

        DEBUG_MSG("Finding pairs ...");

        // find all pairs of spectra with the same precursors
        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          Spectrum* spec1 = &finalClusters[i];
          MZRange nextRange(spec1->parentMass, 0, spec1->parentMassTol);
          if (spec1->size() == 0)
          {
            continue;
          }
          nextPair.setSpec1(spec1);

          // lookup matching precursor masses
          sortedPrecursors.findPeaks(nextRange, &idxCheck);

          for (list<int>::const_iterator idxIt = idxCheck.begin();
              idxIt != idxCheck.end(); idxIt++)
          {
            unsigned int j = floatToInt(sortedPrecursors[*idxIt][1]);

            if (j <= i)
            {
              continue;
            }

            Spectrum* spec2 = &finalClusters[j];

            if (spec2->size() == 0)
            {
              continue;
            }

            if (!MZRange::EqualWithinRange(spec1->parentMass,
                                           spec2->parentMass,
                                           spec1->parentMassTol
                                               + spec2->parentMassTol))
            {
              continue;
            }

            nextPair.setSpec2(spec2);

            //DEBUG_MSG("Considering " << spec1->parentMass << " and " << spec2->parentMass << " w/ tolerance " << spec1->parentMassTol + spec2->parentMassTol);

            pair<int, pair<float, float> > alignScore =
                nextPair.getShiftScore(0, 0, 1);

            float minR = min(alignScore.second.first, alignScore.second.second);

            //DEBUG_VAR(minR);

            if (minR >= minRatioClust)
            {
              tempPair.first = minR;
              tempPair.second[0] = i;
              tempPair.second[1] = j;
              candidatePairs.push_back(tempPair);
            }
          }
        }

        DEBUG_VAR(candidatePairs.size());

        if (candidatePairs.size() == 0)
        {
          break;
        }

        // order pairs in by decreasing matched score ratio
        candidatePairs.sort(SortPairs());

        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          mergedSpecs[i] = false;
        }

        // merge pairs one-by-one, until all are merged
        for (list<pair<float, TwoValues<unsigned int> > >::const_iterator pIt =
            candidatePairs.begin(); pIt != candidatePairs.end(); pIt++)
        {
          unsigned int i = pIt->second[0];
          unsigned int j = pIt->second[1];

          if (mergedSpecs[i] || mergedSpecs[j])
          {
            // this pair has already been merged by transitive pairs
            continue;
          }

          foundPair = true;

          // keep track of clustered indices
          clusterIdxs[i].insert(clusterIdxs[i].end(),
                                clusterIdxs[j].begin(),
                                clusterIdxs[j].end());
          clusterIdxs[j].clear();
          SpecSet pairedSpecs;

          for (list<unsigned int>::iterator cIt = clusterIdxs[i].begin();
              cIt != clusterIdxs[i].end(); cIt++)
          {
            pairedSpecs.push_back((*inOutSpecs)[*cIt]);
          }

          nextCluster.initialize(&pairedSpecs,
                                 enforceDaTolerance,
                                 0,
                                 m_params.getValueFloat("TOLERANCE_PEAK"));

          nextCluster.mergePRMs();
          nextCluster.getMergedSpectrum(&finalClusters[i]);
          finalClusters[i].rankFilterPeaks(rankFiltK);

          finalClusters[j].resize(0);
          numClusters--;
          mergedSpecs[i] = true;
          mergedSpecs[j] = true;
        }
        DEBUG_VAR(numClusters);
      }
    }

    outClusteredScans.resize(clusterIdxs.size());

    unsigned int idxUse = 0;
    DEBUG_VAR(minClustSize);

    SpecSet outSpecs(inOutSpecs->size());
    DEBUG_VAR(inOutSpecs->size());

    float numDiffCharges = 0;
    float totalClusters = 0;

    // Compute final set of condensed clusters
    for (unsigned int i = 0; i < clusterIdxs.size(); i++)
    {
      if (clusterIdxs[i].size() < minClustSize || finalClusters[i].size() == 0)
      {
        continue;
      }

      //DEBUG_VAR(idxUse);
      outSpecs[idxUse] = finalClusters[i];

      outClusteredScans[idxUse].resize(0);
      set<short> parentCharges;
      for (list<unsigned int>::iterator cIt = clusterIdxs[i].begin();
          cIt != clusterIdxs[i].end(); cIt++)
      {
        unsigned int specIdx = *cIt;

        outClusteredScans[idxUse].push_back((*inOutSpecs)[specIdx].scan);
        parentCharges.insert((*inOutSpecs)[specIdx].parentCharge);
      }

      if (clusterIdxs[i].size() > 1)
      {
        totalClusters += 1.0;
        numDiffCharges += (parentCharges.size() > 1) ? 1.0 : 0;
      }

      idxUse++;
    }
    DEBUG_MSG(parseFloat(numDiffCharges * 100.0 / totalClusters, 1) << "\% of clusters contain at least 2 spectra with different precursor charge states");

    outSpecs.resize(idxUse);
    inOutSpecs->operator =(outSpecs);
    outClusteredScans.resize(idxUse);

    DEBUG_VAR(inOutSpecs->size());

    delete initialClusters;
  }

  void ExecPrmClustering::checkSpecset(SpecSet &inOutSpecs)
  {
    DEBUG_MSG("Checking sepecset. Size: " << inOutSpecs.size());

    for (int i = 0; i < inOutSpecs.size(); i++)
    {
      //DEBUG_VAR(inOutSpecs[i].size());
      if (inOutSpecs[i].size() <= 8)
      {
        DEBUG_MSG("Spectrum #" << i << " has " << inOutSpecs[i].size() << " peaks. Resizing to 0.");
        inOutSpecs[i].resize(0);
      }
    }

    DEBUG_MSG("Checking sepecset done.");
  }

} // namespace specnets
