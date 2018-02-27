/*
 * ExecStarGF.cpp
 *
 *  Created on: Nov 15, 2013
 *      Author: aguthals
 */

#include "ExecStarGF.h"
#include "utils.h"
#include "ExecMergeConvert.h"
#include "OutputTable.h"
#include "FdrPeptide.h"
#include "alignment_scoring.h"

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

  const int ExecStarGF::DEBUG_SCAN1 = -1; //194509;
  const int ExecStarGF::DEBUG_SCAN2 = -1; //140656;

  bool ExecStarGF::SavePeptideStarMatches(const string &filename,
                                          const vector<PeptideStarMatch> &outputPSMs)
  {
    OutputTable table;
    table.values.resize(outputPSMs.size() + 1);

    unsigned int r = 0, c = 0;
    const unsigned int numCols = 27;

    /*
     int specIdx;
     int useYEndpts;
     int pepScore;
     int bestIntersectPVIdx;
     double specProb;
     double starProb;
     string peptide;
     string peptideOvlp;
     string protein;
     */

    pair<string, bool> defaultHeader("", true);
    table.values[r].assign(numCols, defaultHeader);
    table.values[r][c++].first = "SpecFile";
    table.values[r][c++].first = "SpecIdx";
    table.values[r][c++].first = "Scan";
    table.values[r][c++].first = "Parent Mass";
    table.values[r][c++].first = "Charge";
    table.values[r][c++].first = "Peptide";
    table.values[r][c++].first = "PeptideLen";
    table.values[r][c++].first = "Y-endpts";
    table.values[r][c++].first = "MSGFProb";
    table.values[r][c++].first = "PepScore";
    table.values[r][c++].first = "SpecProb";
    table.values[r][c++].first = "StarProb";
    table.values[r][c++].first = "AlignGFProb";
    table.values[r][c++].first = "FDR";
    table.values[r][c++].first = "PepFDR";
    table.values[r][c++].first = "Rank";
    table.values[r][c++].first = "PairRank";
    table.values[r][c++].first = "IsDecoy";
    table.values[r][c++].first = "PairIsDecoy";
    table.values[r][c++].first = "PairMP";
    table.values[r][c++].first = "PairSpecProb";
    table.values[r][c++].first = "PairStarProb";
    table.values[r][c++].first = "PairPep";
    table.values[r][c++].first = "PairIdx";
    table.values[r][c++].first = "PairScan";
    table.values[r][c++].first = "PairFile";
    table.values[r][c++].first = "Protein";

    c = 0;
    r++;
    pair<string, bool> defaultStat("", false);
    for (unsigned int i = 0; i < outputPSMs.size(); i++)
    {
      if (outputPSMs[i].specIdx < 0)
      {
        continue;
      }

      table.values[r].assign(numCols, defaultStat);
      table.values[r][c++].first = outputPSMs[i].filename;
      table.values[r][c++].first = parseInt(outputPSMs[i].specIdx);
      table.values[r][c++].first = parseInt(outputPSMs[i].scan);
      table.values[r][c++].first = parseFloat(outputPSMs[i].precursor, 3);
      table.values[r][c++].first = parseInt(outputPSMs[i].charge);
      table.values[r][c++].first = outputPSMs[i].peptide;
      table.values[r][c++].first =
          parseInt(AAJumps::getNumJumps(outputPSMs[i].peptide));
      table.values[r][c++].first = parseInt((int)outputPSMs[i].useYEndpts);
      table.values[r][c++].first = parseDoubleSci(outputPSMs[i].msgfSpecProb,
                                                  5);
      table.values[r][c++].first = parseInt(outputPSMs[i].pepScore);
      table.values[r][c++].first = parseDoubleSci(outputPSMs[i].specProb, 5);
      table.values[r][c++].first = parseDoubleSci(outputPSMs[i].starProb, 5);
      table.values[r][c++].first =
          parseDoubleSci((double)outputPSMs[i].alignGFProb, 3);
      table.values[r][c++].first = parseDoubleSci((double)outputPSMs[i].fdr, 3);
      table.values[r][c++].first = parseDoubleSci((double)outputPSMs[i].pepFdr,
                                                  3);
      table.values[r][c++].first = parseInt(outputPSMs[i].rank);
      table.values[r][c++].first = parseInt(outputPSMs[i].pairRank);
      table.values[r][c++].first = parseInt(outputPSMs[i].isDecoy);
      table.values[r][c++].first = parseInt(outputPSMs[i].pairIsDecoy);
      table.values[r][c++].first = parseInt(outputPSMs[i].pairMP);
      table.values[r][c++].first = parseDoubleSci(outputPSMs[i].pairSpecProb,
                                                  5);
      table.values[r][c++].first = parseDoubleSci(outputPSMs[i].pairStarProb,
                                                  5);
      table.values[r][c++].first = outputPSMs[i].bestPairPeptide;
      table.values[r][c++].first = parseInt(outputPSMs[i].bestIntersectPVIdx);
      table.values[r][c++].first = parseInt(outputPSMs[i].bestPairScan);
      table.values[r][c++].first = outputPSMs[i].bestPairFile;
      table.values[r][c++].first = outputPSMs[i].protein;
      c = 0;
      r++;

    }
    table.values.resize(r);
    if (!table.printToCSV(filename.c_str(), "\t"))
    {
      ERROR_MSG("Failed to write to file \'" << filename << "\'!");
      return false;
    }
    return true;
  }

  ExecStarGF::ExecStarGF(void) :
      ExecBase(), ownInput(true), ownOutput(true), m_inputSpectra(0x0),
          m_inputPSMs(0x0), m_inputPairs(0x0), m_inputAlignedPairs(0x0),
          m_inputClusters(0x0), m_inputProbs(0x0), m_outputSpectra(0x0),
          m_outputPSMs(0x0), m_inputDB(0x0)
  {
    m_name = "ExecStarGF";
    m_type = "ExecStarGF";
    m_inputSpectra = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputPairs = new SpectrumPairSet;
    m_inputAlignedPairs = new SpectrumPairSet;
    m_inputClusters = new ClusterSet;
    m_inputProbs = new vector<vector<double> >;
    m_inputDB = new DB_fasta;
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new vector<PeptideStarMatch>;
  }

  ExecStarGF::ExecStarGF(const ParameterList & inputParams) :
      ExecBase(inputParams), ownInput(true), ownOutput(true),
          m_inputSpectra(0x0), m_inputPSMs(0x0), m_inputPairs(0x0),
          m_inputAlignedPairs(0x0), m_inputClusters(0x0), m_inputProbs(0x0),
          m_outputSpectra(0x0), m_outputPSMs(0x0), m_inputDB(0x0)
  {
    m_name = "ExecStarGF";
    m_type = "ExecStarGF";
    m_inputSpectra = new SpecSet;
    m_inputPSMs = new PeptideSpectrumMatchSet;
    m_inputPairs = new SpectrumPairSet;
    m_inputAlignedPairs = new SpectrumPairSet;
    m_inputClusters = new ClusterSet;
    m_inputProbs = new vector<vector<double> >;
    m_inputDB = new DB_fasta;
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new vector<PeptideStarMatch>;
  }

  ExecStarGF::ExecStarGF(const ParameterList & inputParams,
                         SpecSet * inputSpectra,
                         PeptideSpectrumMatchSet * inputPSMs,
                         ClusterSet *inputClusters,
                         SpectrumPairSet * inputPairs,
                         SpectrumPairSet * inputAlignedPairs,
                         vector<vector<double> > *inputProbs,
                         DB_fasta *inputDB) :
      ExecBase(inputParams), ownInput(false), ownOutput(true),
          m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs),
          m_inputPairs(inputPairs), m_inputAlignedPairs(inputAlignedPairs),
          m_inputClusters(inputClusters), m_inputProbs(inputProbs),
          m_inputDB(inputDB), m_outputSpectra(0x0), m_outputPSMs(0x0)
  {
    m_name = "ExecStarGF";
    m_type = "ExecStarGF";
    m_outputSpectra = new SpecSet;
    m_outputPSMs = new vector<PeptideStarMatch>;
  }

  ExecStarGF::ExecStarGF(const ParameterList & inputParams,
                         SpecSet * inputSpectra,
                         PeptideSpectrumMatchSet * inputPSMs,
                         ClusterSet *inputClusters,
                         SpectrumPairSet * inputPairs,
                         SpectrumPairSet * inputAlignedPairs,
                         vector<vector<double> > *inputProbs,
                         DB_fasta *inputDB,
                         SpecSet *outputSpectra,
                         vector<PeptideStarMatch> *outputPSMs) :
      ExecBase(inputParams), ownInput(false), ownOutput(false),
          m_inputSpectra(inputSpectra), m_inputPSMs(inputPSMs),
          m_inputPairs(inputPairs), m_inputAlignedPairs(inputAlignedPairs),
          m_inputClusters(inputClusters), m_inputProbs(inputProbs),
          m_inputDB(inputDB), m_outputSpectra(outputSpectra),
          m_outputPSMs(outputPSMs)
  {
    m_name = "ExecStarGF";
    m_type = "ExecStarGF";
  }

  ExecStarGF::~ExecStarGF(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
      delete m_inputPairs;
      delete m_inputAlignedPairs;
      delete m_inputPSMs;
      delete m_inputClusters;
      delete m_inputProbs;
      delete m_inputDB;
    }
    if (ownOutput)
    {
      delete m_outputSpectra;
      delete m_outputPSMs;
    }
  }

  ExecBase * ExecStarGF::clone(const ParameterList & inputParams) const
  {
    return new ExecStarGF(inputParams);
  }

  bool compareStarProbs(ExecStarGF::PeptideStarMatch i,
                        ExecStarGF::PeptideStarMatch j)
  {
    return i.starProb < j.starProb;
  }

  bool compareMSGFProbs(ExecStarGF::PeptideStarMatch i,
                        ExecStarGF::PeptideStarMatch j)
  {
    return i.msgfSpecProb < j.msgfSpecProb;
  }

  bool compareSpecProbs(ExecStarGF::PeptideStarMatch i,
                        ExecStarGF::PeptideStarMatch j)
  {
    return i.specProb < j.specProb;
  }

  bool ExecStarGF::invoke(void)
  {
    AAJumps jumps(1, 0.01, -1, AAJumps::NO_MODS, false, true);
    jumps.multiplyMasses(AA_ROUNDING);

    float peakTol = 0.5; //m_params.getValueFloat("TOLERANCE_PEAK");
    const float pmTol = m_params.getValueFloat("TOLERANCE_PM");

    string decoyID = m_params.getValue("DECOY_PROTEIN_ID", "");

    const bool networkOnly = m_params.getValueBool("NETWORK_ONLY", false);

    const int maxRank = m_params.getValueInt("MAX_RANK", 1);

    const int minPeptideLength = m_params.getValueInt("MIN_PEPTIDE_LEN", 7);

    const bool allowSameChargePairs = m_params.getValueBool("SAME_CHARGE_PAIRS",
                                                            true);
    const bool allowDiffChargePairs = m_params.getValueBool("DIFF_CHARGE_PAIRS",
                                                            true);
    const bool allowPartialPairs =
        m_params.getValueBool("PARTIAL_OVERLAP_PAIRS", true);

    const float minRatio = m_params.getValueFloat("MIN_RATIO", 0);

    const int numMatchedPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS", 3);

    const double maxAlignGFPval = m_params.getValueDouble("MAX_ALIGNGF_PVAL",
                                                          1);

    DEBUG_VAR(numMatchedPeaks);
    DEBUG_VAR(minRatio);
    DEBUG_VAR(maxAlignGFPval);

    string fdrMetric = m_params.getValue("FDR_METRIC", "StarProb");
    std::transform(fdrMetric.begin(),
                   fdrMetric.end(),
                   fdrMetric.begin(),
                   ::tolower);

    DEBUG_VAR(fdrMetric);
    DEBUG_VAR(maxRank);

    if (m_inputSpectra->size() == 0)
    {
      WARN_MSG("No input spectra!!");
      return true;
    }

    DEBUG_VAR(m_inputPairs->size());
    DEBUG_VAR(m_inputProbs->size());
    if (m_inputPairs->size() != m_inputProbs->size())
    {
      ERROR_MSG("Input pairs do not align with input intersecting probabilities");
      return false;
    }

    DEBUG_VAR(m_inputClusters->size());
    if (m_inputClusters->size() == 0)
    {
      WARN_MSG("No input clusters!!");
      return true;
    }

    vector<string> peptidesUse(m_inputSpectra->size(), "");
    m_outputPSMs->resize(m_inputSpectra->size());

    vector<int> specToCluster(m_inputSpectra->size(), -1);
    vector<bool> psmsAllowed(m_inputSpectra->size(), false);
    int numAdded = 0;
    set<int> uniqueSpecs;
    for (int i = 0; i < m_inputClusters->size(); i++)
    {
      PeptideSpectrumMatchSet clustPSMs;
      for (int j = 0; j < (*m_inputClusters)[i].size(); j++)
      {
        specToCluster[(*m_inputClusters)[i][j].m_index] =
            (*m_inputClusters)[i].m_index;

        /*
         DEBUG_VAR(i);
         DEBUG_VAR(j);
         DEBUG_VAR((*m_inputClusters)[i][j].m_index);
         DEBUG_VAR((*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].scan);
         DEBUG_VAR((*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].fileName);
         DEBUG_VAR((*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].psmList.size());
         (*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].output(cerr);
         */
        clustPSMs.push_back((*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].psmList.front());
        clustPSMs.m_psmSet.back()->m_dbIndex = (*m_inputClusters)[i][j].m_index;

        if (uniqueSpecs.count((*m_inputClusters)[i][j].m_index) > 0)
        {
          DEBUG_VAR((*m_inputClusters)[i][j].m_index);
        }

        uniqueSpecs.insert((*m_inputClusters)[i][j].m_index);

        if (fdrMetric == "msgf")
        {
          clustPSMs.m_psmSet.back()->m_pValue =
              (*m_inputSpectra)[(*m_inputClusters)[i][j].m_index].psmList.front()->m_strict_envelope_score;
        }
      }
      //DEBUG_VAR(clustPSMs.size());
      sort(clustPSMs.m_psmSet.begin(),
           clustPSMs.m_psmSet.end(),
           FdrPeptide::compareFDRPVal);
      for (int psmIdx = 0; psmIdx < clustPSMs.size() && psmIdx < maxRank;
          psmIdx++)
      {
        numAdded++;
        psmsAllowed[clustPSMs[psmIdx]->m_dbIndex] = true;
        (*m_inputSpectra)[clustPSMs[psmIdx]->m_dbIndex].psmList.front()->m_shared_peaks =
            psmIdx + 1;
      }
    }
    DEBUG_VAR(numAdded);
    DEBUG_VAR(uniqueSpecs.size());
    vector<pair<int, double> > bestSpecPerCluster(m_inputClusters->size(),
                                                  pair<int, double>(-1, 1.0));

    vector<vector<bool> > matchedPRMs(m_inputSpectra->size());

    DEBUG_VAR(m_inputSpectra->size());

    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      Spectrum &spec = (*m_inputSpectra)[i];
      if (spec.psmList.size() == 0)
      {
        continue;
      }
      int clust = specToCluster[i];
      if (clust < 0)
      {
        continue;
      }
      const string &peptide = spec.psmList.front()->m_annotation;

      if (peptide.length() < minPeptideLength)
      {
        continue;
      }

      bool debug = false;
      if ((*m_inputSpectra)[i].scan == DEBUG_SCAN1
          || (*m_inputSpectra)[i].scan == DEBUG_SCAN2)
      {
        DEBUG_VAR(i);
        DEBUG_VAR((*m_inputSpectra)[i].scan);
        DEBUG_VAR(peptide);
        debug = true;
      }

      peptidesUse[i] = AAJumps::stripMods(peptide);

      //replaceAll(peptidesUse[i], "I", "L");
      //replaceAll(peptidesUse[i], "K", "Q");

      (*m_outputPSMs)[i].msgfSpecProb =
          spec.psmList.front()->m_strict_envelope_score;
      (*m_outputPSMs)[i].specProb = spec.psmList.front()->m_pValue;
      (*m_outputPSMs)[i].starProb = 1.0;
      (*m_outputPSMs)[i].useYEndpts = spec.psmList.front()->m_useYendPts == 1;

      (*m_outputPSMs)[i].pepScore = (int)(round(spec.psmList.front()->m_score)
          + 0.1);
      (*m_outputPSMs)[i].protein = spec.psmList.front()->m_protein;
      (*m_outputPSMs)[i].rank = spec.psmList.front()->m_shared_peaks;
      (*m_outputPSMs)[i].peptide = peptide;
      (*m_outputPSMs)[i].peptideOvlp = peptidesUse[i];
      (*m_outputPSMs)[i].filename = (*m_inputClusters)[clust].m_filename;
      (*m_outputPSMs)[i].scan = (*m_inputClusters)[clust].m_scan;
      if (debug)
        DEBUG_VAR((*m_inputClusters)[clust].m_scan);

      (*m_outputPSMs)[i].specIdx = (*m_inputClusters)[clust].m_index;
      (*m_outputPSMs)[i].precursor = spec.psmList.front()->m_mz;
      (*m_outputPSMs)[i].charge = spec.psmList.front()->m_charge;

      Spectrum prmMasses;
      jumps.getRoundedPRMMasses(peptide, prmMasses);
      matchedPRMs[i].assign(prmMasses.size(), false);

      vector<int> idxMatched1;
      vector<int> idxMatched2;
      FindMatchPeaksAll2((*m_inputSpectra)[i],
                         prmMasses,
                         0,
                         0.5,
                         idxMatched1,
                         idxMatched2);

      for (int mIdx = 0; mIdx < idxMatched2.size(); mIdx++)
      {
        matchedPRMs[i][idxMatched2[mIdx]] = true;
      }

      if ((*m_inputClusters)[clust].m_scan == DEBUG_SCAN1
          || (*m_inputClusters)[clust].m_scan == DEBUG_SCAN2)
      {
        DEBUG_VAR(clust);
        DEBUG_VAR((*m_inputClusters)[clust].m_index);
        DEBUG_VAR((*m_inputClusters)[clust].m_scan);
        DEBUG_VAR(i);
        DEBUG_VAR((*m_outputPSMs)[i].peptide);
        DEBUG_VAR((*m_outputPSMs)[i].specProb);
        DEBUG_VAR((*m_outputPSMs)[i].starProb);
        DEBUG_VAR(idxMatched1.size());
      }

      if (psmsAllowed[i]
          && (*m_outputPSMs)[i].specProb < bestSpecPerCluster[clust].second)
      {
        bestSpecPerCluster[clust].first = i;
        bestSpecPerCluster[clust].second = (*m_outputPSMs)[i].specProb;
      }
    }

    DEBUG_VAR(m_outputPSMs->size());

    DEBUG_TRACE;
    m_outputSpectra->resize(m_inputClusters->size());
    vector<map<int, int> > inputAligns(m_inputSpectra->size());
    vector<int> matchingPeaks(m_inputPairs->size(), -1);

    bool filterAlignedPairs = false;
    vector<set<int> > alginedPairIdxs(m_inputClusters->size());
    if (m_inputAlignedPairs->size() > 0)
    {
      DEBUG_MSG("Filtering aligned pairs from input pairs ...");
      filterAlignedPairs = true;
      for (int i = 0; i < m_inputAlignedPairs->size(); i++)
      {
        int clust1 = (*m_inputAlignedPairs)[i].spec1;
        int clust2 = (*m_inputAlignedPairs)[i].spec2;

        if (clust1 > clust2)
        {
          int temp = clust1;
          clust1 = clust2;
          clust2 = temp;
        }

        alginedPairIdxs[clust1].insert(clust2);
      }
    }

    DEBUG_MSG("Building tag index of FASTA_DATABASE...");
    //m_inputDB->populateSequenceIndex(minPeptideLength);
    DEBUG_MSG("finished");

    DEBUG_VAR(m_inputPairs->size());

    int numSkippedOvlp = 0, numSkippedMatch = 0, numSkippedRatio = 0,
        numSkippedAlignGF = 0;

    for (int i = 0; i < m_inputPairs->size(); i++)
    {
      int spec1 = (*m_inputPairs)[i].spec1;
      int spec2 = (*m_inputPairs)[i].spec2;
      int shift = floatToInt((*m_inputPairs)[i].shift1);

      if (spec1 < 0 || spec2 < 0)
      {
        ERROR_MSG("Invalid index");
        return false;
      }
      const string pep1 = peptidesUse[spec1];
      const string pep2 = peptidesUse[spec2];

      bool debug = false;
      if (((*m_inputSpectra)[spec1].scan == DEBUG_SCAN1
          && (*m_inputSpectra)[spec2].scan == DEBUG_SCAN2)
          || ((*m_inputSpectra)[spec1].scan == DEBUG_SCAN2
              && (*m_inputSpectra)[spec2].scan == DEBUG_SCAN1))
      {
        DEBUG_VAR((*m_inputSpectra)[spec1].scan);
        DEBUG_VAR((*m_inputSpectra)[spec2].scan);
        DEBUG_VAR(pep1);
        DEBUG_VAR(pep2);
        debug = true;
      }

      if ((!allowSameChargePairs) && pep1 == pep2 && shift == 0
          && (*m_inputSpectra)[spec1].parentCharge
              == (*m_inputSpectra)[spec2].parentCharge)
      {
        continue;
      }
      else if ((!allowDiffChargePairs) && pep1 == pep2 && shift == 0
          && (*m_inputSpectra)[spec1].parentCharge
              != (*m_inputSpectra)[spec2].parentCharge)
      {
        continue;
      }
      else if ((!allowPartialPairs) && (pep1 != pep2 || shift != 0))
      {
        continue;
      }

      if ((!psmsAllowed[spec1]) || (!psmsAllowed[spec2]))
      {
        if (debug)
          DEBUG_MSG("Skipping bc psm is not allowed");

        continue;
      }

      if (maxAlignGFPval < 1 && (*m_inputPairs)[i].score1 > maxAlignGFPval)
      {
        numSkippedAlignGF++;
        continue;
      }

      float totalOverlap1 = 0;
      float totalOverlap2 = 0;
      float totalMatched1 = 0;
      float totalMatched2 = 0;
      float floatShift = (*m_inputPairs)[i].shift1;

      for (int pIdx = 0; pIdx < (*m_inputSpectra)[spec2].size(); pIdx++)
      {
        const float offsetMass = (*m_inputSpectra)[spec2][pIdx][0] + floatShift;
        if (offsetMass > 0
            && offsetMass
                < (*m_inputSpectra)[spec1][(*m_inputSpectra)[spec1].size() - 1][0])
        {
          totalOverlap2 += (*m_inputSpectra)[spec2][pIdx][1];

          if ((*m_inputSpectra)[spec1].findPeaks(offsetMass, 0.5) >= 0)
          {
            totalMatched2 += (*m_inputSpectra)[spec2][pIdx][1];
          }
        }
      }

      for (int pIdx = 0; pIdx < (*m_inputSpectra)[spec1].size(); pIdx++)
      {
        const float offsetMass = (*m_inputSpectra)[spec1][pIdx][0] - floatShift;
        if (offsetMass > 0
            && offsetMass
                < (*m_inputSpectra)[spec2][(*m_inputSpectra)[spec2].size() - 1][0])
        {
          totalOverlap1 += (*m_inputSpectra)[spec1][pIdx][1];

          if ((*m_inputSpectra)[spec2].findPeaks(offsetMass, 0.5) >= 0)
          {
            totalMatched1 += (*m_inputSpectra)[spec1][pIdx][1];
          }
        }
      }

      const float ratio1 = totalMatched1 / totalOverlap1;
      const float ratio2 = totalMatched2 / totalOverlap2;

      if (debug)
      {
        DEBUG_VAR(ratio1);
        DEBUG_VAR(ratio2);
      }

      if (min(ratio1, ratio1) < minRatio)
      {
        numSkippedRatio++;
        continue;
      }

      //if ((*m_inputSpectra)[spec1].psmList.front()->m_shared_peaks != 1
      //    || (*m_inputSpectra)[spec2].psmList.front()->m_shared_peaks != 1)
      //{
      //  continue;
      // }

      int clust1 = specToCluster[spec1];
      int clust2 = specToCluster[spec2];

      int smallClustIdx = min(clust1, clust2);
      int bigClustIdx = max(clust1, clust2);

      if (filterAlignedPairs
          && alginedPairIdxs[smallClustIdx].count(bigClustIdx) == 0)
      {
        if (debug)
          DEBUG_MSG("Skipping here");

        continue;
      }

      int bestStart1 = -1, bestStart2 = -1, maxOverlap = -1;

      bool foundOverlap = false;
      int maxStart2 = pep2.length() - minPeptideLength;
      for (int start2 = 0; start2 <= maxStart2; start2++)
      {
        int remain2 = pep2.length() - start2;
        remain2 = min(remain2, (int)pep1.length());

        bool res = true;

        for (int strI = 0; strI < remain2; strI++)
        {
          if (pep1[strI] != pep2[strI + start2])
          {
            res = false;
            break;
          }
        }

        if (res && remain2 > maxOverlap)
        {
          maxOverlap = remain2;
          bestStart1 = 0;
          bestStart2 = start2;
        }
      }
      int maxStart1 = pep1.length() - minPeptideLength;
      for (int start1 = 0; start1 <= maxStart1; start1++)
      {
        int remain1 = pep1.length() - start1;
        remain1 = min(remain1, (int)pep2.length());

        bool res = true;

        for (int strI = 0; strI < remain1; strI++)
        {
          if (pep1[strI + start1] != pep2[strI])
          {
            res = false;
            break;
          }
        }

        if (res && remain1 > maxOverlap)
        {
          maxOverlap = remain1;
          bestStart1 = start1;
          bestStart2 = 0;
        }
      }

      int numPks1 = 0, numPks2 = 0, numMatchedPks = 0;
      for (int aaI = 0; aaI < maxOverlap; aaI++)
      {
        if (matchedPRMs[spec1][aaI + bestStart1])
        {
          numPks1++;
        }
        if (matchedPRMs[spec2][aaI + bestStart2])
        {
          numPks2++;
        }
        if (matchedPRMs[spec1][aaI + bestStart1]
            && matchedPRMs[spec2][aaI + bestStart2])
        {
          numMatchedPks++;
        }
      }

      if (numMatchedPks < numMatchedPeaks)
      {
        numSkippedOvlp++;
        continue;
      }

      /*if (clust1 == 170801 && clust2 == 50836)
       {
       DEBUG_VAR(clust1);
       DEBUG_VAR(clust2);
       (*m_inputSpectra)[spec1].output(cerr);
       DEBUG_VAR(spec1);
       DEBUG_VAR(spec2);
       DEBUG_VAR(pep1);
       DEBUG_VAR(pep2);
       Spectrum prmMasses;
       jumps.getPRMMasses(pep1, prmMasses);
       prmMasses.output(cerr);
       DEBUG_VAR(numPks1);
       DEBUG_VAR(numPks2);
       DEBUG_VAR(numMatchedPks);
       for (int pI = 0; pI < matchedPRMs[spec1].size(); pI++)
       {
       DEBUG_MSG("PRM idx " << pI << " = " << matchedPRMs[spec1][pI]);
       }
       DEBUG_VAR(maxOverlap);
       DEBUG_VAR(bestStart1);
       DEBUG_VAR(bestStart2);

       DEBUG_VAR((*m_inputSpectra)[spec1].psmList.front()->m_strict_envelope_score);
       DEBUG_VAR((*m_inputSpectra)[spec1].psmList.front()->m_pValue);

       DEBUG_VAR((*m_inputSpectra)[spec2].psmList.front()->m_strict_envelope_score);
       DEBUG_VAR((*m_inputSpectra)[spec2].psmList.front()->m_pValue);
       abort();
       }*/

      matchingPeaks[i] = numMatchedPks;

      const string &prot1 = (*m_inputSpectra)[spec1].psmList.front()->m_protein;
      const string &prot2 = (*m_inputSpectra)[spec2].psmList.front()->m_protein;
      if (prot1.find(decoyID) == string::npos
          && prot2.find(decoyID) == string::npos)
      {
        string pairStr("");
        for (int aaI = 0; aaI < bestStart1; aaI++)
        {
          pairStr += pep1.at(aaI);
        }
        for (int aaI = 0; aaI < bestStart2; aaI++)
        {
          pairStr += pep2.at(aaI);
        }
        if (pep1.length() - bestStart1 > pep2.length() - bestStart2)
        {
          for (int aaI = bestStart1; aaI < pep1.length(); aaI++)
          {
            pairStr += pep1.at(aaI);
          }
        }
        else
        {
          for (int aaI = bestStart2; aaI < pep2.length(); aaI++)
          {
            pairStr += pep2.at(aaI);
          }
        }
        //DEBUG_VAR(pairStr);

        /*if (!m_inputDB->hasSequenceFromIndex(pairStr))
         {
         numSkippedMatch++;

         if (debug)
         DEBUG_MSG("Skipping here");

         continue;
         }*/
      }

      inputAligns[spec1][spec2] = i;
      inputAligns[spec2][spec1] = i;

      if (debug)
      {
        DEBUG_MSG("Adding pair " << spec1 << "," << spec2 << " - clusters = " << clust1 << "," << clust2);
        DEBUG_VAR(bestSpecPerCluster[clust1].second);
        DEBUG_VAR((*m_inputProbs)[i][0]);
        DEBUG_VAR(bestSpecPerCluster[clust2].second);
        DEBUG_VAR((*m_inputProbs)[i][1]);
      }

      if ((*m_inputProbs)[i][0] < bestSpecPerCluster[clust1].second)
      {
        bestSpecPerCluster[clust1].first = spec1;
        bestSpecPerCluster[clust1].second = (*m_inputProbs)[i][0];
      }

      if ((*m_inputProbs)[i][1] < bestSpecPerCluster[clust2].second)
      {
        bestSpecPerCluster[clust2].first = spec2;
        bestSpecPerCluster[clust2].second = (*m_inputProbs)[i][1];
      }
    }

    DEBUG_VAR(numSkippedOvlp);
    DEBUG_VAR(numSkippedMatch);
    DEBUG_VAR(numSkippedRatio);
    DEBUG_VAR(numSkippedAlignGF);

    DEBUG_TRACE;

    int numSkippedHere1 = 0, numSkippedHere2 = 0, numSkippedHere3 = 0;

    vector<PeptideStarMatch> networkPSMs;
    vector<PeptideStarMatch> nonNetworkPSMs;
    for (int i = 0; i < m_inputClusters->size(); i++)
    {
      int specIdx = bestSpecPerCluster[i].first;

      if (specIdx < 0)
      {
        (*m_outputSpectra)[i].resize(0);
        numSkippedHere1++;
        continue;
      }

      if ((*m_inputSpectra)[specIdx].psmList.size() == 0)
      {
        (*m_outputSpectra)[i].resize(0);
        numSkippedHere2++;
        continue;
      }

      if ((*m_outputPSMs)[specIdx].specIdx < 0)
      {
        (*m_outputSpectra)[i].resize(0);
        numSkippedHere3++;
        continue;
      }

      (*m_outputSpectra)[i] = (*m_inputSpectra)[specIdx];
      (*m_outputSpectra)[i].psmList.clear();

      PeptideStarMatch nextPSM;
      nextPSM = (*m_outputPSMs)[specIdx];

      bool debug = false;
      if ((*m_inputSpectra)[specIdx].scan == DEBUG_SCAN1
          || (*m_inputSpectra)[specIdx].scan == DEBUG_SCAN2)
      {
        DEBUG_VAR(specIdx);
        DEBUG_VAR((*m_outputPSMs)[specIdx].peptide);
        DEBUG_VAR((*m_outputPSMs)[specIdx].specProb);
        DEBUG_VAR((*m_outputPSMs)[specIdx].starProb);
        debug = true;
      }

      int otherIdx = -1;
      int otherClustIdx = -1;
      double bestPVal = 1.0;
      float bestAlignGFPval = 1.0;
      double otherSpecProb = -1;
      double otherStarProb = -1;
      int otherMP = -1;
      for (map<int, int>::const_iterator mIt = inputAligns[specIdx].begin();
          mIt != inputAligns[specIdx].end(); mIt++)
      {
        int alignIdx = mIt->second;
        int spec2 = mIt->first;
        int clust2 = specToCluster[spec2];

        if (debug)
        {
          DEBUG_MSG("alignIdx = " << alignIdx << ", spec2 = " << spec2 << ", clust2 = " << clust2);
        }
        if (spec2 != bestSpecPerCluster[clust2].first)
        {
          continue;
        }
        double pval =
            ((*m_inputPairs)[alignIdx].spec1 == specIdx) ?
                (*m_inputProbs)[alignIdx][0] : (*m_inputProbs)[alignIdx][1];
        double otherPval =
            ((*m_inputPairs)[alignIdx].spec1 == specIdx) ?
                (*m_inputProbs)[alignIdx][1] : (*m_inputProbs)[alignIdx][0];

        if (debug)
        {
          DEBUG_MSG("pval = " << pval << ", otherPval = " << otherPval << ", bestPVal = " << bestPVal);
        }

        if (pval < bestPVal)
        {
          bestPVal = pval;
          otherIdx = spec2;
          otherClustIdx = clust2;
          otherStarProb = otherPval;
          otherSpecProb = (*m_outputPSMs)[spec2].specProb;
          otherMP = matchingPeaks[alignIdx];
          bestAlignGFPval = (*m_inputPairs)[alignIdx].score1;
        }
      }

      if (otherIdx >= 0)
      {
        if (debug)
        {
          DEBUG_MSG("Pairing with " << otherIdx << " ( scan = " << (*m_inputSpectra)[otherIdx].scan << " ) ");
          DEBUG_VAR(bestPVal);
        }
        nextPSM.starProb = bestPVal;
        nextPSM.bestIntersectPVIdx = otherClustIdx;
        if (decoyID.length() > 0)
        {
          nextPSM.pairIsDecoy =
              (*m_inputSpectra)[otherIdx].psmList.front()->m_protein.find(decoyID)
                  != string::npos;
        }
        nextPSM.pairSpecProb = otherSpecProb;
        nextPSM.pairStarProb = otherStarProb;
        nextPSM.bestPairScan = (*m_inputClusters)[otherClustIdx].m_scan;
        nextPSM.bestPairFile = (*m_inputClusters)[otherClustIdx].m_filename;
        nextPSM.bestPairPeptide =
            (*m_inputSpectra)[otherIdx].psmList.front()->m_annotation;
        nextPSM.pairRank =
            (*m_inputSpectra)[otherIdx].psmList.front()->m_shared_peaks;
        nextPSM.pairMP = otherMP;
        nextPSM.alignGFProb = bestAlignGFPval;
        networkPSMs.push_back(nextPSM);
      }
      else
      {
        if (debug)
        {
          DEBUG_MSG("Found no pair");
        }
        nextPSM.starProb = nextPSM.specProb;
        nonNetworkPSMs.push_back(nextPSM);
      }

    }
    DEBUG_VAR(numSkippedHere1);
    DEBUG_VAR(numSkippedHere2);
    DEBUG_VAR(numSkippedHere3);

    DEBUG_VAR(networkPSMs.size());
    DEBUG_VAR(nonNetworkPSMs.size());

    if (networkOnly)
    {
      m_outputPSMs->resize(networkPSMs.size());
    }
    else
    {
      m_outputPSMs->resize(networkPSMs.size() + nonNetworkPSMs.size());
    }

    for (int i = 0; i < networkPSMs.size(); i++)
    {
      (*m_outputPSMs)[i] = networkPSMs[i];
    }

    const int networkCutoff = networkPSMs.size();
    int idxUse = networkCutoff;
    if (!networkOnly)
    {
      for (int i = 0; i < nonNetworkPSMs.size(); i++)
      {
        (*m_outputPSMs)[idxUse++] = nonNetworkPSMs[i];
      }
    }

    DEBUG_VAR(m_outputPSMs->size());
    set<int> identified_spectra;

    if (decoyID.length() > 0)
    {
      PeptideSpectrumMatchSet netPSMs;
      //PeptideSpectrumMatchSet nonNetPSMs;

      DEBUG_TRACE;
      for (int i = 0; i < m_outputPSMs->size(); i++)
      {

        if (networkOnly && i >= networkCutoff)
        {
          continue;
        }

        psmPtr currMatch(new PeptideSpectrumMatch);
        (*m_outputPSMs)[i].setToPSM(*currMatch);
        currMatch->m_dbIndex = i;
        currMatch->m_isDecoy = currMatch->m_protein.find(decoyID)
            != string::npos;

        if (fdrMetric == "msgf")
        {
          currMatch->m_pValue = (*m_outputPSMs)[i].msgfSpecProb;
        }
        else if (fdrMetric == "specprob")
        {
          currMatch->m_pValue = (*m_outputPSMs)[i].specProb;
        }
        else if ((*m_outputPSMs)[i].specProb < (*m_outputPSMs)[i].starProb)
        {
          currMatch->m_pValue = (*m_outputPSMs)[i].specProb;
        }
        netPSMs.push_back(currMatch);

        /*psmPtr currMatch2(new PeptideSpectrumMatch);
         currMatch2->operator =(*currMatch);
         currMatch2->m_pValue = (*m_outputPSMs)[i].specProb;
         nonNetPSMs.push_back(currMatch2);
         */

        if ((*m_outputPSMs)[i].scan == DEBUG_SCAN1
            || (*m_outputPSMs)[i].scan == DEBUG_SCAN2)
        {
          DEBUG_VAR((*m_outputPSMs)[i].scan);
          DEBUG_VAR((*m_outputPSMs)[i].peptide);
          DEBUG_VAR((*m_outputPSMs)[i].specProb);
          DEBUG_VAR((*m_outputPSMs)[i].starProb);
          DEBUG_VAR((*m_outputPSMs)[i].bestIntersectPVIdx);
        }
      }

      DEBUG_VAR(netPSMs.size());
      // DEBUG_VAR(nonNetPSMs.size());

      PeptideSpectrumMatchSet netPSMsTDA;
      //PeptideSpectrumMatchSet nonNetPSMsTDA;

      DEBUG_TRACE;

      FdrPeptide::calculatePValues(netPSMs,
                                   netPSMsTDA,
                                   1,
                                   FDR_SORT_PVALUE,
                                   true);

      DEBUG_TRACE;

      /*
       FdrPeptide::calculatePValues(nonNetPSMs,
       nonNetPSMsTDA,
       1,
       compareFDRPVal,
       true);

       for (int i = 0; i < netPSMsTDA.size(); i++)
       {
       netPSMs[netPSMsTDA[i]->m_dbIndex]->m_pValue = netPSMsTDA[i]->m_fdr;
       }

       if (!networkOnly)
       {
       for (int i = 0; i < nonNetPSMsTDA.size(); i++)
       {
       if (nonNetPSMsTDA[i]->m_dbIndex >= networkCutoff)
       {
       netPSMs[nonNetPSMsTDA[i]->m_dbIndex]->m_pValue =
       nonNetPSMsTDA[i]->m_fdr;
       }
       }
       }

       DEBUG_TRACE;
       FdrPeptide::calculatePValues(netPSMs,
       netPSMsTDA,
       1,
       compareFDRPVal,
       true);
       */

      DEBUG_TRACE;

      int numFDR = 0, numNet = 0, numNonNet = 0;
      double maxCotuff = 0;
      tr1::unordered_set<string> foundPeptides;
      tr1::unordered_set<string> foundPeptidesNet;
      tr1::unordered_set<string> foundPeptidesNonNet;

      for (int i = 0; i < netPSMsTDA.size(); i++)
      {
        int pStarMIdx = netPSMsTDA[i]->m_dbIndex;

        if (networkOnly && pStarMIdx >= networkCutoff)
        {
          (*m_outputPSMs)[pStarMIdx].specIdx = -1;
          continue;
        }

        (*m_outputPSMs)[pStarMIdx].fdr = netPSMsTDA[i]->m_fdr;
        (*m_outputPSMs)[pStarMIdx].pepFdr = netPSMsTDA[i]->m_pepFdr;
        (*m_outputPSMs)[pStarMIdx].isDecoy = (netPSMsTDA[i]->m_isDecoy) ? 1 : 0;

        if ((*m_outputPSMs)[pStarMIdx].fdr <= 0.01
            && !(*m_outputPSMs)[pStarMIdx].isDecoy)
        {
          numFDR++;
          maxCotuff = max(maxCotuff, netPSMsTDA[i]->m_pValue);
          if (pStarMIdx >= networkCutoff)
          {
            numNonNet++;
          }
          else
          {
            numNet++;
          }
          identified_spectra.insert((*m_outputPSMs)[pStarMIdx].scan);
        }

        if ((*m_outputPSMs)[pStarMIdx].pepFdr <= 0.01
            && !(*m_outputPSMs)[pStarMIdx].isDecoy)
        {
          foundPeptides.insert((*m_outputPSMs)[pStarMIdx].peptide);
          if (pStarMIdx >= networkCutoff)
          {
            foundPeptidesNonNet.insert((*m_outputPSMs)[pStarMIdx].peptide);
          }
          else
          {
            foundPeptidesNet.insert((*m_outputPSMs)[pStarMIdx].peptide);
          }
        }
      }
      if (fdrMetric == "msgf")
      {
        sort(m_outputPSMs->begin(), m_outputPSMs->end(), compareMSGFProbs);
      }
      else if (fdrMetric == "specprob")
      {
        sort(m_outputPSMs->begin(), m_outputPSMs->end(), compareSpecProbs);
      }
      else
      {
        sort(m_outputPSMs->begin(), m_outputPSMs->end(), compareStarProbs);
      }

      DEBUG_MSG("Have " << numFDR << " PSMs at 1% spectrum-level FDR (" << numNet << " in network and " << numNonNet << " out of network)");
      DEBUG_MSG("Have " << foundPeptides.size() << " peptides at 1% peptide-level FDR (" << foundPeptidesNet.size() << " in network and " << foundPeptidesNonNet.size() << " out of network)");
      DEBUG_MSG("Maximum probability cutoff is " << maxCotuff << " at 1% spectrum-level FDR");
    }

    if (!saveUnidentifiedMS2(identified_spectra))
    {
      return false;
    }

    return true;
  }

  bool ExecStarGF::loadInputData(void)
  {
    if (m_params.exists("FASTA_DATABASE"))
    {
      DEBUG_MSG("Loading FASTA_DATABASE from \'" << m_params.getValue("FASTA_DATABASE") << "\' ...");
      if (m_inputDB->Load(m_params.getValue("FASTA_DATABASE").c_str()) == 0)
      {
        ERROR_MSG("Failed to db from \'" << m_params.getValue("FASTA_DATABASE") << "\'");
        return false;
      }
    }

    if (m_params.exists("INPUT_SPECTRA"))
    {
      DEBUG_MSG("Loading INPUT_SPECTRA from \'" << m_params.getValue("INPUT_SPECTRA") << "\' ...");
      if (!m_inputSpectra->Load(m_params.getValue("INPUT_SPECTRA").c_str()))
      {
        ERROR_MSG("Failed to load spectra from \'" << m_params.getValue("INPUT_SPECTRA") << "\'");
        return false;
      }
      DEBUG_VAR(m_inputSpectra->size());
      for (int i = 0; i < m_inputSpectra->size(); i++)
      {
        (*m_inputSpectra)[i].psmList.clear();
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!!!");
      return false;
    }

    if (m_params.exists("INPUT_PAIRS"))
    {
      DEBUG_MSG("Loading INPUT_PAIRS from \'" << m_params.getValue("INPUT_PAIRS") << "\' ...");
      if (!m_inputPairs->loadFromBinaryFile(m_params.getValue("INPUT_PAIRS")))
      {
        ERROR_MSG("Failed to load spectrum pairs from \'" << m_params.getValue("INPUT_PAIRS") << "\'");
        return false;
      }
    }

    if (m_params.exists("INPUT_ALIGNED_PAIRS"))
    {
      DEBUG_MSG("Loading INPUT_ALIGNED_PAIRS from \'" << m_params.getValue("INPUT_ALIGNED_PAIRS") << "\' ...");
      if (!m_inputAlignedPairs->loadFromBinaryFile(m_params.getValue("INPUT_ALIGNED_PAIRS")))
      {
        ERROR_MSG("Failed to load spectrum pairs from \'" << m_params.getValue("INPUT_ALIGNED_PAIRS") << "\'");
        return false;
      }
    }

    if (m_params.exists("INPUT_PROBS"))
    {
      string filename = m_params.getValue("INPUT_PROBS");
      DEBUG_MSG("Loading intersecting probabilities from \'" << filename << "\'");
      if (Load_binArray<double>(filename.c_str(), *m_inputProbs) != 1)
      {
        ERROR_MSG("Failed to load intersecting probabilities from \'" << filename << "\'");
        return false;
      }
    }

    if (m_params.exists("INPUT_CLUSTERS"))
    {
      if (!m_inputClusters->loadBinaryFile(m_params.getValue("INPUT_CLUSTERS")))
      {
        ERROR_MSG("Failed to load clusters from \'" << m_params.getValue("INPUT_CLUSTERS") << "\'");
        return false;
      }
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("No input spectra!!!");
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

      int numAdded = m_inputPSMs->addMostSpectra(m_inputSpectra, true, false);

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

  bool ExecStarGF::saveOutputData(void)
  {
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::saveSpecset(m_params.getValue("OUTPUT_SPECTRA"),
              m_outputSpectra,
              true))
      {
        ERROR_MSG("Failed to save spectra to \'"
            << m_params.getValue("OUTPUT_SPECTRA") << "\'");
        return false;
      }
    }
    if (m_params.exists("OUTPUT_STAR_PSMS"))
    {
      string outPSMsFile = m_params.getValue("OUTPUT_STAR_PSMS");
      if (!SavePeptideStarMatches(outPSMsFile, *m_outputPSMs))
      {
        ERROR_MSG("Failed to save intersecting probabilities to \'"
            << outPSMsFile << "\'");
        return false;
      }
    }
    return true;
  }

  bool ExecStarGF::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  bool ExecStarGF::loadOutputData(void)
  {
    return false;
  }

  vector<ExecBase*> const & ExecStarGF::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  bool ExecStarGF::merge(void)
  {
    return false;
  }

  bool ExecStarGF::validateParams(std::string & error)
  {
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
    //VALIDATE_PARAM_EXIST("INPUT_PAIRS");
    //VALIDATE_PARAM_EXIST("INPUT_PROBS");
    VALIDATE_PARAM_EXIST("INPUT_PSMS");
    VALIDATE_PARAM_EXIST("INPUT_CLUSTERS");
    VALIDATE_PARAM_EXIST("FASTA_DATABASE");

    return true;
  }

  bool ExecStarGF::saveUnidentifiedMS2(const std::set<int> & identified_psm_scans)
  {
    SpecSet ms2Specs;
    if (m_params.exists("INPUT_SPECTRA_MS"))
    {
      const string filename = m_params.getValue("INPUT_SPECTRA_MS");
      DEBUG_MSG("Loading \'" << filename << "\' ...");
      if (!ms2Specs.Load(filename.c_str()))
      {
        ERROR_MSG("Failed to load \'" << filename << "\'!");
        return false;
      }
    }
    else
    {
      return true;
    }

    vector<int> msScanToIdx(ms2Specs.size(), -1);
    for (int i = 0; i < ms2Specs.size(); i++)
    {
      const int msScan = ms2Specs[i].scan;

      if (msScan >= msScanToIdx.size())
      {
        msScanToIdx.resize(msScan + 1, -1);
      }

      if (msScanToIdx[msScan] > -1)
      {
        ERROR_MSG("Duplicate scans " << msScanToIdx[msScan]);
        abort();
      }
      msScanToIdx[msScan] = i;
    }

    for (set<int>::const_iterator sIt = identified_psm_scans.begin();
        sIt != identified_psm_scans.end(); sIt++)
    {
      int msIdx = msScanToIdx[*sIt];

      ms2Specs[msIdx].resize(0);
    }

    if (m_params.exists("OUTPUT_UNIDENTIFIED_SPECTRA_MS"))
    {
      const string filename =
          m_params.getValue("OUTPUT_UNIDENTIFIED_SPECTRA_MS");
      DEBUG_MSG("Saving \'" << filename << "\' ...");
      if (!ExecMergeConvert::saveSpecset(filename, &ms2Specs, true))
      {
        ERROR_MSG("Failed to save \'" << filename << "\'!");
        return false;
      }
    }

    return true;
  }
}

