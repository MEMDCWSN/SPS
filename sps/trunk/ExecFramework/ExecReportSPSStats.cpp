/*
 * ExecReportSPSStats.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

// Header Includes
#include "ExecReportSPSStats.h"

using namespace std;
using namespace spsReports;

namespace specnets
{

  ExecReportSPSStats::ExecReportSPSStats(void) :
    ExecBase(), m_contigs(0x0), m_spectraStars(0x0), m_abinfo(0x0),
        m_contigOverlaps(0x0), m_contigProtMatch(0x0), m_fasta(0x0),
        m_model(0x0), m_spectraClustIDs(0x0), m_spectraRawIDs(0x0),
        m_mappedProj(0x0), m_targetProts(0x0), m_proteinSpectra(0x0),
        m_peptides(0x0), ownInput(true), ownOutput(true)
  {
    m_name = "ExecReportSPSStats";
    m_type = "ExecReportSPSStats";
    m_contigs = new SpecSet;
    m_spectraStars = new SpecSet;
    m_abinfo = new abinfo_t;
    m_contigOverlaps = new SpecSet;
    m_contigProtMatch = new vector<vector<int> > ;
    m_fasta = new DB_fasta;
    m_model = new MS2ScoringModel;
    m_spectraClustIDs = new PeptideSpectrumMatchSet;
    m_spectraRawIDs = new PeptideSpectrumMatchSet;
    m_mappedProj = new MappedSpecnets;
  }

  ExecReportSPSStats::ExecReportSPSStats(const ParameterList & inputParams) :

    ExecBase(inputParams), m_contigs(0x0), m_spectraStars(0x0), m_abinfo(0x0),
        m_contigOverlaps(0x0), m_contigProtMatch(0x0), m_fasta(0x0),
        m_model(0x0), m_spectraClustIDs(0x0), m_spectraRawIDs(0x0),
        m_mappedProj(0x0), m_targetProts(0x0), m_proteinSpectra(0x0),
        m_peptides(0x0), ownInput(true), ownOutput(true)
  {

    m_name = "ExecReportSPSStats";
    m_type = "ExecReportSPSStats";
    m_contigs = new SpecSet;
    m_spectraStars = new SpecSet;
    m_abinfo = new abinfo_t;
    m_contigOverlaps = new SpecSet;
    m_contigProtMatch = new vector<vector<int> > ;
    m_fasta = new DB_fasta;
    m_model = new MS2ScoringModel;
    m_spectraClustIDs = new PeptideSpectrumMatchSet;
    m_spectraRawIDs = new PeptideSpectrumMatchSet;
    m_mappedProj = new MappedSpecnets;
  }

  ExecReportSPSStats::~ExecReportSPSStats(void)
  {
    if (ownInput)
    {
      delete m_contigs;
      delete m_spectraStars;
      delete m_abinfo;
      delete m_contigOverlaps;
      delete m_contigProtMatch;
      delete m_fasta;
      delete m_model;
      delete m_spectraClustIDs;
      delete m_spectraRawIDs;
      delete m_targetProts;
      delete m_proteinSpectra;
      delete m_peptides;
    }
    if (ownOutput)
    {
      delete m_mappedProj;
    }
  }

  ExecBase * ExecReportSPSStats::clone(const ParameterList & inputParams) const
  {
    return new ExecReportSPSStats(inputParams);
  }

  bool ExecReportSPSStats::invoke(void)
  {
    DEBUG_TRACE;

    int minContigTag = m_params.getValueInt("MIN_CONTIG_AA_TAG", -1);
    int minDBmp = m_params.getValueInt("MIN_CONTIG_DB_MP", -1);
    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

    AAJumps jumps(1);

    int endsChop = m_params.getValueInt("ENDS_CHOP", 0);

    if (endsChop > 0)
    {
      DEBUG_MSG("Ignoring the first and last " << endsChop << " abruijn vertice(s) of every contig");
      ChopContigEnds(endsChop);
    }

    unsigned int numContigs = 0;
    for (int i = 0; i < m_contigs->size(); i++)
    {
      numContigs += ((*m_contigs)[i].size() > 0) ? 1 : 0;
    }

    DEBUG_VAR(numContigs);

    if (minContigTag > 0)
    {
      DEBUG_MSG("Ignoring any contig with an AA tag less than " << minContigTag << " residues");
      for (int i = 0; i < m_contigs->size(); i++)
      {
        if ((*m_contigs)[i].size() == 0)
        {
          continue;
        }
        vector<MZRange> c_masses((*m_contigs)[i].size());
        for (int j = 0; j < c_masses.size(); j++)
        {
          c_masses[j].set((*m_contigs)[i][j][0], (*m_contigs)[i][j][1], peakTol);
        }
        string tag = getLongestTag(c_masses, jumps);

        //DEBUG_MSG("Contig " << i << " has tag " << tag);
        int tagLen = tag.length();
        if (tagLen < minContigTag)
        {
          (*m_contigs)[i].resize(0);
        }
      }

      numContigs = 0;
      for (int i = 0; i < m_contigs->size(); i++)
      {
        numContigs += ((*m_contigs)[i].size() > 0) ? 1 : 0;
      }

      DEBUG_VAR(numContigs);
    }
    if (minDBmp > 0)
    {
      DEBUG_MSG("Ignoring any contig with less than " << minDBmp << " consecutive residues in the database");
      for (int i = 0; i < m_contigs->size(); i++)
      {
        if ((*m_contigs)[i].size() == 0)
        {
          continue;
        }
        int numConsecMP = 0;
        int lastIdx = -2;
        bool haveMP = false;
        for (unsigned int j = 0; j < (*m_contigOverlaps)[i].size(); j++)
        {
          int nextIdx = floatToInt((*m_contigOverlaps)[i][j][0]);
          if (nextIdx == lastIdx + 1)
          {
            numConsecMP++;
            if (numConsecMP >= minDBmp)
            {
              haveMP = true;
              break;
            }
          }
          else
          {
            numConsecMP = 1;
          }
          lastIdx = nextIdx;
        }
        if (!haveMP)
        {
          (*m_contigs)[i].resize(0);
        }
      }

      numContigs = 0;
      for (int i = 0; i < m_contigs->size(); i++)
      {
        numContigs += ((*m_contigs)[i].size() > 0) ? 1 : 0;
      }

      DEBUG_VAR(numContigs);
    }
    /*
     SpecSet* sps_contigs,
     abinfo_t* sps_components,
     SpecSet* star_spectra,
     SpecSet* matchma_overlaps,
     vector<vector<int> >* matchma_prot_idx,
     vector<string>* _proteins,
     SpecSet* protein_spectra,
     vector<string>* spectrum_ids,
     MS2ScoringModel& input_ion_types,
     float peak_tol
     */
    //DEBUG_MSG("invoking");
    DEBUG_TRACE;
    m_mappedProj->mapProt(m_contigs,
                          m_abinfo,
                          m_spectraStars,
                          m_contigOverlaps,
                          m_contigProtMatch,
                          m_targetProts,
                          m_proteinSpectra,
                          m_peptides,
                          *m_model,
                          peakTol);
    DEBUG_TRACE;
    //DEBUG_MSG("done invoking");
    return true;
  }

  bool ExecReportSPSStats::loadInputData(void)
  {

    string spsDir = m_params.getValue("SPS_PROJECT_DIR");

    if (!m_contigs->Load(m_params.getValue("INPUT_CONTIGS").c_str(), NULL))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIGS"));
      return false;
    }
    DEBUG_VAR(m_contigs->size());

    if (m_contigs->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_CONTIGS or INPUT_CONTIGS_PKLBIN?");
      return false;
    }

    if (!m_spectraStars->Load(m_params.getValue("INPUT_STARS").c_str(), NULL))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_STARS"));
      return false;
    }

    for (unsigned int i = 0; i < m_spectraStars->size(); i++)
    {
      (*m_spectraStars)[i].scan = i + 1;
    }

    if (m_spectraStars->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_STARS?");
      return false;
    }

    DEBUG_MSG("Loading abinfo from <" << m_params.getValue("INPUT_CONTIG_ABINFO") << "> ...");
    if (!Load_abinfo(m_params.getValue("INPUT_CONTIG_ABINFO").c_str(),
                     *m_abinfo))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIG_ABINFO"));
      return false;
    }
    if (m_abinfo->size() == 0)
    {
      ERROR_MSG("Input abinfo size is 0!, did you specify INPUT_CONTIG_ABINFO?");
      return false;
    }

    DEBUG_MSG("Loading matched protein indices from <" << m_params.getValue("INPUT_MP") << "> ...");
    if (!Load_binArray(m_params.getValue("INPUT_MP").c_str(),
                       *m_contigProtMatch))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_MP"));
      return false;
    }
    DEBUG_VAR(m_contigProtMatch->size());

    DEBUG_MSG("Loading matched peak indices from <" << m_params.getValue("INPUT_MIDX") << "> ...");
    if (!m_contigOverlaps->Load(m_params.getValue("INPUT_MIDX").c_str(), NULL))
    {
      ERROR_MSG("Could not load " << m_params.getValue("INPUT_MIDX"));
      return false;
    }
    DEBUG_VAR(m_contigOverlaps->size());
    if (m_contigOverlaps->size() == 0)
    {
      ERROR_MSG("Input matched peak indices size is 0!");
      return false;
    }

    if (m_params.exists("INPUT_REF_INDICES"))
    {
      DEBUG_MSG("Loading mapped contig indices from <" << m_params.getValue("INPUT_REF_INDICES") << "> ...");
      vector<vector<int> > refIdxs;
      if (!Load_binArray(m_params.getValue("INPUT_REF_INDICES").c_str(),
                         refIdxs))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_REF_INDICES"));
        return false;
      }

      SpecSet overlapsAll(m_contigs->size());
      vector<vector<int> > protMatchAll(m_contigs->size());

      for (unsigned int i = 0; i < refIdxs.size(); i++)
      {
        overlapsAll[refIdxs[i][0]] = (*m_contigOverlaps)[i];
        protMatchAll[refIdxs[i][0]] = (*m_contigProtMatch)[i];
      }

      for (unsigned int i = 0; i < m_contigs->size(); i++)
      {
        if (protMatchAll[i].size() == 0)
        {
          overlapsAll[i].resize(0);
          protMatchAll[i].assign(3, -1);
        }
      }
      m_contigOverlaps->operator =(overlapsAll);
      m_contigProtMatch->operator =(protMatchAll);
    }

    if (m_contigProtMatch->size() == 0)
    {
      ERROR_MSG("Input matched protein indices size is 0!, did you specify INPUT_MP?");
      return false;
    }

    if (m_params.exists("INPUT_FASTA"))
    {
      DEBUG_MSG("Loading target proteins from <" << m_params.getValue("INPUT_FASTA") << "> ...");
      if (!m_fasta->Load(m_params.getValue("INPUT_FASTA").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_FASTA"));
        return false;
      }
    }
    if (m_fasta->size() == 0)
    {
      ERROR_MSG("Input fasta proteins size is 0!, did you specify INPUT_FASTA?");
      return false;
    }

    if (m_params.exists("INPUT_ION_TYPES"))
    {
      DEBUG_MSG("Loading ion types model from <" << m_params.getValue("INPUT_ION_TYPES") << "> ...");
      if (!m_model->LoadModel(m_params.getValue("INPUT_ION_TYPES").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_ION_TYPES"));
        return false;
      }
    }
    string specIDFormat = m_params.getValue("SPEC_ID_FORMAT");
    std::transform(specIDFormat.begin(),
                   specIDFormat.end(),
                   specIDFormat.begin(),
                   ::tolower);

    if (m_params.exists("INPUT_SPEC_IDS"))
    {
      bool loadSuccess = true;
      DEBUG_MSG("Loading " << specIDFormat << " results from <" << m_params.getValue("INPUT_SPEC_IDS") << "> ...");
      if (specIDFormat == "inspect")
      {
        if (!m_spectraRawIDs->loadInspectResultsFile(m_params.getValue("INPUT_SPEC_IDS").c_str(),
                                                     true))
        {
          loadSuccess = false;
        }
      }
      else if (specIDFormat == "msgfdb")
      {
        if (!m_spectraRawIDs->loadMSGFDBResultsFile(m_params.getValue("INPUT_SPEC_IDS").c_str(),
                                                    false))
        {
          loadSuccess = false;
        }
      }
      else
      {
        ERROR_MSG("Found unsupported spectrum ID format \'" << specIDFormat << "\'");
        return false;
      }

      if (!loadSuccess)
      {
        ERROR_MSG("Failed to load " << m_params.getValue("INPUT_SPEC_IDS") << " in " << specIDFormat << " format.");
        return false;
      }

      PeptideSpectrumMatchSet nativeIDs;
      nativeIDs = *m_spectraRawIDs;

      vector<string> *inputFileIndex = 0;
      if (m_params.exists("INPUT_FILE_INDEX"))
      {
        // Load the vector of filenames so we can associate file indices with file names
        inputFileIndex = new vector<string> ;

        if (!readFilesFromFile(m_params.getValue("INPUT_FILE_INDEX"),
                               *inputFileIndex))
        {
          ERROR_MSG("Failed to read file name index from \'" << m_params.getValue("INPUT_FILE_INDEX") << "\'");
          delete inputFileIndex;
          return false;
        }

        for (unsigned int i = 0; i < inputFileIndex->size(); i++)
        {
          FilenameManager nameMngr((*inputFileIndex)[i]);
          (*inputFileIndex)[i] = nameMngr.filename;
        }
      }

      ClusterData *clusterInfo = 0;
      if (m_params.exists("INPUT_CLUSTERS_DIR"))
      {
        // Load clusterInfo detailing which MS2 spectra were combined into each cluster
        clusterInfo = new ClusterData;
        DEBUG_MSG("Loading clusters from " << m_params.getValue("INPUT_CLUSTERS_DIR"));
        if (!clusterInfo->loadData(m_params.getValue("INPUT_CLUSTERS_DIR")))
        {
          ERROR_MSG("Failed to load clusters from " << m_params.getValue("INPUT_CLUSTERS_DIR"));
          delete clusterInfo;
          return false;
        }
      }

      vector<vector<int> > *specFileMapping = 0;
      if (m_params.exists("INPUT_FILE_MAPPING") && clusterInfo == 0)
      {
        DEBUG_MSG("Loading file mapping from " << m_params.getValue("INPUT_FILE_MAPPING"));
        // Load mapping of spectrum index to file index and original spectrum index in that file
        specFileMapping = new vector<vector<int> > ;
        if (!Load_binArray(m_params.getValue("INPUT_FILE_MAPPING").c_str(),
                           *specFileMapping))
        {
          ERROR_MSG("Failed to read file mapping from \'" << m_params.getValue("INPUT_FILE_MAPPING") << "\'");
          delete specFileMapping;
          return false;
        }
      }

      vector<vector<int> > *scanRef = 0;
      if (m_params.exists("INPUT_SCAN_REF_FILES"))
      {
        // Load mapping of original spectrum index to scan #
        DEBUG_MSG("Loading bin file list from " << m_params.getValue("INPUT_SCAN_REF_FILES"));
        vector<string> binFiles;
        if (!readFilesFromFile(m_params.getValue("INPUT_SCAN_REF_FILES"),
                               binFiles))
        {
          ERROR_MSG("Failed to read file name index from \'" << m_params.getValue("INPUT_SCAN_REF_FILES") << "\'");
          return false;
        }

        scanRef = new vector<vector<int> > (binFiles.size());
        for (unsigned int fIdx = 0; fIdx < binFiles.size(); fIdx++)
        {
          vector<vector<int> > scans;
          string binFilePath = getPath(spsDir, binFiles[fIdx], false);
          DEBUG_VAR(binFilePath);
          if (!Load_binArray(binFilePath.c_str(), scans))
          {
            ERROR_MSG("Failed to read scan #s from \'" << binFilePath << "\'");
            delete scanRef;
            return false;
          }
          DEBUG_VAR(scans.size());
          (*scanRef)[fIdx].resize(scans.size());

          for (unsigned int i = 0; i < scans.size(); i++)
          {
            (*scanRef)[fIdx][i] = scans[i][0];
          }
        }
      }

      map<int, list<pair<int, string> > > clusterInfoMap;
      if (scanRef != 0 && inputFileIndex != 0)
      {
        // No clustering, use input_mapping.bin to map the scan #s
        clusterInfoMap.clear();

        // scans are always 1-based
        unsigned int newScanIdx = 1;
        for (unsigned int fIdx = 0; fIdx < scanRef->size(); fIdx++)
        {

          for (unsigned int locIdx = 0; locIdx < (*scanRef)[fIdx].size(); locIdx++)
          {
            list<pair<int, string> > scans;
            // Must lookup the scan # from the user input (MS1 spectra are ignored, so they mess up the scan ordering)
            scans.push_back(pair<int, string> ((*scanRef)[fIdx][locIdx],
                                               (*inputFileIndex)[fIdx]));

            clusterInfoMap[newScanIdx++] = scans;
          }
        }
        DEBUG_VAR(m_spectraRawIDs->size());
        DEBUG_MSG("Assigning native scan #s to raw MS/MS IDs ...");
        m_spectraRawIDs->cluster(clusterInfoMap);
        DEBUG_VAR(m_spectraRawIDs->size());
        m_spectraClustIDs->operator =(*m_spectraRawIDs);
        m_rawScanInfo.operator =(clusterInfoMap);
        m_clustScanInfo.operator =(clusterInfoMap);
      }

      if (clusterInfo != 0)
      {
        clusterInfoMap.clear();
        // Use clusterData.bin to map the scan #s
        for (map<unsigned, list<pair<unsigned, unsigned> > >::iterator clustIt =
            clusterInfo->data.begin(); clustIt != clusterInfo->data.end(); clustIt++)
        {
          list<pair<int, string> > scans;

          for (list<pair<unsigned, unsigned> >::iterator scanIt =
              clustIt->second.begin(); scanIt != clustIt->second.end(); scanIt++)
          {
            // Must lookup the scan # from the user input (MS1 spectra are ignored, so they mess up the scan ordering)
            //DEBUG_MSG("File index " << scanIt->first << " spectrum index " << scanIt->second);
            scans.push_back(pair<int, string> ((*scanRef)[scanIt->first][scanIt->second],
                                               (*inputFileIndex)[scanIt->first]));
          }
          // scans are always 1-based, these indices are 0-based
          clusterInfoMap[((int)clustIt->first) + 1] = scans;
        }
        m_spectraClustIDs->operator =(nativeIDs);
        DEBUG_VAR(m_spectraClustIDs->size());
        DEBUG_MSG("Clustering scan raw MS/MS IDs ...");
        m_spectraClustIDs->cluster(clusterInfoMap);
        DEBUG_VAR(m_spectraClustIDs->size());
        m_clustScanInfo.operator =(clusterInfoMap);
      }

      if (m_params.exists("INPUT_CLUST_SPEC_IDS"))
      {
        bool loadSuccess = true;
        DEBUG_MSG("Loading " << specIDFormat << " results from <" << m_params.getValue("INPUT_CLUST_SPEC_IDS") << "> ...");
        if (specIDFormat == "inspect")
        {
          if (!m_spectraClustIDs->loadInspectResultsFile(m_params.getValue("INPUT_CLUST_SPEC_IDS").c_str(),
                                                         true))
          {
            loadSuccess = false;
          }
        }
        else if (specIDFormat == "msgfdb")
        {
          if (!m_spectraClustIDs->loadMSGFDBResultsFile(m_params.getValue("INPUT_CLUST_SPEC_IDS").c_str(),
                                                        false))
          {
            loadSuccess = false;
          }
        }
        else
        {
          ERROR_MSG("Found unsupported spectrum ID format \'" << specIDFormat << "\'");
          return false;
        }

        if (!loadSuccess)
        {
          ERROR_MSG("Failed to load " << m_params.getValue("INPUT_CLUST_SPEC_IDS") << " in " << specIDFormat << " format.");
          return false;
        }
      }

      if (inputFileIndex)
      {
        delete inputFileIndex;
      }
      if (clusterInfo)
      {
        delete clusterInfo;
      }
      if (specFileMapping)
      {
        delete specFileMapping;
      }
      if (scanRef)
      {
        delete scanRef;
      }
    }
    m_targetProts = new vector<string> (m_fasta->size());
    m_proteinSpectra = new SpecSet(m_fasta->size());
    m_peptides = new vector<string> (m_spectraStars->size());

    set<int> target_proteins;
    list<string> strIdxs;

    if (!splitText(m_params.getValue("TARGET_PROTEINS").c_str(), strIdxs, ";"))
    {
      return false;
    }

    for (list<string>::iterator strIt = strIdxs.begin(); strIt != strIdxs.end(); strIt++)
    {
      target_proteins.insert(atoi((*strIt).c_str()));
    }

    FilterFastaProteins(target_proteins, *m_targetProts);

    FilterSpecIds(*m_peptides);

    for (int p = 0; p < m_fasta->size(); p++)
    {
      (*m_proteinSpectra)[p] = m_fasta->getMassesSpec(p);
    }

    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

    m_contigs->setPeakTolerance(peakTol, false);
    m_spectraStars->setPeakTolerance(peakTol, false);

    return true;
  }

  bool ExecReportSPSStats::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  bool ExecReportSPSStats::saveOutputData(void)
  {
    string delim = m_params.getValue("TABLE_DELIM", "\t");
    if (m_params.exists("OUTPUT_SPS_STATS_FILE"))
    {
      MappedSPSStatTable outTab(m_mappedProj);
      outTab.prepareTable();
      DEBUG_MSG("Saving cumulative statistics table to <" << m_params.getValue("OUTPUT_SPS_STATS_FILE") << "> ...");
      if (!outTab.printToCSV(m_params.getValue("OUTPUT_SPS_STATS_FILE").c_str(),
                             delim))
      {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_CONTIG_STATS_FILE"))
    {
      MappedContigSetTable outTab(m_mappedProj);
      outTab.prepareTable();
      DEBUG_MSG("Saving contig statistics table to <" << m_params.getValue("OUTPUT_CONTIG_STATS_FILE") << "> ...");
      if (!outTab.printToCSV(m_params.getValue("OUTPUT_CONTIG_STATS_FILE").c_str(),
                             delim))
      {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_CONTIG_STATS_FILEROOT"))
    {
      string my_fileroot = m_params.getValue("OUTPUT_CONTIG_STATS_FILEROOT");
      string::size_type loc = my_fileroot.find(".tsv");
      if (loc != string::npos)
      {
        my_fileroot.erase(loc);
      }
      string filename;
      string max_idx = parseInt(m_contigs->size() - 1);
      int equalizeLength = max_idx.length();

      string place_holder(equalizeLength, '*');

      DEBUG_MSG("Saving per-contig statistics tables to <" << my_fileroot << place_holder << "> ...");

      MappedContigStatTable outTab(m_mappedProj);
      for (int i = 0; i < m_contigs->size(); i++)
      {
        if ((*m_contigs)[i].size() == 0)
          continue;

        filename = my_fileroot;
        filename.append(parseInt(i, equalizeLength));
        filename.append(".tsv");

        outTab.prepareTable(i);
        if (!outTab.printToCSV(filename.c_str(), delim))
        {
          return false;
        }
      }
    }
    if (m_params.exists("OUTPUT_SEQACC_POS_TABLE"))
    {
      OutputTable seqTab;
      vector<pair<int, int> > values;
      m_mappedProj->getGapAccVsPos(-1, values, &(seqTab.values));
      DEBUG_MSG("Saving Gap Accuracy vs. Position Table to <" << m_params.getValue("OUTPUT_SEQACC_POS_TABLE") << "> ...");
      if (!seqTab.printToCSV(m_params.getValue("OUTPUT_SEQACC_POS_TABLE").c_str(),
                             delim))
      {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PARAMS_FILECOPY"))
    {
      DEBUG_MSG("Saving a copy of the input parameters to <" << m_params.getValue("OUTPUT_PARAMS_FILECOPY") << "> ...");
      if (!m_params.writeToFile(m_params.getValue("OUTPUT_PARAMS_FILECOPY")))
      {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PEPTIDES"))
    {
      DEBUG_MSG("Saving peptides in csv format to <" << m_params.getValue("OUTPUT_PEPTIDES") << "> ...");
      OutputTable csvTable;
      csvTable.values.resize(m_peptides->size());
      for (unsigned int i = 0; i < m_peptides->size(); i++)
      {
        if ((*m_peptides)[i].length() == 0)
        {
          csvTable.values[i].resize(0);
        }
        else
        {
          csvTable.values[i].resize(1);
          csvTable.setValue(i, 0, (*m_peptides)[i]);
        }
      }
      if (!csvTable.printToCSV(m_params.getValue("OUTPUT_PEPTIDES").c_str(),
                               delim))
      {
        return false;
      }
    }
    return true;
  }

  bool ExecReportSPSStats::loadOutputData(void)
  {
    return false;
  }

  vector<ExecBase *>
  const & ExecReportSPSStats::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecReportSPSStats::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecReportSPSStats::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("SPS_PROJECT_DIR");
    VALIDATE_PARAM_EXIST("INPUT_CONTIGS");
    VALIDATE_PARAM_EXIST("INPUT_STARS");
    VALIDATE_PARAM_EXIST("INPUT_CONTIG_ABINFO");
    VALIDATE_PARAM_EXIST("INPUT_MIDX");
    VALIDATE_PARAM_EXIST("INPUT_MP");
    VALIDATE_PARAM_EXIST("INPUT_FASTA");
    VALIDATE_PARAM_EXIST("INPUT_ION_TYPES");
    if (!m_params.exists("INPUT_SPEC_IDS"))
    {
      VALIDATE_PARAM_EXIST("INPUT_CLUST_SPEC_IDS");
    }
    else if (!m_params.exists("INPUT_CLUST_SPEC_IDS"))
    {
      VALIDATE_PARAM_EXIST("INPUT_SPEC_IDS");
    }
    VALIDATE_PARAM_EXIST("SPEC_ID_FORMAT");
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TARGET_PROTEINS");
    VALIDATE_PARAM_EXIST("INPUT_SCAN_REF_FILES");
    VALIDATE_PARAM_EXIST("INPUT_FILE_INDEX");
    VALIDATE_PARAM_EXIST("INPUT_REF_INDICES");

    m_isValid = true;
    return true;
  }

  void ExecReportSPSStats::FilterFastaProteins(set<int>& target_proteins,
                                               vector<string>& put_proteins)
  {
    put_proteins.resize(m_fasta->size());
    for (int i = 0; i < put_proteins.size(); i++)
    {
      if (target_proteins.count(i) > 0 || target_proteins.size() == 0)
      {
        put_proteins[i] = m_fasta->getSequence(i);
        //cout << "protein " << i << " has size " << put_proteins[i].length() << "\n";
      }
      else
      {
        put_proteins[i] = "";
      }
    }
  }

  void ExecReportSPSStats::FilterSpecIds(vector<string>& put_peptides)
  {
    put_peptides.resize(m_spectraStars->size(), "");

    int lastSpecIdx = 0;

    for (int i = 0; i < m_spectraClustIDs->size(); i++)
    {
      if ((*m_spectraClustIDs)[i]->m_annotation.length() == 0)
      {
        continue;
      }
      else
      {
        //(*m_spectraClustIDs)[i]->m_scanNum += 1;
        int specIdx = -1; //(*m_spectraClustIDs)[i]->m_scanNum;
        int j = lastSpecIdx;
        bool init = true;
        while (init || j != lastSpecIdx)
        {
          init = false;

          if ((*m_spectraStars)[j].scan == (*m_spectraClustIDs)[i]->m_scanNum)
          {
            specIdx = j;
            lastSpecIdx = j;
            break;
          }
          j = (j == m_spectraStars->size() - 1) ? 0 : j + 1;
        }
        if (specIdx < 0)
        {
          ERROR_MSG("Failed to locate scan " << (*m_spectraClustIDs)[i]->m_scanNum);
          abort();
        }
        put_peptides[specIdx]
            = makeBracketsMods((*m_spectraClustIDs)[i]->m_annotation);
      }
    }
  }

  void ExecReportSPSStats::ChopContigEnds(int endsChop)
  {
    vector<pair<vector<int> , vector<double> > > *abruijn_verts;
    vector<int>* spectrumIdxs;
    vector<int>* spectrumFlip;

    vector<int> specIdxOut;
    vector<int> specFlipOut;
    vector<pair<vector<int> , vector<double> > > abruijnVertsOut;
    set<int> specIdxUse;
    map<int, int> specFlipUse;

    Spectrum newContig;
    Spectrum newMidx;

    for (int i = 0; i < m_contigs->size(); i++)
    {
      if ((*m_contigs)[i].size() == 0)
      { // ignore empty contigs
        continue;
      }
      if ((*m_contigs)[i].size() <= 2 * endsChop)
      {
        (*m_contigs)[i].resize(0);
        continue;
      }
      abruijn_verts = &(*m_abinfo)[i].second;

      specIdxUse.clear();
      specFlipUse.clear();
      spectrumIdxs = &(*m_abinfo)[i].first.first;
      spectrumFlip = &(*m_abinfo)[i].first.second;

      for (int k = 0; k < spectrumIdxs->size(); k++)
      {
        specFlipUse[(*spectrumIdxs)[k]] = (*spectrumFlip)[k];
      }

      bool foundFirst = false;
      float firstMass = 0;

      bool foundFirstEnd = false;
      float firstEnd = 0;
      float lastEnd = 0;

      // create new contig, abinfo
      newContig = (*m_contigs)[i];
      newMidx = (*m_contigOverlaps)[i];

      newContig.resize((*m_contigs)[i].size() - 2 * endsChop);
      newMidx.resize(newContig.size());
      abruijnVertsOut.resize(newContig.size());
      for (int j = 0; j < (*m_contigs)[i].size(); j++)
      {
        spectrumIdxs = &((*abruijn_verts)[j].first);
        if (j < endsChop)
        { // ignore beginning of contig
          continue;
        }
        else if ((*m_contigs)[i].size() - j <= endsChop)
        { // ignore end of contig
          // remember the cummulative mass removed from the end in order to update parent mass
          if (!foundFirstEnd)
          {
            firstEnd = (*m_contigs)[i][j][0];
            foundFirstEnd = true;
          }
          lastEnd = (*m_contigs)[i][j][0];
        }
        else
        {
          // remember the spectra idxs assembled in remaining abruijn vertices
          for (int k = 0; k < spectrumIdxs->size(); k++)
          {
            specIdxUse.insert((*spectrumIdxs)[k]);
          }
          // remember the cummulative mass removed from the beginning to update heavier peaks and parent mass
          if (!foundFirst)
          {
            firstMass = (*m_contigs)[i][j][0];
            foundFirst = true;
          }

          newContig[j - endsChop] = (*m_contigs)[i][j];
          newContig[j - endsChop][0] -= firstMass;

          // copy assembled peak masses
          abruijnVertsOut[j - endsChop] = (*abruijn_verts)[j];
        }
      }
      // subtract removed mass from parent mass
      newContig.parentMass -= (lastEnd - firstEnd + firstMass);

      specIdxOut.resize(specIdxUse.size());
      specFlipOut.resize(specIdxUse.size());
      int idxUse = 0;
      for (set<int>::iterator specIt = specIdxUse.begin(); specIt
          != specIdxUse.end(); specIt++)
      {
        specIdxOut[idxUse] = *specIt;
        specFlipOut[idxUse] = specFlipUse[*specIt];
        idxUse++;
      }

      idxUse = 0;
      for (int k = 0; k < (*m_contigOverlaps)[i].size(); k++)
      {
        int contigPeakIdx = floatToInt((*m_contigOverlaps)[i][k][0]);
        if (contigPeakIdx >= endsChop && contigPeakIdx < (*m_contigs)[i].size()
            - endsChop)
        {
          newMidx[idxUse] = (*m_contigOverlaps)[i][k];
          newMidx[idxUse][0] -= (float)endsChop;
          idxUse++;
        }
      }
      newMidx.resize(idxUse);

      // update contig SpecSet, abinfo, and overlapped indicies
      (*m_contigs)[i] = newContig;
      (*m_abinfo)[i].second = abruijnVertsOut;
      (*m_abinfo)[i].first.first = specIdxOut;
      (*m_abinfo)[i].first.second = specFlipOut;
      (*m_contigOverlaps)[i] = newMidx;
    }
  }
}
