////////////////////////////////////////////////////////////////////////////////
#include "spsFiles.h"
#include "main_specnets_defs.h"
////////////////////////////////////////////////////////////////////////////////
namespace spsReports
{
  ////////////////////////////////////////////////////////////////////////////////
  SpsFiles::SpsFiles() :
    m_consensusSpecSet(NULL), m_clusterData(NULL), m_fasta(NULL),
        m_abruijn(NULL), m_starSpectra(NULL), m_contigIndices(NULL),
        m_contigNames(NULL), m_input_index(NULL), m_contigsSpectra(NULL),
        m_homglue_ref_midx(NULL), m_homglue_ref_mp(NULL), m_sps_seqs(NULL),
        m_contigs_midx(NULL), m_contigs_mp(NULL), m_homglueMatches(NULL),
        m_homglue_matches_midx(NULL), m_homglue_matches_mp(NULL) //, m_inputMapping(NULL)
  {
  }
  ////////////////////////////////////////////////////////////////////////////////
  SpsFiles::~SpsFiles()
  {
    if (m_consensusSpecSet != NULL)
      delete m_consensusSpecSet;
    if (m_clusterData != NULL)
      delete m_clusterData;
    if (m_fasta != NULL)
      delete m_fasta;
    if (m_abruijn != NULL)
      delete m_abruijn;
    if (m_starSpectra != NULL)
      delete m_starSpectra;
    if (m_contigIndices != NULL)
      delete m_contigIndices;
    //if (m_inputMapping != NULL)
    //  delete m_inputMapping;
    if (m_contigNames != NULL)
      delete m_contigNames;
    if (m_input_index != NULL)
      delete m_input_index;
    if (m_contigsSpectra != NULL)
      delete m_contigsSpectra;
    if (m_homglue_ref_midx != NULL)
      delete m_homglue_ref_midx;
    if (m_homglue_ref_mp != NULL)
      delete m_homglue_ref_mp;
    if (m_sps_seqs != NULL)
      delete m_sps_seqs;
    if (m_contigs_midx != NULL)
      delete m_contigs_midx;
    if (m_contigs_mp != NULL)
      delete m_contigs_mp;
    if (m_homglueMatches != NULL)
      delete m_homglueMatches;
    if (m_homglue_matches_midx != NULL)
      delete m_homglue_matches_midx;
    if (m_homglue_matches_mp != NULL)
      delete m_homglue_matches_mp;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // File loading section.
  ////////////////////////////////////////////////////////////////////////////////
  void *SpsFiles::getData(SpsFileID what, int index)
  {
    switch (what)
    {
    case SPS_FILE_FASTA:
      return m_fasta;
    case SPS_FILE_SPECSET:
      return (index < m_specSet.size() ? &(m_specSet[index]) : NULL);
    case SPS_FILE_MSCLUSTER:
      return m_clusterData;
    case SPS_FILE_STARS:
      return m_starSpectra;
    case SPS_FILE_SEQS:
      return m_sps_seqs;
    case SPS_FILE_ABRUIJN:
      return m_abruijn;
    case SPS_FILE_CONSENSUS_SPECTRA:
      return m_consensusSpecSet;
    default:
      return NULL;
    }

    return NULL;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // File loading section.
  ////////////////////////////////////////////////////////////////////////////////
  // Main load method. Checks for a table to load. If non existent, calls the loadData method, and then saves the table
  int SpsFiles::loadData(ReportGeneratorData &reportGeneratorData)
  {
    int ret = OK;

    if (loadInputSpectraList(reportGeneratorData.projectDir,
                             reportGeneratorData.filenameInputSpectraList)
        == ERROR)
    {
      WARN_MSG("Error loading " << reportGeneratorData.filenameInputSpectraList);
      ret = ERROR;
    }

    if (loadInputSpectraFiles(reportGeneratorData.projectDir,
                              reportGeneratorData.filenameClusterMS) == ERROR)
    {
      WARN_MSG("Error loading " << reportGeneratorData.filenameClusterMS);
      ret = ERROR;
    }

    if (loadScanFiles(reportGeneratorData.projectDir,
                      reportGeneratorData.filenameScanFiles) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameScanFiles);

    if (loadConsensusSpectraFile(reportGeneratorData.projectDir,
                                 reportGeneratorData.filenameConsensusSpectra)
        == ERROR)
    {
      WARN_MSG("Error loading " << reportGeneratorData.filenameConsensusSpectra);
      ret = ERROR;
    }

    string scoredClustFile = SPECTRA_DIR + "/" + SPECS_SCORED_CLUST_FILE;
    string msClustFile = SPECTRA_DIR + "/" + SPECS_MS_CLUST_FILE;

    if (loadMsClusterData(reportGeneratorData.projectDir, scoredClustFile) == ERROR) {
      //                    reportGeneratorData.filenameCluster) == ERROR) {
      if (loadMsClusterData(reportGeneratorData.projectDir, msClustFile) == ERROR) {
        reportGeneratorData.noClusters = true;
        WARN_MSG("Error loading cluster data");
        WARN_MSG("Assuming there is no clusters layer.");
      }
    }

    if (loadProteinsFile(reportGeneratorData.projectDir,
                         reportGeneratorData.filenameProteins) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameProteins);

    if (loadContigIndices(reportGeneratorData.projectDir,
                          reportGeneratorData.filenameContigIndices) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameContigIndices);
    if (loadHomglueRefMidx(reportGeneratorData.projectDir,
                           reportGeneratorData.filenameHomglueRefMidx) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueRefMidx);
    if (loadHomglueRefMp(reportGeneratorData.projectDir,
                         reportGeneratorData.filenameHomglueRefMp) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueRefMp);

    if (loadContigsMidxAll(reportGeneratorData.projectDir,
                           reportGeneratorData.filenameContigsMidx) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameContigsMidx);
    if (loadContigsMpAll(reportGeneratorData.projectDir,
                         reportGeneratorData.filenameContigsMp) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameContigsMp);
    if (loadSpsSeqs(reportGeneratorData.projectDir,
                    reportGeneratorData.filenameSpsSeqs) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameSpsSeqs);

    if (loadStarSpectra(reportGeneratorData.projectDir,
                        reportGeneratorData.filenameStarSpectra) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameStarSpectra);
    if (loadAbruijn(reportGeneratorData.projectDir,
                    reportGeneratorData.filenameAbruijn) == ERROR)
    {
      WARN_MSG("Error loading " << reportGeneratorData.filenameAbruijn);
      ret = ERROR;
    }

    if (loadContigSpectra(reportGeneratorData.projectDir,
                          reportGeneratorData.filenameContigSpectra) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameContigSpectra);

    if (loadCspsMatchesMp(reportGeneratorData.projectDir,
                          reportGeneratorData.filenameHomglueMatchesMp)
        == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatchesMp);
    if (loadCspsMatchesMidx(reportGeneratorData.projectDir,
                            reportGeneratorData.filenameHomglueMatchesMidx)
        == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatchesMidx);
    if (loadCspsSpectra(reportGeneratorData.projectDir,
                        reportGeneratorData.filenameHomglueMatches) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatches);

    if (loadContigNames(reportGeneratorData.projectDir,
                        reportGeneratorData.filenameContigNames) == ERROR)
      WARN_MSG("Error loading " << reportGeneratorData.filenameContigNames);

    // load Input Mapping silently
    //loadInputMapping(reportGeneratorData.projectDir,
    //                 reportGeneratorData.filenameInputMapping);

    return ret;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // File name composition
  string SpsFiles::composeFileName(const string &projectDir,
                                   const string &fileName)
  {
    if (fileName[0] == '/')
      return fileName;

    // Compose output path
    string aux = projectDir;
    if (aux[aux.length() - 1] != '/')
      aux += '/';
    // add filename
    aux += fileName;
    // return composed filename
    return aux;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // File name composition
  int SpsFiles::readPklbin(specnets::SpecSet *&pklbin, const string &fileName)
  {
    // Test if file already loaded
    if (pklbin)
      return OK;
    DEBUG_MSG("Reading pklbin file: " << fileName);
    // if not, create the object
    pklbin = new SpecSet();
    // from <sps>/spectra/specs_ms.pklbin
    if (pklbin->loadPklBin(fileName.c_str()) == -1)
    {
      delete pklbin;
      pklbin = NULL;
      ERROR_MSG("Couldn't open " << fileName);
      return ERROR;
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Read binary array
  int SpsFiles::readBinArray(vector<vector<int> > *&binArray,
                             const string &fileName)
  {
    // Test if file already loaded
    if (binArray)
      return OK;
    DEBUG_MSG("Reading binArray file: " << fileName);
    // if not, create the object
    binArray = new vector<vector<int> > ();
    // load the data
    if (Load_binArray(fileName.c_str(), *binArray) < 0)
    {
      delete binArray;
      binArray = NULL;
      ERROR_MSG("Couldn't open " << fileName);
      return ERROR;
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // read original input spectra filename list
  int SpsFiles::loadInputSpectraList(const string &projectDir,
                                     const string &fileName)
  {
    if (m_input_index)
      return OK;
    // create object
    m_input_index = new vector<string> ();
    // get filename with full path
    string fn = fileName;
    // log file to load
    DEBUG_MSG("Loading input spectra list from file: " << fn);
    //read the file
    //if(!absoluteFilenameInputSpectraList)
    fn = composeFileName(projectDir, fn);
    // read the index file
    if (!readFilesFromFile(fn, *m_input_index))
      return ERROR;
    // remove path from filenames
    for (int i = 0; i < m_input_index->size(); i++)
      (*m_input_index)[i]
          = (*m_input_index)[i].substr((*m_input_index)[i].find_last_of("/\\")
              + 1);
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // read pklbin spectra files
  int SpsFiles::loadInputSpectraFiles(const string &projectDir,
                                      const string &fileName)
  {
    // get filename with full path
    string fn = fileName;
    //if(!absoluteFilenameClusterMS)
    fn = composeFileName(projectDir, fn);
    // log file to load
    DEBUG_MSG("Loading input spectra files names from file: " << fn);
    //read the file
    // read the index file
    if (!readFilesFromFile(fn, m_inputSpectraPklbin))
      return ERROR;
    // read each specset file
    for (int i = 0; i < m_inputSpectraPklbin.size(); i++)
    {
      //fn = composeFileName(projectDir, m_inputSpectraPklbin[i].c_str());
      fn = m_inputSpectraPklbin[i];
      SpecSet currentSpecSet;
      if (currentSpecSet.loadPklBin(fn.c_str()) == -1)
      {
        ERROR_MSG("Couldn't open " << fn);
        return ERROR;
      }
      // store read specset
      m_specSet.push_back(currentSpecSet);
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Load scan numbers files
  int SpsFiles::loadScanFiles(const string &projectDir, const string &fileName)
  {
    // get filename with full path
    string fn = fileName;
    //if(!absoluteFilenameScanFiles)
    fn = composeFileName(projectDir, fn);
    // read the index file
    readFilesFromFile(fn, m_scanNumberFiles);
    // log file to load
    DEBUG_MSG("Loading scan files names from file: " << fn);
    //read the file
    // read each specset file
    for (int i = 0; i < m_scanNumberFiles.size(); i++)
    {
      //fn = composeFileName(projectDir, m_scanNumberFiles[i].c_str());
      fn = m_scanNumberFiles[i];
      // data holder
      vector<vector<int> > data;
      // load the data
      if (Load_binArray(fn.c_str(), data) < 0)
      {
        ERROR_MSG("Couldn't open " << fn);
        return ERROR;
      }
      // store read scan data
      m_specScan.push_back(data);
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Load cluster consensus spectra pklbin file
  int SpsFiles::loadConsensusSpectraFile(const string &projectDir,
                                         const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading consensus spectra from file: " << fn);
    //read the file
    return readPklbin(m_consensusSpecSet, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // load MsCluster data.
  int SpsFiles::loadMsClusterData(const string &projectDir,
                                  const string &fileName)
  {
    // Test if file already loaded
    if (m_clusterData)
      return OK;
    // if not, create the object
    // log file to load
    DEBUG_MSG("Loading clusterData: " << projectDir);
    // first try to load data in binary format

    m_clusterData = new ClusterSet();
    string fileNameComplete(projectDir);
    if(fileNameComplete[fileNameComplete.size()-1] != '/')
      fileNameComplete += '/';
    fileNameComplete += fileName;
    //fileNameComplete += "/spectra/specs_ms.clust";
    int ret = m_clusterData->loadBinaryFile(fileNameComplete);

    //m_clusterData->dump(cerr, false);

    // check for errors
    if (ret == false)
    {
      DEBUG_MSG("Problem opening clusterData file!");
      delete m_clusterData;
      m_clusterData = NULL;
      return ERROR;
    }
    // return
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (4) Load Star Spectra
  int SpsFiles::loadStarSpectra(const string &projectDir,
                                const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading star spectra from file: " << fn);
    //read the file
    return readPklbin(m_starSpectra, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (8) Load Abruijn graph
  int SpsFiles::loadAbruijn(const string &projectDir, const string &fileName)
  {
    // Test if file already loaded
    if (m_abruijn)
      return OK;
    // if not, create the object
    m_abruijn = new abinfo_t();
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Reading abruijn graph from file: " << fn);
    // load the data
    if (Load_abinfo(fn.c_str(), *m_abruijn) == 0)
    {
      delete m_abruijn;
      m_abruijn = NULL;
      ERROR_MSG("Couldn't open " << fn);
      return ERROR;
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Load proteins file (FASTA format).
  int SpsFiles::loadProteinsFile(const string &projectDir,
                                 const string &fileName)
  {
    // Test if file already loaded
    if (m_fasta)
      return OK;
    // if not, create the object
    m_fasta = new DB_fasta();
    // get filename with full path
    string fn = fileName;
    // log file to load
    DEBUG_MSG("Reading proteins from file: " << fn);
    // load the data
    if (m_fasta->Load(fn.c_str()) == 0)
    {
      delete m_fasta;
      m_fasta = NULL;
      ERROR_MSG("Couldn't open " << fn);
      return ERROR;
    }
    return OK;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (5) Load Spectra for contigs that matched a protein
  int SpsFiles::loadContigSpectra(const string &projectDir,
                                  const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading contig spectra from file: " << fn);
    //read the file
    return readPklbin(m_contigsSpectra, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  int SpsFiles::loadContigIndices(const string &projectDir,
                                  const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading contig indices from file: " << fn);
    //read the file
    return readBinArray(m_contigIndices, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  /*int SpsFiles::loadInputMapping(const string &projectDir,
                                 const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading input mapping from file: " << fn);
    //read the file
    return readBinArray(m_inputMapping, fn);
  }*/
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  // from <sps>/homology/homglue_ref_midx.pklbin
  int SpsFiles::loadHomglueRefMidx(const string &projectDir,
                                   const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading homglue ref midx from file: " << fn);
    //read the file
    return readPklbin(m_homglue_ref_midx, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  // from <sps>/homology/homglue_ref_mp.bin
  int SpsFiles::loadHomglueRefMp(const string &projectDir,
                                 const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading homglue ref mp from file: " << fn);
    //read the file
    return readBinArray(m_homglue_ref_mp, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  int SpsFiles::loadContigsMidxAll(const string &projectDir,
                                   const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading contigs midx all from file: " << fn);
    //read the file
    return readPklbin(m_contigs_midx, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  // from <sps>/homology/contigs_mp_all.bin
  int SpsFiles::loadContigsMpAll(const string &projectDir,
                                 const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading contigs mp all from file: " << fn);
    //read the file
    return readBinArray(m_contigs_mp, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // where to load contig indices to
  // from <sps>/assembly/sps_seqs.pklbin
  int SpsFiles::loadSpsSeqs(const string &projectDir, const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading sps_seqs from file: " << fn);
    //read the file
    return readPklbin(m_sps_seqs, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (13) load <sps>/homology/homglue_matches_mp.pklbin
  int SpsFiles::loadCspsMatchesMp(const string &projectDir,
                                  const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading csps matches mp from file: " << fn);
    //read the file
    return readBinArray(m_homglue_matches_mp, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (14) load <sps>/homology/homglue_matches_midx.bin
  int SpsFiles::loadCspsMatchesMidx(const string &projectDir,
                                    const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading csps matches midx from file: " << fn);
    //read the file
    return readPklbin(m_homglue_matches_midx, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (15) load <sps>/homology/homglue_matches.pklbin
  int SpsFiles::loadCspsSpectra(const string &projectDir,
                                const string &fileName)
  {
    // get filename with full path
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading csps spectra from file: " << fn);
    // read the file
    return readPklbin(m_homglueMatches, fn);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // (--) load <sps>/homology/ref_sps_names.txt
  int SpsFiles::loadContigNames(const string &projectDir,
                                const string &fileName)
  {
    // check for object existence
    if (m_contigNames)
      return OK;
    // create object
    m_contigNames = new vector<string> ();
    // get filename with full path
    //string fn = filenameContigNames;
    string fn = composeFileName(projectDir, fileName);
    // log file to load
    DEBUG_MSG("Loading contig names from file: " << fn);
    // read the index file
    readFilesFromFile(fn, *m_contigNames);
    // return ok
    return 1;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Dump the abruijn graph to the screen
  void SpsFiles::dump_abruijn(ostream &sout, abinfo_t *abruijn, bool web)
  {
    if (!abruijn)
      return;

    sout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!";
    if (web)
      sout << "<br>";
    sout << endl;
    sout << "----------------- abruijn graph----------------";
    if (web)
      sout << "<br>";
    sout << endl;

    // DUMP abruijn
    abinfo_t::iterator i0 = abruijn->begin();
    for (; i0 != abruijn->end(); i0++)
    {
      sout << "Contig " << i0->first << endl;

      // output first vector in pair
      for (int i = 0; i < i0->second.first.first.size(); i++)
      {
        sout << "     spectrum: " << i0->second.first.first[i]
            << "  --  flipped: " << i0->second.first.second[i];
        if (web)
          sout << "<br>";
        sout << endl;
      }

      sout << "     ------ ";
      if (web)
        sout << "<br>";
      sout << endl;

      // output second pair -> vector of pairs of vectors
      // output first vector in pair
      std::vector<std::pair<vector<int> , vector<double> > >::iterator i1 =
          i0->second.second.begin();
      for (; i1 != i0->second.second.end(); i1++)
      {
        sout << "     ------ ";
        if (web)
          sout << "<br>";
        sout << endl;

        for (int i = 0; i < i1->first.size(); i++)
        {
          sout << "     spectrum: " << i1->first[i] << "  --  peak mass: "
              << i1->second[i];
          if (web)
            sout << "<br>";
          sout << endl;
        }
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Dump cluster data
#ifdef __USE_CLUSTER_SET__
  void SpsFiles::dump_clusterData(ostream &sout,
                                  ClusterSet *cd,
                                  char *title,
                                  bool web)
#else
  void SpsFiles::dump_clusterData(ostream &sout, ClusterData *cd, char *title, bool web)
#endif
  {
    if (!cd)
      return;

    sout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!";
    if (web)
      sout << "<br>";
    sout << endl;
    sout << "------------''--- cluster data for " << title
        << " ----------------";
    if (web)
      sout << "<br>";
    sout << endl;

    cd->dump(sout, web);
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Dump the bin_array to the screen
  void SpsFiles::dump_binArray(ostream &sout,
                               vector<vector<int> > *a,
                               char *title,
                               bool web)
  {
    if (!a)
      return;

    sout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!";
    if (web)
      sout << "<br>";
    sout << endl;
    sout << "----------------- bin array graph for " << title
        << " ----------------";
    if (web)
      sout << "<br>";
    sout << endl;

    for (int i = 0; i < a->size(); i++)
    {
      sout << i << ":  ";
      for (int j = 0; j < (*a)[i].size(); j++)
      {
        if (j)
          sout << ", ";
        sout << (*a)[i][j];
      }
      if (web)
        sout << "<br>";
      sout << endl;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // Dump the specset to the screen
  void SpsFiles::dump_specset(ostream &sout,
                              specnets::SpecSet *set,
                              char *title,
                              bool web)
  {
    if (!set)
      return;

    sout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!";
    if (web)
      sout << "<br>";
    sout << endl;
    sout << "----------------- specset for " << title << " ----------------";
    if (web)
      sout << "<br>";
    sout << endl;

    for (int i = 0; i < set->size(); i++)
    {
      sout << i << ": size: " << (*set)[i].size() << " : ";
      for (int j = 0; j < (*set)[i].size(); j++)
      {
        sout << " (";
        sout << (*set)[i][j][0] << " ; " << (*set)[i][j][1];
        sout << ")";
      }
      if (web)
        sout << "<br>";
      sout << endl;
    }
  }
////////////////////////////////////////////////////////////////////////////////
}
; // namespace
////////////////////////////////////////////////////////////////////////////////
