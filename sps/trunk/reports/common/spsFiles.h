////////////////////////////////////////////////////////////////////////////////
#ifndef __SPS_FILES_H__
#define __SPS_FILES_H__
////////////////////////////////////////////////////////////////////////////////
// includes section
#include <string>
#include <vector>

#include "spectrum.h"
#include "ClusterData.h"
#include "ClusterSet.h"
#include "db_fasta.h"
#include "abruijn.h"
#include "aminoacid.h"

#include "Defines.h"
#include "ReportDefines.h"

#include "ReportData.h"
////////////////////////////////////////////////////////////////////////////////
// SPS Files
//
// This class contains all the files needed by sps
//
////////////////////////////////////////////////////////////////////////////////
// namespaces section
using namespace std;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Definitions section
typedef enum {
  SPS_FILE_FASTA,
  SPS_FILE_SPECSET,
  SPS_FILE_BIN_ARRAY,
  SPS_FILE_MSCLUSTER,
  SPS_FILE_STARS,
  SPS_FILE_SEQS,
  SPS_FILE_ABRUIJN,
  SPS_FILE_CONSENSUS_SPECTRA
} SpsFileID;
////////////////////////////////////////////////////////////////////////////////

 /*! \brief SPS files container class

   Container for all SPS files used in reports.

   */
class SpsFiles {

 public:

    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Default contructor
     */
  SpsFiles();
    //@}

    //! \name DESTRUCTOR
    //@{
  ~SpsFiles();
    //@}


  /*! \brief
   */
  void *getData(SpsFileID what, int index);

  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data needed for several table fields.

  /*! \brief input spectra filenames (pklbin)
   */
  vector<string> m_inputSpectraPklbin;

  /*! \brief Where to load input spectra data into (pklbin)
  // from <sps>/spectra/specs_ms_*.pklbin
   */
  vector<specnets::SpecSet> m_specSet;

  /*! \brief scan number files (text list)
   */
  vector<string> m_scanNumberFiles;

  /*! \brief Where to load input spectra scan numbers
  // from <sps>/spectra/specs_ms_*.pklbin
   */
  vector<vector<vector<int> > > m_specScan;

  /*! \brief
  // Where to load consensus spectra data into (pklbin)
  // from <sps>/spectra/specs_ms.pklbin
   */
  specnets::SpecSet *m_consensusSpecSet;
  ClusterSet *m_clusterData;

  /*! \brief  where to store clusters information
  // from <sps>/clusterData.bin
   */
  //ClusterData *m_clusterData;

  /*! \brief where to store the protein sequences
   */
  specnets::DB_fasta  *m_fasta;

  /*! \brief where to store the abruijn graph
   */
  specnets::abinfo_t  *m_abruijn;

  /*! \brief where to store the star spectra
   */
  specnets::SpecSet *m_starSpectra;

  /*! \brief where to load contig indices to
  // from <sps>/spectra/contig_indices.bin
   */
  vector<vector<int> >  *m_contigIndices;

  /*! \brief where to load contig names into
   */
  vector<string>        *m_contigNames;

  /*! \brief where to load input spectra file names into
   */
  vector<string>        *m_input_index;

  /*! \brief where to load input spectra mapping information to
  // from <sps>/spectra/input_mapping.bin
   */
  //vector<vector<int> >  *m_inputMapping;


  //////////////////////////////////////////////////////////////////////////////
  // Data needed to generate Reference sequences

  /*! \brief where to store the spectra for contigs
   */
  specnets::SpecSet *m_contigsSpectra;

  /*! \brief where to load contig indices to
  // from <sps>/homology/homglue_ref_midx.pklbin
   */
  specnets::SpecSet *m_homglue_ref_midx;

  /*! \brief where to load contig indices to
  // from <sps>/homology/homglue_ref_mp.bin
   */
  vector<vector<int> > *m_homglue_ref_mp;

  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate Homolog sequences

  /*! \brief where to load contig indices to
  // from <sps>/assembly/sps_seqs.pklbin
   */
  specnets::SpecSet *m_sps_seqs;

  /*! \brief where to load contig indices to
  // from <sps>/homology/contigs_midx_all.pklbin
   */
  specnets::SpecSet *m_contigs_midx;

  /*! \brief where to load contig indices to
  // from <sps>/homology/contigs_mp_all.bin
   */
  vector<vector<int> > *m_contigs_mp;

  //////////////////////////////////////////////////////////////////////////////
  // data related CSPS contigs

  /*! \brief where to store the CSPS spectra
   */
  specnets::SpecSet *m_homglueMatches;

  /*! \brief where to load CSPS contig indices to
  // from <sps>/homology/homglue_matches_midx.pklbin
   */
  specnets::SpecSet *m_homglue_matches_midx;

  /*! \brief where to load CSPS contig indices to
  // from <sps>/homology/homglue_matches_mp.bin
   */
  vector<vector<int> > *m_homglue_matches_mp;


  //////////////////////////////////////////////////////////////////////////////
  // Methods to load data files


  // General load data method
   /*! \brief Top level data loading method.

     Loads all files.
     */
  //virtual int loadData(const ReportGeneratorData &reportGeneratorData);
  virtual int loadData(ReportGeneratorData &reportGeneratorData);


  // pklbin file load
   /*! \brief Load pklbin file

     Loads a generic pklbin

    @param pklbin the pklbin object to load data into
    @param fileName Filename
     */
  virtual int readPklbin(specnets::SpecSet *&pklbin, const string &fileName);

  // binArray file load
   /*! \brief Load binArray file

     Loads a generic binArray

    @param binArray the binArray object to load data into
    @param fileName Filename
     */
  virtual int readBinArray(vector<vector<int> >  *&binArray, const string &fileName);

  // filename composition helper method
   /*! \brief Composes a file name,

     Composes a full filename to open the file, given a project directory and a filename

    @param projectDir project directory
    @param fileName Filename
     */
  virtual string composeFileName(const string &projectDir, const string &fileName);

  // load the input spectra files list file
   /*! \brief Load spectra list file

     Loads the file contiaining the list of spectra files

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadInputSpectraList(const string &projectDir, const string &fileName);

  // (1) Load input spectra
   /*! \brief Load input spectra file

     Loads the file contiaining the input spectra

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadInputSpectraFiles(const string &projectDir, const string &fileName);

  // load scan numbers files
   /*! \brief Load scan numbers file

     Loads the file contiaining the separate scan numbers for input spectra files

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadScanFiles(const string &projectDir, const string &fileName);

  // (3) load MsCluster data
   /*! \brief Load MsCluster file

     Loads the file contiaining the clusters information

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadMsClusterData(const string &projectDir, const string &fileName);

  // (2) Load cluster consensus spectra file
   /*! \brief Load clusters spectra file

     Loads the file contiaining the cluster consensus spectra

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadConsensusSpectraFile(const string &projectDir, const string &fileName);

  // (6) Load contig indices
   /*! \brief Load Contig Indices file

     Loads the file contiaining the contig indices

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadContigIndices(const string &projectDir, const string &fileName);

  // (4) Load Star Spectra
   /*! \brief Load star spectra file

     Loads the file contiaining the stars spectra

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadStarSpectra(const string &projectDir, const string &fileName);

  // (8) Load Abruijn graph
   /*! \brief Load abinfo file

     Loads the file contiaining the Abruijn Graph

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadAbruijn(const string &projectDir, const string &fileName);

  // (5) Load Spectra for contigs that matched a protein
   /*! \brief Load contig spectra file

     Loads the file contiaining the contigs that match a protein

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadContigSpectra(const string &projectDir, const string &fileName);

  // (11) load <sps>/homology/homglue_ref_mp.bin
   /*! \brief Load homglue ref mp file

     Loads the file contiaining the reference homology map

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadHomglueRefMp(const string &projectDir, const string &fileName);

  // (12) load <sps>/homology/homglue_ref_midx.pklbin
   /*! \brief Load homglue ref midx file

     Loads the file contiaining the reference homology index

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadHomglueRefMidx(const string &projectDir, const string &fileName);

  // (7) load <sps>/assembly/sps_seqs.pklbin
  virtual int loadSpsSeqs(const string &projectDir, const string &fileName);

  // (9) load <sps>/homology/contigs_mp_all.bin
   /*! \brief Load contigs mp all file

     Loads the file contiaining the All Contigs map

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadContigsMpAll(const string &projectDir, const string &fileName);

  // (10) load <sps>/homology/contigs_midx_all.pklbin
   /*! \brief Load contigs midx all file

     Loads the file contiaining the All Contigs indexes

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadContigsMidxAll(const string &projectDir, const string &fileName);

  // (13) load <sps>/homology/homglue_matches_mp.pklbin
   /*! \brief Load homglue matches ms file

     Loads the file contiaining the homology matches

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadCspsMatchesMp(const string &projectDir, const string &fileName);

  // (14) load <sps>/homology/homglue_matches_midx.bin
   /*! \brief Load homglue midx file

     Loads the file contiaining the homology indexes

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadCspsMatchesMidx(const string &projectDir, const string &fileName);

  // (15) load <sps>/homology/homglue_matches.pklbin
   /*! \brief Load homglue matches file

     Loads the file contiaining the homology matches

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadCspsSpectra(const string &projectDir, const string &fileName);

  // (16) Load proteins file (FASTA format)
   /*! \brief Load Proteins file

     Loads the file the proteins

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadProteinsFile(const string &projectDir, const string &fileName);

  // (17) Load contig names
   /*! \brief Load ContigNames file

     Loads the file containing the contigs names

    @param projectDir project directory
    @param fileName Filename
     */
  virtual int loadContigNames(const string &projectDir, const string &fileName);

   /*! \brief Load InputMapping file

     Loads the input mapping file

    @param projectDir project directory
    @param fileName Filename
     */
  //virtual int loadInputMapping(const string &projectDir, const string &fileName);




  // debug - dump abruijn graph
  void dump_abruijn(ostream &sout, bool web = false) {dump_abruijn(sout, m_abruijn, web);};

   /*! \brief Dumps an abruijn file.

     Outputs an abruijn file contents.

    @param sout Output stream
    @param abruijn The abruijn graph data structure to output
    @param web true if html tags are to be used
     */
  void dump_abruijn(ostream &sout, specnets::abinfo_t *abruijn, bool web);

   /*! \brief Dumps a clusterData file.

     Outputs a clusterData file contents.

    @param sout Output stream
    @param cd the clusterData object to output
    @param title Title precedign the output
    @param web true if html tags are to be used
     */
  void dump_clusterData(ostream &sout, ClusterSet *cd, char *title, bool web);

   /*! \brief Dumps a specSet file.

     Outputs a specSet file contents.

    @param sout Output stream
    @param set The specset to output
    @param title Title precedign the output
    @param web true if html tags are to be used
     */
  void dump_specset(ostream &sout, specnets::SpecSet *set, char *title, bool web);

   /*! \brief Dumps a binArray file.

     Outputs a binArray file contents.

    @param sout Output stream
    @param a BinArray
    @param title Title precedign the output
    @param web true if html tags are to be used
     */
  void dump_binArray(ostream &sout, vector<vector<int> > *a, char *title, bool web);

   /*! \brief Dumps files contents to the screen.

     Outputs all files' contents to the screen.

     */
  void dump(ostream &sout, bool web = false) {
    dump_clusterData(sout, m_clusterData, "Cluster information", web);
    //dump_binArray(sout, m_inputMapping, "input_mapping", web);
    dump_binArray(sout, m_contigIndices, "contigIndices", web);

    dump_binArray(sout, m_contigs_mp, "contigs_mp", web);
    dump_specset( sout, m_contigs_midx, "contigs_midx_all", web);

    dump_binArray(sout, m_homglue_ref_mp, "homglue_ref_mp", web);
    dump_specset( sout, m_homglue_ref_midx, "homglue_ref_midx", web);

    dump_specset(sout, m_homglueMatches, "homglueMatches", web);
    dump_binArray(sout, m_homglue_matches_mp, "homglue_matches_mp", web);
    dump_specset( sout, m_homglue_matches_midx, "homglue_matches_midx", web);


    dump_specset( sout, m_sps_seqs, "sps_seqs", web);
    dump_specset( sout, m_contigsSpectra, "contig_spectra", web);
    //dump_specset( sout, m_starSpectra, "star spectra", web);
  };

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
