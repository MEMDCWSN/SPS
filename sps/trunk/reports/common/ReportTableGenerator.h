////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_GENERATOR_H__
#define __REPORT_TABLE_GENERATOR_H__
////////////////////////////////////////////////////////////////////////////////
#include <iterator>
#include <vector>

#include "spectrum.h"
#include "ClusterData.h"
#include "db_fasta.h"
#include "abruijn.h"
#include "aminoacid.h"

#include "ReportTableInputSpectra.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
#include "ReportTableProteinCoverage.h"
#include "ReportTableProtein.h"
#include "ReportTableHeader.h"

#include "Defines.h"

#include "ReportData.h"

#include "spsFiles.h"

////////////////////////////////////////////////////////////////////////////////
#define OK      1
#define ERROR  -1

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Helper methods for sequence generation
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 /*! \brief Holds a contig/protein mass index pair unit
   */
struct PairContigProteinMassIdx {

  /*! \brief Index in protein for a mass value
   */
  int proteinMassIdx;

  /*! \brief Index in contig for a mass value
   */
  int contigMassIdx;
};


 /*! \brief holds data for CSPS contigs
   */
struct ContigMatchData {

  /*! \brief List of contig mass values
   */
  std::vector<float>                    contigMass;

  /*! \brief List of contig-protein mass pairs.
   */
  std::vector<PairContigProteinMassIdx> pair;

  /*! \brief Contig index
   */
  int                                   contigIndex;
};

////////////////////////////////////////////////////////////////////////////////
//
 /*! \brief Ordering vector element. Used to estabelish contig order when rendering protein coverage report
   */
struct ContigOrdering {

  /*! \brief Contig index
   */
  int contigIndex;

  /*! \brief Start index, in contig
   */
  int startIndex;

  /*! \brief End index,, in contig
   */
  int endIndex;

  /*! \brief Used by sort method. First order by beggining index, then by ending index.
   */
  bool operator<(const ContigOrdering &o) const
  {return (startIndex == o.startIndex ? endIndex < o.endIndex : startIndex < o.startIndex);};

};
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Amino Acid data cell.Used to annotate a mass interval.
   */
struct aaCell {

  /*! \brief aa sequence for this interval
   */
  std::string   aa;

  /*! \brief mass difference from the expected
   */
  double        delta;

  /*! \brief for how many mass intervals it spans
   */
  int           colspan;

  /*! \brief relative index position
   */
  int           startPosition;

  /*! \brief used when calculating protein sequence coverage
   */
  bool          covered;

  /*! \brief star mass value at this point
   */
  double        massPoint;

  /*! \brief Default constructor.
   */
  aaCell() : delta(0.0), covered(false), colspan(0) {};
};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Structure to hold a contig/protein info in protein details context
   */
struct ProteinDetaisContigInfo {

 /*! \brief Cell start position
   */
  int                       startPosition;

 /*! \brief Cell end position
   */
  int                       endPosition;

 /*! \brief Processed AAs for page generation
   */
  std::vector<aaCell>       processedAA;

 /*! \brief Contig name
   */
  std::string								name;

 /*! \brief base 1 contig index
   */
  int                       base1Idx;
};


// contig information indexed by contig index
typedef std::map<int, ProteinDetaisContigInfo> PdProteinDetail;




 /*! \brief Holds information for an entire protein and it's contigs
   */
struct PdProteinInfo {

  /*! \brief protein sequence
   */
	ProteinDetaisContigInfo proteinDetail;

  /*! \brief Holds information for sps contigs.
   */
	PdProteinDetail spsDetails;

  /*! \brief Holds information for csps contigs.
   */
	PdProteinDetail cspsDetails;

  /*! \brief actual protein index
   */
	int ProteinIndex;

  /*! \brief protein length
   */
	int proteinLength;

  /*! \brief protein name
   */
  string proteinName;

  /*! \brief Structure to hold protein table section
   */
  std::string               proteinSequenceEntryData;

  /*! \brief Structure to hold csps generated table section
   */
  std::string 							cspsEntryData;

  /*! \brief Structure to hold sps generated table section
   */
  std::string 							spsEntryData;
};
////////////////////////////////////////////////////////////////////////////////

  /*! \brief Structure to hold sequence mapping
   */
struct SequenceMapping {

  /*! \brief Cell start position
   */
  int                       startPosition;

  /*! \brief Cell end position
   */
  int                       endPosition;

  /*! \brief Processed AAs for page generation
   */
  std::vector<aaCell>       processedAA;

  /*! \brief Contig name
   */
  std::string								name;

  /*! \brief base 1 contig index
   */
  int                       base1Idx;

  /*! \brief sequence 1st mass point
   */
  double                    startMass;

  /*! \brief Clear data
   */
  void clear() {base1Idx=startPosition=endPosition=-1;name.clear();processedAA.clear();};
};
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Structure to hold differential data
   */
struct DiffData {
  vector<pair<int, double> >  abruijnData;
  vector<double>              contigData;

  vector<double>              abDiff;
  vector<double>              contigDiff;
  vector<double>              diffDiff;
};
////////////////////////////////////////////////////////////////////////////////
// Table data
//

  /*! \brief Structure to hold report internal data
   */
struct ReportInternalData {

  /*! \brief Protein index
   */
  int protein;

  /*! \brief Contig index
   */
  int contig;

  /*! \brief Cluster index
   */
  int cluster;

  /*! \brief Input spectra
   */
  int inputSpectra;


  /*! \brief Protein text
   */
  string proteinText;

  /*! \brief Contig text
   */
  string contigText;
  /*! \brief Cluster text
   */
  string clusterText;
  /*! \brief Input spectra text
   */
  string inputSpectraText;

  /*! \brief Cluster consensus index
   */
  int consensus;

  /*! \brief spectra (consensus or input spectra, depending on the presence of cluster layer) for a contig
   */
  vector<int> spectraOfContig;

  /*! \brief number of contigs that match a protein
   */
  int matchedContigs;

  /*! \brief consensus counter (per protein)
   */
  int consensusAcc;

  /*! \brief Matched contig index
   */
  int matchedContig;
  /*! \brief All contigsContig index
   */
  int allContigsContig;

  /*! \brief Homolog index
   */
  int homolog;

  // protein data

  /*! \brief Protein name
   */
  string proteinName;
  /*! \brief Protein description
   */
  string proteinDesc;

  // contig data

  /*! \brief sequence structure to hold reference and homolog post-processed sequences
   */
  SequenceMapping sequenceMappingReference, sequenceMappingHomolog, sequenceMappingDeNovo;

  /*! \brief cluster and spectra sequences (final string)
   */
  string clusterSequenceReference, clusterSequenceHomolog, clusterSequenceDeNovo;

  /*! \brief Reference sequence to be used in tables
   */
  string clusterSequenceReferenceEfective;

  /*! \brief offset masses
   */
  double offsetHomolog, offsetReference;

  /*! \brief references masses
   */
  vector<float> m_referenceMasses;

  /*! \brief Homolog masses
   */
  vector<float> m_homologMasses;

  /*! \brief DeNovo masses
   */
  vector<float> m_deNovoMasses;

  /*! \brief User masses
   */
  vector<float> m_userMasses;

  /*! \brief true if reference and homolog are equal at cluster level.
   */
  bool auxRefCmp;

  /*! \brief statitics data
   */
  float cluster_B;
  /*! \brief statitics data
   */
  float cluster_Y;
  /*! \brief statitics data
   */
  float cluster_BYint;

  /*! \brief statitics data
   */
  float spectra_B;
  /*! \brief statitics data
   */
  float spectra_Y;
  /*! \brief statitics data
   */
  float spectra_BYint;


  /*! \brief Tool used
   */
  string tool;

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report Table Generator class

   Used to generate the report data tables based on the input files.

   */
class ReportTableGenerator {

  /*! \brief input spectra files
   */
  vector<string> m_inputSpectraFiles;

  /*! \brief Abruijn filename
   */
  string m_fn_abruijn;
  /*! \brief Star spectra filename
   */
  string m_fn_star;
  /*! \brief seqs filename
   */
  string m_fn_seqs;

  /*! \brief Project directory
   */
  string m_projectDir;

  /*! \brief annotation model directory
   */
  string m_annotationModelDirectory;

  /*! \brief annotation model filename
   */
  string m_annotationModel;

  /*! \brief annotation model directory (Prm)
   */
  string m_annotationModelDirectoryPrm;

  /*! \brief annotation model filename (Prm)
   */
  string m_annotationModelPrm;

  /*! \brief mass shift value
   */
  float m_massShift;

  /*! \brief mass shift value (Prm)
   */
  float m_massShiftPrm;

  /*! \brief peak mass tolerance
   */
  float m_peakMassTol;

  /*! \brief parent mass tolerance
   */
  float m_parentMassTol;

  /*! \brief resolution
   */
  float m_resolution;

  /*! \brief job information
   */
  string m_jobName;

  /*! \brief User information
   */
  string m_userName;

  /*! \brief tool used
   */
  int m_tool;

  /*! \brief if cluster layer is used
   */
  bool m_noClusters;

  /*! \brief total number of clusters in clusterData.bin
   */
  int m_clusterCount;

  /*! \brief  total number of spectra mapped in clusterData.bin
   */
  int m_spectraCount;

 protected:


  //////////////////////////////////////////////////////////////////////////////
  // Table pointers.

  /*! \brief Header table
   */
  ReportTableHeader *tableHeader;

  /*! \brief Proteins table
   */
  ReportTableProtein *tableProteins;

  /*! \brief Proteins table
   */
  ReportTableProteinCoverage *tableProteinCoverage;

  /*! \brief Contigs table
   */
  ReportTableContig *tableContigs;

  /*! \brief Cluster consensus table
   */
  ReportTableClusterConsensus *tableClusterConsensus;

  /*! \brief Input spectra table
   */
  ReportTableInputSpectra *tableInputSpectra;


  //////////////////////////////////////////////////////////////////////////////
  // Methods to load data files

  /*! \brief filename subtraction helper method - converts absolute path to relative
   */
  virtual string pathAbsoluteToRelative(const string &projectDir, const string &fileName);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build and edit the table

  /*! \brief Add the table headings
   */
  virtual void buildTableHeadings(void);

  /*! \brief Build the initial page table
   */
  virtual int buildTableHeader(ReportInternalData &data);

  /*! \brief Build proteins table
    The table has the following structure, per row:
     cells[row][0] -> text   --> index of protein
     cells[row][1] -> text   --> protein name
     cells[row][2] -> text   --> protein description
     cells[row][3] -> text   --> number of contigs
     cells[row][4] -> text   --> number of spectra
     cells[row][5] -> text   --> amino acids
     cells[row][6] -> text   --> coverage %
   */
  virtual int buildTableProteins(ReportInternalData &data);

  /*! \brief Build protein detais (civerage) table
    The table has the following structure, per row:
     cells[0][0] -> text   --> protein index (1-based)
     cells[1][0] -> text   --> protein name
     cells[2][0] -> text   --> Protein length
     cells[3][0] -> text   --> Protein sequence
     cells[4][0] -> text   --> csps contig data
     cells[5][0] -> text   --> sps contig data
   */
  virtual int buildTableProteinCoverage(ReportInternalData &data);

  /*! \brief build contigs table - add contigs that map to a protein
   */
  virtual int buildTableContigs(ReportInternalData &data);

  /*! \brief build contigs table - add contigs that do not map to proteins
   */
  virtual int buildTableContigsOrphan(ReportInternalData &data);

  /*! \brief add a single contig to the table
    The table has the following structure, per row:
     cells[row][0] -> text   --> Contig index (1-based)
     cells[row][1] -> text   --> Protein index (1-based)
     cells[row][2] -> text   --> Number of spectra
     cells[row][3] -> text   --> Reference sequence
     cells[row][4] -> text   --> Homolog sequence
     cells[row][5] -> text   --> DeNov sequence
     cells[row][6] -> text   --> User sequence
     cells[row][7] -> text   --> Protein name
     cells[row][8] -> text   --> Protein description
     cells[row][9] -> text   --> spectra (coma separated indexes)
     cells[row][10] -> text  --> Reference intervals
     cells[row][11] -> text  --> Homolog intervals
     cells[row][12] -> text  --> Reference offset
     cells[row][13] -> text  --> Homolog offset
     cells[row][14] -> text  --> reverse flag
     cells[row][15] -> text  --> file names needed: 15 abruijn
     cells[row][16] -> text  --> file names needed: 16 stars
     cells[row][17] -> text  --> file names needed: 17 sps_seqs
     cells[row][18] -> text  --> tool used
     cells[row][19] -> text  --> Grouped homolog sequence
     cells[row][20] -> text  --> Grouped reference sequence
   */
  virtual int buildTableContig(ReportInternalData &data);

  /*! \brief build clustar consensus table
    The table has the following structure, per row:
     cells[row][0] -> text   --> index in specset (consensus spectra index ; 1-based)
     cells[row][1] -> text   --> contig ID (1-based)
     cells[row][2] -> text   --> Protein ID (1-based)
     cells[row][3] -> text   --> Reference sequence --> generateSequence()
     cells[row][4] -> text   --> Homolog sequence --> generateSequence()
     cells[row][5] -> text   --> DeNovo sequence --> generateSequence()
     cells[row][6] -> text   --> User sequence
     cells[row][7] -> text   --> mass value, from specs[i].parentMass
     cells[row][8] -> text   --> charge value, from specs[i].parentCharge
     cells[row][9] -> text   --> B%
     cells[row][10] -> text   --> Y&
     cells[row][11] -> text  --> BY Int %
     cells[row][12] -> text  --> protein name
     cells[row][13] -> text  --> Tool used
     cells[row][14] -> text  --> Fragmentation model
   */
  virtual int buildTableCluster(ReportInternalData &data);

  /*! \brief build input spectra table
    The table has the following structure, per row:
     cells[row][0] -> text   --> table index (used for updates)
     cells[row][1] -> text   --> spectrum index in specset (1-based)
     cells[row][2] -> text   --> scan #
     cells[row][3] -> text   --> cluster index (used when filtered by cluster consensus ; 1-based)
     cells[row][4] -> text   --> Protein name
     cells[row][5] -> text   --> spectrum file index
     cells[row][6] -> text   --> spectrum file name
     cells[row][7] -> text   --> Reference sequence --> generateSequence()
     cells[row][8] -> text   --> Homolog sequence --> generateSequence()
     cells[row][9] -> text   --> DeNovo sequence --> generateSequence()
     cells[row][10] -> text  --> User sequence
     cells[row][11] -> text  --> mass value, from specs[i].parentMass
     cells[row][12] -> text  --> charge value, from specs[i].parentCharge
     cells[row][13] -> text  --> B%
     cells[row][14] -> text  --> Y%
     cells[row][15] -> text  --> BY Int %
     cells[row][16] -> text  --> Original input spectra filename
     cells[row][17] -> text  --> Contig number
     cells[row][18] -> text  --> Protein index
     cells[row][19] -> text  --> Tool used
     cells[row][20] -> text  --> Fragmentation model
   */
  virtual int buildTableInputSpectra(ReportInternalData &data);

  /*! \brief build input spectra table when there is no cluster layer
    The table has the following structure, per row:
     cells[row][0] -> text   --> table index (used for updates)
     cells[row][1] -> text   --> spectrum index in specset (1-based)
     cells[row][2] -> text   --> scan #
     cells[row][3] -> text   --> cluster index (used when filtered by cluster consensus ; 1-based)
     cells[row][4] -> text   --> Protein name
     cells[row][5] -> text   --> spectrum file index
     cells[row][6] -> text   --> spectrum file name
     cells[row][7] -> text   --> Reference sequence --> generateSequence()
     cells[row][8] -> text   --> Homolog sequence --> generateSequence()
     cells[row][9] -> text   --> DeNovo sequence --> generateSequence()
     cells[row][10] -> text  --> User sequence
     cells[row][11] -> text  --> mass value, from specs[i].parentMass
     cells[row][12] -> text  --> charge value, from specs[i].parentCharge
     cells[row][13] -> text  --> B%
     cells[row][14] -> text  --> Y%
     cells[row][15] -> text  --> BY Int %
     cells[row][16] -> text  --> Original input spectra filename
     cells[row][17] -> text  --> Contig number
     cells[row][18] -> text  --> Protein index
     cells[row][19] -> text  --> Tool used
     cells[row][20] -> text  --> Fragmentation model   */
  virtual int buildTableInputSpectra2(ReportInternalData &data);


  //////////////////////////////////////////////////////////////////////////////
  // Methods for data access

  /*! \brief get input spectra original file and spectra indices from combined spectra index. Used when consensus cluster is not present
   */
  void getinputSpectraFromCombinedSpectra(int &fileIndex, int &spectraIndex, int combinedSpectraIndex);

  /*! \brief Get the cluster consensus index, given the file index and input spectra index
   */
  int  getConsensusFromInputSpectra(int fileIndex, int inputSpectra); 

  /*! \brief get the file index and input spectra index pair list given the cluster consensus index
   */
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(int consensus);

  /*! \brief get the contig index given the cluster consensus index
   */
  int  getContigFromConsensus(int consensus);

  /*! \brief get the total number of clusters
   */
  int getClusterCount(void);

  /*! \brief get the total number of spectra
   */
  int getSpectraCount(void);

  /*! \brief Get the cluster consensus index vector given the contig index
   */
  void getConsensusFromContig(int contig, vector<int>&);

  /*! \brief Get contig index given file index and input spectra index pair
   */
  int  getContigFromInputSpectra(int fileIndex, int inputSpectra);

  /*! \brief Get the 'contig that maps to a protein' index given the 'all contigs' index
   */
  int  getContigFromAllContig(int contig);

  /*! \brief Get the 'all contigs' index list given a 'contig that maps to a protein' index list
   */
  void getAllContigFromContig(vector<int> &contigs, vector<int> &allContigs);

  /*! \brief Get the 'all contigs' index given the 'contig that maps to a protein' index
   */
  int  getAllContigFromContig(int contig);

  /*! \brief Get a protein index given a 'contig that maps to a protein' index
   */
  int  getProteinFromHomolog(int homolog);

  /*! \brief Get a 'contig that maps to a protein' index given a protein index
   */
  void getHomologFromProtein(int protein, vector<int> &ret);

  /*! \brief get a protein given a reference index
   */
  int getProteinFromReference(int reference);

  /*! \brief Get a protein name given it's index
   */
  string getProteinName(int protein);

  /*! \brief Get a protein description given it's index
   */
  string getProteinDescription(int protein);

  /*! \brief Get a protein sequence given it's index
   */
  string getProteinSequence(int protein);

  /*! \brief replaces | char by space
   */
  void cleanProteinName(string &proteinName);

  /*! \brief get the model name
   */
  string getModelName(Spectrum *spectrum);


  bool checkContig(int contig);

  //////////////////////////////////////////////////////////////////////////////
  // Protein Coverage related methods

  /*! \brief put protein data in internal structure for processing
   */
  bool processProteinsFile(int protein, PdProteinInfo &proteinData);

  /*! \brief put contig data in internal structure for processing
   */
  bool processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs, bool translate = false, bool csps = false);

  /*! \brief process contig data
   */
  void processContigs(int ProteinIndex, PdProteinInfo &proteinData, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source);

  /*! \brief replace contig IDs by contig names, from contig names file
   */
  void populateContigNames(PdProteinDetail &contig);

  /*! \brief get a specific contig's name
   */
  string getContigName(int i);

  /*! \brief get a specific's contig data
   */
  ContigMatchData *getContigData( std::map<int, std::vector<ContigMatchData> > &contig, int contigIndex);

  /*! \brief generate table field for all contigs
   */
  int generateCoverageOutput(PdProteinInfo &proteinData);

  /*! \brief generate table section for a contig
   */
  int generateOutputContig(PdProteinInfo &proteinData, vector<int> &vectorID, std::string &page, PdProteinDetail &contig, bool link);

  /*! \brief contig sorting
   */
  void getOrder(PdProteinDetail &contig, vector<int> &order);

  /*! \brief get ID from name
   */
  string getIntFromSeqName(string seq);


  //////////////////////////////////////////////////////////////////////////////
  // Sequence and statistics related methods

  /*! \brief abruijn contig reverse
   */
  void getReversedAbruijn(int contig, abContigData_t &nodes);

  /*! \brief abruijn star state (reverse)
   */
  bool getAbruijnStarState(int contig, int star);
  bool getContigState(unsigned contig);
  bool getStarIndex(vector<pair<int, double> > &data, int star, double &value);
  int  getSequenceIndex(SequenceMapping &data, int index);
  bool getCspsContigState(unsigned contig);

  /*! \brief generate deNovo sequence, given a spectrum
   */
  string getSequenceDenovo(specnets::Spectrum &s);

  /*! \brief Generate homolog sequence, given the contig index
   */
  int getSequenceHomolog(ReportInternalData &data);

  /*! \brief Generate Reference sequence, given the contig index
   */
  int getSequenceReference(ReportInternalData &data);

  /*! \brief translate a sequence in SequenceMapping structure into a string sequence
   */
  string translateSequence(SequenceMapping &sequenceMapping, bool group = false);

  /*! \brief translate a sequence form a string to the SequenceMapping structure
   */
  void translateSequenceReverse(SequenceMapping &oSequenceMapping, string &iSequence);

  /*! \brief Remove the () in sequences like MM(MM)MM, where M is any aa
   */
  bool stringExcessParamsFromSequence(string &iSequence);

  /*! \brief get a sequence for clusters consensus spectrum given a contig sequence
   */
  void propagateSequence(SequenceMapping &oSequenceMapping, int contig, int consensus, int homolog, SequenceMapping &iSequenceMapping);

  /*! \brief sequence propagation helper method
   */
  void getSequenceBetween(aaCell &ret, SequenceMapping &data, int start, int end);

  /*! \brief protein coverage calculation
   */
  int getProteinCoverage(int protein, specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::SpecSet &contigSpectra, string &sequence, int &proteinSize, int &covered);

  /*! \brief get mass prefix and suffix for sequences
   */
  void getMassPrefix(int contig, int consensus, double &prefix, double &suffix);

  /*! \brief process mass intervals as a single string
   */
  void processMassIntervals(vector<float> &mass, string &intervals);

  /*! \brief Base sequence generation method
   */
  int getSequence(int contigIdx, specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::Spectrum &contigSpectra, SequenceMapping &sequenceMapping, double &leadMass);

  /*! \brief First step in common sequence generation
   */
  int generateSequenceStep0( SequenceMapping &proteinData, int proteinIndex);

  /*! \brief Second step in common sequence generation
   */
  int generateSequenceStep1(specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::Spectrum &contigSpectra, ContigMatchData  &contigInfo, int contigIdx);

  /*! \brief Third step in common sequence generation
   */
  int generateSequenceStep2(int ProteinIndex, SequenceMapping &proteinData, ContigMatchData &source, SequenceMapping &m_sequenceMapping);

  /*! \brief Third step in protein coverage generation
   */
  int generateSequenceStep2B(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source);

  /*! \brief Generate statistics (B%, Y% and BY intensity %)
   */
  void generateStatistics(specnets::Spectrum &spectrum, string &peptide, float &m_B, float &m_Y, float &m_BYint);
  //void generateStatistics(specnets::Spectrum &spectrum, string &peptide);

  /*! \brief get input spectra based on GenoMS options
   */
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(ReportInternalData &data);

  /*! \brief build mass instervals
   */
  void buildMassIntervals(vector<float> &masses, SequenceMapping &sm);


 public:


  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data needed for several table fields.

  /*! \brief Spectrum data needed for several table fields.
   */
  SpsFiles *spsFiles;

  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Default constructor
   */
  ReportTableGenerator();

    //@}

    //! \name DESTRUCTOR
    //@{

  // Destructor deletes colTypes vector
  /*! \brief Default destructor
   */
  ~ReportTableGenerator();

  /*! \brief set data files
   */
  virtual void setSpsFiles(SpsFiles *f) {spsFiles = f;};

  /*! \brief initialize data. Shoud be invoked after setting sps data files
   */
  virtual void init(const ReportGeneratorData &reportGeneratorData);

  /*! \brief Builds table from input data
   */
  virtual int buildTables(ReportTableHeader *,ReportTableProtein *, ReportTableProteinCoverage *, ReportTableContig *, ReportTableClusterConsensus *, ReportTableInputSpectra *);

  /*! \brief Dumps contig data
   */
  void dump_contigData(ostream &sout, int contig, int star, int homolog, DiffData &diffData, SequenceMapping &iSequenceMapping);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
