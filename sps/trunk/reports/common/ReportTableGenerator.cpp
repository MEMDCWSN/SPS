////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <algorithm>

#include "ReportTableGenerator.h"
#include "utils.h"
#include "Tokenizer.h"

#include "SpectrumAnnotStatistics.h"
#include "spectrum_scoring.h"

#include "Logger.h"
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

using namespace std;
using namespace specnets;

////////////////////////////////////////////////////////////////////////////////
#define TABLE_SEP_L1  '|'
#define TABLE_SEP_L2  '@'
#define TABLE_SEP_L3  '&'
#define TABLE_SEP_L4  '!'

#define TOOL_SPS      "SPS"
#define TOOL_GENOMS   "GenoMS"

////////////////////////////////////////////////////////////////////////////////
// Report table base methods
////////////////////////////////////////////////////////////////////////////////
ReportTableGenerator::ReportTableGenerator() :
  spsFiles(NULL),
  m_spectraCount(-1),
  m_clusterCount(-1)
{
}
////////////////////////////////////////////////////////////////////////////////
ReportTableGenerator::~ReportTableGenerator()
{
}
////////////////////////////////////////////////////////////////////////////////
// File loading section.
////////////////////////////////////////////////////////////////////////////////
// Main load method. Checks for a table to load. If non existent, calls the loadData method, and then saves the table
void ReportTableGenerator::init(const ReportGeneratorData &reportGeneratorData)
{
  m_projectDir = reportGeneratorData.projectDir;

  // annotation model directory
  m_annotationModelDirectory = reportGeneratorData.annotationModelDir;
  // annotation model filename
  m_annotationModel = reportGeneratorData.annotationModel;
  // annotation model directory for PRM model
  m_annotationModelDirectoryPrm = reportGeneratorData.annotationModelDirPrm;
  // Annotation model for PRM spectra
  m_annotationModelPrm = reportGeneratorData.annotationModelPrm;

  // mass shift value
  m_massShift = reportGeneratorData.massShift;
  // mass shift value for PRM
  m_massShiftPrm = reportGeneratorData.massShiftPrm;


  // peak mass tolerance
  m_peakMassTol = reportGeneratorData.peakMassTol;
  // parent mass tolerance
  m_parentMassTol = reportGeneratorData.parentMassTol;
  // resolution
  m_resolution = reportGeneratorData.resolution;

  m_jobName = reportGeneratorData.job;
  m_userName = reportGeneratorData.user;

  m_tool  = reportGeneratorData.tool;

  m_noClusters = reportGeneratorData.noClusters;

  if(m_noClusters) {
    DEBUG_MSG("Cluster layer deactivated");
  } else {
    DEBUG_MSG("Cluster layer activated");
  }

  //get total number of clusters
  m_clusterCount = getClusterCount();
  // get total number of spectra
  m_spectraCount = getSpectraCount();

  m_fn_star = reportGeneratorData.filenameStarSpectra;
  m_fn_abruijn = reportGeneratorData.filenameAbruijn;
  m_fn_seqs = reportGeneratorData.filenameSpsSeqs;
}
////////////////////////////////////////////////////////////////////////////////
// Convert an absolute path to relative.
string ReportTableGenerator::pathAbsoluteToRelative(const string &projectDir, const string &fileName)
{
  size_t found = fileName.find(projectDir);
  if(found == 0) {
    string aux = fileName.substr(projectDir.length());
    return aux;
  }
  return fileName;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTables(ReportTableHeader           *tHeader,
                                      ReportTableProtein          *tProteins,
                                      ReportTableProteinCoverage  *tProteinCoverage,
                                      ReportTableContig           *tContig,
                                      ReportTableClusterConsensus *tClusterConsensus,
                                      ReportTableInputSpectra     *tInputSpectra)
{
  // Header table
  tableHeader           = tHeader;
  // Proteins table
  tableProteins         = tProteins;
  // Protein Coverage table
  tableProteinCoverage  = tProteinCoverage;
  //Contigs table
  tableContigs          = tContig;
  // CLuster consensus table
  tableClusterConsensus = tClusterConsensus;
  // Input spectra table
  tableInputSpectra     = tInputSpectra;

  DEBUG_MSG("Generating tables");

  // add the table headings text
  buildTableHeadings();

  // report table exchange data
  ReportInternalData data;

  // build tables
  if(buildTableProteins(data) == ERROR) {
    ERROR_MSG("Error build tables from protein level.");
    return ERROR;
  }

  // add orphan contigs to tables
  if(buildTableContigsOrphan(data) == ERROR) {
    ERROR_MSG("Error build tables from orphan contigs.");
    return ERROR;
  }
  // buils table header
  if(buildTableHeader(data) == ERROR) {
    ERROR_MSG("Error build header table.");
    return ERROR;
  }


  // sort tables
  DEBUG_MSG("Sorting tables");
  tableContigs->sortTable(0);
  tableClusterConsensus->sortTable(0);
  tableInputSpectra->sortTable(1);

  // populate IDs in input spectra table. This operation should be done after sort() to allow dictomic search when updating the table user sequence
  DEBUG_MSG("Populating spectra tables unique IDs");
  tableInputSpectra->populateIDs();

  // done
  DEBUG_MSG("Done building tables.");

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::buildTableHeadings(void)
{
  // Header table
  if(tableHeader) {
    tableHeader->addHeading("JobName");
    tableHeader->addHeading("UserName");
    tableHeader->addHeading("Status");
    tableHeader->addHeading("ElapsedTime");
    tableHeader->addHeading("Log");
    tableHeader->addHeading("LinkToContigs");
    tableHeader->addHeading("LinkToProteins");
    tableHeader->addHeading("ClusterListFilename");
    tableHeader->addHeading("SpectraFileNames");
    tableHeader->addHeading("SpectraFileLinks");
  }

  // Proteins table
  if(tableProteins) {
    tableProteins->addHeading("ProteinIndex");
    tableProteins->addHeading("ProteinName");
    tableProteins->addHeading("ProteinDescription");
    tableProteins->addHeading("ContigCount");
    tableProteins->addHeading("SpectraCount");
    tableProteins->addHeading("AAsCount");
    tableProteins->addHeading("CoveragePercentage");
    tableProteins->addHeading("ProteinSequence");
  }

  // Protein Coverage table
  if(tableProteinCoverage) {
    tableProteinCoverage->addHeading("ProteinIndex");
    tableProteinCoverage->addHeading("ProteinName");
    tableProteinCoverage->addHeading("ReferenceSequence");
    tableProteinCoverage->addHeading("ProteinSequence");
    tableProteinCoverage->addHeading("CspsContigData");
    tableProteinCoverage->addHeading("SpsContigData");
  }

  //Contigs table
  if(tableContigs) {
    // data for report
    tableContigs->addHeading("ContigIndex");
    tableContigs->addHeading("ProteinIndex");
    tableContigs->addHeading("SpectraCount");
    tableContigs->addHeading("ReferenceSequence");
    tableContigs->addHeading("HomologSequence");
    tableContigs->addHeading("DeNovoSequence");
    tableContigs->addHeading("UserSequence");
    tableContigs->addHeading("ProteinName");
    tableContigs->addHeading("ProteinDescription");
    tableContigs->addHeading("Spectra");
    // data for contplot
    tableContigs->addHeading("ReferenceMassIntervals");
    tableContigs->addHeading("HomologMassIntervals");
    tableContigs->addHeading("ReferenceOffset");
    tableContigs->addHeading("HomologOffset");
    tableContigs->addHeading("ReverseFlag");
    tableContigs->addHeading("AbruijnFilename");
    tableContigs->addHeading("StarsFilename");
    tableContigs->addHeading("SpsSeqsFilename");
    // data for report
    tableContigs->addHeading("Tool");
    // data for contplot
    tableContigs->addHeading("GroupedHomologSequence");
    tableContigs->addHeading("GroupedReferenceSequence");
  }

  // Cluster consensus table
  if(tableClusterConsensus) {
    tableClusterConsensus->addHeading("ClusterIndex");
    tableClusterConsensus->addHeading("ContigIndex");
    tableClusterConsensus->addHeading("ProteinIndex");
    tableClusterConsensus->addHeading("ReferenceSequence");
    tableClusterConsensus->addHeading("HomologSequence");
    tableClusterConsensus->addHeading("DeNovoSequence");
    tableClusterConsensus->addHeading("UserSequence");
    tableClusterConsensus->addHeading("ParentMass");
    tableClusterConsensus->addHeading("Charge");
    tableClusterConsensus->addHeading("B_Per");
    tableClusterConsensus->addHeading("Y_Per");
    tableClusterConsensus->addHeading("BY_Int");
    tableClusterConsensus->addHeading("ProteinName");
    tableClusterConsensus->addHeading("Tool");
    tableClusterConsensus->addHeading("Model");
    //tableClusterConsensus->addHeading("ModelFile");
    //tableClusterConsensus->addHeading("Shift");
  }

  // Input spectra table
  if(tableInputSpectra) {
    tableInputSpectra->addHeading("TableIndex");
    tableInputSpectra->addHeading("SpectrumIndex");
    tableInputSpectra->addHeading("Scan");
    tableInputSpectra->addHeading("ClusterIndex");
    tableInputSpectra->addHeading("ProteinName");
    tableInputSpectra->addHeading("PklbinFileIndex");
    tableInputSpectra->addHeading("PklbinFilename");
    tableInputSpectra->addHeading("ReferenceSequence");
    tableInputSpectra->addHeading("HomologSequence");
    tableInputSpectra->addHeading("DeNovoSequence");
    tableInputSpectra->addHeading("UserSequence");
    tableInputSpectra->addHeading("Mass");
    tableInputSpectra->addHeading("Charge");
    tableInputSpectra->addHeading("B_Per");
    tableInputSpectra->addHeading("Y_Per");
    tableInputSpectra->addHeading("BY_Int");
    tableInputSpectra->addHeading("OriginalFilename");
    tableInputSpectra->addHeading("ContigIndex");
    tableInputSpectra->addHeading("ProteinIndex");
    tableInputSpectra->addHeading("Tool");
    tableInputSpectra->addHeading("FragmentationModel");
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableHeader(ReportInternalData &data)
{
  DEBUG_MSG("Building header table");

  string links, texts;

  if(spsFiles->m_input_index)
    for(int i = 0 ; i < spsFiles->m_input_index->size() ; i++)
      m_inputSpectraFiles.push_back( (*spsFiles->m_input_index)[i] );

  // sort elements for duplicate removal
  //if(m_inputSpectraFiles.size())
  //  sort (m_inputSpectraFiles.begin(), m_inputSpectraFiles.end(), stringSortCompare);

  // duplicate find
  //vector<string>::iterator it;

  // using default comparison
  //it = unique (m_inputSpectraFiles.begin(), m_inputSpectraFiles.end(), stringUniqueCompare);

  // remove extra items
 // m_inputSpectraFiles.resize( it - m_inputSpectraFiles.begin() );

  // cicle thru valid elements
  for(int i = 0 ; i < m_inputSpectraFiles.size() ; i++){

    if(i) {
      texts += TABLE_SEP_L1;
      links += TABLE_SEP_L1;
    }
    texts += m_inputSpectraFiles[i];
    links += parseInt(i);
  }


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;
  // Column 0: Job name
  dataRow.push_back(m_jobName);
  // Column 1: User name
  dataRow.push_back(m_userName);
  // Column 2: Status
  dataRow.push_back("Finished");
  // Column 3: Elapsed time
  dataRow.push_back("");
  // Column 4: Log
  dataRow.push_back("");
  // Column 5: Link to group by contig
  dataRow.push_back("contigs");
  // Column 6: Link to proteins
  dataRow.push_back("proteins");
  // Column 7: File name of cluster list
  dataRow.push_back("cluster.txt");
  // Column 8: File names
  dataRow.push_back(texts);
  // Column 9: File links
  dataRow.push_back(links);

  // store the row
  if(tableHeader)
    tableHeader->addDataRow(dataRow);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableProteins(ReportInternalData &data)
{
  if(!spsFiles->m_fasta) {
    DEBUG_MSG("Proteins file not present. Skipping table generation from protein leve.");
    return OK;
  }

  DEBUG_MSG("Building table 'proteins'");

   // cycle thru all proteins
  for(unsigned i = 0 ; i < spsFiles->m_fasta->IDs.size() ; i++ ) {

    data.protein = i;

    DEBUG_MSG("Protein #" << i);

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table and downstream tables

    // initialize consensus spectra counter for the protein
    data.consensusAcc = 0;
    // get protein name
    data.proteinName = getProteinName(i);
    // get protein desc
    data.proteinDesc = getProteinDescription(i);

    DEBUG_MSG("Protein name: " << data.proteinName);

    string sequence;
    int    proteinSize;
    int    covered;
    double percent;

    //if(spsFiles->m_homglue_ref_midx && spsFiles->m_homglue_ref_mp && spsFiles->m_contigsSpectra)
    //  getProteinCoverage(i, *spsFiles->m_homglue_ref_midx, *spsFiles->m_homglue_ref_mp, *spsFiles->m_contigsSpectra, sequence, proteinSize, covered);

    if(spsFiles->m_contigs_midx && spsFiles->m_contigs_mp && spsFiles->m_contigsSpectra)
      getProteinCoverage(i, *spsFiles->m_contigs_midx, *spsFiles->m_contigs_mp, *spsFiles->m_contigsSpectra, sequence, proteinSize, covered);

    percent = 100.0 * (double)covered / (double)proteinSize;

    ////////////////////////////////////////////////////////////////////////////
    // downstream table rows

    // build contig entries for this protein
    buildTableContigs(data);

    // if the protein has not matched contigs. ignore it
    if(data.matchedContigs < 1)
      continue;

    // protein coverage
    if(buildTableProteinCoverage(data) == ERROR)
      return ERROR;


    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;
    // Column 0: index of protein  (1-based)
    dataRow.push_back(parseInt(i + 1));
    // Column 1: protein name
    dataRow.push_back(data.proteinName);
    // Column 2: protein description
    dataRow.push_back(data.proteinDesc);
    // Column 3: number of contigs
    dataRow.push_back(parseInt(data.matchedContigs));
    // Column 4: number of spectra
    dataRow.push_back(parseInt(data.consensusAcc));
    // Column 5: amino acids
    dataRow.push_back(parseInt(covered));
    // Column 6: coverage %
    dataRow.push_back(parseFloat(percent,1));
    // Column 7: protein sequence
    dataRow.push_back(sequence);

    // store the row
    if(tableProteins)
      tableProteins->addDataRow(dataRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Protein Coverage
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text          --> Protein ID
// cells[row][1] -> text          --> Protein name
// cells[row][2] -> text          --> Protein length (AAs)
// cells[row][3] -> text list     --> Protein sequence, separated by |
// cells[row][4] -> Contig data   --> CSPS Contigs, separated by |
// cells[row][5] -> Contig data   --> SPS Contigs, separated by |
//      Contig data: : items separated by &
//        0 -> Contig ID
//        1 -> Contig name
//        2 -> Contig start
//        3 -> Contig end
//        4 -> Contig items
//           Contig Item: items separated by @
//              0 -> Beginning
//              1 -> Span
//              0 -> Content
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableProteinCoverage(ReportInternalData &data)
{
  DEBUG_MSG("Adding protein coverage entry for protein " << data.protein);

  // where to hold protein structure data
  PdProteinInfo proteinData;

  // Process proteins file
  processProteinsFile(data.protein, proteinData);

  // Structure to hold CSPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_cspsContigInfo;

  // Structure to hold SPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_spsContigInfo;

  // Process csps files
  // CSPS_CONTIG_MATCHES_IDX  /  CSPS_CONTIG_MATCHES  /  CSPS_CONTIG_SPECTRA
  // /homology/homglue_matches_midx.pklbin  /  /homology/homglue_matches_mp.bin  /  /homology/homglue_matches.pklbin
  DEBUG_MSG("Processing cSPS files.");
  processContigFiles(*spsFiles->m_homglue_matches_midx, *spsFiles->m_homglue_matches_mp, *spsFiles->m_homglueMatches, m_cspsContigInfo, false, true);

  // Process sps files
  // SPS_CONTIG_MATCHES_IDX / SPS_CONTIG_MATCHES / SPS_CONTIG_SPECTRA
  // /homology/homglue_ref_midx.pklbin  /  /homology/homglue_ref_mp.bin  /  /spectra/contigs.pklbin
  DEBUG_MSG("Processing SPS files.");
  //processContigFiles(*spsFiles->m_homglue_ref_midx, *spsFiles->m_homglue_ref_mp, *spsFiles->m_contigsSpectra, m_spsContigInfo);

  //Spectrum sp = (*spsFiles->m_sps_seqs)[data.allContigsContig];

  processContigFiles(*spsFiles->m_contigs_midx, *spsFiles->m_contigs_mp, *spsFiles->m_contigsSpectra, m_spsContigInfo);
  //processContigFiles(*spsFiles->m_contigs_midx, *spsFiles->m_contigs_mp, *spsFiles->m_sps_seqs, m_spsContigInfo, true);


  // second step csps processing
  DEBUG_MSG("Processing cSPS contigs.");
  processContigs(data.protein, proteinData, proteinData.cspsDetails, m_cspsContigInfo);

  // second step sps processing
  DEBUG_MSG("Processing SPS contigs.");
  processContigs(data.protein, proteinData, proteinData.spsDetails, m_spsContigInfo);

  // Replace sps contig numbers by their names
  DEBUG_MSG("Populating contig names.");
	populateContigNames(proteinData.spsDetails);

  // Generate table field output data
  DEBUG_MSG("Generating coverage outputs.");
  generateCoverageOutput(proteinData);


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;

  // Column 0: protein index (1-based)
  dataRow.push_back(parseInt(data.protein + 1));
  // column 1: protein name
  dataRow.push_back(data.proteinName);
  // column 2: Protein length
  dataRow.push_back(parseInt(proteinData.proteinLength));
  // column 3: Protein sequence
  dataRow.push_back(proteinData.proteinSequenceEntryData);
  // column 4: csps contig data
  dataRow.push_back(proteinData.cspsEntryData);
  // column 5: sps contig data
  dataRow.push_back(proteinData.spsEntryData);

  // store the row
  if(tableProteinCoverage)
    tableProteinCoverage->addDataRow(dataRow);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Contigs
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> contig index (1-based)
// cells[row][1] -> text   --> number of spectra
// cells[row][2] -> text   --> Reference sequence  --> generateSequence(i, return)
// cells[row][3] -> text   --> Homolog sequence   --> generateSequence(i, return)
// cells[row][4] -> text   --> DeNovo sequence --> generateSequence(i, return)
// cells[row][5] -> text   --> User sequence
// cells[row][6] -> text   --> Protein name
// cells[row][7] -> text   --> Protein description
// cells[row][8] -> text   --> spectra (coma separated indexes)
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContigs(ReportInternalData &data)
{
  DEBUG_MSG("Building table 'contig' mapping to protein #" << data.protein);

  ////////////////////////////////////////////////////////////////////////////
  // Get contigs that match the current protein

  // contig IDs that match a protein
  vector<int> contigsMatchProtein;
  getHomologFromProtein(data.protein, contigsMatchProtein);

  // get number of contigs that match a protein
  data.matchedContigs = contigsMatchProtein.size();

  // get all contigs for this protein
  vector<int> contigsOfProtein;
  getAllContigFromContig(contigsMatchProtein, contigsOfProtein);

  ////////////////////////////////////////////////////////////////////////////
  // cycle thru found contigs

  for(int i = 0 ; i < contigsOfProtein.size() ; i++) {
    data.contig           = i;
    data.matchedContig    = contigsMatchProtein[i];
    data.allContigsContig = contigsOfProtein[i];
    if(buildTableContig(data) == ERROR)
      return ERROR;
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContigsOrphan(ReportInternalData &data)
{
  DEBUG_MSG("Adding orphan contigs");

  if(!spsFiles->m_sps_seqs) {
    ERROR_MSG("sps_seqs data not found.");
    return ERROR;
  }

  data.proteinName = "";
  data.proteinDesc = "";

  data.protein = -1;

  for(int i = 0 ; i < spsFiles->m_sps_seqs->size() ; i++)
    if((*spsFiles->m_sps_seqs)[i].size() != 0) {
      // get homolog
      int homolog = getContigFromAllContig(i);
      if(homolog == -1 || !spsFiles->m_fasta) {
        // add debug message
        DEBUG_MSG("Adding orphan contig #" << i);

        data.contig           = i;
        data.matchedContig    = -1;
        data.allContigsContig = i;
        // add the contig entry to the contig's table
        if(buildTableContig(data) == ERROR)
          return ERROR;
      }
    }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContig(ReportInternalData &data)
{
  DEBUG_MSG("-- Adding contig #" << data.allContigsContig << ", mapped to protein " << data.protein);

  // specify protein ID text
  data.proteinText = (data.protein == -1 ? "-1" : parseInt(data.protein + 1));
  // specify contig ID as text
  data.contigText = parseInt(data.allContigsContig + 1);

  // test for contig already in table
  if(tableContigs->find(data.contigText) != ERROR) {
    DEBUG_MSG("contig #" << data.contig << " already in table! Skipping. ");
    return OK;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Check contig consistency
  if(!checkContig(data.allContigsContig)) {
    WARN_MSG("Contig " << data.allContigsContig << " has empty vertices!!!");
    return OK;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Sequences

  // reference sequence
  string m_sequenceReference, m_sequenceReferenceGroup;
  // homolog sequence
  string m_sequenceHomolog, m_sequenceHomologGroup;
  // deNovo sequence
  string m_sequenceDeNovo;

  // get contig that maps to a protein
  data.homolog = getContigFromAllContig(data.allContigsContig);
  // if ContigIndices not present, use allContigsContig
  if(data.homolog == -2)
    data.homolog = data.allContigsContig;

  // build homolog sequence
  if(spsFiles->m_fasta)
    getSequenceHomolog(data);
  // build reference sequence
  if(spsFiles->m_fasta)
    getSequenceReference(data);

  data.offsetHomolog = -data.offsetHomolog;
  data.offsetReference = -data.offsetReference;

  // build deNovo sequence
  if(spsFiles->m_contigsSpectra) {
    Spectrum sp = (*spsFiles->m_sps_seqs)[data.allContigsContig];
    if(getContigState(data.homolog))
      sp.reverse(0.0 - AAJumps::massH2O);
    m_sequenceDeNovo = getSequenceDenovo(sp);
  }
  // build sequence deNovo structure
  translateSequenceReverse(data.sequenceMappingDeNovo, m_sequenceDeNovo);

  // translate homolog sequence
  m_sequenceHomolog = (spsFiles->m_fasta ? translateSequence(data.sequenceMappingHomolog) : "");
  // translate reference sequence
  m_sequenceReference = (spsFiles->m_fasta ? translateSequence(data.sequenceMappingReference) : "");

  m_sequenceHomologGroup = (spsFiles->m_fasta ? translateSequence(data.sequenceMappingHomolog, true) : "");
  // translate reference sequence
  m_sequenceReferenceGroup = (spsFiles->m_fasta ? translateSequence(data.sequenceMappingReference, true) : "");

  // remove extra () from sequences
  string auxString = m_sequenceHomolog;
  //stringExcessParamsFromSequence(m_sequenceDeNovo);
  //stringExcessParamsFromSequence(m_sequenceHomolog);

  DEBUG_MSG("Homolog (contig): " << auxString << " -> " << m_sequenceHomolog);


  // auxiliary variable to set reference sequence
  string auxRef, auxRef2;
  if(m_sequenceReference.compare(m_sequenceHomolog) && m_sequenceReference.length()) {
    auxRef  = m_sequenceReference;
    auxRef2 = m_sequenceReferenceGroup;
  }


  DEBUG_MSG("Reference: " << m_sequenceReference);
  DEBUG_MSG("Homolog:   " << m_sequenceHomolog);
  DEBUG_MSG("deNovo:    " << m_sequenceDeNovo);
  DEBUG_MSG("Contig reversed? " << getContigState(data.homolog));


  ////////////////////////////////////////////////////////////////////////////
  // Build mass intervals

  buildMassIntervals(data.m_referenceMasses, data.sequenceMappingReference);
  buildMassIntervals(data.m_homologMasses, data.sequenceMappingHomolog);

  // process mass intervals
  string referenceMasses, homologMasses;
  processMassIntervals(data.m_referenceMasses, referenceMasses);
  processMassIntervals(data.m_homologMasses, homologMasses);


  ////////////////////////////////////////////////////////////////////////////
  // Get contig reversed state

  bool contigReversedFlag = getContigState(data.homolog);
  string contigReversedTableEntry;
  if(contigReversedFlag)
    contigReversedTableEntry = '1';

  // clear previous consensus cluster list
  data.spectraOfContig.clear();
  // get spectra (consensus or input spectra) for this contig
  getConsensusFromContig(data.allContigsContig, data.spectraOfContig);


  ////////////////////////////////////////////////////////////////////////////
  // Build the cluster pages

  int ret;
  if(m_noClusters)
    ret = buildTableInputSpectra2(data);
  else
    ret = buildTableCluster(data);

  if(ret == ERROR)
    return ERROR;

  ////////////////////////////////////////////////////////////////////////////
  // Data needed to build this table and downstream tables

  // add concensus counter to accumulator
  data.consensusAcc += data.spectraOfContig.size();

  // build consensus list
  string consensusList;
  for(int j = 0 ; j < data.spectraOfContig.size() ; j++) {
    if(j != 0)
      consensusList += ',';
    consensusList += parseInt(data.spectraOfContig[j]);
  }


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;

  // Column 0: contig index (1-based)
  //dataRow.push_back(parseInt(m_contigs[contig] + 1));
  dataRow.push_back(data.contigText);
  // column 1: protein index (1-based)
  dataRow.push_back(data.proteinText);
  // column 2: number of spectra
  dataRow.push_back(parseInt(data.spectraOfContig.size()));
  // column 3: Reference sequence
  dataRow.push_back(auxRef);
  // column 4: Homolog sequence
  dataRow.push_back(m_sequenceHomolog);
  // column 5: DeNovo sequence
  dataRow.push_back(m_sequenceDeNovo);
  // column 6: User sequence
  dataRow.push_back("");
  // column 7: Protein name
  dataRow.push_back(data.proteinName);
  // column 8: Protein description
  dataRow.push_back(data.proteinDesc);
  // column 9: spectra (coma separated indexes)
  dataRow.push_back(consensusList);

  // data for contplot

  // column 10: Reference intervals
  dataRow.push_back(referenceMasses);
  // column 11: Homolog intervals
  dataRow.push_back(homologMasses);
  // column 12: Reference offset
  dataRow.push_back(parseFloat(data.offsetReference,2));
  // column 13: Homolog offset
  dataRow.push_back(parseFloat(data.offsetHomolog,2));
  // column 14: reverse flag
  dataRow.push_back(contigReversedTableEntry);
  // file names needed: 15 abruijn
  dataRow.push_back(m_fn_abruijn);
  // file names needed: 16 stars
  dataRow.push_back(m_fn_star);
  // file names needed: 17 sps_seqs
  dataRow.push_back(m_fn_seqs);
  // column 18: tool used
  dataRow.push_back(data.tool);
  // column 19: Grouped homolog sequence
  dataRow.push_back(m_sequenceHomologGroup);
  // column 20: Grouped reference sequence
  dataRow.push_back(auxRef2);

  // store the row
  if(tableContigs)
    tableContigs->addDataRow(dataRow);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Consensus cluster table
//////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> index in specset (1-based)
// cells[row][1] -> text   --> cluster index (used when filtered by cluster consensus)
// cells[row][2] -> text   --> contig ID
// cells[row][3] -> text   --> Reference sequence --> generateSequence()
// cells[row][4] -> text   --> Homolog sequence --> generateSequence()
// cells[row][5] -> text   --> DeNovo sequence --> generateSequence()
// cells[row][6] -> text   --> User sequence
// cells[row][7] -> text   --> mass value, from specs[i].parentMass
// cells[row][8] -> text   --> charge value, from specs[i].parentCharge
// cells[row][9] -> text   --> B%
// cells[row][10] -> text  --> Y&
// cells[row][11] -> text  --> BY Int %
// cells[row][12] -> text  --> spectra file name
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableCluster(ReportInternalData &data)
{
  //DEBUG_MSG("Adding clusters that map to contig #" << data.allContigsContig);

  // cycle thru concensus spectra
  for(int i = 0 ; i < data.spectraOfContig.size() ; i++) {

    // current consensus index
    data.consensus = data.spectraOfContig[i];

    DEBUG_MSG("---- cluster #" << data.consensus);

    // spectrum pointer
    Spectrum *spectrum = NULL;

    // if specset exists, set the spectrum
    if(spsFiles->m_consensusSpecSet)
      if(data.consensus < spsFiles->m_consensusSpecSet->size())
        spectrum = &((*spsFiles->m_consensusSpecSet)[data.consensus]);
      else
        WARN_MSG("**** cluster #" << data.consensus << " is greater than specset : " << spsFiles->m_consensusSpecSet->size());

    // specify model to be used
    //string model;
    //if(spectrum) {
    //  if(spectrum->msFragType == Spectrum::FragType_PRM && m_annotationModelPrm.length())
    //    model = m_annotationModelPrm;
    //}

    // specify mass shift to be used
    //float massShift = 0.0;
    //if(spectrum) {
    //  if(spectrum->msFragType == Spectrum::FragType_PRM)
    //    massShift = m_massShiftPrm;
    //}

    string fragmentationModel = getModelName(spectrum);

    ////////////////////////////////////////////////////////////////////////////
    // Sequences

    // sequence structure to hold reference and homolog post-processed sequences at cluster level
    SequenceMapping sequenceMappingReferenceCluster, sequenceMappingHomologCluster, sequenceMappingDeNovoCluster;

    // build reference for cluster
    //sequenceMappingReferenceCluster.clear();
    if(data.homolog != -1 && spsFiles->m_fasta)
      propagateSequence(sequenceMappingReferenceCluster, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingReference);

    // build homolog for cluster
    //sequenceMappingHomologCluster.clear();
    if(data.homolog != -1 && spsFiles->m_fasta)
      propagateSequence(sequenceMappingHomologCluster, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingHomolog);

    // build deNovo for cluster
    propagateSequence(sequenceMappingDeNovoCluster, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingDeNovo);

    // translate reference sequence
    data.clusterSequenceReference = (spsFiles->m_fasta ? translateSequence(sequenceMappingReferenceCluster) : "");
    // translate homolog sequence
    data.clusterSequenceHomolog = (spsFiles->m_fasta ? translateSequence(sequenceMappingHomologCluster) : "");
    // translate deNovo sequence
    data.clusterSequenceDeNovo = translateSequence(sequenceMappingDeNovoCluster);


    // auxiliary variable to set reference sequence
    data.clusterSequenceReferenceEfective = "";
    data.auxRefCmp = false;

    if(spsFiles->m_fasta) {
      data.auxRefCmp = data.clusterSequenceReference.compare(data.clusterSequenceHomolog) && data.clusterSequenceReference.length();
      if(data.auxRefCmp)
        data.clusterSequenceReferenceEfective = data.clusterSequenceReference;
    }

    // remove extra () from sequences
    stringExcessParamsFromSequence(data.clusterSequenceDeNovo);
    //stringExcessParamsFromSequence(data.clusterSequenceHomolog);

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    float spectrumParentMass;
    int   spectrumParentCharge;

    if(spectrum)
      spectrumParentMass = spectrum->parentMass;

    if(spectrum)
      spectrumParentCharge = spectrum->parentCharge;

    // generat statistics data;
    string seq = data.clusterSequenceDeNovo;
    if(data.homolog != -1 && spsFiles->m_fasta)
      seq = (data.auxRefCmp ? data.clusterSequenceReference : data.clusterSequenceHomolog);
    if(spectrum)
      generateStatistics(*spectrum, seq, data.cluster_B, data.cluster_Y, data.cluster_BYint);

    ////////////////////////////////////////////////////////////////////////////
    // generate the input spectra entries for this cluster
    if(buildTableInputSpectra(data) == ERROR)
      return ERROR;

    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;

    // Column 0: consensus index (1-based)
    dataRow.push_back(parseInt(data.consensus + 1));
    // Column 1: contig (1-based)
    dataRow.push_back(data.contigText);
    // Column 2: Protein (1-based)
    dataRow.push_back(data.proteinText);
    // Column 3: sequence - reference
    dataRow.push_back(data.clusterSequenceReferenceEfective);
    // Column 4: sequence - homolog
    dataRow.push_back(data.clusterSequenceHomolog);
    // Column 5: sequence - de novo
    dataRow.push_back(data.clusterSequenceDeNovo);
    // Column 6: sequence - user
    dataRow.push_back("");
    // Column 7: parent mass
    dataRow.push_back(parseFloat(spectrumParentMass, 2));
    // Column 8: charge
    dataRow.push_back(parseInt(spectrumParentCharge));
    // Column 9: B (%)
    dataRow.push_back(parseFloat(data.cluster_B, 2));
    // Column 10: Y (%)
    dataRow.push_back(parseFloat(data.cluster_Y, 2));
    // Column 11: BY Int (%)
    dataRow.push_back(parseFloat(data.cluster_BYint, 2));
    // Column 12: protein name
    dataRow.push_back(data.proteinName);
    // column 13: tool used
    dataRow.push_back(data.tool);
    // column 14: fragmentation model (string)
    dataRow.push_back(fragmentationModel);
    // column 15: model used
    //dataRow.push_back(model);
    // column 16: mass shift used
    //dataRow.push_back(parseFloat(massShift, 2));


    // store the row
    if(tableClusterConsensus)
      tableClusterConsensus->addDataRow(dataRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableInputSpectra(ReportInternalData &data)
{
  DEBUG_MSG("------ Adding input spectra that map to cluster #" << data.consensus);

  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // get input spectra list
  list<pair<unsigned,unsigned> > *inputSpectra = getInputSpectraFromConsensus(data);
  if(!inputSpectra) {
    ERROR_MSG("Got no spectra for cluster " << data.consensus);
    return OK;
  }

  //DEBUG_MSG("Got spectra for cluster: " << inputSpectra->size());

  // build iterator
  list<pair<unsigned,unsigned> >::iterator it = inputSpectra->begin();

  // cycle thru input spectra
  for( ; it != inputSpectra->end() ; it++) {

    // get spectrum
    if(it->first >= spsFiles->m_specSet.size() || it->first < 0) {
      ERROR_MSG("Cluster data file index out of bounds!!!! : " << it->first);
      return ERROR;
    }

    if(it->second >= (spsFiles->m_specSet[it->first]).size() || it->second < 0) {
      ERROR_MSG("Cluster data spectra index out of bounds!!!! : " << it->second);
      return ERROR;
    }

    Spectrum &spectrum = (spsFiles->m_specSet[it->first])[it->second];

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    // get correct sequence for generating statistics
    string  seq = ( data.homolog != -1 && spsFiles->m_fasta ? (data.auxRefCmp ? data.clusterSequenceReference : data.clusterSequenceHomolog) : data.clusterSequenceDeNovo);
    // generat statistics data;
    generateStatistics(spectrum, seq, data.spectra_B, data.spectra_Y, data.spectra_BYint);
    // model name
    string fragmentationModel = getModelName(&spectrum);

    ////////////////////////////////////////////////////////////////////////////
    // Get scan number

    // default is spectrum's scan number
    int scan = spectrum.scan;
    // get scan file for this input spectra file
    //if(spsFiles->m_scanNumberFiles.size() > it->first) {
      // get the file index
    //  const vector<vector<int> > & dataFile = spsFiles->m_specScan[it->first];
      // get the scan number
    //  if(it->second < dataFile.size())
    //    if(dataFile[it->second].size())
     //     scan = dataFile[it->second][0];
    //}

    // compose spectra file filename (make sure it's a full path)
    //string isFileName = composeFileName(m_projectDir, spsFiles->m_inputSpectraPklbin[it->first]);
    string isFileName = pathAbsoluteToRelative(m_projectDir, spsFiles->m_inputSpectraPklbin[it->first]);

    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;

    // Column 0: table index (used for updates)
    dataRow.push_back("");
    // Column 1: spectrum index in specset (1-based)
    dataRow.push_back(parseInt(it->second + 1));
    // Column 2: scan number
    dataRow.push_back(parseInt(scan));
    // column 3: cluster index (1-based)
    dataRow.push_back(parseInt(data.consensus + 1));
    // column 4: Protein name
    dataRow.push_back(data.proteinName);
    // column 5: pklbin file index
    dataRow.push_back(parseInt(it->first));
    // column 6: pklbin file name
    dataRow.push_back(isFileName);
    // column 7: Reference sequence
    dataRow.push_back(data.clusterSequenceReferenceEfective);
    // column 8: Homolog sequence
    dataRow.push_back(data.clusterSequenceHomolog);
    // column 9: DeNovo sequence
    dataRow.push_back(data.clusterSequenceDeNovo);
    // column 10: User sequence
    dataRow.push_back("");
    // column 11: mass value
    dataRow.push_back(parseFloat(spectrum.parentMass, 2));
    // column 12: charge value
    dataRow.push_back(parseInt(spectrum.parentCharge));
    // column 13: B%
    dataRow.push_back(parseFloat(data.spectra_B, 2));
    // column 14: Y&
    dataRow.push_back(parseFloat(data.spectra_Y, 2));
    // column 15: BY Int %
    dataRow.push_back(parseFloat(data.spectra_BYint, 2));
    // column 16: Original filename
    string originalFilename;
    if(spsFiles->m_input_index)
      if( (*spsFiles->m_input_index).size() > it->first)
        originalFilename = (*spsFiles->m_input_index)[it->first];
    dataRow.push_back( originalFilename );
    // column 17: contig number
    dataRow.push_back(data.contigText);
    // column 18: protein index
    dataRow.push_back(data.proteinText);
    // column 19: tool used
    dataRow.push_back(data.tool);
    // column 20: model name
    dataRow.push_back(fragmentationModel);


    // store the row
    if(tableInputSpectra)
      tableInputSpectra->addDataRow(dataRow);

    //if(spsFiles->m_input_index)
      //m_inputSpectraFiles.push_back( (*spsFiles->m_input_index)[it->first] );
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableInputSpectra2(ReportInternalData &data)
{
  DEBUG_MSG("Adding input spectra that map to contig #" << data.allContigsContig);

  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  if(!spsFiles->m_consensusSpecSet) {
    DEBUG_MSG("Combined input spectra file not loaded. Returning.");
    return ERROR;
  }

  // cycle thru input spectra
  for(int i = 0 ; i < data.spectraOfContig.size() ; i++) {

    // current spectra index
    data.consensus = data.spectraOfContig[i];

    // get spectrum index
    int currentSpectrumIndex = data.spectraOfContig[i];
    // test spectrum index
    if(currentSpectrumIndex < 0 || currentSpectrumIndex >= spsFiles->m_consensusSpecSet->size()) {
      DEBUG_MSG("Spectrum index out of bounds: [0;" << spsFiles->m_consensusSpecSet->size() << "] : " << currentSpectrumIndex);
      return ERROR;
    }
    // get spectrum
    Spectrum &spectrum = (*spsFiles->m_consensusSpecSet)[currentSpectrumIndex];

    ////////////////////////////////////////////////////////////////////////////
    // Sequences

    // sequence structure to hold reference and homolog post-processed sequences at cluster level
    SequenceMapping sequenceMappingReferenceSpectra, sequenceMappingHomologSpectra, sequenceMappingDeNovoSpectra;

    // build reference for cluster
    sequenceMappingReferenceSpectra.clear();
    if(data.homolog != -1 && spsFiles->m_fasta)
      propagateSequence(sequenceMappingReferenceSpectra, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingReference);

    // build homolog for cluster
    sequenceMappingHomologSpectra.clear();
    if(data.homolog != -1 && spsFiles->m_fasta)
      propagateSequence(sequenceMappingHomologSpectra, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingHomolog);

    // build deNovo for cluster
    propagateSequence(sequenceMappingDeNovoSpectra, data.allContigsContig, data.consensus, data.homolog, data.sequenceMappingDeNovo);

    // translate reference sequence
    data.clusterSequenceReference = (spsFiles->m_fasta ? translateSequence(sequenceMappingReferenceSpectra) : "");
    // translate homolog sequence
    data.clusterSequenceHomolog = (spsFiles->m_fasta ? translateSequence(sequenceMappingHomologSpectra) : "");
    // translate deNovo sequence
    data.clusterSequenceDeNovo = translateSequence(sequenceMappingDeNovoSpectra);


    // auxiliary variable to set reference sequence
    data.clusterSequenceReferenceEfective = "";
    data.auxRefCmp = false;

    if(spsFiles->m_fasta) {
      data.auxRefCmp = data.clusterSequenceReference.compare(data.clusterSequenceHomolog) && data.clusterSequenceReference.length();
      if(data.auxRefCmp)
        data.clusterSequenceReferenceEfective = data.clusterSequenceReference;
    }

    // remove extra () from sequences
    stringExcessParamsFromSequence(data.clusterSequenceDeNovo);
    //stringExcessParamsFromSequence(data.clusterSequenceHomolog);

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    float spectrumParentMass;
    int   spectrumParentCharge;

    spectrumParentMass = spectrum.parentMass;

    spectrumParentCharge = spectrum.parentCharge;
    // model name
    string fragmentationModel = getModelName(&spectrum);

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    // generat statistics data;
    string seq = ( data.homolog != -1 && spsFiles->m_fasta ? (data.auxRefCmp ? data.clusterSequenceReference : data.clusterSequenceHomolog) : data.clusterSequenceDeNovo);
    generateStatistics(spectrum, seq, data.spectra_B, data.spectra_Y, data.spectra_BYint);

    int fileIndex, spectraIndex;

    getinputSpectraFromCombinedSpectra(fileIndex, spectraIndex, currentSpectrumIndex);


    ////////////////////////////////////////////////////////////////////////////
    // Get scan number

    // default is spectrum's scan number
    int scan = spectrum.scan;

    // compose spectra file filename (make sure it's a full path)
    //string isFileName = ( fileIndex == -1 ? "" : composeFileName(m_projectDir, spsFiles->m_inputSpectraPklbin[fileIndex]));
    string isFileName;
    if(spsFiles)
      if(fileIndex >= 0 && fileIndex <  spsFiles->m_inputSpectraPklbin.size())
        isFileName = pathAbsoluteToRelative(m_projectDir, spsFiles->m_inputSpectraPklbin[fileIndex]);

    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;

    // Column 0: table index (used for updates)
    dataRow.push_back("");
    // Column 1: spectrum index in specset (1-based)
    dataRow.push_back(parseInt(spectraIndex + 1));
    // Column 2: scan number
    dataRow.push_back(parseInt(scan));
    // column 3: cluster index (1-based)
    dataRow.push_back(parseInt(-1));
    // column 4: Protein name
    dataRow.push_back(data.proteinName);
    // column 5: pklbin file index
    dataRow.push_back(parseInt(fileIndex));
    // column 6: pklbin file name
    dataRow.push_back(isFileName);
    // column 7: Reference sequence
    dataRow.push_back(data.clusterSequenceReferenceEfective);
    // column 8: Homolog sequence
    dataRow.push_back(data.clusterSequenceHomolog);
    // column 9: DeNovo sequence
    dataRow.push_back(data.clusterSequenceDeNovo);
    // column 10: User sequence
    dataRow.push_back("");
    // column 11: mass value
    dataRow.push_back(parseFloat(spectrum.parentMass, 2));
    // column 12: charge value
    dataRow.push_back(parseInt(spectrum.parentCharge));
    // column 13: B%
    dataRow.push_back(parseFloat(data.spectra_B, 2));
    // column 14: Y&
    dataRow.push_back(parseFloat(data.spectra_Y, 2));
    // column 15: BY Int %
    dataRow.push_back(parseFloat(data.spectra_BYint, 2));
    // column 16: Original filename
    string originalFilename;
    if(spsFiles->m_input_index)
      if(fileIndex >= 0 && fileIndex < spsFiles->m_input_index->size())
        originalFilename = (*spsFiles->m_input_index)[fileIndex];
    dataRow.push_back( originalFilename );
    // column 17: contig number
    dataRow.push_back(data.contigText);
    // column 18: protein index
    dataRow.push_back(data.proteinText);
    // column 19: tool used
    dataRow.push_back(data.tool);
    // column 20: model name
    dataRow.push_back(fragmentationModel);


    // store the row
    if(tableInputSpectra)
      tableInputSpectra->addDataRow(dataRow);

    //if(spsFiles->m_input_index)
    //  m_inputSpectraFiles.push_back( (*spsFiles->m_input_index)[fileIndex] );
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getinputSpectraFromCombinedSpectra(int &fileIndex, int &spectraIndex, int combinedSpectraIndex)
{
  //cout << "combinedSpectraIndex : " << combinedSpectraIndex << endl;
 
  spectraIndex = -1;
  fileIndex = -1;

  if(spsFiles->m_consensusSpecSet == NULL) return;
    
  if(combinedSpectraIndex < 0 || combinedSpectraIndex >= spsFiles->m_consensusSpecSet->size())
    return;
    
  fileIndex = (*spsFiles->m_consensusSpecSet)[combinedSpectraIndex].fileIndex;
  int scan = (*spsFiles->m_consensusSpecSet)[combinedSpectraIndex].scan;

  if(fileIndex < 0 || fileIndex >= spsFiles->m_specSet.size())
    return;

  for(int i = 0 ; i < spsFiles->m_specSet[fileIndex].size() ; i++) {
    if(spsFiles->m_specSet[fileIndex][i].scan == scan) {
      spectraIndex = i;
      return;
    }
  }

  return;  

/*
  // check for input mapping file
  if(!spsFiles->m_inputMapping)
    return;

  // check for out-of-bounds spectra indexes
  if(combinedSpectraIndex < 0 || combinedSpectraIndex >= (*spsFiles->m_inputMapping).size()) {
    WARN_MSG("combinedSpectraIndex out-of-bounds! : " << combinedSpectraIndex << " ; [" << 0 << ";" << (*spsFiles->m_inputMapping).size() << "]");
    return;
  }

  fileIndex = (*spsFiles->m_inputMapping)[combinedSpectraIndex][0];
  spectraIndex = (*spsFiles->m_inputMapping)[combinedSpectraIndex][1];

  //cout << "fileIndex            : " << fileIndex << endl;
  //cout << "spectraIndex         : " << spectraIndex << endl; 
  */
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::buildMassIntervals(vector<float> &masses, SequenceMapping &sm)
{
  if(spsFiles->m_fasta)
    masses.push_back(sm.startMass);
    for(int j = 0 ; j < sm.processedAA.size() ; j++)
      masses.push_back(sm.processedAA[j].massPoint);
}
////////////////////////////////////////////////////////////////////////////////
// Dump the abruijn graph to the screen
list<pair<unsigned,unsigned> > *ReportTableGenerator::getInputSpectraFromConsensus(ReportInternalData &data)
{
  DEBUG_MSG("Getting input spectra list for consensus " << data.consensus);

  // get input spectra list
  list<pair<unsigned,unsigned> > *inputSpectra = getInputSpectraFromConsensus(data.consensus);

  // set default tool as SPS
  data.tool = TOOL_SPS;

  // check of input spectra file.
  if( (!inputSpectra) && (m_tool == 3) ) {

    // get total number of clusters
    int n = m_clusterCount; // /2;

    DEBUG_MSG("not found. Looking for spectra mapped to consensus " << data.consensus - n);

    // find spectra for consensus 'found' - total
    inputSpectra = getInputSpectraFromConsensus(data.consensus - n);
    // if not found, somethinf is wrong
    if(!inputSpectra) {
      DEBUG_MSG("+++ not found for consensus " << data.consensus - n);
      return NULL;
    }

    DEBUG_MSG("--- found for " << data.consensus - n);

    switch(m_tool) {
    case 1:
      return NULL;
      break;
    case 2:
      data.tool = TOOL_GENOMS;
      break;
    case 3:
      data.tool = TOOL_GENOMS;
      break;
    default:
      return NULL;
    }

  } else {

    if(m_tool == 2)
      data.tool = TOOL_GENOMS;
  }
  return inputSpectra;
}
////////////////////////////////////////////////////////////////////////////////
// Data mapping methods. Return
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::processMassIntervals(vector<float> &mass, string &intervals)
{
  float first = -1.0, second;

  for(int j = 0 ; j < mass.size() ; j++) {

    second = mass[j];

    if(first > second)
      break;

    if(j)
      intervals += TABLE_SEP_L1;

    intervals += parseFloat(second, 1);

    first = second;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Data mapping methods. Return
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getConsensusFromInputSpectra(int fileIndex, int inputSpectra)
{
  if(!spsFiles->m_clusterData) return -1;

  return spsFiles->m_clusterData->getConsensusFromInputSpectra(fileIndex, inputSpectra);
}
////////////////////////////////////////////////////////////////////////////////
list<pair<unsigned,unsigned> > *ReportTableGenerator::getInputSpectraFromConsensus(int consensus)
{
  if(!spsFiles->m_clusterData) //return NULL;
  {
    list<pair<unsigned,unsigned> > *ret = new list<pair<unsigned,unsigned> >();
    
    if(!spsFiles->m_consensusSpecSet) return NULL;
    
    if(consensus < 0 || consensus >= spsFiles->m_consensusSpecSet->size()) return NULL;
    
    // get file index and scan
    int fileIndex = (*spsFiles->m_consensusSpecSet)[consensus].fileIndex;
    int scan      = (*spsFiles->m_consensusSpecSet)[consensus].scan;
    int specIndex = -1;
    
    // get spectrum index
    if(fileIndex < 0 || fileIndex >= spsFiles->m_specSet.size()) return NULL;
    
    for(int i = 0 ; i < spsFiles->m_specSet[fileIndex].size() ; i++)
      if( spsFiles->m_specSet[fileIndex][i].scan == scan) {
        specIndex = i;
        break;
      }
    
      if(specIndex != -1) {
        ret->push_back(make_pair(fileIndex, specIndex));
        return ret;
      }
      return NULL;
  }

  return spsFiles->m_clusterData->getInputSpectraFromConsensus(consensus);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getClusterCount(void)
{
  if(!spsFiles->m_clusterData) return -1;

  return spsFiles->m_clusterData->getClusterCount();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSpectraCount(void)
{
  if(!spsFiles->m_clusterData) return -1;

  return spsFiles->m_clusterData->getSpectraCount();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromConsensus(int consensus)
{
  if(!spsFiles->m_abruijn) return -1;
  if(consensus < 0) return -1;

  abinfo_t::iterator i0 = spsFiles->m_abruijn->begin();
  for( ; i0 != spsFiles->m_abruijn->end() ; i0++)
    // output first vector in pair
    for(int i = 0 ; i < i0->second.first.first.size() ; i++)
      if(i0->second.first.first[i] == consensus) {
        return i0->first;
      }

  return -1;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getConsensusFromContig(int contig, vector<int> &ret)
{
  if(!spsFiles->m_abruijn) return;

  if(contig < 0)
    return;

  abinfo_t::iterator i0 = spsFiles->m_abruijn->find(contig);
  if(i0 == spsFiles->m_abruijn->end())
    return;

  for(int i = 0 ; i < i0->second.first.first.size() ; i++)
    ret.push_back(i0->second.first.first[i]);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromInputSpectra(int fileIndex, int inputSpectra)
{
  // get consensus
  int consensus = getConsensusFromInputSpectra(fileIndex, inputSpectra);
  // get contig
  return getContigFromConsensus(consensus);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getAllContigFromContig(int contig)
{
  if(!spsFiles->m_contigIndices) return contig; // -1;

  if(contig < 0 || contig >= spsFiles->m_contigIndices->size())
    return -1;

  return (*spsFiles->m_contigIndices)[contig][0];
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getAllContigFromContig(vector<int> &contigs, vector<int> &allContigs)
{
  if(!spsFiles->m_contigIndices) {
    allContigs = contigs;
    return;
  }

  for(int i = 0 ; i < contigs.size() ; i++) {
    if((*spsFiles->m_contigIndices)[contigs[i]][0] >= 0)
      allContigs.push_back((*spsFiles->m_contigIndices)[contigs[i]][0]);
    }
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromAllContig(int contig)
{
  if(!spsFiles->m_contigIndices) return -2;

  for(int i = 0 ; i < spsFiles->m_contigIndices->size() ; i++)
    if((*spsFiles->m_contigIndices)[i][0] == contig)
      return i;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinFromHomolog(int homolog)
{
  if(!spsFiles->m_contigs_mp) return -1;
  // determine protein
  int protein = -1;
  if((homolog >= 0) && (homolog < spsFiles->m_contigs_mp->size()))
    protein = (*spsFiles->m_contigs_mp)[homolog][0];
  return protein;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinFromReference(int reference)
{
  if(!spsFiles->m_homglue_ref_mp) return -1;
  // determine protein
  int protein = -1;
  if((reference >= 0) && (reference < spsFiles->m_homglue_ref_mp->size()))
    protein = (*spsFiles->m_homglue_ref_mp)[reference][0];
  return protein;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getHomologFromProtein(int protein, vector<int> &ret)
{
  // spsFiles->m_homglue_ref_mp --> spsFiles->m_contigs_mp
  if(!spsFiles->m_homglue_ref_mp) return;
  // exit empty if protein index invalid
  if(protein < 0) return;

  for(int i = 0 ; i < spsFiles->m_homglue_ref_mp->size() ; i++)  // spsFiles->m_contigs_mp
    if((*spsFiles->m_homglue_ref_mp)[i][0] == protein)
      ret.push_back(i);
}
////////////////////////////////////////////////////////////////////////////////
// Sequence generation methods
////////////////////////////////////////////////////////////////////////////////
// Step 0: get protein sequence into a data structure
int ReportTableGenerator::generateSequenceStep0(SequenceMapping &contigInfo, int proteinIndex)
{
  if(!spsFiles->m_fasta) return ERROR;

  string sequence = spsFiles->m_fasta->getSequence(proteinIndex);
  // cycle thru the sequence
  for(int j = 0 ; j < sequence.size() ; j++) {
    aaCell cell;
    // the aa sequence
    cell.aa = sequence[j];
    // just one
    cell.colspan = 0;
    // its position
    cell.startPosition = j;
    contigInfo.processedAA.push_back(cell);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 1: load MatchesIndex and contigSpectra into the data structure
int ReportTableGenerator::generateSequenceStep1(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, Spectrum &contigSpectra, ContigMatchData &contigInfo, int contigIdx)
{
  // check boundaries
  if(contigIdx < 0 || contigIdx >= contigMatches.size())
    return ERROR;

  // test if there is any data to store
  if(contigMatchesIndex[contigIdx].size() <= 0)
    return OK;

  // if this is the first sequence part, skip empty cells and state sequence beginning
  for(int j = 0 ; j < contigMatchesIndex[contigIdx].size() ; j++) {
    PairContigProteinMassIdx pairItem;
    pairItem.contigMassIdx = (int)contigMatchesIndex[contigIdx][j][0];
    pairItem.proteinMassIdx = (int)contigMatchesIndex[contigIdx][j][1];
    contigInfo.pair.push_back(pairItem);
  }

	// 3rd file: Contig spectra
  if(contigSpectra.size() > 0) {
    DEBUG_MSG("contig: " << contigIdx << "  size: " << contigSpectra.size());
    stringstream str;
    // if this is the first sequence part, skip empty cells and state sequence beginning
 		for(int j = 0 ; j < contigSpectra.size() ; j++) {
      float mass = contigSpectra[j][0];
      contigInfo.contigMass.push_back(mass);
      str << mass << " ; ";
    }
    DEBUG_MSG("masses: " << str.str() );
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 2: Calculate alignment and generate sequence
int ReportTableGenerator::generateSequenceStep2(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source, SequenceMapping &m_sequenceMapping)
{
  // get contig span and location
  int size = source.pair.size();
  if(size) {
    m_sequenceMapping.startPosition = source.pair[0].contigMassIdx;
    m_sequenceMapping.endPosition = source.pair[size-1].contigMassIdx;
  }

/*
  // we may need to add the leading mass
  float leadMass = source.contigMass[0];
  if(fabs(leadMass) > m_peakMassTol ) {
    // cell to hold one data item
    aaCell cell;
    // store mass value
    cell.delta = leadMass;
    // add the cell to the list
    m_sequenceMapping.processedAA.push_back(cell);
  } */

    DEBUG_MSG("**********************");
    DEBUG_MSG("contig: " << contigIdx);
    DEBUG_MSG("size: " << size);

  // Loop thru contig size to generate info cells
  for(int j = 0 ; j < size-1 ; j++) {

    DEBUG_MSG("----------------------");
    DEBUG_MSG("index: " << j);

    // get the corresponding contig indexes
    int contigIndex = source.pair[j].contigMassIdx;
    int contigIndexAfter = source.pair[j+1].contigMassIdx;
    // get the correspondif protein indexes
    int proteinIndex = source.pair[j].proteinMassIdx;
    int proteinIndexAfter = source.pair[j+1].proteinMassIdx;

    // test boundaries
    if((contigIndex >= source.contigMass.size()) || (contigIndexAfter >= source.contigMass.size()))
      break;
    // Get the mass at this point
    float massCsps = source.contigMass[contigIndex];
    // and at the point after
    float massAfter = source.contigMass[contigIndexAfter];

    DEBUG_MSG("protein index Idx: " << proteinIndex);
    DEBUG_MSG("protein index After:  " << proteinIndexAfter);
    DEBUG_MSG("contig index Idx: " << contigIndex);
    DEBUG_MSG("contig index After:  " << contigIndexAfter);
    DEBUG_MSG("contig mass before: " << massCsps);
    DEBUG_MSG("contig mass after:  " << massAfter);


    string sequence;

    // get AA sequence for interval
    //if( (contigIndex >= proteinData.startPosition ) && (contigIndexAfter < proteinData.endPosition) )
      for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.processedAA.size()) ; k++)
        sequence += proteinData.processedAA[k].aa;

    // get masses for sequence
    vector<float> masses;
    getMasses((char*)sequence.c_str(), masses);
    float mass = 0.0;
    for(int k = 0 ; k < masses.size() ; k++)
      mass += masses[k];

    // massAux used to calculate mass difference
    float massAux = (massAfter - massCsps) - mass;

    DEBUG_MSG("mass sum:  " << mass);
    DEBUG_MSG("mass diff: " << massAux);
    DEBUG_MSG("sequence mass:  " << mass);
    DEBUG_MSG("sequence:  " << sequence);

    // cell to hold one data item
    aaCell cell;

    // save it
    cell.aa = sequence; //result;
    // store mass difference
    cell.delta = massAux;
    // calculate the colspan value for this cell
    cell.colspan = source.pair[j+1].proteinMassIdx - source.pair[j].proteinMassIdx - 1;
    // and the start position
    cell.startPosition = source.pair[j].contigMassIdx;
    // mass point at cell end
    cell.massPoint = massAfter;

    if(!j)
      m_sequenceMapping.startMass = massCsps;
    // add the cell to the list
    m_sequenceMapping.processedAA.push_back(cell);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 2: Calculate protein sequence coverage
int ReportTableGenerator::generateSequenceStep2B(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source)
{
  // get contig span and location
  int size = source.pair.size();
  int minimunIndex = source.pair[0].proteinMassIdx;
  int maximunIndex = source.pair[size-1].proteinMassIdx;

  // Loop thru contig size to generate info cells
  for(int j = 0 ; j < size-1 ; j++) {

    // get the corresponding contig indexes
    int contigIndex = source.pair[j].contigMassIdx;
    int contigIndexAfter = source.pair[j+1].contigMassIdx;
    // get the correspondif protein indexes
    int proteinIndex = source.pair[j].proteinMassIdx;
    int proteinIndexAfter = source.pair[j+1].proteinMassIdx;
    // Get the mass at this point
    float massCsps = source.contigMass[contigIndex];
    // and at the point after
    float massAfter = source.contigMass[contigIndexAfter];

    string sequence;

    // mark the protein's aa's covered by the csps contig
    //if( (contigIndex >= proteinData.startPosition ) && (contigIndexAfter < proteinData.endPosition) )
      for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.processedAA.size()) ; k++)
        proteinData.processedAA[k].covered = true;
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getModelName(Spectrum *spectrum)
{
  string fragmentationModel;
  switch(spectrum->msFragType) {
  case Spectrum::FragType_PRM:
    fragmentationModel = "PRM";
    break;
  case Spectrum::FragType_CID:
    fragmentationModel = "CID";
    break;
  case Spectrum::FragType_ETD:
    fragmentationModel = "ETD";
    break;
  case Spectrum::FragType_HCD:
    fragmentationModel = "HCD";
    break;
  default:
    fragmentationModel = " -- ";
    break;
  }
  return fragmentationModel;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getProteinName(int protein)
{
  string result;

  if(!spsFiles->m_fasta) return result;

  if((protein < 0) || (protein >= spsFiles->m_fasta->IDs.size())) {
    DEBUG_MSG("Invalid protein index: " << protein);
  } else {
    result = spsFiles->m_fasta->IDs[protein];
    cleanProteinName(result);
  }
  return result;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getProteinDescription(int protein)
{
  string result;

  if(!spsFiles->m_fasta) return result;

  if((protein < 0) || (protein >= spsFiles->m_fasta->IDs.size())) {
    DEBUG_MSG("Invalid protein index: " << protein);
  } else {
    result = spsFiles->m_fasta->getDesc(protein);
    cleanProteinName(result);
  }
  return result;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getProteinSequence(int protein)
{
  string result;

  if(!spsFiles->m_fasta) return result;

  if((protein < 0) || (protein >= spsFiles->m_fasta->IDs.size())) {
    DEBUG_MSG("Invalid protein index: " << protein);
  } else {
    result = spsFiles->m_fasta->getSequence(protein);
  }
  return result;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getSequenceDenovo(Spectrum &s)
{
  //s.output(cout);
  AAJumps masses(1, m_resolution, m_peakMassTol, AAJumps::NO_MODS, true);
  //AAJumps masses(1, m_resolution, m_peakMassTol, AAJumps::NO_MODS, false);
  //AAJumps masses(1);
  string sequence;
  masses.getPeptideFromSpectrum(s,sequence);
  //cout << sequence << endl;
  return sequence;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::stringExcessParamsFromSequence(string &iSequence)
{
  string ret;
  string aa("ARNDCEQGHILKMFPSTWYV");
  bool found = false;

  // store the position after each tag
  size_t lastPosition = 0;
  int index = 0;

  while(lastPosition < iSequence.length()) {

    string current;
    // in case of an AA
    if(aa.find_first_of(iSequence[lastPosition]) != string::npos) {

      current = iSequence[lastPosition];
      // in case of mass
      lastPosition++;

    } else if(iSequence[lastPosition] == '[') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find end
      while(iSequence[currentPosition] != ']' && currentPosition < iSequence.size())
        currentPosition++;
      // get contents
      current = iSequence.substr(lastPosition, currentPosition - lastPosition + 1);
      // update last position for the next step
      lastPosition = currentPosition + 1;

    } else if(iSequence[lastPosition] == '(') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find comma
      while(iSequence[currentPosition] != ',' && iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
        currentPosition++;

      // found comma, so everything is included
      if(iSequence[currentPosition] == ',') {
        // find end
        while(iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
          currentPosition++;
        // get sub sequence contents
        current = iSequence.substr(lastPosition, currentPosition - lastPosition + 1);
        // update last position
        lastPosition = currentPosition + 1;

      // found ), so the () are eliminated
      } else {
        // get contents
        current = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
        lastPosition = currentPosition + 1;
        found = true;
      }

    } else { // ERROR
    }

    // add the current item to the sequence
    ret += current;
  }

  iSequence = ret;

  return found;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::translateSequence(SequenceMapping &sequenceMapping, bool group)
{
  // return value holder
  string sequence;

  // cycle thru all sequence items
  for(int i = 0 ; i < sequenceMapping.processedAA.size() ; i++) {

    // test delta over peak mass tolerance
    if(fabs(sequenceMapping.processedAA[i].delta) > m_peakMassTol) {
      if(sequenceMapping.processedAA[i].aa.length()) {
        sequence += "(";
        sequence += sequenceMapping.processedAA[i].aa;
        sequence += ",";
        sequence += parseFloat(sequenceMapping.processedAA[i].delta, 2);
        sequence += ")";
      } else {
        sequence += "[";
        sequence += parseFloat(sequenceMapping.processedAA[i].delta, 2);
        sequence += "]";
      }
    } else {
      if((sequenceMapping.processedAA[i].aa.length() > 1) && group)
        sequence += "(";
      sequence += sequenceMapping.processedAA[i].aa;
      if((sequenceMapping.processedAA[i].aa.length() > 1) && group)
        sequence += ")";
    }
  }

  // return the sequence
  return sequence;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::translateSequenceReverse(SequenceMapping &oSequenceMapping, string &iSequence)
{
  string aa("ARNDCEQGHILKMFPSTWYV");
  string aux;

  oSequenceMapping.clear();

  // store the position after each tag
  size_t lastPosition = 0;
  int index = 0;

  while(lastPosition < iSequence.size()) {

    // store the processed sequence
    aaCell cell;
    cell.startPosition = index++;

    // in case of an AA
    if(aa.find_first_of(iSequence[lastPosition]) != string::npos) {

      cell.aa = iSequence[lastPosition];

    // in case of mass
      lastPosition++;

    } else if(iSequence[lastPosition] == '[') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find end
      while(iSequence[currentPosition] != ']' && currentPosition < iSequence.size())
        currentPosition++;
      // get contents
      aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
      // translate to a float
      cell.delta = getFloat(aux.c_str());
      lastPosition = currentPosition + 1;


    } else if(iSequence[lastPosition] == '(') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find comma
      while(iSequence[currentPosition] != ',' && iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
        currentPosition++;
      // get sub sequence contents
      aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);

      lastPosition = currentPosition + 1;

      if(iSequence[currentPosition] == ',') {
        // find end
        while(iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
          currentPosition++;
        // get contents
        aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
        // translate to a float
        cell.delta = getFloat(aux.c_str());

        lastPosition = currentPosition + 1;
      }

    } else { // ERROR
    }

    // add the cell to the list
    oSequenceMapping.processedAA.push_back(cell);
  }

  oSequenceMapping.startPosition = 0;
  oSequenceMapping.endPosition = oSequenceMapping.processedAA.size() - 1;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceReference(ReportInternalData &data)
{
  //DEBUG_MSG("Getting reference sequence for " << data.homolog);
  // clear sequence mapping between protein and contig
  data.sequenceMappingReference.clear();
  // clear reference masses
  data.m_referenceMasses.clear();
  // test for needed files
  if( (!spsFiles->m_homglue_ref_midx) || (!spsFiles->m_homglue_ref_mp) || (!spsFiles->m_contigsSpectra) ) return -1;
  // text for index
  //if((data.allContigsContig < 0) || (data.allContigsContig >= spsFiles->m_sps_seqs->size())) {
  if((data.homolog < 0) || (data.homolog >= spsFiles->m_contigsSpectra->size())) {
    DEBUG_MSG("Attempt to access invalid index in contigs.pklbin:  " << data.homolog);
    return -1;
  }

  // reverse spectrum, if needed
  //Spectrum sp = (*spsFiles->m_sps_seqs)[data.allContigsContig];
  //if(getContigState(data.homolog))
  //  sp.reverse(0.0 - AAJumps::massH2O);
  Spectrum sp = (*spsFiles->m_contigsSpectra)[data.homolog];

  // generate the sequence structure
  return getSequence(data.homolog, *spsFiles->m_homglue_ref_midx, *spsFiles->m_homglue_ref_mp, sp, data.sequenceMappingReference, data.offsetReference);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceHomolog(ReportInternalData &data)
{
  //DEBUG_MSG("Getting homolog sequence for " << data.homolog);
  // clear sequence mapping between protein and contig
  data.sequenceMappingHomolog.clear();
  // clear homolog masses
  data.m_homologMasses.clear();
  // test for needed files
  if( (!spsFiles->m_contigs_midx) || (!spsFiles->m_contigs_mp) || (!spsFiles->m_contigsSpectra) ) return -1;
  // text for index
  //if((data.allContigsContig < 0) || (data.allContigsContig >= spsFiles->m_sps_seqs->size())) {
  if((data.homolog < 0) || (data.homolog >= spsFiles->m_contigsSpectra->size())) {
    DEBUG_MSG("Attempt to access invalid index in contigs.pklbin:  " << data.homolog);
    return -1;
  }
  // reverse spectrum, if needed
  //Spectrum sp = (*spsFiles->m_sps_seqs)[data.allContigsContig];
  //if(getContigState(data.homolog)) {
  //  sp.reverse(0.0 - AAJumps::massH2O);
  //  DEBUG_MSG("Reversing contig spectra (sequences)");
  //}
  Spectrum sp = (*spsFiles->m_contigsSpectra)[data.homolog];

  // generate the sequence structure
  return getSequence(data.homolog, *spsFiles->m_contigs_midx, *spsFiles->m_contigs_mp, sp, data.sequenceMappingHomolog, data.offsetHomolog);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinCoverage(int protein,
    SpecSet &contigMatchesIndex,
    std::vector<vector<int> > &contigMatches,
    SpecSet &contigSpectra,
    string &sequence, int &proteinSize, int &covered)
{
  // Protein sequence
  SequenceMapping proteinData;

  int ret;
  // Get the protein sequence
  ret = generateSequenceStep0(proteinData, protein);

  // get contigs associated with the protein
  vector<int> contigsIdx;
  for(int i = 0 ; i < contigMatches.size() ; i++)
    if(contigMatches[i][0] == protein)
      contigsIdx.push_back(i);

  // cycle thru all contigs
  for(int i = 0 ; i < contigsIdx.size() ; i++) {
    // Structure to hold contig info -- intermidiate step
    ContigMatchData  contigInfo;
    // Middle step
    ret = generateSequenceStep1(contigMatchesIndex, contigMatches, contigSpectra[contigsIdx[i]], contigInfo, contigsIdx[i]);
    // calculate sequence and
    ret = generateSequenceStep2B(contigsIdx[i], proteinData, contigInfo);
  }

  // build sequence for the protein

  // initial state is sequence unmaked
  bool state = false;
  // hold covered AAs count for percentage calculation
  covered = 0;

  // cycle thru all AAs
  for(int i = 0 ; i < proteinData.processedAA.size() ; i++) {

    if(proteinData.processedAA[i].covered)
      covered++;

    // if there is a state change, introduce the mark
    if(proteinData.processedAA[i].covered != state)
      sequence += TABLE_SEP_L1;

    // add the AA
    sequence += proteinData.processedAA[i].aa;

    // set the state as the current state
    state = proteinData.processedAA[i].covered;
  }

  // set the protein size
  proteinSize = proteinData.processedAA.size();

  return ret;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequence(int contigIdx,
    SpecSet &contigMatchesIndex,
    std::vector<vector<int> > &contigMatches,
    Spectrum &contigSpectra,
    SequenceMapping &sequenceMapping,
    double &leadMass)
{
  // Structure to hold contig info -- intermidiate step
  ContigMatchData  contigInfo;
  // Protein sequence
  SequenceMapping proteinData;

  // determine protein
  int protein = -1;
  if(contigIdx < contigMatches.size())
    protein = contigMatches[contigIdx][0];

  if(protein < 0)
    return OK;

  int ret;
  // Get the protein sequence
  ret = generateSequenceStep0(proteinData, protein);
  // Middle step
  ret = generateSequenceStep1(contigMatchesIndex, contigMatches, contigSpectra, contigInfo, contigIdx);
  // calculate sequence and
  ret = generateSequenceStep2(contigIdx, proteinData, contigInfo, sequenceMapping);

  stringstream aux;
  aux << "Protein masses: ";
  for(int i = 0 ; i < contigInfo.pair.size() ; i++)
    aux << contigInfo.pair[i].proteinMassIdx << " ; ";
  aux << endl;

  aux << "Masses: ";
  for(int i = 0 ; i < contigInfo.contigMass.size() ; i++)
    aux << contigInfo.contigMass[i] << " ; ";

  DEBUG_MSG(aux.str());

  // get mass offset
  leadMass = 0.0;
  if(contigInfo.contigMass.size())
    leadMass = contigInfo.contigMass[0];

  DEBUG_VAR(ret);

  return ret;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getContigState(unsigned contig)
{
  // test for needed files
  if(!spsFiles->m_contigs_mp) return false;

  if( contig < spsFiles->m_contigs_mp->size() && contig >= 0)
    return(*spsFiles->m_contigs_mp)[contig][2];
  return false;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getCspsContigState(unsigned contig)
{
  // test for needed files
  if(!spsFiles->m_homglue_matches_mp) return false;

  if( contig < spsFiles->m_homglue_matches_mp->size() && contig >= 0)
    return(*spsFiles->m_homglue_matches_mp)[contig][2];
  return false;
}////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::checkContig(int contig)
{
  // test for needed files
  if(!spsFiles->m_abruijn) return false;

  // contig in abruijn
  abContig_t &ab = (*spsFiles->m_abruijn)[contig];
  // contig data in abruijn
  abContigData_t &cd = ab.second;

  // check abruijn with zero nodes
  if(cd.size() == 0)
    return false;

  // reverse the spectra - travel thru the contig nodes
  for(int i =  cd.size() - 1 ; i >= 0 ; i--) {
    // cycle thru the data in the nodes
    if(!cd[i].first.size())
      return false;
  }

  // check zero length sps_seqs
  Spectrum &sp = (*spsFiles->m_sps_seqs)[contig];
  if(sp.size() == 0)
    return false;

  return true;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getAbruijnStarState(int contig, int star)
{
  // test for needed files
  if(!spsFiles->m_abruijn) return false;

  // contig in abruijn
  abContig_t &ab = (*spsFiles->m_abruijn)[contig];
  // contig data in abruijn
  abContigState_t &cd = ab.first;

  for(int i = 0 ; i < cd.first.size() ; i++) {
    if(cd.first[i] == star)
      return cd.second[i];
  }

  return false;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getReversedAbruijn(int contig, abContigData_t &nodes)
{
  // test for needed files
  if(!spsFiles->m_abruijn) return;

  // contig in abruijn
  abContig_t &ab = (*spsFiles->m_abruijn)[contig];
  // contig data in abruijn
  abContigData_t &cd = ab.second;

  // reverse the spectra - travel thru the contig nodes
  for(int i =  cd.size() - 1 ; i >= 0 ; i--) {
    // star IDs for node
    vector<int> IDs;
    // data for node
    vector<double> data;
    // cycle thru the data in the nodes
    for(int j = 0 ; j < cd[i].first.size() ; j++) {
      int id = cd[i].first[j];
      IDs.push_back(id);
      double parentMass = (*spsFiles->m_starSpectra)[id].parentMass;
      data.push_back(parentMass - AAJumps::massHion - cd[i].second[j]);
    }
    nodes.push_back(make_pair<vector<int>, vector<double> >(IDs, data));
  }

}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getStarIndex(vector<pair<int, double> > &data, int star, double &value)
{
  for(int i = 0 ; i < data.size() ; i++)
    if(data[i].first == star) {
      value = data[i].second;
      return true;
    }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceIndex(SequenceMapping &data, int index)
{
  for(int i = 0 ; i < data.processedAA.size() ; i++)
    if(data.processedAA[i].startPosition == index)
      return i;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getSequenceBetween(aaCell &ret, SequenceMapping &data, int start, int end)
{
  ret.aa = "";
  ret.delta = 0.0;

  for(int i = start ; i < end && i < data.processedAA.size() ; i++) {
    ret.aa += data.processedAA[i].aa;
    ret.delta += data.processedAA[i].delta;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::propagateSequence(SequenceMapping &oSequenceMapping, int contig, int star, int homolog, SequenceMapping &iSequenceMapping)
{
  // test for needed files
  if( (!spsFiles->m_abruijn) || (!spsFiles->m_starSpectra) ) return;

  // contig in abruijn
  abContig_t &ab = (*spsFiles->m_abruijn)[contig];
  // contig data in abruijn
  abContigData_t *cd, nodes;
  // data holder
  DiffData diffData;

  // get abruijn nodes for contig in the right direction
  // check reversal state
  if( (homolog != -1) && getContigState(homolog)) {
    // of reversed, reverse the contig
    getReversedAbruijn(contig, nodes);
    // and use that one
    cd = &nodes;
  } else {
    // otherwise, use the original one
    cd = &(ab.second);
  }

  // gather start spectra data from abruijn
  for(int i = 0 ; i < cd->size() ; i++)
    for(int j = 0 ; j < (*cd)[i].first.size() ; j++)
      if((*cd)[i].first[j] == star)
        diffData.abruijnData.push_back(make_pair<int,double>(i, (*cd)[i].second[j]));

  // gather contig data from contig spectra
  Spectrum sp = (*spsFiles->m_sps_seqs)[contig];
  if(getContigState(homolog))
    sp.reverse(0.0 - AAJumps::massH2O);
  for(int i = 0 ; i < sp.size() ; i++)
    diffData.contigData.push_back(sp[i][0]);

  // clear sequence container
  oSequenceMapping.clear();

  int proteinItemAtContigIndex = 0;

  int lastAaItem = iSequenceMapping.processedAA.size() +  iSequenceMapping.startPosition;
  // get first item
  //int start =  diffData.abruijnData[0].first;
  // variables used for data differences
  double s0 = 0.0, c0 = 0.0, c1, s1, sd;
  int proteinItemAtLastContigIndex = -1;
  bool firstItem = true;
  bool mark = false;

  //cout << contig << " ; " << homolog << " --> ";
  //for(int i = 0 ; i < iSequenceMapping.processedAA.size() ; i++) {
  //  cout << '(';
  //  cout << iSequenceMapping.processedAA[i].aa;
  //  cout << '|' << iSequenceMapping.processedAA[i].startPosition << ')';
  //}
  //cout << " ; ";

  // cycle thru all star data (masses)
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++) {
    // get contig node index
    int i0 = diffData.abruijnData[i].first;
    // get star data (mass)
    s1 = diffData.abruijnData[i].second;
    // get contig data (mass)
    if(i0 < diffData.contigData.size())
      c1 = diffData.contigData[i0];
    else {
      ERROR_MSG("DATA INCONSISTENCY: Detected data inconsistensy between the abruijn graph and sps_seqs!!!");
      return;
    }
    // calc difference of differences
    sd = (s1-s0) - (c1-c0);
    // declare cell object
    aaCell cell;
    // make sure delta value is 0
    cell.delta = 0.0;
    // in case the current star node is:
    proteinItemAtContigIndex = getSequenceIndex(iSequenceMapping, i0);

    // if no sequence item was found for this contig item
    if((proteinItemAtContigIndex == -1) && (proteinItemAtLastContigIndex == -1)) {
      //cout << "0(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      // ignore contig and star peaks with no AA correspondence
      continue;

    // the first:
    } else if(firstItem) { //if(i0 <= start) {
      //cout << "1(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      // just store the mass. This is the mass prefix
      // the aa string goes empty, the mass value is the interval mass
      cell.delta = s1;
      // say prefix item already done
      firstItem = false;

    // in a hole
    } else if(proteinItemAtContigIndex == -1 && proteinItemAtLastContigIndex != -1 && !mark) {
      //cout << "2(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";
      // get the sequence item
      //cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
      // apply mass difference
      //cell.delta += sd;

      mark = true;

      //oSequenceMapping.processedAA.push_back(cell);

      //s0 = s1;
      //c0 = c1;

      continue;

    } else if(proteinItemAtContigIndex == -1) {
      //cout << "3(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";
      // mark as gapped and move to next
      mark = true;
      continue;

    // at the end of a hole
    } else if(mark) { //proteinItemAtLastContigIndex == -1) {
      //cout << "4(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";
      mark = false;
      // get the sequence item
      cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
      // apply mass difference
      cell.delta += sd;
      // get sequence between the 2 points from consensus
      getSequenceBetween(cell, iSequenceMapping, proteinItemAtLastContigIndex, proteinItemAtContigIndex);

    // inside sequence scope and star scope
    } else {
      //cout << "5(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";
      if(proteinItemAtContigIndex - proteinItemAtLastContigIndex > 1)
        getSequenceBetween(cell, iSequenceMapping, proteinItemAtLastContigIndex, proteinItemAtContigIndex);

      else {
        // get the sequence item
        cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
        // apply mass difference
        cell.delta += sd;
      }
    }


    // store the cell
    oSequenceMapping.processedAA.push_back(cell);

    // update lower mass intervalls items for next iteration as this interaction's upper masses
    proteinItemAtLastContigIndex = proteinItemAtContigIndex;
    s0 = s1;
    c0 = c1;
  }

  // if there are unprocessed items:
  if(mark) {
    // declare cell object
    aaCell cell;    // get the sequence item
    cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
    // apply mass difference
    cell.delta += sd;
    // get sequence between the 2 points from consensus
    getSequenceBetween(cell, iSequenceMapping, proteinItemAtLastContigIndex, proteinItemAtContigIndex);
    // store the cell
    oSequenceMapping.processedAA.push_back(cell);
  }

  // suffix mass -- after star end, until contig end

  // get the last item in abruijn graph
  double ddc = s0; //diffData.abruijnData[diffData.abruijnData.size() - 1].second;
  // subtract from parent mass (along with massMH)
  double ddd = (*spsFiles->m_starSpectra)[star].parentMass - ddc - AAJumps::massMH;
  // if greather tham peak mass tolerance
  if(fabs(ddd) > m_parentMassTol) {
    // declare the cell
    aaCell cell2;
    // set the mass value
    cell2.delta = ddd;
    // store the cell
    oSequenceMapping.processedAA.push_back(cell2);
  }

  //cout << " --> ";
  //for(int i = 0 ; i < oSequenceMapping.processedAA.size() ; i++)
  //  cout << oSequenceMapping.processedAA[i].aa << ".";
  //cout << endl;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::dump_contigData(ostream &sout, int contig, int star, int homolog, DiffData &diffData, SequenceMapping &iSequenceMapping)
{
  // get the flipped state
  int flipped = getContigState(homolog);

  // calc diffs (star)
  double s0 = 0.0, s1;
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++) {
    s1 = diffData.abruijnData[i].second;
    double sd = s1 - s0;
    diffData.abDiff.push_back(sd);
    s0 = s1;
  }

  // calc diffs (contig)
  double c0 = 0.0, c1;
  for(int i = 0 ; i < diffData.contigData.size() ; i++) {
    c1 = diffData.contigData[i];
    double cd = c1 - c0;
    diffData.contigDiff.push_back(cd);
    c0 = c1;
  }


  // diference of differences
  for(int i = 0 ; i < diffData.abDiff.size() ; i++) {
    int i1 = diffData.abruijnData[i].first;
    double df = diffData.abDiff[i] - diffData.contigDiff[i1];
    diffData.diffDiff.push_back(df);
  }

  sout << "---- " << contig << " -----" << star << " ------------" << endl;
  sout << "flipped: " << flipped << endl << endl;

  for(int i = 0 ; i < iSequenceMapping.processedAA.size() ; i++) {
    if(iSequenceMapping.processedAA[i].aa.length()) {
       if(fabs(iSequenceMapping.processedAA[i].delta) < m_peakMassTol)
         sout << iSequenceMapping.processedAA[i].aa;
       else
         sout << "(" << iSequenceMapping.processedAA[i].aa << "," << iSequenceMapping.processedAA[i].delta << ")";
    } else
      sout << "[" << iSequenceMapping.processedAA[i].delta << "]";
  }
  sout << endl << endl;

  sout << "abruijn: ";
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++)
    sout << "(" << diffData.abruijnData[i].first << " ; " << diffData.abruijnData[i].second << ")  ";
  sout << endl << endl;

  sout << "contig: ";
  for(int i = 0 ; i < diffData.contigData.size() ; i++)
    sout << diffData.contigData[i] << " ; ";
  sout << endl << endl;

  sout << "abruijn diff: ";
  for(int i = 0 ; i < diffData.abDiff.size() ; i++)
    sout << diffData.abDiff[i] <<  " ; ";
  sout << endl << endl;

  sout << "config diff: ";
  for(int i = 0 ; i < diffData.contigDiff.size() ; i++)
    sout << diffData.contigDiff[i] << " ; ";
  sout << endl << endl;

  sout << "diff of diffs: ";
  for(int i = 0 ; i < diffData.abDiff.size() ; i++)
    sout << diffData.diffDiff[i] << " ; ";
  sout << endl << endl;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getMassPrefix(int contig, int consensus, double &prefix, double &suffix)
{
  // test for needed files
  if(!spsFiles->m_abruijn) return;

  abinfo_t::iterator i0 = spsFiles->m_abruijn->find(contig);

  vector<double> masses;

  abContigData_t::iterator i1 = i0->second.second.begin();
  for(; i1 != i0->second.second.end() ; i1++)
    for(int i = 0 ; i < i1->first.size() ; i++)
      if(i1->first[i] == consensus)
        masses.push_back(i1->second[i]);

  // make sure its ordered
  sort(masses.begin(), masses.end());

  prefix = suffix = 0.0;
  if(!masses.size())
    return;

  // first is mass prefix
  prefix = masses[0];
  // last is suffix
  suffix = masses[masses.size() - 1];
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::generateStatistics(Spectrum &spectrum, string &peptide, float &m_B, float &m_Y, float &m_BYint)
{
  m_Y = m_B = m_BYint = 0.0;

  if(!spectrum.size()) {
    WARN_MSG("Generate statistics: spectrum size is 0. Returning.");
    return;
  }

  psmPtr psm(new PeptideSpectrumMatch);
  //Set spectrum
  psm->m_spectrum = &spectrum;
  //Set PSM on spectrum (probably not necessary, but might make things easier
  spectrum.psmList.push_back(psm);

	MS2ScoringModel model;
	// Define annotation file location
	string annotationFile = m_annotationModelDirectory;
	annotationFile += '/';
	annotationFile += m_annotationModel;
	// Set annotatioln file
	model.LoadModel((char*)annotationFile.c_str());

  // ions
  vector<string> bIons, yIons;
  string b, y;
  model.getBreakIonsNames(bIons, yIons);
  for(int i = 0 ; i < bIons.size() ; i++) {
    if(b.length())
      b += ',';
    b += bIons[i];
  }
  for(int i = 0 ; i < yIons.size() ; i++) {
    if(y.length())
      y += ',';
    y += yIons[i];
  }

/*
  switch(spectrum.msFragType) {
  case Spectrum::FragType_PRM:
    b = "b";
    y = "y";
    break;
  case Spectrum::FragType_CID:
    b = "b,b++";
    y = "y,y++";
    break;
  case Spectrum::FragType_ETD:
    b = "c,c++";
    y = "z,z.++,z.++-iso";
    break;
  case Spectrum::FragType_HCD:
    b = "b,b++";
    y = "y,y++";
    break;
  default:
    b = "b,b++";
    y = "y,y++";
    break;
  }
*/

  DEBUG_MSG("Ions for spectrum scan: " << spectrum.scan);
  DEBUG_MSG("y: " << y);
  DEBUG_MSG("b: " << b);

	string allIons("all");
  float aux = m_massShift - AAJumps::massHion;

  // Annotate the spectrum
  psm->annotate(peptide, allIons, model, aux, aux, m_peakMassTol);

  // get statistics
  SpectrumAnnotStatistics stats;

  SpectrumAnnotParameter allParam;
  allParam.ionNames = allIons;
  SpectrumAnnotParameter bParam;
  bParam.ionNames = b;
  SpectrumAnnotParameter yParam;
  yParam.ionNames = y;

  m_BYint = stats.percentExplainedIntensity(*psm, allParam);
  m_B     = stats.observedBreaks(*psm, bParam);
  m_Y     = stats.observedBreaks(*psm, yParam);
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::cleanProteinName(string &proteinName)
{
  for(int i = 0 ; i < proteinName.length() ; i++)
    if(proteinName[i] == '|' || proteinName[i] == ';' || proteinName[i] == '\t' || proteinName[i] == '_')
      proteinName[i] = ' ';
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::processProteinsFile(int protein, PdProteinInfo &proteinData)
{
  // do we have a proteins' file?
  if(!spsFiles->m_fasta) return false;
  // is the protein index consistent with the proteins' file?
  if(spsFiles->m_fasta->size() <= protein) return false;

  // set protein name
  proteinData.proteinName = getProteinName(protein);
  // get the sequence
  string sequence = getProteinSequence(protein);
  // cycle thru the sequence
  for(int j = 0 ; j < sequence.size() ; j++) {
    aaCell cell;
    // the aa sequence
    cell.aa = sequence[j];
    // just one
    cell.colspan = 0;
    // its position
    cell.startPosition = j;
    proteinData.proteinDetail.processedAA.push_back(cell);
  }
  // set protein size
  proteinData.proteinDetail.startPosition = 0;
  proteinData.proteinDetail.endPosition = proteinData.proteinDetail.processedAA.size();
  proteinData.proteinLength = proteinData.proteinDetail.endPosition - proteinData.proteinDetail.startPosition;
  // set internal protein index value
  proteinData.ProteinIndex = protein;

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs, bool translate, bool csps)
{
  // First file: homglue matches
	for(int i = 0 ; i < contigMatches.size() ; i++)
    // if equal to -1, then has no match
	  if(contigMatches[i][0] != -1) {
      // create strcuture for current contig
      ContigMatchData contigMatchData;
      // Store contig index
      contigMatchData.contigIndex = i;
	    // find current protein
	    std::map<int, std::vector<ContigMatchData> >::iterator it;
	    it = contigs.find(contigMatches[i][0]);
	    // if key exists, insert the new value
	    if(it != contigs.end()) {
	      it->second.push_back(contigMatchData);
	    } else {
	      // if it doesn't exist, create a new entry.
        // create the vector
        vector<ContigMatchData> aux;
        aux.push_back(contigMatchData);
       // Insert CSPS contig (index) info indexed by protein index
        contigs.insert(pair<int, vector<ContigMatchData> > (contigMatches[i][0], aux));
	    }
  	}

  // Second file: CSPS contig matches
	for(int i = 0 ; i < contigMatchesIndex.size() ; i++)
    // test if there is any data to store
	  if(contigMatchesIndex[i].size() > 0) {

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, i);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig)
    		for(int j = 0 ; j < contigMatchesIndex[i].size() ; j++) {
  		    PairContigProteinMassIdx pairItem;
  		    pairItem.contigMassIdx = (int)contigMatchesIndex[i][j][0];
  		    pairItem.proteinMassIdx = (int)contigMatchesIndex[i][j][1];
  		    contig->pair.push_back(pairItem);
  		  }
	  }

	// 3rd file: Contig spectra
	for(int i = 0 ; i < contigSpectra.size() ; i++)
	  if(contigSpectra[i].size() > 0) {

      int ctg = getContigFromAllContig(i);

      int q = (translate ? ctg : i);

      if(q == -1)
        continue;

      Spectrum sp = contigSpectra[i];

      // reverse spectrum, if needed
      if(translate) {
        if(getContigState(ctg)) {
          sp.reverse(0.0 - AAJumps::massH2O);
          DEBUG_MSG("Reversing contig spectra (coverage): " << i);
        }
      }

      // reverse spectrum, if needed
      if(csps) {
        if(getCspsContigState(i)) {
          sp.reverse(0.0 - AAJumps::massH2O);
          DEBUG_MSG("Reversing cSPS contig spectra (coverage): " << i);
        }
      }

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, q);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig) {
        if(csps) {
          DEBUG_MSG("cSPS contig: " << contig->contigIndex << ", idx: " << q << "  size: " << sp.size());
        } else {
          DEBUG_MSG("SPS contig: " << contig->contigIndex << ", idx: " << q << "  size: " << sp.size());
        }
        stringstream str;
    		for(int j = 0 ; j < sp.size() ; j++) {
  		    float mass = sp[j][0];
  		    contig->contigMass.push_back(mass);
  		    str << mass << " ; ";
  		  }
 		    DEBUG_MSG("masses: " << str.str());
      }
	  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::processContigs(int ProteinIndex, PdProteinInfo &proteinData, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source)
{
  // Get csps contig information. If not found for this protein, NULL is returned.
  std::map<int, vector<ContigMatchData> >::iterator it2;
  it2 = source.find(proteinData.ProteinIndex);

  // Add CSPS contig information (if exists)
  if(it2 == source.end())
    return;

  DEBUG_MSG("###############################");
  DEBUG_MSG("Protein: #" << proteinData.ProteinIndex);

  // There may be several contigs for this protein. Must cycle thru them
  for(int i = 0 ; i < (*it2).second.size() ; i++) {

    DEBUG_MSG("###############################");
    DEBUG_MSG("Contig: #" << (*it2).second[i].contigIndex);

    // get csps contig span and location
    int size = (*it2).second[i].pair.size();
    if(!size) {
      WARN_MSG("cSPS contig has o size.");
      continue;
    }
    int minimunIndex = (*it2).second[i].pair[0].proteinMassIdx;
    int maximunIndex = (*it2).second[i].pair[size-1].proteinMassIdx;

    // Define the data structure to hold it and set some initial data
    ProteinDetaisContigInfo contig;
    contig.startPosition = minimunIndex;
    contig.endPosition = maximunIndex;


    // Loop thru contig size to generate info cells
    for(int j = 0 ; j < size-1 ; j++) {

      // get the corresponding contig indexes
      int contigIndex = (*it2).second[i].pair[j].contigMassIdx;
      int contigIndexAfter = (*it2).second[i].pair[j+1].contigMassIdx;

      // get the correspondif protein indexes
      int proteinIndex = (*it2).second[i].pair[j].proteinMassIdx;
      int proteinIndexAfter = (*it2).second[i].pair[j+1].proteinMassIdx;
      // Get the mass at this point

      if(contigIndex >= (*it2).second[i].contigMass.size() || contigIndex >= (*it2).second[i].contigMass.size()) {
        DEBUG_MSG("Contig #" << j << ": Invalid mass index: " << contigIndex << " > " << (*it2).second[i].contigMass.size());
        return;
      }

      float massCsps = (*it2).second[i].contigMass[contigIndex];
      // and at the point after
      float massAfter = (*it2).second[i].contigMass[contigIndexAfter];

      string sequence;

      // get AA sequence for interval
      if( (contigIndex >= proteinData.proteinDetail.startPosition ) && (contigIndexAfter < proteinData.proteinDetail.endPosition) )
        for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.proteinDetail.processedAA.size()) ; k++)
          sequence += proteinData.proteinDetail.processedAA[k].aa;

      // get masses for sequence
      vector<float> masses;
      getMasses((char*)sequence.c_str(), masses);
      float mass = 0.0;
      for(int k = 0 ; k < masses.size() ; k++)
        mass += masses[k];

      string result = sequence;

      // massAux used to calculate mass difference
      float massAux = (massAfter - massCsps) - mass;

      if(fabs(massAux) > m_peakMassTol) {
        result = "(";
        result += sequence;
        result += ", ";
        result += parseFloat(massAux, 2);
        result += ")";
      }


      DEBUG_MSG("---------- " << j << " ----------");
      DEBUG_MSG("contig mass index       " << contigIndex);
      DEBUG_MSG("contig mass index after " << contigIndexAfter);
      DEBUG_MSG("protein mass index      " << proteinIndex);
      DEBUG_MSG("protein mass index after " << proteinIndexAfter);
      DEBUG_MSG("contig mass             " << massCsps);
      DEBUG_MSG("contig mass after       " << massAfter);
      DEBUG_MSG("protein sequence        " << sequence);
      DEBUG_MSG("protein sequence mass   " << mass);
      DEBUG_MSG("contig mass difference  " << (massAfter - massCsps) );
      DEBUG_MSG("mass difference         " << massAux);



      // cell to hold one data item
      aaCell cell;

      // save it
      cell.aa = result;
      // calculate the colspan value for this cell
      cell.colspan = (*it2).second[i].pair[j+1].proteinMassIdx - (*it2).second[i].pair[j].proteinMassIdx - 1;
      // and the start position
      cell.startPosition = (*it2).second[i].pair[j].proteinMassIdx;

      // add the cell to the list
      contig.processedAA.push_back(cell);
    }

    // set the contig name
    contig.name = parseInt((*it2).second[i].contigIndex);

    // set the base 1 index
    contig.base1Idx = (*it2).second[i].contigIndex + 1;

    // add the csps contig info to the csps contig list
    target.insert(pair<int, ProteinDetaisContigInfo > ((*it2).second[i].contigIndex, contig));
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::populateContigNames(PdProteinDetail &contig)
{
  PdProteinDetail::iterator it;
  for(it = contig.begin() ; it != contig.end() ; it++) {
    int contigID = it->first;
   	string name = getContigName(contigID);
    it->second.name = (name.length() > 0 ? name : parseInt(contigID + 1) );
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getContigName(int i)
{
  if(spsFiles->m_contigNames)
  	if( (i>= 0) && i < spsFiles->m_contigNames->size() )
	  	return (*spsFiles->m_contigNames)[i];
	stringstream ret;
  int aux = getAllContigFromContig(i);
  if(aux < 0)
    return string("--");
	ret << "Contig:" << (aux + 1);
	return ret.str();
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
ContigMatchData *ReportTableGenerator::getContigData( std::map<int, std::vector<ContigMatchData> > &contig, int contigIndex)
{
  std::map<int, vector<ContigMatchData> >::iterator it;
  // search for contig linked to given contig index
  for( it = contig.begin() ; it != contig.end() ; it++)
    // Cycle thru contig vector (for each protein)
    for(int i = 0 ; i < (*it).second.size() ; i++)
      if( (*it).second[i].contigIndex == contigIndex)
        // return pointer to struct
        return &((*it).second[i]);

  return NULL;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::generateCoverageOutput(PdProteinInfo &proteinData)
{
  // test for proteins with no contig mapped to them, and do not generate them
  if((proteinData.cspsDetails.size() == 0) && (proteinData.spsDetails.size() == 0))
    return 1;

  // Get the protein sequence. From this, we know the lenght
  ProteinDetaisContigInfo & proteinSequence = proteinData.proteinDetail;

  // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
  vector<int> spsID, cspsID;
  // get the csps contig indexes
  getOrder(proteinData.cspsDetails, cspsID);
  // get the sps contig indexes
  getOrder(proteinData.spsDetails, spsID);

  // generate protein sequence
  for(int j = 0 ;  j < proteinSequence.processedAA.size() ; j++) {
    if(j)
      proteinData.proteinSequenceEntryData += TABLE_SEP_L1;
    proteinData.proteinSequenceEntryData += proteinSequence.processedAA[j].aa;
  }

  // Add CSPS contig information (if exists)
  generateOutputContig(proteinData, cspsID, proteinData.cspsEntryData, proteinData.cspsDetails, false);

  // Add SPS contig information (if exists)
  generateOutputContig(proteinData, spsID, proteinData.spsEntryData, proteinData.spsDetails, true);

	return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::generateOutputContig(PdProteinInfo &proteinData, vector<int> &vectorID, std::string &page, PdProteinDetail &contig, bool link)
{
  // Add CSPS contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = contig[vectorID[j]];

    // get protein size
    int contigSize =  contigSequence.processedAA.size();

    // add contig separator
    if(j) page += TABLE_SEP_L1;

    // Write the contig id and link
    if(link) {
      // add ID
      page += getIntFromSeqName(contigSequence.name);
      page += TABLE_SEP_L2;
      // add name
      page += contigSequence.name;
      page += TABLE_SEP_L2;
    } else {
      // add id
      page += parseInt(contigSequence.base1Idx);
      page += TABLE_SEP_L2;
      // add name
      page += parseInt(contigSequence.base1Idx);
      page += TABLE_SEP_L2;
    }

    // start position
    page += parseInt(contigSequence.startPosition);
    page += TABLE_SEP_L2;
    // end position
    page += parseInt(contigSequence.endPosition);
    page += TABLE_SEP_L2;

    // cycle thru
    for(int k = 0 ; k < contigSize ; k++) {

      // add contig item separator
      if(k) page += TABLE_SEP_L3;

      // contig item start position
      page += parseInt(contigSequence.processedAA[k].startPosition);
      page += TABLE_SEP_L4;
      // contig item colspan
      page += parseInt(contigSequence.processedAA[k].colspan);
      page += TABLE_SEP_L4;
      // contig item content
      page += contigSequence.processedAA[k].aa;
      //page += TABLE_SEP_L4;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getOrder(PdProteinDetail &contig, vector<int> &order)
{
  // contig iterator
  PdProteinDetail::iterator it;
  // cycle thru all contigs
  for( it = contig.begin() ; it != contig.end() ; it++ ) {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = it->second;
    // Check for empty contigs
    if(!contigSequence.processedAA.size()) continue;
    // add item
    order.push_back(it->first);
  }

  // Sort
  sort(order.begin(), order.end());
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getIntFromSeqName(string seq)
{
  vector<string> aux;
  stringSplit(seq, aux, ":");
  if(aux.size() > 1)
    return aux[1];
  return "";
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
