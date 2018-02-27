#ifndef __EXEC_REPORT_PROTEIN_COVERAGE_H__
#define __EXEC_REPORT_PROTEIN_COVERAGE_H__


// External Includes
#include "ExecBase.h"
#include "ParameterList.h"
#include "db_fasta.h"
#include "spectrum.h"

// System Includes
#include <string>
#include <vector>
#include <map>


namespace specnets {


////////////////////////////////////////////////////////////////////////////////
// Ordering vector element. Used to estabelish contig order when rendering protein coverage report
struct ContigOrdering {
  int contigIndex;
  int startIndex;
  int endIndex;

  // Used by sort method. First order by beggining index, then by ending index.
  bool operator<(const ContigOrdering &o) const
  {return (startIndex == o.startIndex ? endIndex < o.endIndex : startIndex < o.startIndex);};

};
////////////////////////////////////////////////////////////////////////////////
// Holds a contig/protein mass index pair unit
struct PairContigProteinMassIdx {
  int proteinMassIdx;
  int contigMassIdx;
};


// holds data for CSPS contigs
struct ContigMatchData {
  std::vector<float>                    contigMass;
  std::vector<PairContigProteinMassIdx> pair;
  int                                   contigIndex;
};

////////////////////////////////////////////////////////////////////////////////
// holds aa sequence for a cell and the number of columns it ocupies
struct aaCell {
  std::string   aa;
  int           colspan;
  int           startPosition;
};


// Structure to hold a contig/protein info in protein details context
struct ProteinDetaisContigInfo {
	// Cell start position
  int                       startPosition;
  // Cell end position
  int                       endPosition;
  // Processed AAs for page generation
  std::vector<aaCell>       processedAA;
	// Contig name
  std::string								name;
  //base 1 contig index
  int                       base1Idx;
};


// contig information indexed by contig index
typedef std::map<int, ProteinDetaisContigInfo> PdProteinDetail;



// Holds information for an entire protein and it's contigs
struct PdProteinInfo {

	// protein sequence
	ProteinDetaisContigInfo proteinDetail;

	// Holds information for sps contigs.
	PdProteinDetail spsDetails;

	// Holds information for csps contigs.
	PdProteinDetail cspsDetails;


	// Input and output file names
	string oFileName;
	string oFileNameCsv;
	string iFileName;
	// actual protein index
	int ProteinIndex;
  // protein name
  string proteinName;
  // Structure to hold page data
  std::string 							page;
  // Structure to hold csv file
  std::string               csvFile;
};



// Structure to hold all proteins info. Mapped by protein index.
typedef std::map<int, PdProteinInfo> PdProtein;


////////////////////////////////////////////////////////////////////////////////
class ExecReportProteinCoverage  : public specnets::ExecBase {

	// FASTA file loaded
  //std::vector<std::string> m_references;
  DB_fasta  m_fasta;

  // Contig Names
  std::vector<std::string> m_contigNames;

  int indexBase;

  // SPS Contig Matches Indexes
  SpecSet m_spsContigMatchesIndex;

  // SPS Contig Matches;
  std::vector<vector<int> > m_spsContigMatches;

  // SPS Contig Spectra
  SpecSet m_spsContigSpectra;




  // CSPS Contig Matches Indexes
  SpecSet m_cspsContigMatchesIndex;

  // CSPS Contig Matches;
  std::vector<vector<int> > m_cspsContigMatches;

  // CSPS Contig Spectra
  SpecSet m_cspsContigSpectra;



	// Structure to hold information for all proteins
	PdProtein 	m_protein;

  // Structure to hold CSPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_cspsContigInfo;

  // Structure to hold SPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_spsContigInfo;



	// function prototypes

	// Generate HTML page header
	void genheader(std::string & shtml);

	// Generate page footer
  void genFooter(std::string & page);
  void genTableFooter(std::string & page);

	// Reads a single file into a buffer
	char *readFile(std::string &iFileName, int &) ;

	// Parses an overlap file into the m_protein structure
	int parseFile(int, char *buffer);

	// Condenses parsed overlap information into cell info, for page generation
	int processData(int);

	// Put loaded contig files in a structure usable for fast page generation
  bool processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs);


	// 2nd processing stage, where mass values are compared against AAs and mass differences calculated
	int replaceMassesWithAaReferences(int);

	// Replaces contig IDs with contig names
  void populateContigNames(PdProteinDetail &contig);

	// Generates the output HTML file based on the processed overlap information
	int generateOutputFile2(int, std::string &page, int cellsPerLine);

	// Generates the output CSV file based on the processed overlap information
  int generateOutputFileCsv(int, std::string &page, int cellsPerLine);

	// Writes a string to a file. Used to write HTML files.
	int outputFile(std::string &oFileName, std::string &page);


  // Gets a protein name from the protein reference list
//	int getProteinReferences(std::string &refIndex);

  // Reads the contig names file
	int getContigNames(std::string &contigNamesFileName);

	// Gets a contig name from the contig name list
	std::string getContigName(int i);

	ContigMatchData *getContigData(std::map<int, std::vector<ContigMatchData> > &contig, int proteinIndex);

  bool getContigMassAtProteinIndex(ContigMatchData &contig, int proteinIndex, float &contigMassAtIndex);

  void processContigs(int ProteinIndex, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source);

  // Generate contig HTML sequence for report
  int generateOutputContig(int &i, int proteinIdx, vector<int> &vectorID, std::string &page, int cellPerLine, PdProteinDetail &contig, bool link);

  // Generate contig CSV sequence for report
  int generateOutputContigCsv(int &i, int proteinIdx, vector<int> &vectorID, std::string &page, int cellPerLine, PdProteinDetail &contig, bool link);

  void processProteinsFile(void);

  // Debug purposes.
  void dumpContig(std::map<int, std::vector<ContigMatchData> > &contigs);

  // Order contig to output based on rendering index.
  void getOrder(PdProteinDetail &contig, int i, int, vector<int> &order);

  //
  string getIntFromSeqName(string seq);


 public:

  //! \name CONSTRUCTORS
  //@{
  /*! \brief The exemplar constructor.

   Generally this constructor should not be used. It is only used by the
   module execution factory in order to create an exemplar object (without
   valid parameters) which is then used to create a real (valid) object
   using the clone() method.
   @sa clone()
   */
  ExecReportProteinCoverage();

  /*! \brief The default constructor

   This is the default constructor. A valid set of parameters must be
   supplied in the input parameter. The parameters can then be verified
   using the validateParams() method.
   @sa validateParams()
   @param inputParams structure containing all input parameters necessary for execution
   */

	ExecReportProteinCoverage(const specnets::ParameterList &params);
  //@}

  //! \name DESTRUCTOR
  //@{
	~ExecReportProteinCoverage() {};
  //@}

  //! \name ACCESSORS
  //@{
  /*! \brief Creates a new module of the virtual class with the given params.

   The parameters should be sufficient for the derived class to be invoked properly.

   @return A pointer to the newly created object
   */
  virtual ExecBase * clone(const ParameterList & input_params) const;


  /*! \brief Returns validity of the module.

   This is a non-virtual method that simply returns the value of the internal
   validity flag. This flag should have been set before calling this method
   by calling the validate() method.

   @return The value of the internal object validity flag.
   @sa validateParams(std::string & error)
   */
  bool isValid(void) const;
  //@}

  //! \name MODIFIERS
  //@{
  /*! \brief Executes the module.

   In order to call this method succesfully all the necessary data for
   execution must already be loaded into memory (data members). This can
   be accomplished using the loadInputData() method.

   @return True if execution finished successfully, false otherwise.
   @sa loadInputData()
   */
  virtual bool invoke(void);

  /*! \brief Loads the input data from the files specified in the params.

   Loads all the data from files into the data members. This method is
   primarily used by the execution module to load necessary data when
   executing in a separate process..

   @return True if data was loaded successfully, false otherwise.
   @sa ExecBase(const ParameterList & input_params), saveOutputData()
   */
  virtual bool loadInputData(void);

  /*! \brief Saves all the result data to files specified in the params.

   Saves all the result data generated by the module to files specified
   by values in the params. This is used to either save the data permanantly,
   or to be loaded back in after remote execution is finished and results
   need to be merged by the merge() method.

   @param filenames A list of file names that contain the data necessary to run the module
   @return True if data was saved successfully, false otherwise.
   @sa ExecBase(const ParameterList & input_params), loadInputData(), merge()
   */
  virtual bool saveOutputData(void);

  /*! \brief Saves all the internal data into files specified in the params.

   Saves all the data required for an external process to execute the
   module. The external process would call loadInputData() to reload the
   data into the members before calling invoke(). The user passes a vector
   that will contain the names of all files necessary to run the module as
   a separate process. The first file in this list will always be the main
   parameter file.

   @param filenames A list of file names that contain the data necessary to run the module
   @return True if data was saved successfully, false otherwise.
   @sa ExecBase(const ParameterList & input_params), loadInputData()
   */
  virtual bool saveInputData(std::vector<std::string> & filenames);

  /*! \brief Loads the output data from the files specified in the params.

   Loads all the data from output files into the data members. The purpose
   of this is to ready "child" modules to be merged back together after
   being executed separately.

   @return True if data was loaded successfully, false otherwise.
   @sa ExecBase(const ParameterList & input_params), saveOutputData()
   */
  virtual bool loadOutputData(void);

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.

     @param numSplit Number of separate modules to split into
     @return The set of sub-modules (children) that the original module has been split into.
     An empty vector implies an error.
     @sa merge()
     */
    virtual std::vector<ExecBase *> const & split(int numSplit);

    /*! \brief Merges the child modules back into a complete result

     This method is only called when split is used for parallel execution
     and after the execution of all children has been completed. This method
     will merge the results generated by each of the child modules into
     one cohesive result as if the module had been run as a single entity.

     @return True if merge could be performed successfully, false otherwise.
     @sa split(int numSplit)
     */
    virtual bool merge(void);

  /*! \brief Performs validation of the input parameters.

   Checks the parameters structure provided in the constructor to see if they
   are sufficient and correct to invoke the module. Also sets the internal
   validity flag so that isValid() will return the correct result.

   @param error A description of the error (if any occurs)
   @return True if the parameters for the module are valid, false otherwise.
   @sa isValid()
   */
  virtual bool validateParams(std::string & error);
  //@}


  void indexBase0(void) {indexBase=0;};

  void indexBase1(void) {indexBase=1;};


};

}

#endif
