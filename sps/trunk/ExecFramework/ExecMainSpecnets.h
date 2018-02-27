#ifndef __ExecMainSpecnets_H__
#define __ExecMainSpecnets_H__

// Module Includes
#include "ExecBase.h"

// SpecNets Includes
#include "spectrum.h"
#include "SpectrumPairSet.h"
#include "inspect_parse.h"
#include "PeptideSpectrumMatch.h"
#include "db_fasta.h"
#include "AlignmentPenaltyBased.h"
#include "ExecSvmStatistics.h"



// System Includes
#include <string>
#include <vector>

namespace specnets
{
	using namespace std;

	/*! \brief Execution class for Spectrum Networks projects after ExecFilterStarPairs

   Projects peptide identifications from identified to unidentified spectra, filters out
   bad identifications using FDR thresholds on PSM SVM scores and aggregates annotated spectra
   into final spectral networks using spectral alignments and overlapping PSMs (ExecHomologyAssembly).

   */
  class ExecMainSpecnets : public ExecBase
  {
  protected:
    /*! \brief Determines annotations for specnets nodes

     Used from invoke() to get the initial+propagated specnets
     node annotations,
     @sa ExecSpecNetsPropagation
     @param db Database used to search for the seed annotations
     @param psms Output set of Peptide Spectrum Matches
     @param psmSpectra PSM-filtered star spectra (containing only annotated PRMs for identified spectra)
     @param svm Input svm model
     @param storeMatchInfo If true propagation sets m_psmSpectra, m_psms_midx and m_psms_mp
     @param aminoacids Amino acid masses used for node annotations
     @return True if no errors, false otherwise
     */
  	bool annotateNodes(DB_fasta *db,
  	                   PeptideSpectrumMatchSet * psms,
  			           SpecSet *psmSpectra,
                       ExecSvmStatistics * svm,
  			           bool storeMatchInfo,
  			           AAJumps &aminoacids);

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
  	ExecMainSpecnets(void);

    /*! \brief Constructor which takes only the set of parameters

     Constructor which takes only the set of parameters required to execute the
     module. The parameters specify all data files for input and output, and requires
     that the user call loadInputData() in order to read all data from files into
     the members of the object. The parameters can be verified using the
     validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
  	ExecMainSpecnets(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     executing the alignment. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param spectra Input spectra to be projected from/to
     @param pairs Input set of pairs used for projections
     @param projSpectra Output spectra with projected IDs
     @param projPeaks Output indices of pairs of peaks matched by projection
     @param projPairs Output set of pairs used for projections
     @param projInfo Per-spectrum propagation statistics
     */
	  ExecMainSpecnets(const ParameterList & inputParams,
        SpecSet * msSpectra,
        SpecSet * scoredSpectra,
        SpecSet * starSpectra,
        SpectrumPairSet * pairs,
        MS2ScoringModel * model,
        DB_fasta * db,
        PenaltyMatrix * penaltyMatrixBlosum,
        PenaltyMatrix * penaltyMatrixMods,
        PeptideSpectrumMatchSet * psms,
        PeptideSpectrumMatchSet * origPsms,
        SpecSet * psms_spectra,
        SpecSet * psms_midx,
        vector<vector<int> > * psms_mp,
        SpecSet * snets_contigs,
        SpecSet * snets_midx,
        vector<vector<int> > * snets_mp);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecMainSpecnets(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & input_params) const;

    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module.

     In order to call this method succesfully all the necessary data for
     execution must already be loaded into memory (data members). This can
     be accomplished using the loadInputData() method or by calling the
     appropriate constructor.

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
     data into the members before calling invoke().

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
     @sa ExecBase(const ParameterList & input_params), saveInputData()
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

    private:
  	// Input data
    SpecSet * m_specsMsSpectra; //! Input MS spectra
    SpecSet * m_specsScoredSpectra; //! Input scored spectra
    SpecSet * m_starSpectra;    //! Input star spectra
    PeptideSpectrumMatchSet * m_origPsms; //! Input PSMs from MSGFDB, Inspect, etc.
    SpectrumPairSet * m_pairs;  //! Input set of pairs usable for propagations/networks
    MS2ScoringModel * m_model;  //! ion fragmentation model (for PSM SVM scores / FDR)
    DB_fasta * m_db;
    bool ownInput; //! Does this object "own" the input data pointers (and hence has to free them)

  	// Output data
  	PeptideSpectrumMatchSet * m_psms; //! (--> m_psms) FDR-filtered PSMs (modified and unmodified)
    SpecSet * m_psmSpectra; //! PSM-filtered star spectra (containing only annotated PRMs for identified spectra)
    SpecSet * m_psms_midx;            //! Per-PSM indices of matched spectrum/protein masses
    vector<vector<int> > * m_psms_mp; //! Per-PSM: best-matched protein (-1 if none), # mods, prefix/suffix match (cols 0-2, respectively)

    SpecSet * m_snets_contigs;  //! Spectral network de novo sequences as determined by ExecHomologyAssembly
  	SpecSet * m_snets_midx;     //! Per-spectral-network indices of matched spectrum/protein masses
  	vector<vector<int> > * m_snets_mp; //! Per-spectral-network: best-matched protein (-1 if none), # mods, prefix/suffix match (cols 0-2, respectively)
  	bool ownOutput; //! Does this object "own" the output data pointers (and hence has to free them)

  	//Alignment parameters
  	PenaltyMatrix * m_penaltyMatrixBlosum;
  	PenaltyMatrix * m_penaltyMatrixMods;

  	// Internal data
//    PeptideSpectrumMatchSet * m_inspectResults; //! ( --> m_psms_in) Initial set of PSMs (usually Inspect annotations used to seed specnets)
//	  PeptideSpectrumMatchSet m_psms_fdr; 	//! PSMs available to project from
//    SpecSet m_pepsAsSpecs;          //! Per-PSM: matched peptide sequence represented as a spectrum
//    SpecSet m_propagation_spectra;        //! Per-PSM indices of matched spectrum/protein masses
//    PeptideSpectrumMatchSet * m_specnetsResults; //! ( --> m_propagation_psms) Unfiltered set of spectral networks PSMs
  	SpectrumPairSet m_propagation_pairs;  //! Set of pairs used for propagations
    SpecSet         m_propagation_midx;   //! Per-propagation indices of matched spectrum/spectrum masses (in m_pepsAsSpecs)
  	vector<vector<int> > m_propagation_info; //! Per-spectrum propagation/PSM info (see ExecSpecNetsPropagation->m_projInfo)
  };

} // namespace specnets

#endif // __ExecMainSpecnets_H__
