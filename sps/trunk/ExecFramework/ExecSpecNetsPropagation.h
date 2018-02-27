#ifndef __ExecSpecNetsPropagation_H__
#define __ExecSpecNetsPropagation_H__

// Module Includes
#include "ExecBase.h"

// SpecNets Includes
#include "spectrum.h"
#include "SpectrumPairSet.h"
#include "PeptideSpectrumMatch.h"
#include "db_fasta.h"
#include "ExecStatistics.h"
#include "ExecSvmStatistics.h"
#include "AlignmentPenaltyBased.h"



// System Includes
#include <string>
#include <vector>

namespace specnets
{
  using namespace std;

  /*! \brief Execution class for propagation of peptide identifications over spectrum networks

   Finds new peptide identifications for unidentified spectrum network nodes by using
   spectrum alignments to propagate identifications from annotated spectrum nodes.

   */
  class ExecSpecNetsPropagation : public ExecBase
  {
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
    ExecSpecNetsPropagation(void);

    /*! \brief Constructor which takes only the set of parameters

     Constructor which takes only the set of parameters required to execute the
     module. The parameters specify all data files for input and output, and requires
     that the user call loadInputData() in order to read all data from files into
     the members of the object. The parameters can be verified using the
     validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ExecSpecNetsPropagation(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     executing the alignment. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param msms_spectra Input MS/MS spectra to be projected from/to
     @param prm_spectra Input PRM spectra to be projected from/to
     @param star_spectra Input star spectra to be projected from/to
     @param svm SVM used for PSM scores
     @param svm_model Ion fragmentation model (for PSM SVM scores)
     @param pairs Input set of pairs used for propagations
     @param db Fasta sequence database
     @param psms Previous Peptide Spectrum Matches (e.g., InsPecT search)
     @param propagation_pairs Output set of pairs used for propagations
     @param propagation_info Per-spectrum propagation statistics
     @param psms_spectra Spectrum peaks matched by output PSMs
     */
    ExecSpecNetsPropagation(const ParameterList & inputParams,
                            SpecSet * msms_spectra,
                            SpecSet * prm_spectra,
                            SpecSet * star_spectra,
                            AAJumps * aminoacids,
                            ExecSvmStatistics * svm,
                            MS2ScoringModel * svm_model,
                            SpectrumPairSet * pairs,
                            DB_fasta * db,
                            PenaltyMatrix * penaltyMatrixMods,
                            PenaltyMatrix * penaltyMatrixBlosum,
                            PeptideSpectrumMatchSet * psms,
                            SpectrumPairSet * propagation_pairs,
                            vector<vector<int> > * propagation_info,
                            SpecSet * psms_spectra);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecSpecNetsPropagation(void);
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

     In order to call this method successfully all the necessary data for
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

    //! \name MODIFIERS
    //@{
    /*! \brief (DEPRECATED) Executes the module.

     Old implementation that iterates over annotated nodes and propagates
     identifications to non-annotated neighbors.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    virtual bool invoke_deprecated(void);

    SpecSet * m_msms_spectra; //! Input MS/MS spectra
    SpecSet * m_prm_spectra; //! Input PRM spectra
    SpecSet * m_star_spectra; //! Input star spectra
    AAJumps * m_aminoacids; //! Input amino acid masses used for propagation
    SpectrumPairSet * m_pairs; //! Input set of pairs used for propagations
    DB_fasta * m_db; //! FASTA database of protein sequences
    PenaltyMatrix * m_penaltyMatrixBlosum; //! Penalty matrix for alignment
    PenaltyMatrix * m_penaltyMatrixMods; //! Penalty matrix for alignment
    PeptideSpectrumMatchSet * m_psms; //! Input PSMs (if available)
    ExecSvmStatistics * m_svm; //! SVM used for PSM scores
    MS2ScoringModel * m_svm_model; //! Ion fragmentation model (for PSM SVM scores)
    bool ownInput; //! Does this object "own" the input data pointers (and hence has to free them)

    SpecSet * m_psms_spectra; //! Output PSM spectra
    //               SpecSet * m_psms_midx;     //! Output indices of PSM/Protein matches
    //  vector<vector<int> > * m_psms_mp;       //! Output matched protein per spectrum

    //               SpecSet * m_propagation_peaks; //! Output indices of pairs of peaks matched by propagation
    SpectrumPairSet * m_propagation_pairs; //! Output set of pairs used for propagations
    vector<vector<int> > * m_propagation_info; //! Per-spectrum:
    //  proximate level of annotation (col.1),
    //  index of the spectrum where the annotation came from (col.2),
    //  index of the first peak to the right of the annotation (col.3),
    //  100*percentage of explained score for the selected annotation (col.4)
    //  100*percentage of matched peptide breakpoints (either b or y but not both) (col.5)
    bool ownOutput; //! Does this object "own" the output data pointers (and hence has to free them)
  };

} // namespace specnets

#endif // __ExecSpecNetsPropagation_H__
