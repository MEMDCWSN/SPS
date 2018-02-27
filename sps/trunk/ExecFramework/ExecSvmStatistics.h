#ifndef __ExecSvmStatistics_H__
#define __ExecSvmStatistics_H__

// External Includes
#include "ExecBase.h"
#include "ExecStatistics.h"
#include "ParameterList.h"
#include "SvmScaleParameterList.h"
#include "SpectrumAnnotParameterList.h"
#include "DelimitedTextReader.h"
#include "SvmModel.h"
#include "PeptideSpectrumMatch.h"
#include "spectrum.h"
#include "utils.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"

// System Includes
#include <stdio.h>
#include <string>
#include <vector>

namespace specnets
{
  /*! \brief Class to generate svm scores on a per spectrum basis

   */
  class ExecSvmStatistics : public ExecBase
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
    ExecSvmStatistics(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ExecSvmStatistics(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     generating svm scores. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     */
    ExecSvmStatistics(const ParameterList & inputParams,
                      std::vector<string> * starSpectraStatsHeader,
                      std::vector<vector<float> > * starSpectraStats,
                      std::vector<string> * specsMsSpectraStatsHeader,
                      std::vector<vector<float> > * specsMsSpectraStats,
                      std::vector<string> * specsScoredSpectraStatsHeader,
                      std::vector<vector<float> > * specsScoredSpectraStats,
                      PeptideSpectrumMatchSet * peptideResults,
                      PeptideSpectrumMatchSet * peptideOutputResults);
    /*! \brief Constructor which takes both the parameters and the peptide results

     This constructor is intended to be used with invokeWithStatistics.
     @sa validateParams()
     */
    ExecSvmStatistics(const ParameterList & inputParams,
                      PeptideSpectrumMatchSet * peptideResults,
                      PeptideSpectrumMatchSet * peptideOutputResults);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecSvmStatistics(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & inputParams) const;
    //! \name MODIFIERS
    //@{
    /*! \brief Generates svm stats and also invokes Statistics module.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    bool invokeWithStatistics(SpecSet * specsMsSpectra,
                              SpecSet * specsScoredSpectra,
                              SpecSet * starSpectra,
                              MS2ScoringModel * model);

    //@{
    /*! \brief Generates svm stats and also invokes Statistics module for single psm.

     Sets psm.m_score on original object based on SVM score.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    bool getSvmScore(Spectrum * msSpectra,
                     Spectrum * scoredSpectra,
                     Spectrum * starSpectra,
                     PeptideSpectrumMatch &psm,
                     MS2ScoringModel * model);

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
     @sa ExecSvmStatistics(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void);

    /*! \brief Loads the svm models
     *
     * Loads the charge 1, 2 and 3 models and scaling parameters.
     * Called by invoke.
     */
    bool loadModels(void);

    /*! \brief Loads the statistics parameter files
     *
     * Loads ms, star and specsScored annotation parameter files
     * Called by invokeWithStatistics and getSvmScore
     */
    bool loadAnnotStatisticsParams(void);

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecSvmStatistics(const ParameterList & input_params), loadInputData(), merge()
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
     @sa ExecSvmStatistics(const ParameterList & input_params), loadInputData()
     */
    virtual bool saveInputData(std::vector<std::string> & filenames);

    /*! \brief Loads the output data from the files specified in the params.

     Loads all the data from output files into the data members. The purpose
     of this is to ready "child" modules to be merged back together after
     being executed separately.

     @return True if data was loaded successfully, false otherwise.
     @sa ExecSvmStatistics(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadOutputData(void);

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.
     The split method should divide the work into "nodes" and sub-divide into
     "cpus" for each node. Because the precise implementation of parallel execution
     is not known apriori, nodes may be thought of as large batches, while the
     number of CPUs can be thought of as modules which will be run simultaneously
     within those batches. For example, split may be called with numNodes = 1
     and numCpus = 4 on a quad-processor machine to prep for a threaded implementation.
     Or it may be called with numNodes = 4 and numCpus = 1 to prepare for
     distributed execution on 4 separate machines on a network.

     @param numNodes Number of separate nodes
     @param numCpus Number of CPUs per node
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
     @sa split(int numNodes, int numCpus)
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
    /*! \brief Performs validation of the input parameters if we're generating statistics in module.

     Checks the parameters structure provided in the constructor to see if they
     are sufficient and correct to invoke the module. Also sets the internal
     validity flag so that isValid() will return the correct result.

     @param error A description of the error (if any occurs)
     @return True if the parameters for the module are valid, false otherwise.
     @sa isValid()
     */
    bool validateStatisticsParams(std::string & error);
    //@}

  private:
    bool ownInput; //! Does this object "own" the input data pointers (and hence have to free them)
    std::vector<string> * m_starSpectraStatsHeader; //! Header to contain statistics names
    std::vector<vector<float> > * m_starSpectraStats; //! Vector of vectors containing peptide statistics
    std::vector<string> * m_specsMsSpectraStatsHeader; //! The statistics that we are considering
    std::vector<vector<float> > * m_specsMsSpectraStats; //! Vector of vectors containing peptide statistics
    std::vector<string> * m_specsScoredSpectraStatsHeader; //! The statistics that we are considering
    std::vector<vector<float> > * m_specsScoredSpectraStats; //! Vector of vectors containing peptide statistics

    SvmScaleParameterList * m_svmScaleParamsCharge1;
    SvmScaleParameterList * m_svmScaleParamsCharge2;
    SvmScaleParameterList * m_svmScaleParamsCharge3;
    SvmModel * m_svmModelCharge1;
    SvmModel * m_svmModelCharge2;
    SvmModel * m_svmModelCharge3;

    //statistics params
    SpectrumAnnotParameterList * m_msAnnotParams;
    SpectrumAnnotParameterList * m_starAnnotParams;
    SpectrumAnnotParameterList * m_specsScoredAnnotParams;

    bool ownModels; //! Does this object "own" the svm models or statistics parameters.
    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
    bool ownStats; //! Does this class own the statistics vectors. Intended for use with invokeWithStatistics.

    PeptideSpectrumMatchSet * m_peptideResults; //! The peptide results we're considering.

    PeptideSpectrumMatchSet * m_peptideOutputResults;

    void buildKeyMap(vector<vector<unsigned int> > &keyMapping,
                     SvmModel * currModel);

    void addKeys(vector<vector<unsigned int> > &keyMapping,
                 vector<string> * currHeader,
                 SvmModel * currModel,
                 int fileEnum);

    bool buildStatsVector(vector<vector<unsigned int> > &keyMapping, vector<
        float> &outputVector, int lineIndex);
  };

} // namespace specnets

#endif // __ExecSvmStatistics_H__
