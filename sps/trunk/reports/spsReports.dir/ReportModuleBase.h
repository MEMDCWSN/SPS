////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_MODULE_BASE_H__
#define __REPORT_MODULE_BASE_H__
////////////////////////////////////////////////////////////////////////////////

// External Includes
#include "ParameterList.h"

// System Includes
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

#define VALIDATE_PARAM_EXIST(value)  \
    if (!m_params.exists(value))                \
    {                                         \
      error = value;                          \
      error += " not present";                \
      return false;                           \
    }


////////////////////////////////////////////////////////////////////////////////
using namespace specnets;
using namespace std;
////////////////////////////////////////////////////////////////////////////////

namespace spsReports {
  /*! \brief Base class for all execution framework modules.

   Declares a set of virtual functions required to execute an Specnets module
   within the specnets execution framework. Implements some very basic
   functions such as the isValid() method and holds the parameters structure
   that all modules will require as well as the member data necessary for
   splitting into sub-modules and re-merging the results.

   An execution module of this type may be used in two basic ways:<br>
   1) The module is created using the constructor that takes an InputParams
   parameter and then run using the invoke() method.<br>
   2) The module is created using the constructor that takes an InputParams
   parameter and then one of the classes derived from ParallelExecution is
   created to split the module into sub-modules which will be run in parallel.
   */
  class ReportModuleBase
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
    ReportModuleBase(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ReportModuleBase(const ParameterList & inputParams);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ReportModuleBase(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ReportModuleBase * clone(const ParameterList & input_params) const = 0;

    virtual ReportModuleBase * clone(void) const = 0;

    /*! \brief Returns the class name of the module.

     By default the name of the object is same as the class itself, however it may
     be set to a different value using the setName() method. Or in the case of
     modules that are children or other modules, the name may be a variation of
     the parent's name.

     @sa setName()
     @return The name of the class
     */
    std::string getName(void) const;

    /*! \brief Returns the type name of the module.

     The type name is automatically set by the constructor to be the same name
     as the class. This method may be used to determine the true class from
     a pointer to the base (this) class.

     @return The type name of the class
     */
    std::string getType(void) const;

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
    virtual bool invoke(void) = 0;


    virtual bool invoke(int argc, char ** argv) {};


    virtual bool getData(string &data) {};


    virtual bool setData(int type, void *data) {};


    /*! \brief Loads the input data from the files specified in the params.

     Loads all the data from files into the data members. This method is
     primarily used by the execution module to load necessary data when
     executing in a separate process..

     @return True if data was loaded successfully, false otherwise.
     @sa ReportModuleBase(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void) = 0;

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ReportModuleBase(const ParameterList & input_params), loadInputData(), merge()
     */
    virtual bool saveOutputData(void) = 0;


    /*! \brief Sets the name of the module.

     By default the name of the object is same as the class itself, however this
     method may be used to set the name to a different value.

     @sa getName()
     @return The name of the class
     */
    void setName(std::string name);

    /*! \brief Performs validation of the input parameters.

     Checks the parameters structure provided in the constructor to see if they
     are sufficient and correct to invoke the module. Also sets the internal
     validity flag so that isValid() will return the correct result.

     @param error A description of the error (if any occurs)
     @return True if the parameters for the module are valid, false otherwise.
     @sa isValid()
     */
    virtual bool validateParams(std::string & error) = 0;
    //@}

    ParameterList m_params;

  protected:
    /*! \brief Creates a unique file name for a sub_module

     Convenience method that creates a filename based on the base name given,
     along with the number of the split.

     @param baseFilename the base name of the file on which will be appended the
     numSplit
     @return The full filename
     */
    std::string makeName(std::string baseFilename, int numSplit) const;

    bool m_isValid;
    std::string m_name;
    std::string m_type;

    int m_numNodes;
    int m_numCpus;
    std::vector<ReportModuleBase *> m_subModules;

  private:
    //! \name Constructor
    //@{
    ReportModuleBase(const ReportModuleBase & that);
    //@}
  };

////////////////////////////////////////////////////////////////////////////////
} // namespace specnets
////////////////////////////////////////////////////////////////////////////////
#endif // __ReportModuleBase_H__
////////////////////////////////////////////////////////////////////////////////
