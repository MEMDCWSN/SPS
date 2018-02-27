#ifndef _CommandLineParser_H_
#define _CommandLineParser_H_

// Module Includes
#include "ParameterList.h"

// External Module Includes

// System Includes
#include <vector>

namespace specnets
{
  class CommandLineParser
  {
  public:
    /*! \brief Structure to hold description of parameters.

     This structure is used to hold the description of a single option
     that would be present on the command line.
     Example #1:
     If one of the command line options was "-f filename", the
     structure representation would be:
     option.flag = "f";            The option flag is "f"
     option.name = "filename";     The name of the option is "filename"
     option.numValues = 1;         This option requires 1 value (the filename)
     */
    struct Option
    {
      Option(const std::string & flag,
             const std::string & name,
             unsigned int numValues) :
        m_flag(flag), m_name(name), m_numValues(numValues)
      {
        // EMPTY
      }

      std::string m_flag; // The command line flag
      std::string m_name; // The option name
      unsigned int m_numValues; // How many parameters does this option require?
    };

    /*! \brief This is the default constructor.

     The construct parses the command line parameters given in the first two
     parameters according to allowable options given in the listOptions. The
     first N parameters are skipped according to the numNonOptions argument.

     @param argc Number of arguments
     @param argv Array of argument strings
     @param numNonOptions Number of arguments that are not options (does not include command)
     @param listOptions  vector containing 
     */
    CommandLineParser(int argc,
                      char ** argv,
                      int numNonOptions,
                      const std::vector<Option> & listOptions);
    virtual ~CommandLineParser();

    /*! \brief Get the values of the parameters as a ParameterList

     Returns the list of actual parameters found as a ParameterList structure.
     Options that are found will be set as parameters with their values
     (if they ave values) will be be set as the value of the parameter. For
     options that do not require values, the value will be set as an empty
     string but can be checked using the exists() method of ParameterList. 
     
     @param params Returned list of option values
     */
    void getOptionsAsParameterList(ParameterList & params) const;

    bool validate(std::string & error) const;

  private:
    std::string m_error;
    std::vector<std::pair<std::string, std::string> > m_listOptions;
  };

} //namespace specnets

#endif // _CommandLineParser_H_
