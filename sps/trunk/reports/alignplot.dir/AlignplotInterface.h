///////////////////////////////////////////////////////////////////////////////
#ifndef __ALIGNPLOT_INTERFACE_H__
#define __ALIGNPLOT_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "spectrum.h"
#include "SpecSet.h"
#include "PlotContig.h"
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
  /*! \brief Interface class for Alignplot tool

   Provides command line user interface 

   */
class AlignplotInterface {

  // abruijn graph
  abinfo_t  abinfo;
  // star spectra
  SpecSet   star;
  // seqs file
  SpecSet   seqs;


  // Spectrum draw object -- used to draw image
  PlotContig  plotContig;

  // Spectrum index.
  int m_spectrumIndex;

  // spectrum scan
  string m_spectrumScan;


 public:


  // Constructors and destructor
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Generally this constructor should not be used. It is only used by the
     module execution factory in order to create an exemplar object (without
     valid parameters) which is then used to create a real (valid) object
     using the clone() method.
     @sa clone()
     */
  AlignplotInterface();
    //@}

    //! \name DESTRUCTOR
    //@{
  ~AlignplotInterface();
    //@}

  // Option parsing
    /*! \brief Processes the command line params.

     Initializes the default parameter values and parses the command line against the declared parameters.

     @return 0 if successful; -1 if there was an error
     */
  int processOptions(int argc, char **argv);

  // Execution based on options
    /*! \brief Processes the options.

     Top level program flow. Program segments and functions are called based on input options.

     @return 0 if successful; -1 if there was an error
     */
  int processOptions(ParameterList & commandLineParams);


  // Output help
    /*! \brief Help().

     Outputs program options.

     @return 0 if successful; -1 if there was an error
     */
  int help(ostream &);

  // Output version information
    /*! \brief version().

     Outputs program version information.

     @return 0 if successful; -1 if there was an error
     */
  int version(ostream &);

  // Output error messages
    /*! \brief display error message.

     Outputs an error message and exits program.

     @return 0 if successful; -1 if there was an error
     */
  int error(const string &);

  // Spectrum drawing method
    /*! \brief Plot main routine.

     Call the drawing top level routine.

     @return 0 if successful; -1 if there was an error
     */
  int plot(void);


  void getData(string & data) {data = plotContig.getImage();};


};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
