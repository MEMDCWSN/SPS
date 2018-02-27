///////////////////////////////////////////////////////////////////////////////
#ifndef __CONTPLOT_INTERFACE_H__
#define __CONTPLOT_INTERFACE_H__
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
  /*! \brief Interface class for Contplot tool

   Provides command line user interface 

   */
class ContplotInterface {

  // abruijn graph
  abinfo_t  *abinfo;
  // star spectra
  SpecSet   *star;
  // seqs file
  SpecSet   *seqs;

  bool abinfo_own, star_own, seqs_own;


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

     Default contructor
     */
  ContplotInterface();
    //@}

    //! \name DESTRUCTOR
    //@{
  ~ContplotInterface();
    //@}

  // Option parsing
    /*! \brief Processes the command line params.

     Initializes the default parameter values and parses the command line against the declared parameters.

     @return 0 if successful; -1 if there was an error
     */
  int processOptions(int argc, char **argv);

  // Execution based on options
    /*! \brief Processes the params from the params object.

     Initializes internal objects based on parameters.

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
  int plot(void);


  void getData(string & data) {data = plotContig.getImage();};

  void setData(int type, void *data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
