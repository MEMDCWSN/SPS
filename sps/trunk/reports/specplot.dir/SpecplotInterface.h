///////////////////////////////////////////////////////////////////////////////
#ifndef __SPECPLOT_INTERFACE_H__
#define __SPECPLOT_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "spectrum.h"
#include "SpecSet.h"
#include "PlotSpectrum.h"
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
  /*! \brief Specplot interface class
   */
class SpecplotInterface {

 protected:

  // Spectrum draw object -- used to draw image
  /*! \brief Spectrum draw object -- used to draw image
   */
  PlotSpectrum  plotSpectrum;

  // Specset object - used to load specset from file
  /*! \brief Specset object - used to load specset from file
   */
  SpecSet *specSet;

  // Is the specset object owned by specplot?
  /*! \brief Is the specset object owned by specplot?
   */
  bool specSet_own;

  // Spectrum index.
  /*! \brief Spectrum index.
   */
  int m_spectrumIndex;
  
  int m_origspectrumIndex;
  
  bool m_loadFull;
  
  bool m_verbose;
  
  int m_useindex;
  
  /*! \brief Spectrum ID.
   */
  string m_spectrumID;

  // spectrum scan
  /*! \brief spectrum scan
   */
  string m_spectrumScan;
  
  
  /*! \brief Filename of input spectra
   */
  string m_inputSpectraFilename;

  /*! \brief Filename of input spectra
   */
  string m_outputImageFilename;


 public:


  // Constructors and destructor
  //! \name CONSTRUCTORS
  //@{
  SpecplotInterface();
  //@}

  //! \name DESTRUCTOR
  //@{
  ~SpecplotInterface();
  //@}

  // Option parsing
  /*! \brief Option parsing
   */
  int processOptions(int argc, char **argv);

  // Execution based on options
  /*! \brief Execution based on options
   */
  int processOptions(ParameterList & commandLineParams);

  // load files
  /*! \brief load files
   */
  virtual int load(string &fn, string &ext, bool spectrumScanPresent, bool spectrumIdPresent, bool usePwizFirst);

  // Output help
  /*! \brief Output help
   */
  int help(ostream &);

  // Output version information
  /*! \brief Output version information
   */
  int version(ostream &);

  // Output error messages
  /*! \brief Output error messages
   */
  int error(const string &);

  // Spectrum drawing method
  /*! \brief Spectrum drawing method
   */
  int plot(void);


  /*! \brief Gets the image into internal structures
   */
  void getData(string & data) {data = plotSpectrum.getImage();};

  void setData(int type, void *data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
