///////////////////////////////////////////////////////////////////////////////
#ifndef __PAIRS_INFO_INTERFACE_H__
#define __PAIRS_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"
#include "SpectrumPairSet.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
  /*! \brief PairsInfo tool interface class
   */
class PairsInfoInterface {


  // cluster data object
  /*! \brief cluster data object
   */
  SpectrumPairSet pairSet;

  // Input data
  /*! \brief Input data
   */
  string inputFilename, projectdir;

  // output file and directory
  /*! \brief output file and directory
   */
  string outdir, outFileName;


  float m_minScore1, m_minScore2, m_edgeTopKBoth, m_edgeTopKOne;
  int m_component_size;
  bool  m_exists_minScore1, m_exists_minScore2, m_exists_edgeTopKBoth, m_exists_edgeTopKOne, m_exists_component_size;


  /*! \brief Writes the output file
   */
  int writeOutFile(void);
  
  


 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  PairsInfoInterface();
  /*! \brief destructor
   */
  ~PairsInfoInterface();

  // Option parsing
  /*! \brief Option parsing
   */
  int processOptions(int argc, char **argv);

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

};

///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
