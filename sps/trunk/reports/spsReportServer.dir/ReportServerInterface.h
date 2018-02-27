///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_SERVER_INTERFACE_H__
#define __REPORT_SERVER_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "ParameterList.h"
#include "ReportTableBase.h"

#include "ReportData.h"

///////////////////////////////////////////////////////////////////////////////
using namespace std;
///////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
// Defines
#define STR(s) #s
#define XSTR(s) STR(s)

#define DEFAULT_ROWS_PER_TABLE  20
///////////////////////////////////////////////////////////////////////////////
class ReportServerInterface {

  // tables names
  /*! \brief tables names
   */
  string m_tableNameHeader;
  /*! \brief tables names
   */
  string m_tableNameProtein;
  /*! \brief tables names
   */
  string m_tableNameProteinCoverage;
  /*! \brief tables names
   */
  string m_tableNameContig;
  /*! \brief tables names
   */
  string m_tableNameCluster;
  /*! \brief tables names
   */
  string m_tableNameSpectra;
  /*! \brief tables names
   */
  string m_tableConfig;

  // build directory path by concatenation with given path
  /*! \brief build directory path by concatenation with given path
   */
  string composeFileName(const string &projectDir, const string &fileName);

  // get a table object
  /*! \brief get a table object
   */
  ReportTableBase *getTableObject(ReportData &data);


 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportServerInterface();
  /*! \brief destructor
   */
  ~ReportServerInterface();

  // Option parsing
  /*! \brief Option parsing
   */
  int parseOptions(int argc, char **argv);
  // Option processing
  /*! \brief Option processing
   */
  int processOptions(specnets::ParameterList &commandLineParams);

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
  int updateTable(ReportData &data);

  // Spectrum drawing method
  /*! \brief Spectrum drawing method
   */
  int getTableData(ReportData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
