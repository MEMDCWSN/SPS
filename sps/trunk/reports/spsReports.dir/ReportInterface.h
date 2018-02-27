///////////////////////////////////////////////////////////////////////////////
#ifndef __SPS_PLOT_H__
#define __SPS_PLOT_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "ReportTableGenerator.h"
#include "ParameterList.h"
#include "ReportBase.h"
#include "ReportRendererBase.h"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
///////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
// Defines
#define STR(s) #s
#define XSTR(s) STR(s)

///////////////////////////////////////////////////////////////////////////////
class ReportInterface {

  // verbose flag
  /*! \brief verbose flag
   */
  bool m_verbose;

  // tables names
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
  
  // configuratin file name
  /*! \brief configuratin filename
   */
  string m_tableConfig;

  // build directory path by concatenation with given path
  /*! \brief build directory path by concatenation with given path
   */
  string composeFileName(const string &projectDir, const string &fileName);

  // build directory path by concatenation with path
  /*! \brief build directory path by concatenation with path
   */
  int buildDirectoryPath(const string &dir);

  void dump_abruijn(ReportGeneratorData &data);
  void dump_binArray(ReportGeneratorData &data);


 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportInterface();
  /*! \brief destructor
   */
  ~ReportInterface();

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
  /*! \brief Tables building method
   */
  int buildTables(ReportGeneratorData &data);

  // Generate HTML report pages
  /*! \brief Generate HTML report pages
   */
  int generateReportHtml(ReportGeneratorData &data);

  // Generate HTML report pages with pagination on client
  /*! \brief Generate HTML report pages with pagination on client
   */
  int generateReportHtmlClient(ReportGeneratorData &data);

  // Generate HTML report in a single page
  /*! \brief Generate HTML report in a single page
   */
  int generateReportHtmlSingle(ReportGeneratorData &data);

  // Generate dynamic HTML report entry page
  /*! \brief Generate dynamic HTML report entry page
   */
  int generateReportHtmlDynamic(ReportGeneratorData &data);

  // Generate PDF report pages
  /*! \brief Generate PDF report pages
   */
  int generateReportPdf(ReportGeneratorData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
