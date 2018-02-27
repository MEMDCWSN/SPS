////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_SINGLE_H__
#define __REPORT_RENDERER_HTML_SINGLE_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtml.h"



////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererHtmlSingle
//
//
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Class to render a report in a single HTML page
   */
class ReportRendererHtmlSingle : public ReportRendererHtml {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Reports pages rendering methods (specific pages)

  /*! \brief Renders a table
   */
  virtual int renderTable(ReportTableBase *table, ostream &outstream);

  /*! \brief Renders a table row
   */
  virtual int renderTableRow(vector<string> &row, ostream &outstream);

  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  /*! \brief Renders page prolog
   */
  virtual int renderProlog(ostream &outstream);
  /*! \brief Renders page prolog
   */
  virtual int renderProlog2(ostream &outstream);
  // epilog
  /*! \brief Renders page epilog
   */
  virtual int renderEpilog(ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // main page
  /*! \brief Entry point for main page rendering
   */
  virtual int renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream);

  // table cell renderer
  /*! \brief Renders table cell
   */
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss);

  // Builders for the diferent cell types
  /*! \brief Renders image on demand cell
   */
  virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, stringstream &ss);

  // write the table
  /*! \brief Writes the table to file
   */
  virtual int writeTable(ReportTableBase *table, ostream &outstream, int idx);


 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportRendererHtmlSingle()  {};
  /*! \brief destructor
   */
  ~ReportRendererHtmlSingle() {};

  /*! \brief Report generation entry point
   */
  virtual int  generateReport(ReportGeneratorData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
