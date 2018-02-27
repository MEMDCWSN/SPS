////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_CLIENT_H__
#define __REPORT_RENDERER_HTML_CLIENT_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtml.h"



////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererHtmlClient
//
//
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Class for rendering an HTML report with client functionality (pagination, filtering)
   */
class ReportRendererHtmlClient : public ReportRendererHtml {

 protected:

  int curr_table, curr_row, curr_col;


  //////////////////////////////////////////////////////////////////////////////
  // Reports pages rendering methods (specific pages)

  // Header page generation
  //virtual int  generateReportPageHeader(ReportGeneratorData &data);

  /*! \brief Render protein list page
   */
  virtual int  generateReportPageProteins(ReportGeneratorData &data);

  /*! \brief Render single protein page
   */
  virtual int  generateReportPageProtein(ReportGeneratorData &data);

  /*! \brief ender contig list page
   */
  virtual int  generateReportPageContigs(ReportGeneratorData &data);

  /*! \brief Render single contig page
   */
  virtual int  generateReportPageContig(ReportGeneratorData &data);

  /*! \brief Render single cluster page
   */
  virtual int  generateReportPageCluster(ReportGeneratorData &data);

  /*! \brief render spectra page
   */
  virtual int  generateReportPageSpectra(ReportGeneratorData &data);


  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  /*! \brief Render page prolog
   */
  virtual int renderProlog(ReportBase *table, ostream &outstream);
  // epilog
  /*! \brief Render page epilog
   */
  virtual int renderEpilog(ReportBase *table, ostream &outstream);
  // between tables
  /*! \brief Render section between tables
   */
  virtual int renderInterTable(ReportBase *table, ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table content renderers

  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
  /*! \brief render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
   */
  virtual int renderTable(ReportTableBase *table, ostream &outstream);
  // render a header row method. Renders a table header row.
  /*! \brief render a header row method. Renders a table header row.
   */
  virtual int renderTableHeaderRow(vector<string> &row, ostream &outstream);
  // render header cell method. Renders a header cell
  //virtual int renderTableHeaderCell(string &cell, ostream &outstream);
  // build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
  /*! \brief build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
   */
  virtual int renderTableRow(vector<string> &row, ostream &outstream);
  // build cell comon method. builds a specific cell based on ColumnType specifications and a row of data
  //virtual int renderTableCell(string &cell, ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  // render a header row comon method. Cycles thu all columnType cells and renders each one.
  //virtual int buildTableHeaderRow(ReportTableBase *table, vector<string> &renderedRow);
  // render header cell comon method. Renders all ColumnTypes
  /*! \brief render header cell comon method. Renders all ColumnTypes
   */
  virtual int buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &ss);


  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // build table common method. Cycles thru all the table rows and invokes renderTableRow() method to render a row
  //virtual int buildTable(ReportTableBase *table);
  // build a row comon method. Cycles thu all columnType cells and renderCell for each ColumnType item
  //virtual int buildTableRow(ReportTableBase *table, vector<string> *row, vector<string> &renderedRow);
  // build cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
  /*! \brief build cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
   */
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss);

  // Builders for the diferent cell types
  //virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, stringstream &ss);
  //virtual int buildCellString(ReportColumnTypeString *, vector<string> *row, stringstream &ss);
  //virtual int buildCellBox(ReportColumnTypeBox *, vector<string> *row, stringstream &ss);
  // Method to process 'multiple' table cells, to have several lines of data in the same table cell
  //virtual int buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, stringstream &ss);


  //////////////////////////////////////////////////////////////////////////////
  // Table exception renderers

  // exception for rendering main page
  //virtual int renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream);



 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportRendererHtmlClient()  {};
  /*! \brief destructor
   */
  ~ReportRendererHtmlClient() {};

  //virtual int  generateReport(ReportGeneratorData &data);


};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
