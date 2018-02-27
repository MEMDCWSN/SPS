///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_BASE_H__
#define __REPORT_RENDERER_BASE_H__
///////////////////////////////////////////////////////////////////////////////
#include "ReportBase.h"
#include "ReportData.h"
#include "ReportModuleFactory.h"

#include "ReportTableHeader.h"
#include "ReportTableProtein.h"
#include "ReportTableProteinCoverage.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
#include "ReportTableInputSpectra.h"

#include "ReportHeader.h"
#include "ReportProtein.h"
#include "ReportProteinCoverage.h"
#include "ReportProteinCoverageCsv.h"
#include "ReportInputSpectra.h"
#include "ReportContig.h"
#include "ReportCluster.h"

#include "spsFiles.h"
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
#define PAGE_INITIAL               10
#define PAGE_PROTEINS              20
#define PAGE_PROTEIN               30
#define PAGE_PROTEIN_COVERAGE      40
#define PAGE_CONTIGS               50
#define PAGE_CONTIG                60
#define PAGE_CLUSTER               70
#define PAGE_SPECTRA               80
#define PAGE_PROTEIN_COVERAGE_CSV  90
///////////////////////////////////////////////////////////////////////////////
// class ReportRendererBase
///////////////////////////////////////////////////////////////////////////////
  /*! \brief Report rendering base class
   */
class ReportRendererBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // SPS data files

  /*! \brief SPS data files structure
   */
  SpsFiles *spsFiles;


  //////////////////////////////////////////////////////////////////////////////
  // tables names
  /*! \brief table filename
   */
  string m_tableNameHeader;
  /*! \brief table filename
   */
  string m_tableNameProtein;
  /*! \brief table filename
   */
  string m_tableNameProteinCoverage;
  /*! \brief table filename
   */
  string m_tableNameContig;
  /*! \brief table filename
   */
  string m_tableNameCluster;
  /*! \brief table filename
   */
  string m_tableNameSpectra;

  // current row and col being parsed. This is used to process <row> and <col> tags.
  /*! \brief urrent row and col being parsed. This is used to process <row> and <col> tags.
   */
  int currentRow, currentCol;

  // verbose flag
  /*! \brief verbose flag
   */
  bool m_verbose;

  // rendered image
  /*! \brief rendered image
   */
  string m_image;

  // number of available CPUs. Defaults to 1
  /*! \brief number of available CPUs. Defaults to 1
   */
  int m_ncpu;

  // executables directory ; needed to call image generation modules
  string m_exeDir;
  string m_tablesDir;
  string m_reportDir;
  string m_projectDir;
  string m_projectDirRel;
  string m_serverLocation;
  string m_fastaFilename;
  string m_msFilename;
  string m_aasFilename;
  string m_aasFile;
  string m_aasFileDir;

  /*! \brief Number of cells per line in protein coverage pages
   */
  int    m_cellPerLine;

  /*! \brief Current level of display.
  This will be used to filter out element that fall below the current level.
   */
  int    m_displayLevel;

  /*! \brief True if the clustering layer doesn't exist in the report
   */
  bool   m_noClusters;

  /*! \brief If the report will allow to realign
   */
  bool   m_allowRealign;

  /*! \brief If the report is dynamic
   */
  bool   m_dynamic;

  /*! \brief If the report has a password
   */
  bool   m_pwd;

  // where to render table cells to
  /*! \brief where to render table cells to
   */
  vector<vector<string> > reportTableCells;

  // navigation bar used to render the pages
  /*! \brief navigation bar used to render the pages
   */
  string m_navBar;


  // build directory path by concatenation with given path
  /*! \brief build directory path by concatenation with given path
   */
  string composeFileName(const string &projectDir, const string &fileName);

  //////////////////////////////////////////////////////////////////////////////
  // verbose output

  void verboseOutput(ostream &os, const double &v, bool nl = true);
  void verboseOutput(ostream &os, const char *str, bool nl = true);
  void verboseOutput(ostream &os, const char *str, const char *val, bool nl = true);
  void verboseOutput(ostream &os, const char *str, const char *val, const char *term, bool nl = true);


  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog, epilogue and between tables section renderers

  // prolog
  /*! \brief page prolog
   */
  virtual int renderProlog(ReportBase *table, ostream &outstream)     {return OK;};
  // epilog
  /*! \brief page epilog
   */
  virtual int renderEpilog(ReportBase *table, ostream &outstream)     {return OK;};
  // between tables
  /*! \brief section between tables
   */
  virtual int renderInterTable(ReportBase *table, ostream &outstream) {return OK;};


  //////////////////////////////////////////////////////////////////////////////
  // Table content renderers

  // render table method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
  /*! \brief render table method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
   */
  virtual int renderTable(ReportTableBase *table, ostream &outstream);
  // render a header row method. Renders a table header row.
  /*! \brief render a header row method. Renders a table header row.
   */
  virtual int renderTableHeaderRow(vector<string> &row, ostream &outstream);
  // render header cell method. Renders a header cell
  /*! \brief render header cell method. Renders a header cell
   */
  virtual int renderTableHeaderCell(string &cell, ostream &outstream);
  // render a row method. Renders a table row
  /*! \brief render a row method. Renders a table row
   */
  virtual int renderTableRow(vector<string> &row, ostream &outstream);
  // render cell method. Renders a cell
  /*! \brief render cell method. Renders a cell
   */
  virtual int renderTableCell(string &cell, ostream &outstream);

  // transpose a row for vertical rendering
  /*! \brief transpose a row for vertical rendering
   */
  virtual void transposeRow(int col, vector<string> &row);
  // render table vertically
  /*! \brief render table vertically
   */
  virtual int renderTableVertically(ostream &outstream);


  ////////////////////////////////////////////////////////////////////////////////
  // Header section rendering (top of page)

  // Header renderer entry point
  /*! \brief Header renderer entry point
   */
  virtual int buildElement(ReportElementsBase *elem, ostream &outstream);
  // builders for the different element types
  /*! \brief builders for the different element types: DIV
   */
  virtual int buildElementDiv(ReportElementsDiv *div, ostream &outstream)        {return OK;};
  /*! \brief builders for the different element types: table
   */
  virtual int buildElementTable(ReportElementsTable *tab, ostream &outstream)    {return OK;};
  /*! \brief builders for the different element types: row
   */
  virtual int buildElementRow(ReportElementsRow *row, ostream &outstream)        {return OK;};
  /*! \brief builders for the different element types: cell
   */
  virtual int buildElementCell(ReportElementsCell *cell, ostream &outstream)     {return OK;};
  /*! \brief builders for the different element types: string
   */
  virtual int buildElementString(ReportElementsString *str, ostream &outstream)  {return OK;};


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  // build a header row comon method. Cycles thu all columnType cells and builds each one.
  /*! \brief build a header row comon method. Cycles thu all columnType cells and builds each one.
   */
  virtual int buildTableHeaderRow(ReportTableBase *table, vector<string> &);
  // build header cell comon method. builds all ColumnTypes
  /*! \brief build header cell comon method. builds all ColumnTypes
   */
  virtual int buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &) {return OK;};


  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // build table comon method. Cycles thru all the table rows and invokes buildTableRow() method to build a row
  /*! \brief table comon method. Cycles thru all the table rows and invokes buildTableRow() method to build a row
   */
  virtual int buildTable(ReportTableBase *table, int, int);
  // build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
  /*! \brief build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
   */
  virtual int buildTableRow(ReportTableBase *table, vector<string> *row, vector<string> &);
  // build cell comon method. builds a specific cell based on ColumnType specifications and a row of data
  /*! \brief build cell comon method. builds a specific cell based on ColumnType specifications and a row of data
   */
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &);

  // builders for the diferent cell types
  /*! \brief builders for the diferent cell types: Image on demand
   */
  virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, stringstream &)   {return OK;};
  /*! \brief builders for the diferent cell types: string
   */
  virtual int buildCellString(ReportColumnTypeString *, vector<string> *row, stringstream &)                 {return OK;};
  /*! \brief builders for the diferent cell types: column box
   */
  virtual int buildCellBox(ReportColumnTypeBox *, vector<string> *row, stringstream &)                       {return OK;};
  /*! \brief builders for the diferent cell types: multiple strings
   */
  virtual int buildCellStringMultiple(ReportColumnTypeStringMultiple *, vector<string> *row, stringstream &) {return OK;};

  // parse the params vector
  /*! \brief parse the params vector
   */
  virtual void parseParamsVector(vector<string> &params, vector<ReportParamsOption> &options, vector<string> *row);
  // parse a template string and return the parsed string (wrapper)
  /*! \brief parse a template string and return the parsed string (wrapper)
   */
  virtual string parseTemplatesAll(const string &str, vector<string> *row,
          char primaryOpen = '<', char primaryClose = '>',
          char secundaryOpen = '«', char secundaryClose = '»',
          int add = 1);
  // parse a template string and return the parsed string (specify delimiters)
  /*! \brief parse a template string and return the parsed string (specify delimiters)
   */
  virtual string parseTemplates(const string &str, vector<string> *row, char tagOpen = '<' ,char tagClose = '>', int add = 0);
  // Helper method to translate the tag. If 'add' is different from 0, it tries to convert to an integer and add the 'add' value.
  /*! \brief Helper method to isolate the tag contents.
   */
  virtual string getTag(const string &tag, vector<string> *row, char tagOpen = '<', char tagClose = '>', int add = 0);
  /*! \brief Helper method to translate the tag. If 'add' is different from 0, it tries to convert to an integer and add the 'add' value.
   */
  virtual string translateTag(const string &tag, vector<string> *row, int add = 0, bool preserveEmpty = false);
  // Helper method to split text into chunks
  /*! \brief method to split text into chunks
   */
  virtual string splitText(const string &sequence, const char *sep);
  // Helper method to split sequences into chunks
  /*! \brief Helper method to split sequences into chunks
   */
  virtual string splitSequence(const string &sequence, const char *sep);


  //////////////////////////////////////////////////////////////////////////////
  // Image renderers

  // method to generate images
  /*! \brief method to generate images
   */
  virtual int renderImage(const string &object, vector<string> &params, vector<ReportParamsFiles> &files) {return OK;};


 public:


  // Constructors and destructor
  ReportRendererBase();
  ~ReportRendererBase();

  // set data files
  /*! \brief set data files
   */
  virtual void setSpsFiles(SpsFiles *f) {spsFiles = f;};

 // render report. Cycles thru all the report tables and renders them, by invoing renderTable()
  //virtual int renderReport(void);
  /*! \brief render report. Cycles thru all the report tables and renders them, by invoing renderTable()
   */
  virtual int renderReport(ReportBase *, ostream &outstream, string &navBar, bool renderDecoration = true, int start = 0, int count = -1);
  // render report with a pagination ruler
  /*! \brief render report with a pagination ruler
   */
  virtual void setNavigationBar(string &navBar) {m_navBar = navBar;};
  /*! \brief Clear the navigation bar temp var
   */
  virtual void clearNavigationBar(void)         {m_navBar.clear();};


  /*! \brief Set exe location
   */
  virtual void setExeDir(const string &ed)          {m_exeDir         = ed;};
  /*! \brief Set tables location
   */
  virtual void setTablesDir(const string &pd)       {m_tablesDir      = pd;};
  /*! \brief set report location
   */
  virtual void setReportDir(const string &pd)       {m_reportDir      = pd;};
  /*! \brief set project location
   */
  virtual void setProjectDir(const string &pd)      {m_projectDir     = pd;};
  /*! \brief Set ptoject location to relative
   */
  virtual void setProjectDirRel(const string &pd)   {m_projectDirRel  = pd;};
  /*! \brief Set the number of cells per line in the protein coverage page
   */
  virtual void setCellsPerLine(const int cpl)       {m_cellPerLine    = cpl;};
  /*! \brief Set server location string. Server location refers to the server URL used in dynamic reports
   */
  virtual void setServerLocation(const string &sl)  {m_serverLocation = sl;};
  /*! \brief Fasta file name
   */
  virtual void setFastaFilename(const string &s)    {m_fastaFilename  = s;};
  /*! \brief Set display level
   */
  virtual void setDisplayLevel(const int dl)        {m_displayLevel   = dl;};
  /*! \brief Set number of CPUs to be used. This translates to simulataneous parallel processes when generating static reports.
   */
  virtual void setCPUs(const int cpu)               {m_ncpu           = cpu;};
  /*! \brief Enable or disable the clustering layer
   */
  virtual void setClusterLayerState(const bool nc)  {m_noClusters     = nc;};
  /*! \brief Set MS file name
   */
  virtual void setMsFilename(const string &ms)      {m_msFilename     = ms;};
  /*! \brief If realignment is allowed from the protein coverage page. Only relevant for dynamic reports.
   */
  virtual void setAllowRealign(const bool &ar)      {m_allowRealign   = ar;};
  /*! \brief Set reports as dynamic
   */
  virtual void setDynamic(const bool &d)            {m_dynamic        = d;};
  /*! \brief Set AAs file name
   */
  virtual void setAasFile(const string &aa)         {m_aasFile        = aa;};
  /*! \brief Set AAs file location
   */
  virtual void setAasFileDir(const string &ad)      {m_aasFileDir     = ad;};
  /*! \brief Set both file name and location of AAs file.
   */
  virtual void setAasFilename(void)                 {m_aasFilename=m_aasFileDir;if(m_aasFilename.length()) if(m_aasFilename[m_aasFilename.length()-1] !=  '\\' &&  m_aasFilename[m_aasFilename.length()-1] !=  '/') m_aasFilename+='/';m_aasFilename+=m_aasFile;};
  /*! \brief Set reports password. Only relevant for dynamic reports
   */
  virtual void setPwd(bool p)                       {m_pwd = p;};

  /*! \brief Builds the navigation bar
   */
  virtual void buildNavigationBar(string &navBar, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames) {};

  /*! \brief Report generation entry point
   */
  virtual int generateReport(ReportGeneratorData &data) {};

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
