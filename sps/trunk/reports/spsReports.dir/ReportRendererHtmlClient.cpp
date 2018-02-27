////////////////////////////////////////////////////////////////////////////////

#include "ReportRendererHtmlClient.h"
#include "Tokenizer.h"
#include "copyright.h"
#include "Timer.h"
#include "ReportDefines.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageProteins(ReportGeneratorData &data)
{
  // define report protein list object, and load tables
  ReportProtein pl(data.tablesDir, m_tableNameProtein);
  // load table
  pl.load();
  // generate pages with pagination
  paginate(data, pl, -1, "proteins.", "proteins", 10);
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageProtein(ReportGeneratorData &data)
{
  // define specific protein report object
  ReportProtein p2(data.tablesDir, m_tableNameProtein, m_tableNameContig);
  // load table
  p2.load();
  // get the proteins column from the proteins table
  vector<string> aux = p2.getTableColumn(0, TABLE_PROTEIN_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    p2.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "protein.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, p2, -1, fnamePrefix, "protein", 11);
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageContigs(ReportGeneratorData &data)
{
  // define report object for contig list
  ReportContig cl(data.tablesDir, m_tableNameContig);
  // load table
  cl.load();
  // generate pages with pagination
  paginate(data, cl, -1, "contigs.", "contigs", 20);
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageContig(ReportGeneratorData &data)
{
  // Child (specific) contig report pages
  clearNavigationBar();
  // define specific contig report object
  ReportContig rc(data.tablesDir, m_tableNameContig, m_tableNameCluster, m_tableNameSpectra, data.noClusters);
  // load table
  rc.load();
  // get the contig id column from the contigs table
  vector<string> aux = rc.getTableColumn(0, TABLE_CONTIG_FILTER_COL_CONTIG);
  // cycle tru all contigs
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster/spectra table
    rc.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "contig.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rc, -1, fnamePrefix, "contig", 21);
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageCluster(ReportGeneratorData &data)
{
  if(data.noClusters)
    return OK;
  // cluster pages are generated depending on flag
  if(!data.noClusters) {
    clearNavigationBar();
    // define report cluster object, and load tables
    ReportCluster clust(data.tablesDir, m_tableNameCluster, m_tableNameSpectra);
    // load tables
    clust.load();
    // get the cluster id column from the cluster table
    vector<string> aux = clust.getTableColumn(0, TABLE_CLUSTER_FILTER_COL_CLUSTER);
    // cycle tru all clusters
    for(int i = 0 ; i < aux.size() ; i++) {
      // set the filter for a specific contig in cluster table
      clust.applyFilter(aux[i]);
      // set filename
      string fnamePrefix = "cluster.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
      // generate pages with pagination
      paginate(data, clust, -1, fnamePrefix, "cluster", 31);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlClient::generateReportPageSpectra(ReportGeneratorData &data)
{
  clearNavigationBar();
  // define report cluster object, and load tables
  ReportInputSpectra is(data.tablesDir, m_tableNameSpectra, data.noClusters);
  // load tables
  is.load();
  // get the file index column from the input spectra table
  vector<string> aux = is.getTableColumn(0, TABLE_SPECTRA_FILTER_COL_FILE);
  // sort elements for duplicate removal
  sort(aux.begin(), aux.end(), stringSortCompare);
  // duplicate find
  vector<string>::iterator it;
  // using default comparison
  it = unique (aux.begin(), aux.end(), stringUniqueCompare);
  // remove extra items
  aux.resize( it - aux.begin() );
  // cycle tru all input files
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    is.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "spectra.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, is, -1, fnamePrefix, "spectra", 40);
  }
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererHtmlClient::renderProlog(ReportBase *table, ostream &outstream)
{
  outstream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n";
  outstream << "  \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
  outstream << "<HTML xmlns=\"http://www.w3.org/1999/xhtml\">\n";

  outstream << "<head>\n";
  outstream << "  <meta http-equiv=\"Content-Type\" content=\"text/shtml; charset=ISO-8859-1\" />\n";
  outstream << "  <title>" << REPORT_TITLE << "</title>\n";

  // sps reports Client
  outstream << "<script type=\"text/javascript\" language=\"javascript\" src=\"js/spsReportClient.js\"></script>";
  // sps stylles
  outstream << "<link href=\"styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n";
  // icon
  outstream << "<link rel=\"shortcut icon\" href=\"images/favicon.ico\" type=\"image/icon\" />\n";

  // lightbox items
  outstream << "<script type='text/javascript' src='js/prototype.js'></script>\n";
  outstream << "<script type='text/javascript' src='js/scriptaculous.js?load=effects,builder'></script>\n";
  outstream << "<script type='text/javascript' src='js/lightbox.js'></script>\n";
  outstream << "<link rel='stylesheet' href='css/lightbox.css' type='text/css' media='screen' />\n";

  // sorttable items
  //outstream << "<script src=\"js/sorttable.js\"></script>\n";

  outstream << "</head>\n";

  // body section
  outstream << "<body>\n";

  outstream << "<div class='wrapper' align='center'>";

  // Logo image
  //outstream << "<div><h3 align=center><img src='data:image/jpg;base64," << ICON_LOGO <<" ' /></h3></div><br /><br />";
  outstream << "<div><h3 align=center><img src='images/logo.jpg' /></h3></div><br /><br />";

  // data containers
  outstream << "<div id='bodyWrapper'>";
  outstream << "<div id='textWrapper'>";

  // language initializer section
  outstream << "<script language='javascript'>";
  outstream << "var D = [];var DH = [];";
  outstream << "var T=0;var R=0;var C=0;var RH=0;var CH=0;";
  outstream << "</script>";

  curr_table = 0;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtmlClient::renderEpilog(ReportBase *table, ostream &outstream)
{
  // language section - epilog - call for code initializers
  outstream << "<script language='javascript'>";
  outstream << "if(window.init) init(D, DH);";
  outstream << "</script>";

  outstream << "</div></div>";

  // page footers
  outstream << "<div class='push'></div>";
  outstream << "</div>";

  outstream << "<div class='footer'><table width='100%'><tr><td class='VHSep'></td></tr><tr><td class='HSep'></td><td class='ln'></td><td class='HSep'></td></tr><tr><td class='HSep'></td><td class='Footer'>";
  outstream << REPORT_FOOTER_TEXT;
  outstream << "</td></tr></table></div>";

  // body end
  outstream << "</body>\n";
  // html end
  outstream << "</html>\n";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// between tables
int ReportRendererHtmlClient::renderInterTable(ReportBase *table, ostream &outstream)
{
  // output a line break in a table row
  outstream << "<div> </div>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererHtmlClient::renderTable(ReportTableBase *table, ostream &outstream)
{
  switch(table->getRenderingException()) {

  case RENDERING_EXCEPTION_NONE:
    break;

  case RENDERING_EXCEPTION_PROTEIN_HEADER:
    return renderTableExceptionProteinHeader(table, outstream);
    break;

  case RENDERING_EXCEPTION_MAIN_PAGE:
    return renderTableExceptionMainPage(table, outstream);

  case RENDERING_EXCEPTION_PROTEIN_COVERAGE_PAGE:
    return renderTableExceptionProteinCoverage(table, outstream);

  case RENDERING_EXCEPTION_PROTEIN_COVERAGE_CSV_PAGE:
    return renderTableExceptionProteinCoverageCSV(table, outstream);

  default:
    break;
  }

  // HTML class definition
  string aux = "class='result sortable'";
  // boder definition
  if(!table->doBorders())
    aux = "border='0'";

  // holder for top navigation bar
  outstream << "<div id='bar_top_" << curr_table << "'></div>";
  // HTML table start
  outstream << "<table " << aux << " align='center' id='tab_" << curr_table << "'></table>";
  // holder for bottom navigation bar
  outstream << "<div id='bar_bot_" << curr_table << "'></div>";
  // language section
  outstream << "<script language='javascript'>";
  // initialize vars
  outstream << "R=0;C=0;RH=0;CH=0;";
  outstream << "D[T]=[];DH[T]=[];";
  // render contents
  ReportRendererBase::renderTable(table, outstream);
  // HTML table terminator
  outstream << "T++;";
  // end language section
  outstream << "</script>";

  curr_table++;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header row method. Renders a table header row.
int ReportRendererHtmlClient::renderTableHeaderRow(vector<string> &row, ostream &outstream)
{

  // render the header
  outstream << "DH[T][RH]=[];";
  ReportRendererBase::renderTableHeaderRow(row, outstream);
  outstream << "RH++;";

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererHtmlClient::renderTableRow(vector<string> &row, ostream &outstream)
{
  outstream << "D[T][R]=[];";
  ReportRendererBase::renderTableRow(row, outstream);
  outstream << "R++;C=0;";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// build header cell comon method. Renders all ColumnTypes
int ReportRendererHtmlClient::buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &ss)
{
  // cell begin
  ss << "DH[T][RH][CH++]=\"";
  // cell content
  ss << base->columnLabel;
  // cell end HTML tag
  ss << "\";";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
int ReportRendererHtmlClient::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss)
{
  // don't render cells with dynamic content
  if(base->dynamic)
    return OK;

  // check if this cell is a StringMultiple. If it is, we shouldn't issue <a> tags here.
  ReportColumnTypeStringMultiple *auxM = dynamic_cast<ReportColumnTypeStringMultiple*>(base);

  // auxiliary variables
  stringstream cls, link;

  // gather needed attributes
  if(base->link.size() && (auxM == NULL))
    link << parseTemplates(base->link, row); // parseTemplatesAll

  if(base->cssClass.size())
    cls << " class='" << base->cssClass << "'";


  // cell begin
  ss << "D[T][R][C++]=\"" << cls.str();
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "<a href='" << link.str() << "'>";
  // process base class cell renderer
  ReportRendererBase::buildTableCell(base, row, ss);
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "</a>";
  // cell end
  ss << "\";";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
